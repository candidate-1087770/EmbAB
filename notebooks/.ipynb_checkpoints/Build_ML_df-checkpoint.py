

import pandas as pd

import sys
sys.path.insert(0, '../notebooks')



from thermo_features import *
from Structural_attributes import *
from Amino_acid_attributes import *



class df_for_ML:
    '''create df ready for ML study from clinical data + structural features'''
    
    def __init__(self):
        '''housekeeping and import clinical + structural data'''
        
        #clinical data
        mutation_data = pd.read_csv('../base_data/EMB_MUTATIONS.csv.gz', low_memory=False)
        self.mutation_data = mutation_data.drop_duplicates()
        phenotype_data = pd.read_csv('../base_data/EMB_PHENOTYPES.csv.gz', low_memory=False)
        self.phenotype_data = phenotype_data.drop_duplicates()

                
    def filter_mutations(self):
        '''filters data for non-synonymous solos'''
        
        #coding region + remove synonymous
        mutation_df = self.mutation_data[(self.mutation_data['IN_CDS'])&(~self.mutation_data['IS_SYNONYMOUS'])]
        #solos
        mutation_df = mutation_df.drop_duplicates(subset=['UNIQUEID'], keep=False)
        #no indels and filter pass
        mutation_df = mutation_df[(~mutation_df['IS_INDEL'])&(mutation_df['IS_FILTER_PASS'])]
        #no unknowns
        List = [i for i in mutation_df['MUTATION'] if i[-1] == 'X']
        mutation_df = mutation_df[~mutation_df['MUTATION'].isin(List)]
        #only EmbAB and EmbC (not EmbR)
        mutation_df = mutation_df[mutation_df['GENE']!='embR'].copy()
        
        return mutation_df
    
    
    def filter_phenotypes(self):
        '''filters data for high quality data'''
        
        self.phenotype_df = self.phenotype_data[(self.phenotype_data['PHENOTYPE_QUALITY']!='LOW')&
                                                (self.phenotype_data['BINARY_PHENOTYPE'].notna())]
        return self.phenotype_df
    
    
    def merge_clinical(self):
        '''merge phenotype and mutation data + filter for necessary data only + clarify numbering'''
        
        self.mutation_df = self.filter_mutations()
        self.phenotype_df = self.filter_phenotypes()
        
        #merge mutation and phenotype datasets
        merged = pd.merge(self.mutation_df, self.phenotype_df, on=['UNIQUEID'], how='inner')
        #filter for necessary columns only
        Merged = merged[['UNIQUEID', 'GENE', 'MUTATION', 
                         'BINARY_PHENOTYPE', 'MIC', 'AMINO_ACID_NUMBER']].copy()
        #insert pdb mutation numbering into df
        Merged['residue_ID(pdb)'] = [(Merged['MUTATION'][i][:1] + 
                                      str(int(Merged['MUTATION'][i][1:-1]))) for i in Merged.index]
        #insert genbank mutaiton number into df
        Merged['residue_ID(genbank)'] = [(str(i[:1]) + str(int(i[1:])-6)) for i in Merged['residue_ID(pdb)']]
        #create aa mutaiton codes
        Merged['dAA']  = [(i[0]+i[-1]) for i in Merged['MUTATION']]
        #insert column for 'chain' as in the pdb file       
        List = ['A' if i == 'embA' else 'B' if i == 'embB' else ['A','B'] for i in Merged['GENE']]
        Merged['pdb_chain'] = List

        self.Merged = Merged
        return self.Merged
        
    
    def merged_structural(self):
        '''merge clinical data with strcutural features'''
        
        #split clinical data into EmbAB and EmbC
        clinical = self.merge_clinical()
        EmbAB = clinical[(clinical['GENE']=='embA') | (clinical['GENE']=='embB')].reset_index(0, drop=True)
        EmbC = clinical[(clinical['GENE']=='embC')].reset_index(0, drop=True)
        
        #remove EmbAB mutations at residues below 24 in chain B (not in pdb structure)
        EmbAB = EmbAB[EmbAB['AMINO_ACID_NUMBER']>=24].copy()
        #run distance feature pulldown for tb EmbAB
        EmbAB_distances = distances().output('EmbAB_tb')
        EmbAB_depth = depth().output('7bvf')
        EmbAB_ss = secondary_structure('../pdb_files/7bvf').binary_output()
        merged = pd.merge(EmbAB, EmbAB_distances, on=['AMINO_ACID_NUMBER', 'GENE'], how='inner')
        merged = pd.merge(merged, EmbAB_depth, on=['AMINO_ACID_NUMBER', 'GENE'], how='inner')
        merged = pd.merge(merged, EmbAB_ss, on=['AMINO_ACID_NUMBER', 'GENE'], how='inner')
        #build foldx model
        dG_df = stability(EmbAB, '7bvf').dG()
        merged = pd.merge(merged, dG_df, on=['MUTATION', 'GENE'], how='inner')
        #load in aa data
        d_aa = biochem_features().biochem_change()
        merged = pd.merge(merged, d_aa, on=['dAA'], how='inner')
        MAPP = biochem_features().MAPP()
        merged = pd.merge(merged, MAPP, on=['GENE', 'MUTATION'], how='inner')
        
        return merged
        

        
        
                                      
    
        
        
        
    

