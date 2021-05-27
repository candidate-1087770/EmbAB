
import pandas as pd
import re
import os.path
import sys


class stability:
    '''If individual list files are not present, creates files.
       If individual list files are present, but Dif files are not, tells user to run foldx.
       If Dif files are present, pulls out dG'''

    def __init__(self, df, PDB):
        '''housekeeping'''
        
        self.ML_df = df
        self.pdb = PDB 
        
        
    def individual_list(self, pdb):
        '''creates individual list files'''
        
        f = open(f'../intermediate_data_files/foldx_thermo/individual_list_{self.pdb}.txt', 'w')
        List = []
        if len(self.ML_df['pdb_chain'][0])==1:
            for i in self.ML_df.index:
                foldx_mutation = (self.ML_df['MUTATION'][i][:1] + self.ML_df['pdb_chain'][i] +
                                  str(int(int(self.ML_df['AMINO_ACID_NUMBER'][i]))) + 
                                  self.ML_df['MUTATION'][i][-1])
                if foldx_mutation not in List:
                    f.write(foldx_mutation + ';\n')
                List.append(foldx_mutation)
        elif len(self.ML_df['pdb_chain'][0])>1:
            #there is more than one identical chain - therefore multimer (e.g EmbCC)
            for i in self.ML_df.index:
                for j in range(0, len(self.ML_df['pdb_chain'][i])):
                    foldx_mutation = (self.ML_df['MUTATION'][i][:1] + self.ML_df['pdb_chain'][i][j] +
                                      str(int(int(self.ML_df['AMINO_ACID_NUMBER'][i]))) + 
                                      self.ML_df['MUTATION'][i][-1])
                    if (foldx_mutation not in List) & (j<len(self.ML_df['pdb_chain'][0])-1):
                        f.write(foldx_mutation + ',')
                        List.append(foldx_mutation)
                    elif (foldx_mutation not in List) & (j==len(self.ML_df['pdb_chain'][0])-1):
                        f.write(foldx_mutation + ';\n')
                        List.append(foldx_mutation)
            
            
    def dG(self):
        '''pulls dG stability for each mutation from Dif files and checks if foldx has been run'''
        
        if (os.path.isfile(f'../intermediate_data_files/foldx_thermo/Dif_{self.pdb}_ResiduesOnly_Repair.fxout')) == False: 
            self.individual_list(f'../pdb_files/{self.pdb}_Repair.pdb')
            print ('need to run foldx on individual_list files')
            sys.exit()
        
        Dict = {}
        with open(f'../intermediate_data_files/foldx_thermo/individual_list_{self.pdb}.txt', 'r') as f:
            with open(f'../intermediate_data_files/foldx_thermo/Dif_{self.pdb}_ResiduesOnly_Repair.fxout', 'r') as g:
                F, G = f.readlines(), g.readlines()
                for i in range(0, len(F)):
                    #after line number 108 (7bvf_ResiduesOnly_Repair_99), columns shift right by several spaces
                    if (9+i) < 108:
                        #create dictionary with mutations
                        Dict[F[i][:-2]] = re.sub(r'[\n\t\s]*', '', (G[9+i][30:38].strip('b').strip('.').strip('\t')))
                    else:
                        #create dictionary with mutations
                        Dict[F[i][:-2]] = re.sub(r'[\n\t\s]*', '', (G[9+i][32:40].strip('b').strip('.').strip('\t')))
                        
                    
        #create df from Dict
        df = (pd.DataFrame.from_dict(Dict, orient='index').reset_index(0).
            rename(columns={'index':'fx_mutation', 0:'dG_stability'}))    
        
        chain_list, mutation_list = [], []
        for i in df['fx_mutation']:
            chain_list.append(i[1])
            if i[1] == 'C':
                mutation_list.append(str(i[0]) + str(int(i[2:-1])-6) + i[-1])
            else:
                mutation_list.append(str(i[0]) + str(i[2:]))
        df['pdb_chain_fx'] = chain_list
        df['MUTATION'] = mutation_list 
        df['GENE'] = [('emb' + i) for i in df['pdb_chain_fx']]
        return df
        
            
    
            
        
    

