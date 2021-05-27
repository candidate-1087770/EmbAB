
import pandas as pd
from Bio import SeqIO


class biochem_features():
    '''returns dataframes for biochemical changes and MAPP score'''
    
    def __init__(self, *args):
        '''housekeeping'''
        
        self.gene = ['A', 'B', 'C']
        
        self.aa_list = ['A','C','D','E','F','G','H','I','K','L',
                       'M','N','P','Q','R','S','T','V','W','Y']
        
        self.aa_volumes = {'A': 88.6, 'R': 173.4, 'N': 114.1, 'D': 111.1, 'C': 108.5,
                           'Q': 143.8, 'E': 138.4, 'G': 60.1, 'H': 153.2, 'I': 166.7,
                           'L': 166.7, 'K': 168.6, 'M': 162.9, 'F': 189.9, 'P': 112.7,
                           'S': 89.0, 'T': 116.1, 'W': 227.8, 'Y': 193.6, 'V': 140.0}
        
        self.aa_MW = {'A': 89.1, 'R': 174.2, 'N': 132.1, 'D': 133.1, 'C': 121.1,
                      'E': 147.1, 'Q': 146.2, 'G':75.1, 'H': 155.2, 'I': 131.2,
                      'L': 131.2, 'K': 146.2, 'M': 149.2, 'F': 165.2, 'P': 115.1,
                      'S': 105.1, 'T': 119.1, 'W': 204.2, 'Y': 181.2, 'V': 117.1}
        
        self.hydropathy_index = {'A': 1.8, 'R': -4.5, 'N': -3.5, 'D': -3.5,
                                 'C': 2.5, 'E': -3.5, 'Q': -3.5, 'G': -0.4,
                                 'H': -3.2, 'I': 4.5, 'L': 3.8, 'K': -3.9,
                                 'M': 1.9, 'F': 2.8, 'P': -1.6, 'S': -0.8, 
                                 'T': -0.7, 'W': -0.9, 'Y': -1.3, 'V': 4.2}
        
        self.Pi = {'A': 6, 'R': 10.76, 'N': 5.41, 'D': 2.77, 'C': 5.07, 'E': 3.22,
                   'Q': 5.65, 'G': 5.97, 'H': 7.59, 'I': 6.02, 'L': 5.98, 'K': 9.74,
                   'M': 5.74, 'F': 5.48, 'P': 6.3, 'S': 5.68, 'T': 5.6, 'W': 5.89,
                   'Y': 5.66, 'V': 5.96}
    
    
    def biochem_change(self):
        '''calculates change in biochemical feature for every possible mutation'''
        
        #create df from attribute function data
        df = pd.DataFrame.from_dict([self.aa_volumes, self.aa_MW, 
                                     self.hydropathy_index, self.Pi]).T.reset_index(0)
        DF = df.rename(columns = {'index':'Amino_acid', 0:'volume', 1:'MW', 
                                  2:'hydropathy',3:'Pi'})
        
        #create dicts for delta attributes from df
        d_volume, d_Pi, d_MW, d_hydropathy = {},{},{},{}
        for i in DF.index:
            for j in DF.index:
                mutation = DF['Amino_acid'][i] + DF['Amino_acid'][j]
                d_volume[mutation] = DF['volume'][i] - DF['volume'][j]
                d_hydropathy[mutation] = DF['hydropathy'][i] - DF['hydropathy'][j]
                d_Pi[mutation] = DF['Pi'][i] - DF['Pi'][j]
                d_MW[mutation] = DF['MW'][i] - DF['MW'][j]
                
        #creates df from delta dictionaries
        d_chem_df = pd.DataFrame.from_dict([d_volume, d_MW, 
                                            d_hydropathy, d_Pi]).T.reset_index(0)
        d_chem_DF = d_chem_df.rename(columns={'index':'dAA', 0:'d_volume', 
                                              1:'d_MW', 2:'d_hydropathy',3:'Pi'})
        
        return d_chem_DF
    
    
    def MAPP(self):
        '''generates df of MAPP scores for every aa combination'''
        
        #pulls data from each fasta and MAPP output files for each gene and
        #stores in dicts/df
        record_dict, file_dict = {},{}
        for i in self.gene:
            record_dict['records{0}'.format(i)] = (
                list(SeqIO.parse('../intermediate_data_files/MAPP/Emb{0}_align.fasta'.
                                 format(i), 'fasta')))
            file_dict['Emb{0}'.format(i), 'MAPP'] = (
                pd.read_excel('../intermediate_data_files/MAPP/Emb{0}_output.xlsx'.
                              format(i)))
            
        #creates df for each gene containing MAPP for every possible mutation for the first
        #sequence in the fasta alignment file
        gene_list = []
        for i in self.gene:
            Dict, count = {}, 0
            sequence = record_dict['records{0}'.format(i)][0].seq
            output = file_dict['Emb{0}'.format(i), 'MAPP']
            for g in range(0, len(sequence)):
                chcter = sequence[g]
                if chcter != '-':
                    count += 1
                    for j in self.aa_list:
                        Dict[chcter+str(count)+j] = output[j][g]
                        gene_list.append('emb{0}'.format(i))
                        
            if i == 'A':
                df_A = pd.DataFrame.from_dict(Dict, orient='index').reset_index()
            elif i == 'B':
                df_B = pd.DataFrame.from_dict(Dict, orient='index').reset_index()
            elif i == 'C':
                df_C = pd.DataFrame.from_dict(Dict, orient='index').reset_index()
            
        #adds each df to df_A, then creates master df from df_A
        df_A = df_A.append(df_B.append(df_C, ignore_index=True), ignore_index=True)
        df_A['GENE'] = gene_list
        df = df_A.rename(columns={'index':'MUTATION', 0:'MAPP'})
        
        return df

