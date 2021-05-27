
import pandas as pd
import numpy
import MDAnalysis
from MDAnalysis.analysis.distances import distance_array
from MDAnalysis import AtomGroup
import freesasa




class distances():
    '''Calculates distance from different reference points for each residue'''
    
    def __init__(self):
        '''housekeeping'''
        
        self.maxasa = {'ALA':129, 'ARG':274, 'ASN':195, 'ASP':193, 'CYS':167, 
                       'GLU':223, 'GLN':225, 'GLY':104, 'HIS':224, 'ILE':197, 
                       'LEU':201, 'LYS':236, 'MET':224, 'PHE':240, 'PRO':159, 
                       'SER':155, 'THR':172, 'TRP':285, 'TYR':263, 'VAL':174}
        
        
        #tb
        self.EmbAB_tb_chain_codes = {'A':'EmbA', 'B':'EmbB', 'P':'AcpM2'}
        
        self.EmbAB_tb_Universe = MDAnalysis.Universe('../pdb_files/7bvf.pdb')
        self.CA_EmbAB_tb = self.EmbAB_tb_Universe.select_atoms('name CA')
        
        #smeg:
        self.EmbAB_chain_codes = {'A':'EmbA', 'B':'EmbB', 'C':'AcpM2'}
        
        self.EmbAB_eth_Universe = MDAnalysis.Universe('../pdb_files/7bvc.pdb')
        self.EmbAB_Universe = MDAnalysis.Universe('../pdb_files/7bvg.pdb')
        self.CA_EmbAB = self.EmbAB_Universe.select_atoms('name CA')
        self.CA_EmbAB_eth =self.EmbAB_eth_Universe.select_atoms('name CA')
    
    
    def ligand_distance(self, Universe, CA):
        '''creates df with distances from CA to ethambutol in each protein'''
        eth = Universe.select_atoms('resname 95E').residues
        no = 0
        for i in eth:
            distances = distance_array(i.atoms.center_of_mass(), CA.positions)
            try:
                df['Ligand{0}_Distance'.format(no)] = distances[0]
            except NameError:
                df = pd.DataFrame(distances[no], columns = ['Ligand0_Distance'])
            no += 1
        return df

    
    def Ca_distance(self, Universe, CA):
        '''creates df with distances from CA to Ca in each protein'''
        #Ca provide strucutral stability
        Ca = Universe.select_atoms('resname Ca') # = resname calcium
        distances = distance_array(Ca.positions, CA.positions)
        df = pd.DataFrame(distances[0], columns = ['Ca_Distance'])
        return df
    
    
    def cardiolipin_distance(self, Universe, CA):
        '''creates df with distances from CA to cardiolipin'''
        #the distance to the nearest cariolipin atom should roughly
        #equal the distance to the EmbA-B interface
        cardiolipin = Universe.select_atoms('resname CDL')
        distances = distance_array(cardiolipin.positions, CA.positions)
        smallest_distance = []
        for i in distances.T:
            smallest_distance.append(min(i))
        df = pd.DataFrame(distances[0], columns = ['cardiolipin_Distance'])
        return df
        
    
    def DPA_head_distance(self, Universe, CA):
        '''creates df with distances from CA to DPA head groups'''
        #DPA = donor substrate
        #DPA head group binds in the active site
        DPA = Universe.select_atoms('resname F8L and type O')
        distances = distance_array(DPA.positions, CA.positions)
        #using min distances as oppose to com due to DPA's length
        smallest_distance = []
        for i in distances.T:
            smallest_distance.append(min(i))
        df = pd.DataFrame(smallest_distance, columns = ['DPA_Distance'])
        return df
    
    
    def diaribinose_distance(self, Universe, CA):
        '''creates df with distances from CA to diaribinose'''
        #diaribinose = acceptor substrate in EmbAB
        diaribinose = Universe.select_atoms('resname BXY')
        distances = distance_array(diaribinose.center_of_mass(), CA.positions)
        df = pd.DataFrame(distances[0], columns = ['Diaribinose_Distance'])
        return df
    
    def aribinose_distance(self, Universe, CA):
        ''' creates df with distances from CA to specific aribinose'''
        #EmbC diaribinose = represents 1 terminal residue position of acceptor 
            # + 1 binding site for donor
        List = [Universe.select_atoms('segid E and resname BXY').residues,
        Universe.select_atoms('segid J and resname BXY').residues]
        no = 0
        for i in List:
            for j in i:
                distances = distance_array(j.atoms.center_of_mass(), CA.positions)
                try:
                    df['Arabinose{0}_Distance'.format(no)] = distances[0]
                except NameError:
                    df = pd.DataFrame(distances[no], columns = ['Arabinose0_Distance'])
                no += 1
        return df

    
    def phosphate_distance(self, Universe, CA):
        '''creates df with distances from CA to each phosphorous'''
        #in EmbC the phosphate is trapped in the active site and IS the DPA phosphate
        P = Universe.select_atoms('resname PO4 and name P')
        distances = distance_array(P.positions, CA.positions)
        for i in range(0, len(P)):
            try:
                df['Phosphate{0}_Distance'.format(i)] = distances[i]
            except NameError:
                df = pd.DataFrame(distances[i], columns = ['Phosphate0_Distance'])
        return df
    
    
    def compatible_df(self, df):
        '''edits and organises df for merging with CRyPTIC data'''
        
        df['residue_no.(pdb)'] = df['pdb_residue']
        df['residue_ID(pdb)'] = df['rescode'] + df['residue_no.(pdb)'].astype(str)
              
        #sort out chain numbering to match CRyPTIC dataset pdb numbering
        gene_list = []
        df = df[df['pdb_chain'].isin(self.chain_codes.keys())].copy()
        for i in df['pdb_chain']:
            for k, v in self.chain_codes.items():
                if i == k:
                    gene_list.append(v)
        df['GENE'] = gene_list
        
        #insert genbank residue numbering
        genbank_ID_list = []
        for i in df['residue_ID(pdb)']:
            genbank_ID_list.append(i[0] + str((int(i[1:]) - 6)))
        df['residue_ID(genbank)'] = genbank_ID_list
        return df
    
    
    def compile(self, dimer):
        '''compiles methods to return list of 2 dataframes; 
        1 df for EmbAB and 1 df for EmbC'''
        
        if dimer == 'EmbAB':
            df = pd.merge(pd.merge(pd.merge(pd.merge(pd.merge(
                self.ligand_distance(self.EmbAB_eth_Universe, self.CA_EmbAB_eth),
                self.Ca_distance(self.EmbAB_Universe, self.CA_EmbAB), 
                    left_index=True, right_index=True, how='inner'),
                self.cardiolipin_distance(self.EmbAB_Universe, self.CA_EmbAB), 
                    left_index=True, right_index=True, how='inner'),
                self.DPA_head_distance(self.EmbAB_Universe, self.CA_EmbAB), 
                    left_index=True, right_index=True, how='inner'),
                self.diaribinose_distance(self.EmbAB_Universe, self.CA_EmbAB), 
                    left_index=True, right_index=True, how='inner'),
                self.phosphate_distance(self.EmbAB_Universe, self.CA_EmbAB), 
                    left_index=True, right_index=True, how='inner')
            #call in the form distances().compile('EmbAB')
            return df
        
        
        elif dimer =='EmbAB_tb':
            df = pd.merge(pd.merge(self.ligand_distance(self.EmbAB_tb_Universe,
                self.CA_EmbAB_tb), self.Ca_distance(self.EmbAB_tb_Universe, 
                self.CA_EmbAB_tb), left_index=True, right_index=True, how='outer'),
                self.cardiolipin_distance(self.EmbAB_tb_Universe, self.CA_EmbAB_tb),
                    left_index=True, right_index=True, how='outer')
            return df
        
        
    def output(self, dimer):
        '''creates df with residues and associated distances'''
        
        if dimer == 'EmbAB_tb':
            df = self.compile(dimer)
            df['AMINO_ACID_NUMBER'] = [i for i in self.CA_EmbAB_tb.residues.resnums]
            df['GENE'] = [('emb' + i) for i in self.CA_EmbAB_tb.residues.segids]
        elif dimer == 'EmbAB':
            df = self.compile(dimer)
            df['AMINO_ACID_NUMBER'] = [i for i in self.CA_EmbAB.residues.resnums]
        
        return df
        
        
    


# In[ ]:



class depth():
    '''Calculates surface depth, depth from lipid head groups, and depth from acyl tails'''
    
    def __init__(self):
        '''housekeeping'''
        
        self.res_codes = {'ALA':'A', 'ARG':'R', 'ASN':'N', 'ASP':'D', 'CYS':'C', 
                          'GLU':'E', 'GLN':'Q', 'GLY':'G', 'HIS':'H', 'ILE':'I', 
                          'LEU':'L', 'LYS':'K', 'MET':'M', 'PHE':'F', 'PRO':'P', 
                          'SER':'S', 'THR':'T', 'TRP':'W', 'TYR':'Y', 'VAL':'V'}
        
        #define max total sasa for each residue (also useful to define aa)
        self.maxasa_dict = {'ALA':129, 'ARG':274, 'ASN':195, 'ASP':193, 'CYS':167, 'GLU':223,
                            'GLN':225, 'GLY':104, 'HIS':224, 'ILE':197, 'LEU':201, 'LYS':236,
                            'MET':224, 'PHE':240, 'PRO':159, 'SER':155, 'THR':172, 'TRP':285,
                            'TYR':263, 'VAL':174}
        
        #tb
        self.EmbAB_tb_sasa_structure = freesasa.Structure('../pdb_files/7bvf.pdb')
        self.EmbAB_tb_Universe = MDAnalysis.Universe('../pdb_files/7bvf.pdb')
        self.EmbAB_tb_head_Universe = MDAnalysis.Universe('../pdb_files/7bvf_dppc-head-contacts.pdb')
        self.EmbAB_tb_tail_Universe = MDAnalysis.Universe('../pdb_files/7bvf_dppc-tail-contacts.pdb')
        self.CA_EmbAB_tb = self.EmbAB_tb_Universe.select_atoms('name CA')
        self.CA_EmbAB_tb_tail = self.EmbAB_tb_tail_Universe.select_atoms('name CA')
        self.CA_EmbAB_tb_head = self.EmbAB_tb_head_Universe.select_atoms('name CA')
        
        #smeg
        self.EmbAB_sasa_structure = freesasa.Structure('../pdb_files/7bvg.pdb')
        self.EmbAB_Universe = MDAnalysis.Universe('../pdb_files/7bvg.pdb')
        self.EmbAB_head_Universe = MDAnalysis.Universe('../pdb_files/7bvg_dppc-head-contacts.pdb')
        self.EmbAB_tail_Universe = MDAnalysis.Universe('../pdb_files/7bvg_dppc-tail-contacts.pdb')
        self.CA_EmbAB = self.EmbAB_Universe.select_atoms('name CA')
        self.CA_EmbAB_tail = self.EmbAB_tail_Universe.select_atoms('name CA')
        self.CA_EmbAB_head = self.EmbAB_head_Universe.select_atoms('name CA')

        
    def sasa_df(self, Universe, sasa_structure):
        '''calculate SASA for each residue and create df'''

        sasa_result = freesasa.calc(sasa_structure)
        SA_dict_atom = {}
        for i in range(0, sasa_structure.nAtoms()):
            SA_dict_atom[i] = sasa_result.atomArea(i)
        SASA_df = pd.DataFrame.from_dict(SA_dict_atom, orient='index', columns=['atom_sasa'])
        
        return SASA_df
        
        
    def fill_df(self, Universe, sasa_structure):
        '''label atoms and residues'''
        
        #label
        SASA_df = self.sasa_df(Universe, sasa_structure)
        SASA_df['atom'] = [sasa_structure.atomName(i).strip() for i in range(0,sasa_structure.nAtoms())]
        SASA_df['residue'] = [sasa_structure.residueName(i).strip() for i in range(0,sasa_structure.nAtoms())]
        SASA_df['resnum'] = [sasa_structure.residueNumber(i).strip() for i in range(0,sasa_structure.nAtoms())]
        SASA_df['chain'] = [sasa_structure.chainLabel(i).strip() for i in range(0,sasa_structure.nAtoms())]
        SASA_df['combined_ID'] = SASA_df['chain'] + '_' + SASA_df['resnum'] + '_' + SASA_df['residue']
        
        #ensure atoms exists within amino acids only
        for i in SASA_df.index:
            if SASA_df['residue'][i] not in self.res_codes.keys():
                SASA_df.drop(i, inplace=True)
        SASA_df.reset_index(0, inplace=True)
        
        return SASA_df
    
    
    def sum_sasa(self, Universe, sasa_structure):
        '''sum atom sasa for each residue and insert index for each residue'''
        
        SASA_df = self.fill_df(Universe, sasa_structure)
        residue_list, residue_ID_list, CA_atom_ID_list = [], [], []
        Dict = {}
        SASA_sum = 0    
        for i in SASA_df.index:
            try:
                if SASA_df['combined_ID'][i] == SASA_df['combined_ID'][i+1]:
                    #atoms in same residue, so add current sasa and keep iterating
                    SASA_sum += SASA_df['atom_sasa'][i]
                else:
                    #atoms not in same residue, so add sasa and reset SASA_sum
                    SASA_sum += SASA_df['atom_sasa'][i]
                    Dict[SASA_df['combined_ID'][i]] = SASA_sum
                    SASA_sum = 0
                    residue_list.append(SASA_df['residue'][i])
                    residue_ID_list.append(SASA_df['resnum'][i])
            except KeyError:
                #for last index of df, so add sasa and reset SASA_sum
                SASA_sum += SASA_df['atom_sasa'][i]
                Dict[SASA_df['combined_ID'][i]] = SASA_sum
                SASA_sum = 0
                residue_list.append(SASA_df['residue'][i])
                residue_ID_list.append(SASA_df['resnum'][i])
            if SASA_df['atom'][i].strip() == 'CA':
                #list indexes of CA 
                CA_atom_ID_list.append(SASA_df['index'][i])
        
        sasa_sum_df = pd.DataFrame.from_dict(Dict, orient='index', 
                columns=['SASA_residue']).reset_index(0).rename(columns={'index':'combined_ID'})
        sasa_sum_df['residue_ID'] = residue_ID_list
        sasa_sum_df['CA_ID'] = CA_atom_ID_list
        
        return sasa_sum_df
        
        
    def relative_sasa(self, Universe, sasa_structure):
        '''calculate relative SASA using max sasa'''
        
        df = self.sum_sasa(Universe, sasa_structure)
        #insert predefined max sasa for each residue
        maxasa_list = []
        for i in df['combined_ID']:
            for key, value in self.maxasa_dict.items():
                if i[-3:] == key:
                    maxasa_list.append(value)
        df['maxASA'] = maxasa_list
        #calculate and insert relative SASA for each residue
        df['RSA'] = df['SASA_residue']/df['maxASA']
        #filter for residues with SASA > 0.25 -- NEED TO FIND PAPER THAT RECOMMENDS GREATER THAT 0.25
        surface_residues_df = df[df['RSA'] >= 0.26].reset_index(0)
        
        return surface_residues_df
    
    
    def distances(self, Universe, CA_Universe, sasa_structure):
        '''calculate distance array and construct distance df'''
        
        df = self.relative_sasa(Universe, sasa_structure)
        sa_CA = AtomGroup(df['CA_ID'], Universe)
        distances = distance_array(CA_Universe.positions, sa_CA.positions)
        distance_df = pd.DataFrame([min(i) for i in distances], columns=['Depth'])
        distance_df['pdb_residue'] = [i.resid for i in CA_Universe]
        distance_df['pdb_chain'] = [str(i.segment)[9] for i in CA_Universe]
        
        #convert resnamaes to rescode IDs
        resnames, rescodes = [i.resname for i in self.CA_EmbAB], {}
        for i in resnames:
            for k,v in rescodes.items():
                if i == k:
                    rescodes.append(v)
        
        return distance_df

    
    def lipid_distance(self, h_or_t, pdb, universe):
        '''returns list of distances from nearest lipid head or tail (h_or_t) group;
           substitute 'head' or 'tail' for 'h_or_t' '''
        
        with open('../pdb_files/{0}_dppc-{1}-contacts.txt'.format(pdb, h_or_t), 'r') as file:
            ht_list = [i[61:66] for i in file.readlines()[5:] if i[13:16]=='CA ']
        #convert head contacts into binary if there is sufficient contact (0.2 = arbitrary)
        binary = [(1 if float(i.strip()) >= 0.2 else 0) for i in ht_list]
        
        atom_list = [i for i in universe.atoms.names]
        df_1 = pd.DataFrame(atom_list, index=[i for i in range(1, len(atom_list)+1)])
        df_1['resnum'] = universe.atoms.resnums
        df_2 = pd.DataFrame(binary, index=[i for i in range(1, len(binary)+1)])
        df = pd.merge(df_1, df_2, how='outer', right_index=True, left_on=['resnum'])
        df = df[(df['0_y']==1) & (df['0_x']=='CA')]
        
        all_CA = universe.select_atoms('name CA')
        ht_CA = AtomGroup(df.index, universe)
        distances = distance_array(all_CA.positions, ht_CA.positions)
        distance_list = [min(i) for i in distances]
        
        return distance_list
        
    
    def compile(self, pdb):
        '''compiles and returns depth df'''
        
        if pdb == '7bvg':
            df = self.distances(self.EmbAB_Universe, self.CA_EmbAB, self.EmbAB_sasa_structure)
            df['lipid_head_dis'] = self.lipid_distance('head', pdb, self.EmbAB_head_Universe)
            df['lipid_tail_dis'] = self.lipid_distance('tail', pdb, self.EmbAB_tail_Universe)
        elif pdb == '7bvf':
            df = self.distances(self.EmbAB_tb_Universe, self.CA_EmbAB_tb, self.EmbAB_tb_sasa_structure)
            df['lipid_head_dis'] = self.lipid_distance('head', pdb, self.EmbAB_tb_head_Universe)
            df['lipid_tail_dis'] = self.lipid_distance('tail', pdb, self.EmbAB_tb_tail_Universe)
        return df
    
    
    def output(self, pdb):
        '''creates df with residues and associated depths'''
        
        df = self.compile(pdb)
        if pdb == '7bvg':
            df['AMINO_ACID_NUMBER'] = [i for i in self.CA_EmbAB.residues.resnums]
        elif pdb == '7bvf':
            df['AMINO_ACID_NUMBER'] = [i for i in self.CA_EmbAB_tb.residues.resnums]
            df['GENE'] = [('emb' + i) for i in self.CA_EmbAB_tb.residues.segids]
        return df
            

    
    
    


# In[ ]:


import subprocess

class secondary_structure():
    
    def __init__(self, pdb_file):
        '''run mkdssp and housekeeping'''
        
        subprocess.run(['mkdssp', '-i', pdb_file, '-o', '../intermediate_data_files/EmbAB_ss.txt'])
       
        self.text_file = '../intermediate_data_files/EmbAB_ss.txt'
        
        
    def pull_data(self):
        '''pull data from text_file and create df'''
        
        structure_list, pdb_residue, chain_list, rescode_list = [],[],[],[]
        with open(self.text_file) as filename:
            lines = filename.readlines()
            for i in range(0, len(lines)):
                if lines[i][2] == '#':
                    for j in lines[i+1:]:
                        structure_list.append(j[16])
                        pdb_residue.append(j[6:10].strip())
                        chain_list.append(j[11])
                        rescode_list.append(j[13])
                    break
                              
        df = pd.DataFrame({'Structure':structure_list, 'AMINO_ACID_NUMBER':pdb_residue, 
                           'GENE':[('emb'+i) for i in chain_list], 'rescode':rescode_list})
        df = df[df['AMINO_ACID_NUMBER'] != ''].copy()
        df['AMINO_ACID_NUMBER'] = [int(i) for i in df['AMINO_ACID_NUMBER']]
        
        return df
        
        
    def binary_output(self):
        '''convert secondary structure into binary readouts and return df'''
        
        df = self.pull_data()
        List = ['H', 'B', 'E', 'G', 'I', 'T', 'S', ' ']
        Dict = {i: [] for i in List}
            
        for i in df['Structure']:
            for j in List:
                if i == j:
                    Dict[j].append('1')
                else:
                    Dict[j].append('0')
                    
        for k,v in Dict.items():
            df[k] = v
        df.rename(columns={' ':'NaN'}, inplace=True)

        return df        
    
        

