import pandas as pd
import numpy as np
import os

os.chdir('E:\Documents\VIII Sem (Project)\File Processing Order')

A_gene_df = pd.read_csv('(1.4) Gene List.csv')

A_molecule_df = pd.read_csv('(5) Gene_drug_list.csv')
A_molecule_df.set_index('Gene', inplace = True)
A_molecule_df.drop(columns = ['InChI'], inplace = True)
A_molecule_df.rename(columns = {'Molecular Formula' : 'Formula'}, inplace = True)

A_similarity_df = pd.read_csv('(9) Similar_Molecule_List.csv')
A_similarity_df.rename(columns = {'Primary Molecule' : 'Molecule Name'}, inplace = True)

B_set_sim = set(A_similarity_df['Molecule Name'])

final_df = pd.DataFrame()

x = 0
for gene in A_gene_df['Gene']:
    if gene in A_molecule_df.index:
        pass
    else:
        continue
    
    B_set_molecule = set(A_molecule_df.loc[gene , 'Molecule Name'])
    
    C_intersect = list(B_set_molecule.intersection(B_set_sim))
    C_diff = list(B_set_molecule.difference(B_set_sim))

    
    D_sim_df = A_similarity_df[A_similarity_df['Molecule Name'].isin(C_intersect)]
    D_sim_df['Gene'] = gene
    
    D_mol_df = A_molecule_df[A_molecule_df['Molecule Name'].isin(C_diff)]
    D_mol_df.reset_index(inplace = True)
    D_mol_df = D_mol_df[D_mol_df['Gene'] == gene]
    
    E_temp_df = pd.concat([D_sim_df, D_mol_df])
    
    final_df = pd.concat([final_df, E_temp_df])
    

final_df.set_index('Gene', inplace = True)


final_df.to_csv('(10) Complete_Molecule_List.csv')
