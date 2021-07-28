from chembl_webresource_client.new_client import new_client
import pandas as pd
import numpy as np
import os

os.chdir('E:\Documents\VIII Sem (Project)\File Processing Order')

df = pd.read_csv('(2) Gene_Drug_List_DGiDB.csv')

for i in range(len(df)):

    try:
        molecule = new_client.molecule
        m1 = molecule.get(df['Drug ID'][i])
        df.loc[i , 'SMILES'] = m1['molecule_structures']['canonical_smiles']
        df.loc[i , 'InChI'] = m1['molecule_structures']['standard_inchi']
        df.loc[i , 'Molecular Formula'] = m1['molecule_properties']['full_molformula']
        df.loc[i , 'Mass'] = m1['molecule_properties']['full_mwt']
    except:
        pass


df.to_csv('(2) Gene_Drug_List_DGiDB.csv')
    
    


