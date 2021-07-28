import pubchempy as pcp
import pandas as pd
import numpy as np
import os

os.chdir('E:\Documents\VIII Sem (Project)\File Processing Order')

df = pd.read_csv('(4) Gene_Drug_list_PubChem.csv')

for i in range(len(df)):
    try:
        cid = int(df['cid'][i])
        compound = pcp.get_compounds(cid, 'cid')[0]
    except:
        continue
    
    df.loc[i , 'SMILES'] = compound.canonical_smiles
    df.loc[i , 'InChI'] = compound.inchi
    df.loc[i , 'Charge'] = compound.charge
    df.loc[i , 'Mass'] = compound.exact_mass
    df.loc[i , 'Formula'] = compound.molecular_formula
    
    
#df.to_csv('(4) Gene_Drug_list_PubChem.csv')