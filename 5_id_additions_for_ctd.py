import pubchempy as pcp
from chembl_webresource_client.new_client import new_client
molecule = new_client.molecule
import pandas as pd
import numpy as np
import os

os.chdir('E:\Documents\VIII Sem (Project)\File Processing Order')

df = pd.read_csv('(3) Gene_Drug_List_CTD.csv')
x = 0
for i in range(len(df)):
    name = df['ChemicalName'][i]
    try:
        compound = pcp.get_compounds(name, 'name')[0]
        df.loc[i , 'SMILES'] = compound.canonical_smiles
        df.loc[i , 'InChI'] = compound.inchi
        df.loc[i , 'Mass'] = compound.exact_mass
        df.loc[i , 'Molecular Formula'] = compound.molecular_formula
        df.loc[i , 'Charge'] = compound.charge
    except:
        try:
            compound = molecule.get(name)
            df.loc[i , 'SMILES'] = m1['molecule_structures']['canonical_smiles']
            df.loc[i , 'InChI'] = m1['molecule_structures']['standard_inchi']
            df.loc[i , 'Molecular Formula'] = m1['molecule_properties']['full_molformula']
            df.loc[i , 'Mass'] = m1['molecule_properties']['full_mwt']
        except:
            pass

df.to_csv('(3) Gene_Drug_List_CTD.csv')