import pandas as pd
import numpy as np
import os

os.chdir('E:\Documents\VIII Sem (Project)\File Processing Order')


chebi = pd.read_csv('(6) Similar_Molecules_ChEBI.csv')
chebi.drop(columns = ['Unnamed: 0', 'CHEBI', 'InChI'], axis = 1, inplace = True)
chebi.drop_duplicates(subset = ['Primary Molecule', 'SMILES'], inplace = True)
chebi.dropna(subset = ['SMILES'], inplace = True)

pubchem = pd.read_csv('(7) Similar_Molecules_Pubchem.csv')
pubchem.drop(columns = ['cid'], axis = 1, inplace = True)
pubchem.drop_duplicates(subset = ['Primary Molecule', 'SMILES'], inplace = True)
pubchem.dropna(subset = ['SMILES'], inplace = True)

chembl = pd.read_csv('(8) Similar_Molecules_ChEMBL.csv')
chembl.drop(columns = ['Unnamed: 0', 'InChI'], inplace = True)
chembl.drop_duplicates(subset = ['Primary Molecule', 'SMILES'], inplace = True)
chembl.rename(columns = {'Molecular Formula' : 'Formula'}, inplace = True)
chembl.dropna(subset = ['SMILES'], inplace = True)
chembl['Tanimoto Score'] = 0.95

df = pd.concat([chebi, chembl, pubchem])
df.drop_duplicates(subset = ['Primary Molecule', 'SMILES'], inplace = True)
#df.set_index('Primary Molecule', inplace = True)
#df.to_csv('(9) Similar_Molecule_List.csv')