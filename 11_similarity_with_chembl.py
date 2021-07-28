from chembl_webresource_client.new_client import new_client
import pandas as pd
import numpy as np
import os

os.chdir('E:\Documents\VIII Sem (Project)\File Processing Order')

df = pd.read_csv('(5) Gene_drug_list.csv')
similarity = new_client.similarity

final_list = []
for i in range(len(df)):
    smiles_1 = df['SMILES'][i]
    
    try:
        compounds_list = similarity.filter(smiles = smiles_1, similarity = 95)
        for molecule in compounds_list:
            try:
                temp_list = []
                temp_list.append(df['Molecule Name'][i])
                temp_list.append(molecule['molecule_chembl_id'])
                temp_list.append(molecule['molecule_structures']['canonical_smiles'])
                temp_list.append(molecule['molecule_structures']['standard_inchi'])
                temp_list.append(molecule['molecule_properties']['full_molformula'])
                temp_list.append(molecule['molecule_properties']['full_mwt'])
                final_list.append(temp_list)
            except:
                pass
                
    except:
        pass
    
temp_df = pd.DataFrame(data = final_list, columns = ['Primary Molecule',
                                                    'Similar Molecule',
                                                    'SMILES',
                                                    'InChI',
                                                    'Molecular Formula',
                                                    'Mass'])

temp_df = temp_df.astype({'Mass' : 'float32'})

temp_df = temp_df[temp_df['Mass'] <= 500.0]

temp_df.to_csv('(8) Similar_Molecules_ChEMBL.csv')
