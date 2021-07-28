import pubchempy as pcp
import pandas as pd
import os

os.chdir('E:\Documents\VIII Sem (Project)\File Processing Order')


df = pd.read_csv('(5) Gene_drug_list.csv')
final_list = []

for i in range(len(df)):
    smiles = df['SMILES'][i]
    try:
        compounds_list = pcp.get_compounds(smiles, 'smiles',
                                           searchtype = 'similarity',
                                           listkey_count = 50,
                                           Threshold = 95)
    except:
        pass
    
    if type(compounds_list) == list:
        for compounds in compounds_list:
            if float(compounds.exact_mass) <= 500.0:
                if int(compounds.charge) == 0:
                    temp_list = []
                    temp_list.append(df['Molecule Name'][i])
                    temp_list.append(compounds.iupac_name)
                    temp_list.append(compounds.cid)
                    temp_list.append(compounds.canonical_smiles)
                    temp_list.append(0.95)
                    temp_list.append(compounds.molecular_formula)
                    temp_list.append(compounds.exact_mass)
                    temp_list.append(compounds.charge)
                    final_list.append(temp_list)
                else:
                    continue
            else:
                continue
    
    elif type(compounds_list) != list:
        continue
    

    
similarity_df = pd.DataFrame(data = final_list, columns = ['Primary Molecule',
                                                           'Similar Molecule',
                                                           'cid',
                                                           'SMILES',
                                                           'Tanimoto Score',
                                                           'Formula',
                                                           'Mass',
                                                           'Charge'])

similarity_df.to_csv('(7) Similar_Molecules_Pubchem.csv')

