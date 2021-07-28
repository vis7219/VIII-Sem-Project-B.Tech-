import os
import pandas as pd
import numpy as np

os.chdir('E:\Documents\VIII Sem (Project)\File Processing Order\(10.2) SwissADME Files')
file_list = os.listdir()
file_list.sort(key=os.path.getmtime)

file_list.remove('(10.1) Complete_Molecule_list_start_removed .csv')

i = 0
final_df = pd.DataFrame()
for file in file_list:
    file_df = pd.read_csv(file)
    final_df = pd.concat([final_df , file_df])
    #print(final_df.shape)
    i = i + 1


df = pd.read_csv('(10.1) Complete_Molecule_list_start_removed .csv')
final_df.reset_index(inplace = True)

for i in range(0 , 4368):
    molecule = 'Molecule' + str(i + 1)
    df_molecule = str(final_df.loc[i , 'Molecule'])
    if df_molecule != molecule:
        break

temp_df = df.loc[100 : 119 , 'SMILES']
temp_df = list(temp_df)

#final_df.to_csv('(11) SwissADME.csv')
    

