import pandas as pd
import numpy as np
import os

os.chdir('E:/Documents/VIII Sem (Project)/Results/Drug_Gene_csvfiles')

ctd_df = pd.DataFrame()
gene_list = pd.read_csv('Gene List.csv')
for gene in gene_list['Gene']:
    temp_df = pd.read_csv(gene + '_drug_list_ctd.csv')
    ctd_df = pd.concat([ctd_df , temp_df])

ctd_df.to_csv('Gene_Drug_List_CTD.csv')