import os
import pandas as pd
import requests
import numpy as np

os.chdir('E:\Documents\VIII Sem (Project)\Results')

gene_list = pd.read_csv('Gene List.csv')

drug_gene_list = {'Gene' : [],
                  'Drug ID' : [],
                  'Drug Name' : [],
                  'Interaction ID' : [],
                  'Interaction Type' : [],
                  'pmID' : [],
                  'Score' : [],
                  'Source' : []
                  }

for gene in gene_list['Gene']:
    DGIdb_url = 'https://dgidb.org/api/v2/interactions.json?genes=' + gene
    url_request = requests.get(DGIdb_url)
    DGI_json = url_request.json()
    DGI_list = DGI_json['matchedTerms'][0]['interactions']
    if len(DGI_list) == 0:
        drug_gene_list['Gene'].append(gene)
        drug_gene_list['Drug ID'].append(np.nan)
        drug_gene_list['Drug Name'].append(np.nan)
        drug_gene_list['Interaction ID'].append(np.nan)
        drug_gene_list['Interaction Type'].append(np.nan)
        drug_gene_list['pmID'].append(np.nan)
        drug_gene_list['Score'].append(np.nan)
        drug_gene_list['Source'].append(np.nan)
    else:
        for i in range(len(DGI_list)):
            if DGI_list[i]['score'] > 0.1:
                drug_gene_list['Gene'].append(gene)
                drug_gene_list['Drug ID'].append(DGI_list[i]['drugConceptId'][7:])
                drug_gene_list['Drug Name'].append(DGI_list[i]['drugName'])
                drug_gene_list['Interaction ID'].append(DGI_list[i]['interactionId'])
                drug_gene_list['Interaction Type'].append(DGI_list[i]['interactionTypes'])
                drug_gene_list['pmID'].append(DGI_list[i]['pmids'])
                drug_gene_list['Score'].append(DGI_list[i]['score'])
                drug_gene_list['Source'].append(DGI_list[i]['sources'])
            
        
drug_gene_df = pd.DataFrame.from_dict(drug_gene_list)

drug_list = drug_gene_df['Drug Name'].drop_duplicates()
drug_list.to_csv('Drug_List_Total_DGiDB.csv')

drug_gene_df = drug_gene_df.set_index('Drug Name')
drug_gene_df.to_csv('Gene_Drug_List_DGiDB_Total.csv')


    