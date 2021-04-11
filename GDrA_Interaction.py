import pandas as pd
import requests
import os
import numpy as np

os.chdir('E:\Documents\VIII Sem (Project)\Gene-Disease Association Network')


def DGIdb_gene_drug(Gene_name):
    DGIdb_url = 'https://dgidb.org/api/v2/interactions.json?genes=' + Gene_name
    url_request = requests.get(DGIdb_url)
    DGI_json = url_request.json()
    try:
        Drug_gene_list = DGI_json['matchedTerms'][0]['interactions']
        #Drug_gene_df = pd.DataFrame(Drug_gene_list)
        return(Drug_gene_list)
    except:
        Drug_gene_df = np.nan
        return(Drug_gene_df)

gene_metadata = pd.read_csv('Gene Metadata.csv')
gene_metadata = gene_metadata['Gene']
gene_metadata = list(gene_metadata.drop_duplicates())

gene_drug_dict = {'Gene': [], 'Drug Name': [], 'Score': [],
                  'Interaction ID': [], 'Interaction Type': [], 
                  'Source': [],'PMIDs': []}

for gene in gene_metadata:
    gene_drug_metadata = DGIdb_gene_drug(gene)
    if gene_drug_metadata != np.nan:
        for drug in gene_drug_metadata:
            if drug['score'] >= 0.1:
                gene_drug_dict['Gene'].append(gene)
                gene_drug_dict['Interaction ID'].append(drug['interactionId'])
                gene_drug_dict['Interaction Type'].append(drug['interactionTypes'])
                gene_drug_dict['Drug Name'].append(drug['drugName'])
                gene_drug_dict['Score'].append(drug['score'])
                gene_drug_dict['Source'].append(drug['sources'])
                gene_drug_dict['PMIDs'].append(drug['pmids'])
            
gene_drug_pd = pd.DataFrame(gene_drug_dict)
gene_drug_pd.to_csv('Gene Drug Association.csv')
        
        
        
        
        
        
        
        
        