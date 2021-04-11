import pandas as pd
import numpy as np
import requests
import ast
import os
import time
import json

os.chdir('E:\Documents\VIII Sem (Project)\Gene-Disease Association Network')


pathway_list = ['Base Excision Repair', 'Nucleotide Excision Repair', 'DNA Damage Bypass',
                'DNA Damage Reversal', 'HDR through Homologous Recombination (HRR)',
                'HDR through Single Strand Annealing (SSA)', 'HDR through MMEJ (alt-NHEJ)',
                'Nonhomologous End-Joining (NHEJ)', 'Fanconi Anemia Pathway', 'Mismatch Repair']


STRING_dict = {'Pathway' : [],'Gene A' : [], 'Gene B' : [], 'Gene A ID' : [], 'Gene B ID' : [],
               'Score' : []}

for pathway in pathway_list:
    Gene_data = pd.read_csv('Gene List.csv')
    for gene in Gene_data['Gene']:
        STRING_API = 'https://string-db.org/api/json/network?identifiers=' + gene
        STRING_request = requests.get(STRING_API)
        
        if STRING_request.status_code == 200:
            STRING_request = STRING_request.content
            PPI_data = ast.literal_eval(STRING_request.decode('UTF-8'))
            time.sleep(1)
            for interaction in PPI_data:
                if interaction['score'] >= 0.700:
                    STRING_dict['Pathway'].append(pathway)
                    STRING_dict['Gene A'].append(interaction['preferredName_A'])
                    STRING_dict['Gene B'].append(interaction['preferredName_B'])
                    STRING_dict['Gene A ID'].append(interaction['stringId_A'])
                    STRING_dict['Gene B ID'].append(interaction['stringId_B'])
                    STRING_dict['Score'].append(interaction['score'])

PPI_pd = pd.DataFrame(STRING_dict)
PPI_pd.to_csv('DNA Damage Repair PPI_STRING.csv')            
            