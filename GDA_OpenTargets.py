import os 
import requests
from opentargets import OpenTargetsClient
import pandas as pd


pathway_list = ['Base Excision Repair', 'Nucleotide Excision Repair', 'DNA Damage Bypass',
                'DNA Damage Reversal', 'HDR through Homologous Recombination (HRR)',
                'HDR through Single Strand Annealing (SSA)', 'HDR through MMEJ (alt-NHEJ)',
                'Nonhomologous End-Joining (NHEJ)', 'Fanconi Anemia Pathway', 'Mismatch Repair']

def open_targets(Gene_Ensembl_ID):
    ot = OpenTargetsClient()
    gene_disease_list = []
    Gene_Disease_associations = ot.get_associations_for_target(Gene_Ensembl_ID)
    for Disease in Gene_Disease_associations:
        temp_list = []
        if Disease['association_score']['overall'] >= 0.5:
            temp_list.append(Disease['id'])
            temp_list.append(Disease['association_score']['overall'])
            temp_list.append(Disease['disease']['id'])
            temp_list.append(Disease['disease']['efo_info']['label'])
            temp_list.append(Disease['disease']['efo_info']['therapeutic_area']['labels'])
            gene_disease_list.append(temp_list)
    gene_disease_df = pd.DataFrame(gene_disease_list, columns = ['ID','Score','Disease ID',
                                                      'Disease Name','Therapeutic Area'])
    return(gene_disease_df)

os.chdir('E:\Documents\VIII Sem (Project)\Gene-Disease Association Network')


GDA_dict = {'Pathway' : [], 'Gene' : [], 'ID' : [], 'Disease ID' : [], 'Disease Name' : [],
                                          'Therapeutic Area' : [], 'Score' : []}

for pathway in pathway_list:
    gene_data = pd.read_csv(pathway + '.csv')
    for i in range(len(gene_data['Ensembl ID'])):
        temp_GDA = open_targets(gene_data['Ensembl ID'][i])
        for j in range(len(temp_GDA)):
            for k in range(len(temp_GDA['Therapeutic Area'][j])):
                GDA_dict['Pathway'].append(pathway)
                GDA_dict['Gene'].append(gene_data['Gene'][i])
                GDA_dict['ID'].append(temp_GDA['ID'][j])
                GDA_dict['Disease ID'].append(temp_GDA['Disease ID'][j])
                GDA_dict['Disease Name'].append(temp_GDA['Disease Name'][j])
                GDA_dict['Therapeutic Area'].append(temp_GDA['Therapeutic Area'][j][k])
                GDA_dict['Score'].append(temp_GDA['Score'][j])

open_targets_GDA = pd.DataFrame(GDA_dict)
 
open_targets_GDA.to_csv('DNA Damage Repair GDA_OpenTargets.csv')