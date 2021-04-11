# -*- coding: utf-8 -*-
"""
Created on Sun Apr  4 19:15:17 2021

@author: VISHAKMADHWARAJKADAM
"""
import os
import pandas as pd
import numpy as np

os.chdir('E:\Documents\VIII Sem (Project)\Gene-Disease Association Network')

pathway_list = ['Nucleotide Excision Repair', 'HDR through Homologous Recombination (HRR)',
                'HDR through Single Strand Annealing (SSA)', 'Fanconi Anemia Pathway',
                'Mismatch Repair']

disease_class = ['Congenital, Hereditary, and Neonatal Diseases and Abnormalities',
                 'Nutritional and Metabolic Diseases', 'Skin and Connective Tissue Diseases',
                 'Neoplasms']

######### REMOVING PATHWAYS AND IT'S ASSOCIATED GENE AND DISEASES 
######### WHICH ARE NOT TO BE CONSIDERED ##########

GDA_disgenet = pd.DataFrame(columns = ['Pathway', 'Gene', 'Disease Name', 
                                       'Disease Class','Disease ID', 'Score'])

for pathway in pathway_list:
    disgenet = pd.read_csv(pathway + ' GDA_DisGeNet.csv')
    disgenet = disgenet[disgenet['Disease Class'].isin(disease_class)]
    disgenet = disgenet.drop_duplicates(subset = 'Disease Name')

    disgenet = disgenet.drop('Unnamed: 0', axis = 1)
    GDA_disgenet = GDA_disgenet.append(disgenet)
    
GDA_disgenet = GDA_disgenet.set_index('Pathway')

GDA_disgenet.to_csv('DNA Damage Repair GDA_DisGeNet.csv')

######### CREATING A TABLE CONTAINING METADATA ABOUT THE CHOSEN PATHWAYS ##########

GDA_disgenet = pd.read_csv('DNA Damage Repair GDA_DisGeNet.csv')
GDA_disgenet = GDA_disgenet.set_index('Pathway')
gene_metadata = pd.DataFrame(columns = ['Pathway', 'UniProt ID', 'Gene',
                                        'Ensembl ID', 'Entrez ID', 'PharmGKB ID', 'Alias'
                                        , 'PDB'])

for pathway in pathway_list:
    gene = pd.read_csv(pathway + '.csv')
    GDA_pathway = GDA_disgenet.loc[pathway]
    pathway_gene = GDA_pathway['Gene']
    pathway_gene = list(pathway_gene.drop_duplicates())
    gene = gene[gene['Gene'].isin(pathway_gene)]
    gene.insert(0, 'Pathway', pathway)
    gene_metadata = gene_metadata.append(gene)

gene_metadata = gene_metadata.drop('Unnamed: 0', axis = 1)
gene_metadata.to_csv('Gene Metadata.csv')

########## FINDING GENES PRESENT IN ALL PATHWAYS ##########

gene = pd.read_csv('Gene Metadata.csv')
gene = gene['Gene'].drop_duplicates()
gene.to_csv('Gene List.csv')

########## FINDING THE LIST OF ALL DRUGS RELATED TO THE GENES ##########

gene_drug_data = pd.read_csv('Gene Drug Association.csv')
drug_list = gene_drug_data['Drug Name'].drop_duplicates()
drug_list.to_csv('Gene Drug List.csv')


temp_tsv = pd.read_csv('DNA Damage Repair PPI.tsv', delimiter = '\t')
temp_tsv.to_csv('DNA Damage Repair PPI.csv')
