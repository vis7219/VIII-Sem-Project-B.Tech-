import pypdb
import requests
from IPython.display import HTML
from pypdb import *
import pandas as pd
import numpy as np 
from opentargets import OpenTargetsClient
import time
import os


def Pathway_Protein_Extraction(input_pathway_code):
    pathway = str(input_pathway_code)
    pathway_participants = 'https://reactome.org/ContentService/data/participants/' + pathway
    pathway_proteins = requests.get(pathway_participants)
    pathway_proteins = pathway_proteins.json()
    proteins = []
    count = -1
    for i in range (len(pathway_proteins)):
        for j in range(len(pathway_proteins[i]['refEntities'])):
            count = count + 1
            if pathway_proteins[i]['refEntities'][j]['schemaClass'] == 'ReferenceMolecule':
                count = count - 1
                continue
            proteins.append(pathway_proteins[i]['refEntities'][j]['displayName'])
            proteins[count] = proteins[count].replace('UniProt:','')
    proteins = list(set(proteins))
    for i in range(len(proteins)):
        proteins[i] = proteins[i].split(' ')
        for j in range(len(proteins[i][0])-1):
            if proteins[i][0][j] == '-':
                isoform_acc_rem = slice(j)
                temp_prot = proteins[i][0]
                temp_prot = temp_prot[isoform_acc_rem]
                proteins[i][0] = temp_prot
                continue
    return (proteins)

def Annotation_Extraction(proteins):
    # This function finds the Ensemble ID, Entrez ID, Pharmgkb ID, different names of the gene
    # and all its possible protein structure.
    for i in range(len(proteins)):
        annote_prot = proteins[i][0]
        # A query is made to mygene.info to get the protein's entrez ID using UniProt ID
        mygeneinfo_query = 'https://mygene.info/v3/query?q='+ annote_prot +'&species=human'
        mygeneinfo_query_metadata = requests.get(mygeneinfo_query)
        mygeneinfo_query_metadata = mygeneinfo_query_metadata.json()
        entrez_id = mygeneinfo_query_metadata['hits'][0]['entrezgene']
        # A query is made to find the other IDs using the entre ID
        mygeneinfo_annotation = 'https://mygene.info/v3/gene/' + entrez_id
        annotation_metadata = requests.get(mygeneinfo_annotation)
        annotation_metadata = annotation_metadata.json()
        if type(annotation_metadata['ensembl']) == list:
            ensemble_id = annotation_metadata['ensembl'][0]['gene']
        else:
            ensemble_id = annotation_metadata['ensembl']['gene']
        proteins[i].append(ensemble_id)
        proteins[i].append(entrez_id)
        pharmgkb_id = annotation_metadata['pharmgkb']
        proteins[i].append(pharmgkb_id)
        if 'alias' in annotation_metadata.keys():
            protein_alias = annotation_metadata['alias']
            proteins[i].append(protein_alias)
        if 'pdb' in annotation_metadata.keys():
            protein_structure = annotation_metadata['pdb']
            proteins[i].append(protein_structure)
        time.sleep(0.2)
    return(proteins)

def PDB_Extraction_Reactome(proteins):
    for i in range(len(proteins)):
        if len(proteins[i]) == 6:
            continue
        if type(proteins[i][6])== list:
            for j in proteins[i][6]:
                PDB_metadata = get_info(j)
                if 'rcsb_entry_info' in PDB_metadata:
                    if 'resolution_combined' in PDB_metadata['rcsb_entry_info']:
                        if PDB_metadata['rcsb_entry_info']['resolution_combined'][0] > 2.5:
                            proteins[i][6].remove(j)
                elif 'rcsb_accession_info' in PDB_metadata:
                    if int(PDB_metadata['rcsb_accession_info']['initial_release_date'][0:4]) < 2010:
                        proteins[i][6].remove(j)
        if type(proteins[i][6]) == str:
            PDB_metadata = get_info(proteins[i][6])
            if 'rcsb_entry_info' in PDB_metadata:
                if 'resolution_combined' in PDB_metadata['rcsb_entry_info']:
                    if PDB_metadata['rcsb_entry_info']['resolution_combined'][0] > 2.5:
                        proteins[i][6] = np.nan
            elif 'rcsb_accession_info' in PDB_metadata:
                if int(PDB_metadata['rcsb_accession_info']['initial_release_date'][0:4]) < 2010:
                    proteins[i][6] = np.nan
    time.sleep(0.2)
    return (proteins)

def disgenet_gene_disease(Gene_uniprot_ID, Gene_name):
    Gene_ID = Gene_uniprot_ID
    disgenet_Gene_disease_association = 'https://www.disgenet.org/api/gda/gene/uniprot/' + Gene_ID
    GDA_request = requests.get(disgenet_Gene_disease_association)
    disease_list = []
    
    if GDA_request.status_code == 200:
        GDA_json = GDA_request.json()
        
        check = 0
        for disease in GDA_json:
            if disease['score'] >= 0.5:
                try:
                    disease_class_name = disease['disease_class_name'].split(';')
                    for i in range(len(disease_class_name)):
                        disease_class_name[i] = disease_class_name[i].strip()
                except AttributeError:
                    disease_class_name = ['No Disease Class']
                finally:
                    for i in range(len(disease_class_name)):
                        temp_list = []
                        temp_list.append(Gene_name)
                        temp_list.append(disease['disease_name'])
                        temp_list.append(disease_class_name[i])
                        temp_list.append(disease['diseaseid'])
                        temp_list.append(disease['score'])
                        disease_list.append(temp_list)
                    check = check + 1

            elif check == 0:
                disease_dict = [[Gene_name,np.nan,np.nan,np.nan,np.nan]]
                disease_df = pd.DataFrame(disease_dict, columns = ['Gene','Disease Name','Disease Class','Disease ID','Score'])
                #print(disease_df)
                return(disease_df)
            
            elif check != 0:
                disease_df = pd.DataFrame(disease_list, columns = ['Gene','Disease Name','Disease Class','Disease ID','Score'])
                #print(disease_df)
                return(disease_df)
    else:
        disease_dict = [[Gene_name,np.nan,np.nan,np.nan,np.nan]]
        disease_df = pd.DataFrame(disease_dict, columns = ['Gene','Disease Name','Disease Class','Disease ID','Score'] )
        #print(disease_df)
        return(disease_df)


#startTime = time.time()

REACTOME_Pathway = ['R-HSA-73884.2','R-HSA-5696398.2','R-HSA-73893.1','R-HSA-73942.1',
                         'R-HSA-5685942.1','R-HSA-5685938.1','R-HSA-5685939.1','R-HSA-5693571.1',
                         'R-HSA-6783310.2','R-HSA-5358508.1']

for pathway in REACTOME_Pathway:
    
    pathway_name_request = requests.get('https://reactome.org/ContentService/data/query/' + pathway)
    pathway_name_json = pathway_name_request.json()
    pathway_name = pathway_name_json["displayName"]
    
    Gene_list = Pathway_Protein_Extraction(pathway)
    Gene_list = Annotation_Extraction(Gene_list)
    Gene_list = pd.DataFrame(Gene_list, columns = ['UniProt ID','Gene','Ensembl ID',
                                                   'Entrez ID','PharmGKB ID','Alias','PDB'])
    Gene_list.to_csv(pathway_name + '.csv')
    Gene_list['Pathway'] = pathway_name
    Gene_list['Score'] = 1
    Gene_list.drop(labels = ['UniProt ID','Ensembl ID','Entrez ID','PharmGKB ID','Alias','PDB'], axis = 1)
    Gene_list.to_csv(pathway_name + ' Pathway-Gene.csv')
    
    for i in range(len(Gene_list)):
        if i ==0:
            final_disease_df = disgenet_gene_disease(Gene_list['UniProt ID'][0],Gene_list['Gene'][0])
        if i > 0:
            disease_list_final = disgenet_gene_disease(Gene_list['UniProt ID'][i],Gene_list['Gene'][i])
            final_disease_df = pd.concat([final_disease_df,disease_list_final], axis = 0)
    final_disease_df.to_csv(pathway_name + ' GDA.csv')
    final_disease_df.drop(labels = ['Disease Class','Disease ID'], axis = 1)
    final_disease_df.to_csv(pathway_name + ' GDA Network.csv')

#executionTime = (time.time() - startTime)
#print('Execution time in seconds: ' + str(executionTime))


########## NOT IMPORTANT ##########
os.chdir('E:\Documents\VIII Sem (Project)\Gene-Disease Association Network')

for pathway in REACTOME_Pathway:
    pathway_name_request = requests.get('https://reactome.org/ContentService/data/query/' + pathway)
    pathway_name_json = pathway_name_request.json()
    pathway_name = pathway_name_json["displayName"]
    
    temp_df = pd.read_csv(pathway_name + ' GDA Network.csv')
    temp_df = temp_df.drop(labels = ['Disease Class','Disease ID'], axis = 1)
    temp_df.to_csv(pathway_name + ' GDA Network.csv')
    
    temp_df_two = pd.read_csv(pathway_name + ' Pathway-Gene.csv')
    temp_df_two = temp_df_two.drop(labels = ['UniProt ID','Ensembl ID','Entrez ID',
                                     'PharmGKB ID','Alias','PDB'], axis = 1)
    temp_df_two = temp_df_two.reindex(columns = ['Pathway','Gene','Score'])
    temp_df_two.to_csv(pathway_name + ' Pathway-Gene.csv')

for pathway in REACTOME_Pathway:
    pathway_name_request = requests.get('https://reactome.org/ContentService/data/query/' + pathway)
    pathway_name_json = pathway_name_request.json()
    pathway_name = pathway_name_json["displayName"]
    
    temp_df = pd.read_csv(pathway_name + ' GDA.csv')
    temp_df['Pathway'] = pathway_name
    temp_df = temp_df.reindex(columns = ['Pathway','Gene','Disease Name','Disease Class','Disease ID','Score'])
    temp_df.to_csv(pathway_name + ' GDA.csv')
    
for i in range(len(REACTOME_Pathway)):
    pathway_name_request = requests.get('https://reactome.org/ContentService/data/query/' + REACTOME_Pathway[i])
    pathway_name_json = pathway_name_request.json()
    pathway_name = pathway_name_json["displayName"]
    
    if i == 0:
        initial_GDA_df = pd.read_csv(pathway_name + ' GDA.csv')
       
    if i != 0:
        final_GDA_df = pd.read_csv(pathway_name + ' GDA.csv')
        initial_GDA_df = pd.concat([initial_GDA_df,final_GDA_df], axis = 0)
    

pathway_name_list = []
for i in range(len(REACTOME_Pathway)):
    pathway_name_request = requests.get('https://reactome.org/ContentService/data/query/' + REACTOME_Pathway[i])
    pathway_name_json = pathway_name_request.json()
    pathway_name = pathway_name_json["displayName"]
    pathway_name_list.append(pathway_name)

Disease_class = initial_GDA_df['Disease Class'].dropna()
Disease_class = Disease_class.unique()
Disease_class = list(Disease_class)
Pathway_class = ['Pathway']
Disease_column = Pathway_class + Disease_class
Gene_Disease_df = pd.DataFrame(columns = Disease_column, index = [0,1,2,3,4,5,6,7,8,9] )

for i in range(len(Gene_Disease_df)):
    Gene_Disease_df['Pathway'][i] = pathway_name_list[i]
Gene_Disease_df = Gene_Disease_df.fillna(0)
Gene_Disease_df = Gene_Disease_df.set_index ('Pathway')  

for j in range(len(initial_GDA_df)):
    pathway = initial_GDA_df['Pathway'][j]
    try:
        Disease_class_name = initial_GDA_df['Disease Class'][j]
        Gene_Disease_df[Disease_class_name][pathway] = Gene_Disease_df[Disease_class_name][pathway] + 1
    except KeyError:
        continue

os.chdir('E:\Documents\VIII Sem (Project)')
Gene_Disease_df.to_csv('Final GDA Association.csv')







