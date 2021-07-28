import os
import pandas as pd

os.chdir('E:\Documents\VIII Sem (Project)\Genes\Extra files\Pubchem_Predecessors')

gene_list = pd.read_csv('Gene List.csv')
x = 0
final_df = pd.DataFrame()
for gene in gene_list['Gene']:
    
    gene_df = pd.DataFrame(columns = ['Gene',
                                      'Primary Molecule',
                                      'ID',
                                      'Type of Action',
                                      'cid',
                                      'Molecule Name',
                                      'Gene Uniprot ID',
                                      'Database',
                                      'Drug Group',
                                      'Units',
                                      'Affinity'
                                      ]
                           )
    
    try:
        chembl = pd.concat([pd.read_csv('Gene_' + gene + '_chembldrugtargets.csv'),
                            pd.read_csv('Gene_' + gene + '_chembldrugtargets (1).csv')])
        chembl.drop(['mecid', 'moa', 'targetchemblid', 'targetname',
                              'geneids', 'pmids', 'dois'], axis = 1 ,inplace = True)
        chembl.rename(columns = {'drugname' : 'Primary Molecule',
                                 'action' : 'Type of Action',
                                 'protacxns' : 'Gene Uniprot ID',
                                 'cmpdname' : 'Molecule Name',
                                 'chemblid' : 'ID'},
                      inplace = True)
        chembl['Gene'] = gene
        columns = ['Gene', 'Primary Molecule', 'ID', 'Type of Action', 'cid', 'Molecule Name', 'Gene Uniprot ID']
        chembl = chembl[columns]
        chembl['Database'] = 'ChEMBL'
    except:
        pass
    
    try:
        ctd = pd.read_csv('Gene_' + gene + '_chembldrugtargets.csv')
        ctd.drop(['mecid', 'moa', 'targetchemblid', 'targetname',
                              'geneids', 'pmids', 'dois'], axis = 1 ,inplace = True)
        ctd.rename(columns = {'drugname' : 'Primary Molecule',
                                 'action' : 'Type of Action',
                                 'protacxns' : 'Gene Uniprot ID',
                                 'cmpdname' : 'Molecule Name',
                                 'chemblid' : 'ID'},
                      inplace = True)
        ctd['Gene'] = gene
        columns = ['Gene', 'Primary Molecule', 'ID', 'Type of Action', 'cid', 'Molecule Name', 'Gene Uniprot ID']
        ctd = ctd[columns]
        ctd['Database'] = 'CTD'
    except:
        pass
    
    try:
        dgidb = pd.read_csv('Gene_' + gene + '_dgidb.csv')
        dgidb.drop(['geneid', 'geneclaimname','sid', 'pmids', 'dois',
                    'drugclaimname', 'drugclaimprimaryname','interactionclaimsource'], axis = 1, inplace = True)
        dgidb.rename(columns = {'genename' : 'Gene',
                                'interactiontypes' : 'Type of Action',
                                'drugname' : 'Primary Molecule',
                                'drugchemblid' : 'ID',
                                'cmpdname' : 'Molecule Name'},
                     inplace = True)
        columns = ['Gene', 'Primary Molecule', 'ID', 'Type of Action', 'cid', 'Molecule Name']
        dgidb = dgidb[columns]
        dgidb['Database'] = 'DGIdb'
    except:
        pass
    
    try:
        drugbank = pd.read_csv('Gene_'+ gene + '_drugbank.csv')
        drugbank.drop(['geneid', 'drugtype', 'drugdetail', 'targettype', 'targetid', 'targetname', 'targetcomponent',
                       'targetcomponentname', 'generalfunc', 'specificfunc', 'pmids', 'dois'], axis = 1, inplace = True)
        drugbank.rename(columns = {'protacxn' : 'Gene Uniprot ID',
                                   'genesymbol' : 'Gene',
                                   'drugname' : 'Primary Molecule',
                                   'druggroup' : 'Drug Group',
                                   'drugaction' : 'Type of Action',
                                   'cmpdname' : 'Molecule Name'},
                        inplace= True)
        columns = ['Gene', 'Primary Molecule', 'Type of Action', 'cid', 'Molecule Name', 'Drug Group', 'Gene Uniprot ID']
        drugbank = drugbank[columns]
        drugbank['Database'] = 'DrugBank'
    except:
        pass
    
    try:
        gtp = pd.read_csv('Gene_' + gene + '_gtopdb.csv')
        gtp.drop(['geneid', 'ligandid', 'primarytarget', 'type', 'pmids', 'targetid',
                  'targetname', 'targetspecies', 'dois'], axis = 1, inplace = True)
        gtp.rename(columns = {'genesymbol' : 'Gene',
                              'protacxn' : 'Gene Uniprot ID',
                              'ligand' : 'Primary Molecule',
                              'action' : 'Type of Action',
                              'cmpdname' : 'Molecule Name',
                              'units' : 'Units',
                              'affinity' : 'Affinity'},
                   inplace = True)
        columns = ['Gene', 'Primary Molecule', 'Type of Action', 'cid', 'Molecule Name', 'Gene Uniprot ID',
                   'Units', 'Affinity']
        gtp = gtp[columns]
        gtp['Database'] = 'Guide To Pharmacology'
    except:
        pass
    
    for database in [chembl,drugbank,dgidb,ctd,gtp]:
        try:
            gene_df = pd.concat([gene_df, database])
        except:
            continue
    gene_df.convert_dtypes()        
    gene_df = gene_df.set_index('Gene')
    final_df = pd.concat([final_df , gene_df])

#final_df.to_csv('Gene_Drug_list_PubChem.csv')

    