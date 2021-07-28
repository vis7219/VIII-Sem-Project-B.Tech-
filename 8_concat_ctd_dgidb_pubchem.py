import pandas as pd
import numpy as np
import os

os.chdir('E:\Documents\VIII Sem (Project)\File Processing Order')

ctd = pd.read_csv('(3) Gene_Drug_List_CTD.csv')
ctd.drop_duplicates(subset = ['SMILES', 'GeneSymbol'],
                    inplace = True)

ctd.drop(columns = ['ChemicalID', 'CasRN', 'PubMedIDs',
                    'No. Of Interactions'],
         inplace = True)

ctd.rename({'ChemicalName' : 'Molecule Name',
            'GeneSymbol' : 'Gene'},
           axis = 1,
           inplace = True)


dgidb = pd.read_csv('(2) Gene_Drug_List_DGiDB.csv')
dgidb.drop_duplicates(subset = ['SMILES', 'Gene'], inplace = True)
dgidb.drop(columns = ['Drug ID', 'Interaction ID', 'Interaction Type',
                      'pmID', 'Score', 'Source'],
           inplace = True)
dgidb.rename({'Chemical Names' : 'Molecule Name'},
             inplace = True,
             axis = 1)


pubchem = pd.read_csv('(4) Gene_Drug_list_PubChem.csv')
pubchem.drop_duplicates(subset = ['SMILES', 'Gene'], inplace = True)
pubchem.drop(columns = ['ID', 'Type of Action', 'cid', 'Gene Uniprot ID',
                        'Database', 'Drug Group', 'Units', 'Affinity',
                        'Molecule Name'],
             inplace = True)
pubchem.rename({'Formula' : 'Molecular Formula',
                'Primary Molecule' : 'Molecule Name'},
               inplace = True,
               axis = 1)

final_df = pd.concat([ctd, dgidb, pubchem])
final_df = final_df.drop_duplicates(subset = ['SMILES', 'Gene'])

final_df = final_df.dropna(subset = ['SMILES'])
#final_df.set_index(['Molecule Name'], inplace = True)
#final_df.to_csv('(5) Gene_drug_list.csv')
