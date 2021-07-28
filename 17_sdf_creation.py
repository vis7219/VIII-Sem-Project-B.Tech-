import pubchempy as pcp
import pandas as pd
import os

x = 'XRCC3'
os.chdir('E:\Documents\VIII Sem (Project)\File Processing Order\(12.1) Molecule_with_swissADME_files\{}'.format(x))

df = pd.read_csv('{}_swissADME.csv'.format(x))
count = 0
for i in range(len(df)):
    smiles = df['SMILES'][i]
    try:
        cid = pcp.get_compounds(identifier = smiles,
                                namespace = 'smiles')[0].cid
        if type(df['Similar Molecule'][i]) == float:
            name = str(smiles) + '.sdf'
        else:
            name = str(df['Similar Molecule'][i]) + '.sdf'
        print(i , name)
    except:
        print('error: {}'.format(str(df['Similar Molecule'][i])))
        continue
    
    try:
        pcp.download('SDF' , name , cid , 'cid')
    except:
        print('error - Duplicate: {}'.format(df['Similar Molecule'][i]))
        continue

#x = pcp.get_compounds(identifier = 'Nc1nc2n(COCCO)cnc2c(=O)[nH]1',
#                      namespace = 'smiles')[0]
#cid = x.cid

#y = pcp.get_sdf(identifier = str(cid))

#pcp.download('SDF' , 'Acyclovir.sdf' , cid , 'cid')
