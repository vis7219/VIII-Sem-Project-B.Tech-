import os
import pandas as pd

# This script extracts information related to miRNA-Disease Interaction
# The databases from which information is extracted -
# 1. Tarbase v8
# 2. HMDD v3.2
# 3. dbDEMC
# 4. miRcancer


os.chdir('E:\Documents\VIII Sem (Project)\TarBase_v8_download')
chunk = pd.read_csv('TarBase_v8_download.txt', delimiter = '\t',chunksize = 100000)
df = pd.concat(chunk)
def tarbase(Gene_id):
    x = 0