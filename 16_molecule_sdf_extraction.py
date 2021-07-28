import pandas as pd
import numpy as np
from selenium import webdriver
import time

df = pd.read_csv('ABL1_swissADME.csv')
df['CHEBI ID'] = np.nan
driver = webdriver.Chrome()

for i in range(len(df)):
    smiles = df['SMILES'][i]
    driver.get("https://www.ebi.ac.uk/chebi/advancedSearchFT.do")
    
    try:
        search_box = driver.find_element_by_name('searchString')
        search_box.send_keys(smiles)
        submit_box = driver.find_element_by_name('submit').click()
    
        df.loc[i , 'CHEBI ID'] = driver.find_element_by_xpath('//*[@id="content"]/table[2]/tbody/tr/td/table[1]/tbody/tr[1]/td/table/tbody/tr[1]/td[2]/table/tbody/tr[3]/td[2]/b').text
        time.sleep(7)
    
    except:
        pass

df1 = df.set_index('CHEBI ID') 
df1.to_csv('ABL1_swissADME.csv')
    