import os
import pandas as pd
import numpy as np
from selenium import webdriver

os.chdir('E:\Documents\VIII Sem (Project)')

df = pd.read_csv('Cytarabine_list.csv')
#df = df[df['Mass'] <= 250]
#df.reset_index(inplace = True)
#df.drop(columns = ['index', 'Unnamed: 0'], inplace = True)

driver = webdriver.Chrome()

l = 161
for i in range(160 , len(df), 10):
    k = 9 + i
    A_df = list(df.loc[i : k , 'SMILES'])

    driver.get("http://www.swissadme.ch/")

    search_box = driver.find_element_by_xpath('//*[@id="smiles"]')

    for j in range(len(A_df)):
        search_box.send_keys(A_df[j])
        search_box.send_keys(' ')
        search_box.send_keys('Molecule' + str(l))
        search_box.send_keys('\n')
        l = l + 1
    
    run_button = driver.find_element_by_xpath('//*[@id="submitButton"]').click()
    csv_button = driver.find_element_by_xpath('//*[@id="sib_body"]/div[7]/a[1]').click()
    

#df.to_csv('(10.1) Complete_Molecule_list_start_removed .csv')


    


