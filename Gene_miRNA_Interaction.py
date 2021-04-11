import pandas as pd
import requests
import numpy as np
import os
from bs4 import BeautifulSoup

import selenium
from selenium.webdriver.common.keys import Keys
from selenium.webdriver.chrome.options import Options
from seleniumwire import webdriver




def miRNA_MicroT_CDS(Gene_Ensembl):
    Gene_ID = str(Gene_Ensembl)
    miRNA_name = []
    miRNA = []
    #There is a threshold value in the url. If needed, variable to be created to manage threshold
    MicroT_CDS_csv_url = 'http://diana.imis.athena-innovation.gr/DianaTools/index.php?r=download/microT_CDS&keywords=' + Gene_ID + '&genes=' + Gene_ID + '%20&mirnas=&descr=&threshold=0.7'
    MicroT_CDS_response = requests.get(MicroT_CDS_csv_url)
    if MicroT_CDS_response.status_code == 200:
        open('temporary data.csv', 'wb').write(MicroT_CDS_response.content)
        MicroT_miRNA = pd.read_csv('temporary data.csv')
        for i in range(len(MicroT_miRNA)):
            temp_miRNA = []
            if MicroT_miRNA['Mirna Name'][i][0:3] == 'hsa':
                miRNA_name.append(MicroT_miRNA['Mirna Name'][i])
                temp_mitg = MicroT_miRNA['miTG score'][i]
                temp_name = MicroT_miRNA['Mirna Name'][i]
                continue
            else:
                temp_miRNA.append(temp_name)
                temp_miRNA.append(temp_mitg)
                temp_miRNA.append(MicroT_miRNA['Transcript Id'][i])
                temp_miRNA.append(MicroT_miRNA['Gene Id(name)'][i])
                miRNA.append(temp_miRNA)
        MicroT_df = pd.DataFrame(miRNA, columns = ['miRNA Name','miTG Score',
                                                   'Transcript Region','Chr Location'])
        return(MicroT_df)
    else:
        miRNA = np.nan
        return(miRNA)

def miRNA_miRcode(Gene_name):
    Gene_name = str(Gene_name)
    miRNA = []
    miRcode_url = 'http://mircode.org/?gene=' + Gene_name + '&mirfam=&class=&cons=&trregion='
    miRcode_response = requests.get(miRcode_url)
    if miRcode_response.status_code == 200:
        miRcode_list = pd.read_html(miRcode_url)
        if len(miRcode_list) == 7:
            miRNA_df = miRcode_list[6]
            for i in range(5,len(miRNA_df)):
                if miRNA_df[0][i].count('/') > 0:
                    miRNA_name = miRNA_df[0][i][4:len(miRNA_df[0][i])].split('/')
                    species_name = miRNA_df[0][i][0:4]
                    for j in range(len(miRNA_name)):
                        temp_list = []
                        miRNA_name[j] = 'hsa-' + species_name + miRNA_name[j]
                        temp_list.append(miRNA_name[j])
                        temp_list.append(miRNA_df[1][i])
                        temp_list.append(miRNA_df[2][i])
                        temp_list.append(miRNA_df[3][i])
                        miRNA.append(temp_list)
                else:
                    temp_list = []
                    miRNA_name = 'hsa-' + miRNA_df[0][i]
                    temp_list.append(miRNA_name)
                    temp_list.append(miRNA_df[1][i])
                    temp_list.append(miRNA_df[2][i])
                    temp_list.append(miRNA_df[3][i])
                    miRNA.append(temp_list)
            miRNA_df = pd.DataFrame(miRNA, columns = ['miRNA Name','Seed Position','Seed Type','Transcript Region'])
            return(miRNA_df)
        else:
            miRNA = np.nan
            return(miRNA)
    else:
        miRNA = np.nan
        return(miRNA)

def miRNA_Comirnet(Gene_name):
    try:
        comirnet_url = 'http://comirnet.di.uniba.it:8080/interactionsRetrieves?gene=' + Gene_name + '&mirna=&min=&max=&numTuplesForPage=&showInputAndWeight=on'
        miRNA_df = pd.read_html(comirnet_url)
        miRNA_df = miRNA_df[0]
        return(miRNA_df)
    except ValueError:
        miRNA = np.nan
        return(miRNA)

def miRNA_miRDB(Gene_name):
    chrome_options = Options()
    chrome_options.add_argument("--headless")
    driver = webdriver.Chrome(options = chrome_options)
    miRDB_url = 'http://mirdb.org/index.html'
    driver.get(miRDB_url)
    gene_search_box = driver.find_element_by_xpath('//*[@id="table2"]/tbody/tr[4]/td/form/p/input[1]')
    gene_search_box.send_keys(Gene_name)
    search_button = driver.find_element_by_xpath('//*[@id="table2"]/tbody/tr[4]/td/form/p/input[2]')
    search_button.click()
    page_source = driver.page_source
    soup = BeautifulSoup(page_source,'lxml')
    table = soup.find_all('table')[1]
    rows = table.find_all('tr')
    temp_dict = {'miRNA Name':[],'Target Score':[],'Target Rank':[]}
    #temp_dict = {'Target Rank': [],'Target Score': [],'miRNA Name': []}
    k = 0
    for i in rows:
        data = i.find_all('td')
        for j in data:
            data_temp = j.text
            if k == 1:
                temp_dict['Target Rank'].append(data_temp)
            if k == 2:
                temp_dict['Target Score'].append(data_temp)
            if k == 3:
                temp_dict['miRNA Name'].append(data_temp)
            k = k + 1
            if k ==6:
                k = 0
    try:
        temp_dict['Target Rank'].remove('Target Rank')
        temp_dict['Target Score'].remove('Target Score')
        temp_dict['miRNA Name'].remove('miRNA Name')
        miRNA_df = pd.DataFrame(temp_dict)
        return(miRNA_df)
    except ValueError:
        miRNA = np.nan
        return(miRNA)
    finally:
        driver.quit()

def miRNA_miRabel(Gene_name):
    driver = webdriver.Chrome()
    miRabel_url = 'http://bioinfo.univ-rouen.fr/mirabel/index.php?page=gene'
    driver.get(miRabel_url)
    driver.implicitly_wait(1)
    gene_box = driver.find_element_by_xpath('//*[@id="center"]/div/form/div/p/input')
    gene_box.send_keys(Gene_name)
    driver.implicitly_wait(1)
    search_button = driver.find_element_by_xpath('//*[@id="center"]/div/form/p[2]/input[3]')
    search_button.click()
    driver.implicitly_wait(2)
    try:
        search_limit = driver.find_element_by_xpath('//*[@id="' + Gene_name + '_length"' + ']/label/select/option[6]')
        search_limit.click()
        driver.implicitly_wait(5)
        network_connections_list = []
        for request in driver.requests:
            temp_list = []
            if request.response:
                temp_list.append(request.url)
                temp_list.append(request.response.status_code)
                network_connections_list.append(temp_list)
        for network_connections in network_connections_list:
            if 'http://bioinfo.univ-rouen.fr/mirabel/model/gene.ajax.php?ncbiGeneId=' in network_connections[0]:
                table_request_url = network_connections[0]
        table_request = requests.get(table_request_url)
        table_json = table_request.json()
        table = table_json['aaData']
        df = pd.DataFrame(table, columns = ['miR','miRabel Score','PITA',',MiRanda',
                                          'SVMircO','TargetScan','ExpVal',"5'UTR",
                                          'CDS'] )
        return(df)
    except:
        x = 'No such gene Exists in the database'
        return(x)








