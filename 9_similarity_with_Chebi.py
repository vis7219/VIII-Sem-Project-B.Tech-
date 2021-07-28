import pandas as pd

from selenium import webdriver
from selenium.webdriver.chrome.options import Options
chrome_options = Options()
chrome_options.add_argument("--headless")

#struct_sim_dict = {'Primary Molecule' : [],
#                      'Similar Molecule' : [],
#                      'CHEBI' : [],
#                      'Tanimoto Score' : [],
#                      'Formula' : [],
#                      'Mass' : [],
#                      'Charge' : []
#                      }


final_list = []

drug_df = pd.read_csv('(5) Gene_drug_list.csv')

driver = webdriver.Chrome()
driver.get("https://www.ebi.ac.uk/chebi/advancedSearchFT.do")

xxxx = 0
for i in range(0, len(drug_df)):
    if type(drug_df['SMILES'][i]) == float:
        continue
    
    
    search_box = driver.find_element_by_name('searchString')
    search_box.send_keys(drug_df['SMILES'][i])
    submit_box = driver.find_element_by_name('submit').click()
    
    
    try:
        similarity_link = driver.find_element_by_xpath('/html/body/div/div[1]/table[2]/tbody/tr/td/table[1]/tbody/tr/td/table/tbody/tr[2]/td/ul/form/li[2]/a')
        driver.execute_script("arguments[0].click();", similarity_link)
    except:
        continue
    
    
    for l in range(2,5000):
        
        try:
            no_of_rows = driver.find_elements_by_xpath('//*[@id="content"]/table/tbody/tr[3]/td/table/tbody/tr[2]/td/table/tbody/tr')
            if len(no_of_rows) == 1:
                no_of_rows.append(str(2))
            for j in range(1, len(no_of_rows)):
                no_of_columns = driver.find_elements_by_xpath('//*[@id="content"]/table/tbody/tr[3]/td/table/tbody/tr[2]/td/table/tbody/tr[' + str(j) + ']/td')
            
                for k in range(1,len(no_of_columns) + 1):
                    temp_list = []
                    temp_list.append(drug_df['Molecule Name'][i])
#                    struct_sim_dict['Primary Molecule'].append(drug_df['Chemical Names'][i])
                    
#                    struct_sim_dict['Similar Molecule'].append(driver.find_element_by_xpath('/html/body/div/div[1]/table/tbody/tr[3]/td/table/tbody/tr[2]/td/table/tbody/tr[' + str(j) + ']/td[' + str(k) + ']/table/tbody/tr[1]/td/a').text)
                    temp_list.append(driver.find_element_by_xpath('/html/body/div/div[1]/table/tbody/tr[3]/td/table/tbody/tr[2]/td/table/tbody/tr[' + str(j) + ']/td[' + str(k) + ']/table/tbody/tr[1]/td/a').text)
                    
                    temp_attributes = driver.find_element_by_xpath('//*[@id="content"]/table/tbody/tr[3]/td/table/tbody/tr[2]/td/table/tbody/tr[' + str(j) + ']/td[' + str(k) + ']/table/tbody/tr[2]/td[2]').text
                    temp_attributes = temp_attributes.split('\n')
                    temp_attributes.remove('Stars:')
                    for attributes in temp_attributes:
                        x,y = attributes.split(':')
                        temp_list.append(y)
#                        struct_sim_dict[x].append(y)
                    final_list.append(temp_list)
            next_page = driver.find_element_by_link_text(str(l))
            next_page.click()
        except:
            break



    
driver.quit()
similarity_df = pd.DataFrame(final_list, columns = ['Primary Molecule',
                                          'Similar Molecule',
                                          'CHEBI',
                                          'Tanimoto Score',
                                          'Formula',
                                          'Mass',
                                          'Charge'])

#similarity_df.to_csv('(6) Similar_Molecules_ChEBI.csv')
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    