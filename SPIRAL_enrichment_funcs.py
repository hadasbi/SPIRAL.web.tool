#!/home/yaellab/SPIRAL/bin/python3.8

import os.path
import pandas as pd
import numpy as np
from scipy.stats import kendalltau, spearmanr, pearsonr, norm
from time import time, sleep
from datetime import datetime
import pickle
import mechanicalsoup
import gseapy
from gseapy.parser import Biomart
from SPIRAL_basic_funcs import *
import urllib.request
import urllib.error
import urllib.parse


####################################################################################################################
def add_GO_terms(sigfile_GO, sigfile_GO_temp, genetable_file, data_path, species,
                 sigtable=None, sigtable_file=None,
                 impute_method='agg_wald', start_from_scratch=False,
                 pvals=[0.000001, 0.0000001], real_samp_name='sample'
                 # log10_pavg_of_genes_thr=-10
                 ):
    print('add_GO_terms!')

    gorilla_url = 'http://cbl-gorilla.cs.technion.ac.il/'

    if species is None:
        species = open(os.path.join(data_path, 'species.txt'), "r").read()

    if start_from_scratch or (not os.path.exists(sigfile_GO_temp)):
        if sigtable is None:
            if sigtable_file is None:
                print('No sigtable path and no sigtable argument.')
            else:
                sigtable = pd.read_excel(sigtable_file, index_col=0, engine='openpyxl')
        sigtable['Gorilla_access_time'] = 99999
    else:
        sigtable = load_excel_with_openpyxl_and_convert_to_pd_DataFrame(sigfile_GO_temp)

    if impute_method == 'no_imputation':
        samp_name = real_samp_name
    else:
        samp_name = 'repcell'

    if len(sigtable) > 0:
        num_samples = int(sigtable.iloc[0, :]['num_' + samp_name + 's'])
        num_genes = int(sigtable.iloc[0, :]['num_genes'])
        # print('num_samples:', num_samples, ' num_genes:', num_genes)

        # Load genes_table
        with open(genetable_file, 'rb') as fp:
            genes_table = pickle.load(fp)
        background_list = str(list(genes_table.index)).replace("'", "").replace('[', '').replace(']', '').replace(', ',
                                                                                                                  '\n')
    pvals.sort(reverse=True)

    structs_lst = sigtable.index[sigtable['Gorilla_access_time'] == 99999]

    for i in structs_lst:
        # for i in list(sigtable.index[sigtable['Gorilla_access_time']==99999]):

        if len(sigtable.loc[i, 'genes']) == 32767:  # The genelist was too long to be written in one excel cell
            structs_file = structs_filename(data_path=data_path, impute_method=impute_method,
                                            num_stds_thresh=sigtable.loc[i, 'num_stds_thresh'],
                                            mu=sigtable.loc[i, 'mu'],
                                            path_len=sigtable.loc[i, 'path_len'],
                                            num_iters=sigtable.loc[i, 'num_iters'])
            structs = read_struct_list_from_file(structs_file)
            geneslist = list(genes_table[structs[sigtable.loc[i, 'old_struct_num']][0]].index)
        else:  # read the gene list from the excel file
            geneslist = read_gene_list(sigtable.loc[i, 'genes'])

        target_list = str(geneslist).replace("'", "").replace('[', '').replace(']', '').replace(', ', '\n')
        # print('target_list:', target_list)

        # text_file = open("target_list.txt", "w")
        # text_file.write(target_list)
        # text_file.close()

        # text_file = open("background_list.txt", "w")
        # text_file.write(background_list)
        # text_file.close()

        browser = mechanicalsoup.StatefulBrowser()
        browser.open(gorilla_url)

        # print(browser.url)
        # print(browser.page)

        browser.select_form()

        browser["species"] = species

        # <option value="ARABIDOPSIS_THALIANA">Arabidopsis thaliana</option>
        # <option value="SACCHAROMYCES_CEREVISIAE">Saccharomyces cerevisiae</option>
        # <option value="CAENORHABDITIS_ELEGANS">Caenorhabditis elegans</option>
        # <option value="DROSOPHILA_MELANOGASTER">Drosophila melanogaster</option>
        # <option value="DANIO_RERIO">Danio rerio (Zebrafish)</option>
        # <option value="HOMO_SAPIENS" selected="selected">Homo sapiens</option>
        # <option value="MUS_MUSCULUS">Mus musculus</option>
        # <option value="RATTUS_NORVEGICUS">Rattus norvegicus</option>

        browser["run_mode"] = "hg"
        # "hg" for Two unranked lists of genes (target and background lists) 
        # "mhg" for Single ranked list of genes

        browser["target_set"] = target_list
        browser["background_set"] = background_list

        browser["db"] = "all"
        # "all" OR "proc" (process) OR "func" (function) OR "comp" (component)

        browser["pvalue_thresh"] = "0.001"
        browser["analysis_name"] = ""
        browser["user_email"] = ""

        browser["output_excel"] = True

        # Uncomment to launch a real web browser on the current page.
        # browser.launch_browser()

        # Uncomment to display a summary of the filled-in form
        # browser.form.print_summary()

        response = browser.submit_selected()
        # time.sleep(30)
        # browser.get_current_page()
        # print(response.text)
        new_link = browser.get_url()
        # print('new_link:', new_link)

        sleep(10)
        browser = mechanicalsoup.StatefulBrowser()
        browser.open(new_link)
        new_link2 = browser.get_url()

        c = 0
        while ("GOResults" not in new_link2) and (c < 10):
            sleep(3)
            browser = mechanicalsoup.StatefulBrowser()
            browser.open(new_link)
            new_link2 = browser.get_url()
            c += 1
        if ("GOResults" not in new_link2):
            print('PROBLEM: "GOResults" not in new_link2')

        # print('new_link2:', new_link2)

        ind = [m.start() for m in re.finditer('/', new_link2)][-1]
        for ontology_init, ontology_html in zip(['proc', 'func', 'comp'],
                                                ['GOResultsPROCESS.html', 'GOResultsFUNCTION.html',
                                                 'GOResultsCOMPONENT.html']):
            new_link3 = new_link2[:ind + 1] + ontology_html
            # print('new_link3', ontology_init, ':', new_link3)
            """
            # To save the links locally, please uncomment the following section. All links will be saved in the 'links' 
            # directory within the same data path. You'll then need to manually move those files to the 
            # 'templates/go/{datapath}/' directory."
            if not os.path.exists(data_path + "/links/"):
                os.makedirs(data_path + "/links/")
            if not os.path.exists(data_path + "/links/" + str(i)):
                os.makedirs(data_path + "/links/" + str(i))
            save_to_file(url=new_link3, file_path=os.path.join(data_path, "links", str(i)), file_name=ontology_html, 
            browser=browser)
            new_link3 = 'https://spiral.technion.ac.il/results/{}/{}/{}'.format(data_n_to_url(data_n=data_path[-1]), 
                                                                                str(i), 
                                                                                ontology_html.replace('.html', ''))
            """
            hyperlink_txt = '=hyperlink("' + new_link3 + '","link")'
            # print(hyperlink_txt)
            sigtable.loc[i, ontology_init + '_link'] = hyperlink_txt

            try:
                dfs = pd.read_html(new_link3)
                df = dfs[1]
                df.columns = df.iloc[0, :]
                df = df.drop([0])
                # print(df)
                for pval in pvals:
                    df = df[df['P-value'].astype(float) <= pval]
                    if len(df) > 0:
                        sigtable.loc[i, ontology_init + '_GOterms_below_' + str(pval)] = str(
                            [(a + ':' + b + ' (qval' + str(c) + ')') for a, b, c in
                             zip(df['GO term'], df['Description'], df['FDR q-value'])])
                    else:
                        sigtable.loc[i, ontology_init + '_GOterms_below_' + str(pval)] = 'NO TERMS'
            except ValueError as error:
                if str(error) == 'No tables found':
                    for pval in pvals:
                        sigtable.loc[i, ontology_init + '_GOterms_below_' + str(pval)] = 'NO TERMS'
                else:
                    for pval in pvals:
                        sigtable.loc[i, ontology_init + '_GOterms_below_' + str(pval)] = str(error)
        sigtable.loc[i, 'Gorilla_access_time'] = datetime.now().strftime("%d/%m/%Y, %H:%M:%S")

        sigtable.to_excel(sigfile_GO_temp)

    sigtable.to_excel(sigfile_GO)
    try:
        os.remove(sigfile_GO_temp)
    except:
        pass

    return sigtable

####################################################################################################################


def save_to_file(url, file_path, file_name, browser):
    """
    Saves a webpage and its images to a file on the local file system.

    Parameters:
    - url (str): URL of the webpage to be saved.
    - file_path (str): directory path where the file will be saved.
    - file_name (str): name of the file to be saved.
    - browser (WebBrowser object): web browser object to fetch the webpage and its images.

    """
    browser.open(url)
    for img in browser.page.find_all('img'):
        imgSrc = img['src']
        imgUrl = re.sub(r'[a-zA-Z]*.html$', imgSrc, url)
        response = urllib.request.urlopen(imgUrl)
        content = response.read()
        with open(os.path.join(file_path, imgSrc), "wb") as f:
            f.write(content)
            f.close()
    response = urllib.request.urlopen(url)
    web_content = response.read().decode('UTF-8')
    with open(os.path.join(file_path, file_name), "w") as f:
        f.write(web_content)
        f.close()
