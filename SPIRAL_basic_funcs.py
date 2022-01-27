import pandas as pd
import numpy as np
from sklearn.cluster import AgglomerativeClustering
from scipy.cluster.hierarchy import dendrogram
from sklearn.decomposition import PCA
from sklearn.preprocessing import normalize
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import scipy.io
from scipy import special
import csv
import scipy.sparse as sp_sparse
import tables as tb
import re
import os
import openpyxl
from gseapy.parser import Biomart
import glob
from time import sleep, time


#### Data loading and pre-processing #################################################################################

def get_MTgenes(genelist, species):
    # returns a sublist of mitochondrial genes
    print('get_MTgenes!')

    error = None

    if species in ['synthetic', 'other']:
        return [], error

    # gene symbols
    if 'ACTB' in genelist or 'GAPDH' in genelist or 'PGK1' in genelist:
        print('Found gene symbols')
        return [gene for gene in genelist if gene[:3] == 'MT-'], error
    else:
        if species == "ARABIDOPSIS_THALIANA":
            bm = Biomart(host="plants.ensembl.org")
            dataset = "athaliana_eg_gene"
        elif species == "SACCHAROMYCES_CEREVISIAE":
            bm = Biomart()
            dataset = 'scerevisiae_gene_ensembl'
        elif species == "CAENORHABDITIS_ELEGANS":
            bm = Biomart()
            dataset = 'celegans_gene_ensembl'
        elif species == "DROSOPHILA_MELANOGASTER":
            bm = Biomart()
            dataset = 'dmelanogaster_gene_ensembl'
        elif species == "DANIO_RERIO":
            bm = Biomart()
            dataset = 'drerio_gene_ensembl'
        elif species == "HOMO_SAPIENS":
            bm = Biomart()
            dataset = 'hsapiens_gene_ensembl'
        elif species == "MUS_MUSCULUS":
            bm = Biomart()
            dataset = 'mmusculus_gene_ensembl'
        elif species == "RATTUS_NORVEGICUS":
            bm = Biomart()
            dataset = 'rnorvegicus_gene_ensembl'

        # ENSEMBL IDs
        if genelist[0][:3] == 'ENS' and genelist[1][:3] == 'ENS':
            print('Found ENSEMBL IDs')
            ENSG_symbols_table = None
            count = 0
            while ENSG_symbols_table is None and count < 3:
                ENSG_symbols_table = bm.query(dataset=dataset,
                                              attributes=['ensembl_gene_id', 'external_gene_name'],
                                              filters={'ensembl_gene_id': genelist})
                sleep(10)
                count += 1
            if ENSG_symbols_table is None:
                error = "Due to a connection error with Biomart, SPIRAL could not get a list of mitochondrial genes."
                print(error)
                return [], error
            ENSG_symbols_table = ENSG_symbols_table.dropna()
            ENSG_symbols_table = ENSG_symbols_table.loc[
                [i for i in ENSG_symbols_table.index if ENSG_symbols_table.loc[i, 'ensembl_gene_id'] in genelist]
            ]
            ENSG_symbols_table_MT = ENSG_symbols_table.loc[
                [i for i in ENSG_symbols_table.index if 'MT-' in ENSG_symbols_table.loc[i, 'external_gene_name']]]
            return list(set(ENSG_symbols_table_MT.loc[:, 'ensembl_gene_id'])), error

        # unknownk gene ids
        else:
            ENSG_symbols_table = None
            count = 0
            while ENSG_symbols_table is None and count < 3:
                ENSG_symbols_table = bm.query(dataset=dataset,
                                              attributes=['ensembl_gene_id', 'hgnc_id', 'uniprot_gn_id',
                                                          'entrezgene_id', 'external_gene_name'])
                sleep(10)
                count += 1
            if ENSG_symbols_table is None:
                error = "Due to a connection error with Biomart, SPIRAL could not get a list of mitochondrial genes."
                print(error)
                return [], error
            ENSG_symbols_table = ENSG_symbols_table.dropna()
            s = int(np.ceil(len(genelist) / 2))
            correct_col = None
            for col in ['ensembl_gene_id', 'hgnc_id', 'uniprot_gn_id', 'entrezgene_id']:
                new_s = len(set(ENSG_symbols_table.loc[:, col]).intersection(set(genelist)))
                if new_s > s:
                    s = new_s
                    correct_col = col
            if correct_col is not None:
                print('Found ', correct_col)
                ENSG_symbols_table = ENSG_symbols_table.loc[
                    [i for i in ENSG_symbols_table.index if ENSG_symbols_table.loc[i, correct_col] in genelist]
                ]
                ENSG_symbols_table_MT = ENSG_symbols_table.loc[
                    [i for i in ENSG_symbols_table.index if 'MT-' in ENSG_symbols_table.loc[i, 'external_gene_name']]]
                return list(set(ENSG_symbols_table_MT.loc[:, correct_col])), error

    error = 'Gene type was not identified.'
    print(error)
    return [], error


####################################################################################################################
def data_path_name(analysis_folder, data_n):
    data_path = os.path.join(analysis_folder, 'data' + str(data_n))
    return data_path


####################################################################################################################
def data_norm_loc_name(data_path):
    data_loc = os.path.join(data_path, 'counts_norm.csv')
    return data_loc


####################################################################################################################
def data_norm_filt_loc_name(data_path):
    return os.path.join(data_path, 'counts_norm_filt.csv')


####################################################################################################################
def spatial_norm_filt_loc_name(data_path):
    return os.path.join(data_path, 'spatial_coors_norm_filt.csv')


####################################################################################################################
def vln_plot_filename_(static_path, data_n):
    return os.path.join(static_path, 'vln_data' + str(data_n) + '.jpg')


####################################################################################################################
def vln_plot_mt_filename_(static_path, data_n):
    return os.path.join(static_path, 'vln_mt_data' + str(data_n) + '.jpg')


####################################################################################################################
def structs_filename(data_path, impute_method, num_stds_thresh, mu, path_len, num_iters):
    return os.path.join(data_path, impute_method + '_std_' + str(num_stds_thresh) + '_mu_' + str(
        mu) + '_path_len_' + str(path_len) + '_num_iters_' + str(num_iters) + '_structs.txt')


####################################################################################################################
def sig_filename(data_path, impute_method, num_stds_thresh, mu, path_len, num_iters):
    return os.path.join(data_path, impute_method + '_sigtable_std_' + str(num_stds_thresh) + '_mu_' + str(
        mu) + '_path_len_' + str(path_len) + '_num_iters_' + str(num_iters) + '.xlsx')


####################################################################################################################
def final_sig_filename(data_path, impute_method):
    return os.path.join(data_path, impute_method + '_sigtable_filt_GO_vis.xlsx')


####################################################################################################################
def pic_folder(data_path):
    return os.path.join(data_path, 'structure_layouts')


####################################################################################################################
def spatial_norm_loc_name(data_path):
    data_loc = os.path.join(data_path, 'spatial_coors_norm.csv')
    return data_loc


####################################################################################################################
def delete_genes_and_cells_with_zero_reads(table):
    print(table.shape)
    table = table.loc[(table.sum(axis=1) != 0), (table.sum(axis=0) != 0)]
    print(table.shape)
    return table


def delete_genes_with_zero_reads(table):
    print(table.shape)
    table = table.loc[(table.sum(axis=1) != 0), :]
    print(table.shape)
    return table


def median_count_normalization(table, with_log=False, mode='median'):
    # Normalization to library size factor
    cells_num_of_reads = table.sum().astype(int)
    if mode == 'median':
        factor = cells_num_of_reads.median()
        print('median_total_counts_across_cells =', factor)
    elif mode == 'mean':
        factor = cells_num_of_reads.mean()
        print('mean_total_counts_across_cells =', factor)
    else:
        print('ERROR: unknown mode')

    table = table.astype('float64').divide(cells_num_of_reads.values, axis=1).multiply(factor)

    # log(x+1)
    if with_log:
        table = np.log10(table + 1)

    return table


def norm_0_max(data, maxx):
    # The function normalizes the range of each row of the data to 0-maxx
    # adjust the mins of all genes to zero
    data_norm = data.subtract(data.min(axis=1), axis='index')
    # adjust the range to 0-1
    data_norm = data_norm.divide(data_norm.max(axis=1) - data_norm.min(axis=1), axis='index')
    data_norm = maxx * data_norm
    return data_norm


def switch_10x_to_txt(matrix_mtx_file, features_tsv_file, barcodes_tsv_file, new_txt_file):
    the_matrix = scipy.io.mmread(matrix_mtx_file).todense()

    with open(features_tsv_file) as fd:
        rd = csv.reader(fd, delimiter="\t", quotechar='"')
        genes = []
        for row in rd:
            # print(row)
            genes.append(row[1])
        # print(genes)

    with open(barcodes_tsv_file) as fd:
        rd = csv.reader(fd, delimiter="\t", quotechar='"')
        cells = []
        for row in rd:
            # print(row)
            cells.append(row[0])
        # print(cells)

    orig_scRNA_data = pd.DataFrame(the_matrix, columns=cells, index=genes)
    validate_folder(new_txt_file)
    orig_scRNA_data.to_csv(new_txt_file, sep='\t')


def switch_10x_to_txt_in_sparse(matrix_mtx_file, features_tsv_file, barcodes_tsv_file, new_txt_file):
    print('reading matrix')
    the_matrix = scipy.io.mmread(matrix_mtx_file)

    print('reading features/genes')
    with open(features_tsv_file) as fd:
        rd = csv.reader(fd, delimiter="\t", quotechar='"')
        genes = []
        for row in rd:
            # print(row)
            genes.append(row[1])
        # print(genes)

    print('reading barcodes of cells')
    with open(barcodes_tsv_file) as fd:
        rd = csv.reader(fd, delimiter="\t", quotechar='"')
        cells = []
        for row in rd:
            # print(row)
            cells.append(row[0])
        # print(cells)

    time = time()
    print('constructing sparse dataframe')
    orig_scRNA_data = pd.DataFrame.sparse.from_spmatrix(the_matrix, columns=cells, index=genes)
    print('finished in {0} seconds'.format(time() - time))
    time = time()
    print('writing sparse dataframe to file')
    validate_folder(new_txt_file)
    orig_scRNA_data.to_csv(new_txt_file, sep='\t')
    print('finished in {0} seconds'.format(time() - time))


def switch_h5_to_txt(filename, new_file):
    with tb.open_file(filename, 'r') as f:
        mat_group = f.get_node(f.root, 'matrix')
        barcodes = f.get_node(mat_group, 'barcodes').read()
        data = getattr(mat_group, 'data').read()
        indices = getattr(mat_group, 'indices').read()
        indptr = getattr(mat_group, 'indptr').read()
        shape = getattr(mat_group, 'shape').read()
        matrix = sp_sparse.csc_matrix((data, indices, indptr), shape=shape)

        feature_ref = {}
        feature_group = f.get_node(mat_group, 'features')
        feature_ids = getattr(feature_group, 'id').read()
        feature_names = getattr(feature_group, 'name').read()
        feature_types = getattr(feature_group, 'feature_type').read()
        feature_ref['id'] = feature_ids
        feature_ref['name'] = feature_names
        feature_ref['feature_type'] = feature_types
        tag_keys = getattr(feature_group, '_all_tag_keys').read()
        for key in tag_keys:
            feature_ref[key] = getattr(feature_group, key.decode("utf-8")).read()
        orig_scRNA_data = pd.DataFrame(matrix.todense(), columns=[x.decode("utf-8") for x in barcodes],
                                       index=[x.decode("utf-8") for x in feature_ref['name']])
        validate_folder(new_file)
        orig_scRNA_data.to_csv(new_file)


def write_struct_list_to_file(filename, structlist):
    with open(filename, 'w') as f:
        for item in structlist:
            # f.write("%s\n" % item)
            f.write('{0}\n'.format(item))


def read_struct_list_from_file(filename):
    structlist = []
    with open(filename, 'r') as f:
        for line in f:
            # remove linebreak which is the last character of the string
            currentPlace = line[:-1]
            # add item to the list
            structlist.append(currentPlace)
    for i in range(len(structlist)):
        a = [s.strip('(').strip(')').strip('[').strip(']') for s in structlist[i].split('], [')]
        a = [s.split(', ') for s in a]
        a[0] = [int(s) for s in a[0]]
        a[1] = [int(s) for s in a[1]]
        structlist[i] = a
    return structlist


def write_gene_list(genelist):
    return str(genelist).replace("['", "").replace("']", "").replace("', '", ",")


def read_gene_list(genelist):
    return genelist.split(",")


def validate_folder(file):
    # The function checks if the folder of the input file exist, it it doesn't- it creates it
    ind = [m.start() for m in re.finditer('/', file)][-1]
    if not os.path.exists(file[:ind]):
        os.makedirs(file[:ind])


def load_excel_with_openpyxl_and_convert_to_pd_DataFrame(xlsx_file):
    # this method preserves the hyperlinks in the excel file (unlike pd.read_excel)
    book = openpyxl.load_workbook(xlsx_file)
    sheet = book.active
    data = sheet.values
    columns = next(data)[0:]
    # Create a DataFrame based on the second and subsequent lines of data
    df = pd.DataFrame(data, columns=columns)
    df.index = df.iloc[:, 0]
    df = df.drop([None], axis=1)
    return df


def chooseln(N, k):
    return (special.gammaln(N + 1) - special.gammaln(N - k + 1) - special.gammaln(k + 1))


def correct_for_multiple_testing(num_genes, num_rows, num_sets, num_cols):
    # return np.log10(comb(num_genes,num_rows))+np.log10(comb(num_sets,num_cols))
    return ((chooseln(num_genes, num_rows) + chooseln(num_sets, num_cols)) / np.log(10))


#### data visualization functions ##########################################################################################

def plot_pca(datamatrix, n_components=2, annotate=False):
    pca = PCA(n_components=n_components)
    datamatrix = normalize(datamatrix)
    principalComponents = pca.fit_transform(datamatrix)
    print(principalComponents.shape)

    fig = plt.figure(figsize=(8, 8))
    ax = fig.add_subplot(1, 1, 1)
    ax.set_xlabel('Principal Component 1', fontsize=15)
    ax.set_ylabel('Principal Component 2', fontsize=15)
    ax.set_title('2 component PCA', fontsize=20)
    ax.scatter(principalComponents[:, 0], principalComponents[:, 1], c="g")
    if annotate:
        for i, txt in enumerate(np.array(range(principalComponents.shape[0])).astype('str')):
            ax.annotate(txt, (principalComponents[i, 0], principalComponents[i, 1]))
    ax.grid()


def compare_pca(datamatrix1, datamatrix2, datamatrix3, n_components=2, annotate1=False, annotate2=False,
                annotate3=False):
    fig = plt.figure(figsize=(20, 10))
    j = 1
    for (datamatrix, annotate) in zip([datamatrix1, datamatrix2, datamatrix3], [annotate1, annotate2, annotate3]):
        pca = PCA(n_components=n_components)
        datamatrix = normalize(datamatrix)
        principalComponents = pca.fit_transform(datamatrix)
        print(principalComponents.shape)

        ax = fig.add_subplot(1, 3, j)
        ax.set_xlabel('Principal Component 1', fontsize=15)
        ax.set_ylabel('Principal Component 2', fontsize=15)
        ax.set_title('2 component PCA - datamatrix' + str(j), fontsize=20)
        ax.scatter(principalComponents[:, 0], principalComponents[:, 1], c="g")
        if annotate:
            for i, txt in enumerate(np.array(range(principalComponents.shape[0])).astype('str')):
                ax.annotate(txt, (principalComponents[i, 0], principalComponents[i, 1]))
        ax.grid()
        j += 1


def plot_clustering(datamatrix, labels, plt_colors=list(colors._colors_full_map.values()), title=None):
    pca = PCA(n_components=2)
    principalComponents = pca.fit_transform(datamatrix)

    fig = plt.figure(figsize=(17, 17))
    ax = fig.add_subplot(1, 1, 1)
    ax.set_xlabel('Principal Component 1', fontsize=15)
    ax.set_ylabel('Principal Component 2', fontsize=15)
    ax.set_title('2 component PCA', fontsize=20)
    targets = range(np.max(labels) + 1)
    colors = plt_colors[:np.max(labels) + 1]
    for target, color in zip(targets, colors):
        indicesToKeep = (labels == target)
        ax.scatter(principalComponents[indicesToKeep, 0]
                   , principalComponents[indicesToKeep, 1]
                   , c=color
                   , s=50)
    ax.legend(targets)
    ax.grid()


def plot_dendrogram(datamatrix, **kwargs):
    # create model
    model = AgglomerativeClustering(linkage='ward', distance_threshold=0, n_clusters=None)
    model = model.fit(datamatrix)

    # create the counts of samples under each node
    counts = np.zeros(model.children_.shape[0])
    n_samples = len(model.labels_)
    for i, merge in enumerate(model.children_):
        current_count = 0
        for child_idx in merge:
            if child_idx < n_samples:
                current_count += 1  # leaf node
            else:
                current_count += counts[child_idx - n_samples]
        counts[i] = current_count

    linkage_matrix = np.column_stack([model.children_, model.distances_, counts]).astype(float)

    # Plot the corresponding dendrogram
    plt.figure(figsize=(20, 20))
    plt.title('Hierarchical Clustering Dendrogram')
    dendrogram(linkage_matrix, **kwargs)
    plt.xlabel("Number of points in node (or index of point if no parenthesis).")
    a, b = plt.yticks()
    plt.yticks(np.arange(0, a[-1], 10))
    plt.xticks(fontsize=14, rotation=90)
    plt.show()


#### general functions #####################################################################################################

def print_full(x):
    pd.set_option('display.max_rows', len(x))
    print(x)
    pd.reset_option('display.max_rows')
    return np.nan
