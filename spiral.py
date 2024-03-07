#!/var/www/html/SPIRAL.web.tool/spiral_venv/bin/python3.10

from numba.core.errors import NumbaDeprecationWarning, NumbaPendingDeprecationWarning
import warnings
warnings.simplefilter('ignore', category=NumbaDeprecationWarning)
warnings.simplefilter('ignore', category=NumbaPendingDeprecationWarning)

import platform
import matplotlib
if platform.system() == 'Windows':
    matplotlib.use('TkAgg')

import shutil
import pathlib
import gzip
from SPIRAL_pipeline_funcs import *

def dataset_number(path):
    existing_folders = [d for d in os.listdir(path) if os.path.isdir(os.path.join(path, d))]
    existing_folders = [d for d in existing_folders if d[:4] == 'data' and len(d) > 4]
    if existing_folders:
        new_n = max([int(d[4:]) for d in existing_folders]) + 1
    else:
        new_n = 1
    return new_n

ANALYSIS_FOLDER = './static/analysis'
STATIC_FOLDER = './static'

if __name__ == '__main__':
    # On Windows calling this function is necessary.
    mp.freeze_support()

    print("#############################################################################")
    print("########                                                             ########")
    print("########                           SPIRAL                            ########")
    print("########                                                             ########")
    print("#############################################################################")

    # create static and analysis folders if they do not exist
    if not os.path.exists(STATIC_FOLDER):
        os.mkdir(STATIC_FOLDER)
    if not os.path.exists(ANALYSIS_FOLDER):
        os.mkdir(ANALYSIS_FOLDER)

    # move the ensembl folder if needed
    if not os.path.exists('./ensembl') and os.path.exists('./_internal/ensembl'):
        shutil.move('./_internal/ensembl', './ensembl')

    # create a folder for the new dataset
    data_n = dataset_number(ANALYSIS_FOLDER)
    data_path = os.path.join(ANALYSIS_FOLDER, 'data' + str(data_n))
    os.mkdir(data_path)

    # save dataset name
    print("\n" + color.BLUE + "Please provide a name for your data set:" + color.END)

    dataset_name = input()
    with open(os.path.join(data_path, 'dataset_name.txt'), 'w') as text_file:
        text_file.write(dataset_name)

    # save species
    print("\n" + color.BLUE + "Please select a species:" + color.END)

    species_dict = {
        "ARABIDOPSIS_THALIANA": "Arabidopsis thaliana",
        "SACCHAROMYCES_CEREVISIAE": 'Saccharomyces cerevisiae',
        "CAENORHABDITIS_ELEGANS": 'Caenorhabditis elegans',
        "DROSOPHILA_MELANOGASTER": 'Drosophila melanogaster',
        "DANIO_RERIO": 'Danio rerio (Zebrafish)',
        "HOMO_SAPIENS": 'Homo sapiens',
        "MUS_MUSCULUS": 'Mus musculus',
        "RATTUS_NORVEGICUS": 'Rattus norvegicus',
        'other': 'other (GO enrichment analysis will not be performed)',
        'synthetic': 'synthetic (GO enrichment analysis will not be performed)',
    }

    species = None
    while species is None:
        print("Enter the number next to the correct species:")
        for index, val in enumerate(list(species_dict.values())):
            print("{} for {}".format(index + 1, val))
        try:
            chosen_index = int(input()) - 1
            if 0 <= chosen_index < len(species_dict):
                species = list(species_dict)[chosen_index]
        except:
            pass

    with open(os.path.join(data_path, 'species.txt'), 'w') as text_file:
        text_file.write(species)

    # save sample name
    samples_names = ["samples", "cells", "spots"]
    chosen_samples_names = None
    while chosen_samples_names is None:
        print("\n" + color.BLUE + "How are your samples called?" + color.END)
        for index, val in enumerate(samples_names):
            print("{} for {}".format(index + 1, val))
        try:
            chosen_index = int(input()) - 1
            if 0 <= chosen_index < len(samples_names):
                chosen_samples_names = samples_names[chosen_index]
        except:
            pass


    with open(os.path.join(data_path, 'samp_name.txt'), 'w') as text_file:
        text_file.write(chosen_samples_names)

    # choose type of data file
    data_types = ["csv or txt file", "feature-barcode matrix (as downloaded from the Visium website)"]
    chosen_data_type = None
    while chosen_data_type is None:
        print("\n" + color.BLUE + "What is the file type of your gene expression data set?" + color.END)
        for index, val in enumerate(data_types):
            print("{} for {}".format(index + 1, val))
        try:
            chosen_index = int(input()) - 1
            if 0 <= chosen_index < len(data_types):
                chosen_data_type = data_types[chosen_index]
        except:
            pass


    # save counts
    if chosen_data_type == "csv or txt file":
        print("\n" + color.BLUE + "Please upload the count matrix file to this folder: {}, and name the file 'counts'".format(data_path) + color.END)
        print("The file type should be either csv or txt.")
        print("Gene names should appear in the first column. To be compatible with GOrilla (for the enrichment step), "
              "the supported formats are: gene symbol (preferred), gene and protein RefSeq, Uniprot, Unigene and Ensembl.")
        print("Sample names should appear in the first row. If your dataset has labels (for example: 'day0','day1'...), "
              "the labels should appear at the beginning of each sample name, before an underscore. For example: "
              "'day1_ACGTGTGA'.")
        #print("If your data has less than 200 cells\\spots\\samples, SPIRAL will run on the original data. If your data has"
        #      " more than that, SPIRAL will run in a repcell-mode.")
        input('when done press enter')

        while not (os.path.isfile(os.path.join(data_path, 'counts.txt')) or
                   os.path.isfile(os.path.join(data_path, 'counts.csv'))):
            print("Something went wrong. please make sure you uploaded the file to the correct folder, and named it "
                  "correctly.")
            print(
                "Folder name should be: {}. The file name should be 'counts'. The file type should be either csv or txt.".
                format(data_path))
            input('when done press enter')

        # labels checkbox
        labels_checkbox = None
        while labels_checkbox is None:
            print("\n" + color.BLUE + "Does your data set has labels? Type yes/no" + color.END)
            labels_ans = input()
            if labels_ans.lower() in ['yes', 'no']:
                labels_checkbox = labels_ans.lower() == 'no'
                with open(os.path.join(data_path, 'labels_checkbox.txt'), 'w') as text_file:
                    text_file.write(str(labels_checkbox))

        # spatial coordinates
        spatial_flag = None
        while spatial_flag is None:
            print("\n" + color.BLUE + "Do you have a spatial coordinates file (for spatial datasets)? Type yes/no" + color.END)
            spatial_ans = input()
            if spatial_ans.lower() in ['yes', 'no']:
                spatial_flag = spatial_ans.lower() == 'yes'

        if spatial_flag:
            print("\n" + color.BLUE + "Please upload the spatial coordinates file to this folder: {}, and name the file 'spatial_coors'"
                  .format(data_path) + color.END)
            print("The file type should be either csv or txt.")
            print("The number of rows should be the same as the number of columns in the count matrix (a row for every "
                  "spot\\cell\\sample).")
            print("There should be two or three columns:")
            print("     Two column version (the order of rows should match the order of columns in the count matrix):")
            print("         First column: X coordinates.")
            print("         Second column: Y coordinates.")
            print("     Three column version:")
            print("         First column: Sample names. If the sample names fit the column names in the count matrix, we will "
                  "use that. Otherwise, we will assume that the order of rows match the order of columns in the count matrix.")
            print("         Second column: X coordinates.")
            print("         Third column: Y coordinates.")
            input('when done press enter')

            while not (os.path.isfile(os.path.join(data_path, 'spatial_coors.txt')) or
                       os.path.isfile(os.path.join(data_path, 'spatial_coors.csv'))):
                print(
                    "Something went wrong. please make sure you uploaded the files to the correct folder, and named it "
                    "correctly.")
                print("Folder name should be: {}. The file name should be 'spatial_coors.txt' or "
                      "'spatial_coors.csv'".format(data_path))
                input('when done press enter')

    elif chosen_data_type == "feature-barcode matrix (as downloaded from the Visium website)":
        print("\n" + color.BLUE + "Please upload three files to this folder: {}: 'matrix.mtx.gz', "
                           "'features.tsv.gz' and 'barcodes.tsv.gz'".format(data_path) + color.END)
        input('when done press enter')
        print("processing... This might take a few minutes")

        while not (os.path.isfile(os.path.join(data_path, 'matrix.mtx.gz')) and os.path.isfile(
                os.path.join(data_path, 'features.tsv.gz')) and os.path.isfile(
                os.path.join(data_path, 'barcodes.tsv.gz'))):
            print("Something went wrong. please make sure you uploaded the files to the correct folder, and named it "
                  "correctly.")
            print("Folder name should be: {}. The file names should be 'matrix.mtx.gz', "
                  "'features.tsv.gz' and 'barcodes.tsv.gz'".format(data_path))
            input('when done press enter')
            print("processing... This might take a few minutes")

        # transform into "counts.csv"
        # read in MEX format matrix as table
        mat_filtered = scipy.io.mmread(os.path.join(data_path, "matrix.mtx.gz"))
        # list of gene names, e.g. 'SAMD11'
        gene_names = [row[1] for row in csv.reader(gzip.open(
            os.path.join(data_path, "features.tsv.gz"), mode="rt"), delimiter="\t")]
        # list of barcodes, e.g. 'AAACATACAAAACG-1'
        barcodes = [row[0] for row in csv.reader(gzip.open(
            os.path.join(data_path, "barcodes.tsv.gz"), mode="rt"), delimiter="\t")]
        # transform table to pandas dataframe and label rows and columns
        matrix_df = pd.DataFrame.sparse.from_spmatrix(mat_filtered)
        matrix_df.columns = barcodes
        matrix_df.insert(loc=0, column="gene", value=gene_names)
        # save the table as a CSV (note the CSV will be a very large file)
        matrix_df.to_csv(os.path.join(data_path, "counts.csv"), index=False)

        # labels checkbox
        # there are no labels in this type of data
        with open(os.path.join(data_path, 'labels_checkbox.txt'), 'w') as text_file:
            text_file.write(str(True))

        # spatial coordinates
        spatial_flag = None
        while spatial_flag is None:
            print("\n" + color.BLUE + "Do you have a spatial coordinates file (for spatial datasets)? Type yes/no" + color.END)
            spatial_ans = input()
            if spatial_ans.lower() in ["yes", "no"]:
                spatial_flag = spatial_ans.lower() == "yes"

        if spatial_flag:
            spatial_types = ["Visium: columns are 'barcode', 'in_tissue', 'array_row', 'array_col', "
                             "'pxl_row_in_fullres', 'pxl_col_in_fullres'",
                             "other format "]
            chosen_spatial_type = None
            while chosen_spatial_type is None:
                print(
                    color.BLUE + "Is your spatial coordinates file designed as in Visium or has a different format?" +
                    color.END)
                for index, val in enumerate(spatial_types):
                    print("{} for {}".format(index + 1, val))
                try:
                    chosen_index = int(input()) - 1
                    if 0 <= chosen_index < len(spatial_types):
                        chosen_spatial_type = spatial_types[chosen_index][:spatial_types[chosen_index].find(":")]
                except:
                    pass

            if chosen_spatial_type == "Visium":
                print("\n" +
                    color.BLUE + "Please upload Visium's 'tissue_positions.csv' file to this folder: {}"
                    .format(data_path) + color.END)
                print("This file is usually located inside the 'Spatial imaging data' folder, that can be downloaded "
                      "for each dataset from the Visium website.")
                input('when done press enter')

                while not os.path.isfile(os.path.join(data_path, 'tissue_positions.csv')):
                    print(
                        "Something went wrong. please make sure you uploaded the file to the correct folder, and named it "
                        "correctly.")
                    print(
                        "Folder name should be: {}. The file name should be 'tissue_positions.csv'.".format(data_path))
                    input('when done press enter')

                spatial_coors = pd.read_csv(os.path.join(data_path, 'tissue_positions.csv'))
                spatial_coors.drop(spatial_coors.iloc[:, 1:4], axis=1, inplace=True)
                spatial_coors.to_csv(os.path.join(data_path, "spatial_coors.csv"), index=False)

            elif chosen_spatial_type == "other format":
                print("\n" +
                    color.BLUE + "Please upload the spatial coordinates file to this folder: {}, and name the file 'spatial_coors'"
                    .format(data_path) + color.END)
                print("The file type should be either csv or txt.")
                print(
                    "The number of rows should be the same as the number of columns in the count matrix (a row for every "
                    "spot\\cell\\sample).")
                print("There should be two or three columns:")
                print(
                    "     Two column version (the order of rows should match the order of columns in the count matrix):")
                print("         First column: X coordinates.")
                print("         Second column: Y coordinates.")
                print("     Three column version:")
                print(
                    "         First column: Sample names. If the sample names fit the column names in the count matrix, we will "
                    "use that. Otherwise, we will assume that the order of rows match the order of columns in the count matrix.")
                print("         Second column: X coordinates.")
                print("         Third column: Y coordinates.")
                input('when done press enter')

                while not (os.path.isfile(os.path.join(data_path, 'spatial_coors.csv')) or
                           os.path.isfile(os.path.join(data_path, 'spatial_coors.txt'))):
                    print(
                        "Something went wrong. please make sure you uploaded the file to the correct folder, and named it "
                        "correctly.")
                    print(
                        "Folder name should be: {}. The file name should be 'spatial_coors.csv' or 'spatial_coors.txt'.".format(data_path))
                    input('when done press enter')

    # offer the user to use non-default parameters
    num_stds_thresh_lst_default = [0.5, 0.75, 1.0]
    num_stds_thresh_lst = num_stds_thresh_lst_default
    mu_lst_default = [0.9, 0.95]
    mu_lst = mu_lst_default
    print("\n" + color.BLUE + "We normally run SPIRAL with different sets of parameters to achieve maximal performance. "
                       "The recommended sets of parameters are:" + color.END)
    print('\N{GREEK SMALL LETTER ALPHA}', ': ', num_stds_thresh_lst_default)
    print('\N{GREEK SMALL LETTER MU}', ': ', mu_lst_default)
    change_params_flag = None
    while change_params_flag is None:
        print(color.BLUE + "Would you like to change these default parameters? Type yes/no" + color.END)
        change_params_ans = input()
        if change_params_ans.lower() in ["yes", "no"]:
            change_params_flag = change_params_ans.lower() == "yes"

    if change_params_flag:
        # change alpha
        change_alpha_flag = None
        while change_alpha_flag is None:
            print("\n" + color.BLUE + "Would you like to change the values of " + '\N{GREEK SMALL LETTER ALPHA}' + "? Type yes/no"
                  + color.END)
            change_alpha_ans = input()
            if change_alpha_ans.lower() in ["yes", "no"]:
                change_alpha_flag = change_alpha_ans.lower() == "yes"

        if change_alpha_flag:
            print("\n" + color.BLUE + "Type decimal numbers between 0 and 5 to be used as " + '\N{GREEK SMALL LETTER ALPHA}'
                  + ". After every number press enter. When done type 'done'."
                  + color.END)
            print("Note that every additional " + '\N{GREEK SMALL LETTER ALPHA}' + " extends the running time.")
            done = False
            num_stds_thresh_lst = []
            while not done:
                inp = input()
                try:
                    if inp.lower() == "done":
                        done = True
                    else:
                        inp = float(inp)
                        if inp >= 0 and inp <= 5:
                            if inp not in num_stds_thresh_lst:
                                num_stds_thresh_lst += [inp]
                        else:
                            print("Please type a decimal number between 0 and 5. When done type 'done'.")
                except:
                    print("We could not read the last number. Please type a decimal number between 0 and 5. "
                          "When done type 'done'.")

            if not num_stds_thresh_lst: # if the list is empty
                print("\n" + color.BLUE + "We could not capture any legal values for " + '\N{GREEK SMALL LETTER ALPHA}'
                      + ". SPIRAL will run with the default values."
                      + color.END)
                num_stds_thresh_lst = num_stds_thresh_lst_default
            else:
                print("\n" + color.BLUE + "SPIRAL will run with the following set of values for "
                      + '\N{GREEK SMALL LETTER ALPHA}' + ":"
                      + color.END)
                num_stds_thresh_lst.sort()
                print('\N{GREEK SMALL LETTER ALPHA}', ': ', num_stds_thresh_lst)

        # change mu
        change_mu_flag = None
        while change_mu_flag is None:
            print("\n" +
                color.BLUE + "Would you like to change the values of " + '\N{GREEK SMALL LETTER MU}' + "? Type yes/no"
                + color.END)
            change_mu_flag = input()
            if change_mu_flag.lower() in ["yes", "no"]:
                change_mu_flag = change_mu_flag.lower() == "yes"

        if change_mu_flag:
            print("\n" + color.BLUE + "Type decimal numbers between 0 and 1 to be used as " + '\N{GREEK SMALL LETTER MU}'
                  + ". After every number press enter. When done type 'done'."
                  + color.END)
            print("Note that every additional " + '\N{GREEK SMALL LETTER MU}' + " extends the running time.")
            done = False
            mu_lst = []
            while not done:
                inp = input()
                try:
                    if inp.lower() == "done":
                        done = True
                    else:
                        inp = float(inp)
                        if inp >= 0 and inp <= 1:
                            if inp not in mu_lst:
                                mu_lst += [inp]
                        else:
                            print("Please type a decimal number between 0 and 1. When done type 'done'.")
                except:
                    print("We could not read the last number. Please type a decimal number between 0 and 1. "
                          "When done type 'done'.")

            if not mu_lst: # if the list is empty
                print("\n" + color.BLUE + "We could not capture any legal values for " + '\N{GREEK SMALL LETTER MU}'
                      + ". SPIRAL will run with the default values."
                      + color.END)
                mu_lst = mu_lst_default
            else:
                print("\n" + color.BLUE + "SPIRAL will run with the following set of values for "
                      + '\N{GREEK SMALL LETTER MU}' + ":"
                      + color.END)
                mu_lst.sort()
                print('\N{GREEK SMALL LETTER MU}', ': ', mu_lst)

    ensembl_folder = ensembl_loc(production=False, hostname=None)
    load_data_first_time(analysis_folder=ANALYSIS_FOLDER, data_n=data_n, median_count_normalization_flag=True,
                         with_log=False)
    compute_violin_plots(analysis_folder=ANALYSIS_FOLDER, data_n=data_n, static_path=STATIC_FOLDER, species=species,
                         local_run=True, ensembl_folder=ensembl_folder)
    outcome = run_SPIRAL_pipeline(analysis_folder=ANALYSIS_FOLDER, data_n=data_n, species=species,
                                  min_nFeatures=None, max_nFeatures=None, max_mtpercent=None,
                                  numerical_shapes=None, local_run=True,
                                  num_stds_thresh_lst=num_stds_thresh_lst, mu_lst=mu_lst,
                                  num_iters_lst=[10000], path_len_lst=[3], ensembl_folder=ensembl_folder)
    with open(outcome_path(data_path), 'w') as text_file:
        text_file.write(str(outcome))

    ###################################################################################################
    # design the final folder
    files_in_data_path = set([f for f in os.listdir(data_path) if os.path.isfile(os.path.join(data_path, f))])
    impute_method = open(os.path.join(data_path, 'imputation_method.txt'), "r").read()

    ###### spiral_results.zip ########
    # choose files for spiral_results.zip
    files_to_zip = set(
        # needed files for the results panel:
        # results excel file
        # Jaccard matrix for the slider
        # gene table and structures, in case gene lists were too long for the excel file
        [f for f in files_in_data_path if f in
         ["dataset_name.txt", "GO_flag.txt", "imputation_method.txt", "Jaccard_mat_genes.npy",
          "Jaccard_thr_genes.txt", "outcome.txt", os.path.basename(final_sig_filename(data_path)),
          os.path.basename(genetable_file_name(data_path, impute_method))] or f[-12:] == "_structs.txt"])

    # create "spiral_results" folder and move needed files and folders
    os.mkdir(os.path.join(data_path, "spiral_results"))
    for f in files_to_zip:
        shutil.move(os.path.join(data_path, f), os.path.join(data_path, "spiral_results", f))
    if os.path.exists(os.path.join(data_path, "structure_layouts")):
        shutil.move(os.path.join(data_path, "structure_layouts"),
                    os.path.join(data_path, "spiral_results", "structure_layouts"))

    # zip the "spiral_results" folder into "spiral_results.zip"
    directory = pathlib.Path(os.path.join(data_path, "spiral_results"))
    with ZipFile(os.path.join(data_path, "spiral_results.zip"), mode="w") as archive:
        for file_path in directory.rglob("*"):
            archive.write(file_path, arcname=file_path.relative_to(directory))

    # delete the "spiral_results" folder
    shutil.rmtree(os.path.join(data_path, "spiral_results"))

    ###### spiral_aux_files.zip ########
    # choose files for spiral_aux_files.zip
    aux_files_to_zip = set([f for f in files_in_data_path if "PCA" in f or "UMAP" in f or "counts_norm" in f
                            or f in ["labels_checkbox.txt", "max_nFeatures.txt", "min_nFeatures.txt",
                                     "MT_list.txt", os.path.basename(with_mt_filename_(data_path)),
                                     "max_mtpercent.txt", "samp_name.txt", "species.txt"]])

    # create "spiral_aux_files" folder and move needed files and folders
    os.mkdir(os.path.join(data_path, "spiral_aux_files"))
    for f in aux_files_to_zip:
        shutil.move(os.path.join(data_path, f), os.path.join(data_path, "spiral_aux_files", f))
    if os.path.exists(os.path.join(data_path, "repcell_partition")):
        shutil.move(os.path.join(data_path, "repcell_partition"),
                    os.path.join(data_path, "spiral_aux_files", "repcell_partition"))

    # zip the "spiral_aux_files" folder into "spiral_aux_files.zip"
    directory = pathlib.Path(os.path.join(data_path, "spiral_aux_files"))
    with ZipFile(os.path.join(data_path, "spiral_aux_files.zip"), mode="w") as archive:
        for file_path in directory.rglob("*"):
            archive.write(file_path, arcname=file_path.relative_to(directory))

    # delete the "spiral_aux_files" folder
    shutil.rmtree(os.path.join(data_path, "spiral_aux_files"))

    ###### delete unnecessary files ########
    data_files_to_keep = set([f for f in files_in_data_path if "mtx.gz" in f or "tsv.gz" in f or f in
                              ["counts.csv", "counts.txt", "spatial_coors.csv", "spatial_coors.txt",
                               "tissue_positions.csv"]])
    files_to_delete = [
        f for f in files_in_data_path if f not in data_files_to_keep.union(aux_files_to_zip).union(files_to_zip)]
    for f in files_to_delete:
        os.remove(os.path.join(data_path, f))

    print(color.BLUE + "SPIRAL results are available in this folder: {}".format(data_path) + color.END)
    exit()
