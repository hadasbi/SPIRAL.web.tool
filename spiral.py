#!/var/www/html/SPIRAL.web.tool/spiral_venv/bin/python3.10

from numba.core.errors import NumbaDeprecationWarning, NumbaPendingDeprecationWarning
import warnings
warnings.simplefilter('ignore', category=NumbaDeprecationWarning)
warnings.simplefilter('ignore', category=NumbaPendingDeprecationWarning)

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

    # create a folder for the new dataset
    data_n = dataset_number(ANALYSIS_FOLDER)
    data_path = os.path.join(ANALYSIS_FOLDER, 'data' + str(data_n))
    os.mkdir(data_path)

    # save dataset name
    print(color.BLUE + "Please provide a name for your data set:" + color.END)

    dataset_name = input()
    with open(os.path.join(data_path, 'dataset_name.txt'), 'w') as text_file:
        text_file.write(dataset_name)

    # save species
    print(color.BLUE + "Please select a species:" + color.END)

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
        chosen_index = int(input()) - 1
        if 0 <= chosen_index < len(species_dict):
            species = list(species_dict)[chosen_index]

    with open(os.path.join(data_path, 'species.txt'), 'w') as text_file:
        text_file.write(species)

    # save sample name
    samples_names = ["samples", "cells", "spots"]
    chosen_samples_names = None
    while chosen_samples_names is None:
        print(color.BLUE + "How are your samples called?" + color.END)
        for index, val in enumerate(samples_names):
            print("{} for {}".format(index + 1, val))
        chosen_index = int(input()) - 1
        if 0 <= chosen_index < len(samples_names):
            chosen_samples_names = samples_names[chosen_index]

    with open(os.path.join(data_path, 'samp_name.txt'), 'w') as text_file:
        text_file.write(chosen_samples_names)

    # save counts
    print(color.BLUE + "Please upload the count matrix file to this folder: {}, and name the file 'counts'".format(data_path) + color.END)
    print("The file type should be either csv or txt.")
    print("Gene names should appear in the first column. To be compatible with GOrilla (for the enrichment step), "
          "the supported formats are: gene symbol (preferred), gene and protein RefSeq, Uniprot, Unigene and Ensembl.")
    print("Sample names should appear in the first row. If your dataset has labels (for example: 'day0','day1'...), "
          "the labels should appear at the beginning of each sample name, before an underscore. For example: "
          "'day1_ACGTGTGA'.")
    print("If your data has less than 200 cells\\spots\\samples, SPIRAL will run on the original data. If your data has"
          " more than that, SPIRAL will run in a repcell-mode.")

    input('when done press enter')

    print("processing...")

    while not (os.path.isfile(os.path.join(data_path, 'counts.txt')) or os.path.isfile(
            os.path.join(data_path, 'counts.csv'))):
        print("Something went wrong. please make sure you uploaded the file to the correct folder, and named it "
              "correctly")
        print(
            "Folder name should be: {}. The file name should be 'counts'. The file type should be either csv or txt.".
            format(data_path))
        input('when done press enter')

    # labels checkbox
    print(color.BLUE + "Does your data set has labels? Type yes/no" + color.END)
    labels_ans = input()
    labels_checkbox = labels_ans not in ["yes", "Yes", "YES", "Y", "y"]
    #print("Does your data set has labels?")
    #labels_checkbox = input('press 1 if it has labels, anything else otherwise.') != '1'
    with open(os.path.join(data_path, 'labels_checkbox.txt'), 'w') as text_file:
        text_file.write(str(labels_checkbox))

    # spatial coordinates
    print(color.BLUE + "Do you have a spatial coordinates file (for spatial datasets)? Type yes/no" + color.END)
    spatial_ans = input()
    if spatial_ans in ["yes", "Yes", "YES", "Y", "y"]:
        print(color.BLUE + "Please upload the spatial coordinates file to this folder: {}, and name the file 'spatial_coors'"
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

    ensembl_folder = ensembl_loc(production=False, hostname=None)

    load_data_first_time(analysis_folder=ANALYSIS_FOLDER, data_n=data_n, median_count_normalization_flag=True,
                         with_log=False)
    compute_violin_plots(analysis_folder=ANALYSIS_FOLDER, data_n=data_n, static_path=STATIC_FOLDER, species=species,
                         local_run=True, ensembl_folder=ensembl_folder)
    print(color.BLUE + "SPIRAL is now runnning. This might take up to a few hours, depending on the computer power and "
                       "your data. When done, the results would appear in this folder: {}".format(data_path) + color.END)
    run_SPIRAL_pipeline(analysis_folder=ANALYSIS_FOLDER, data_n=data_n, species=species,
                        min_nFeatures=None, max_nFeatures=None, max_mtpercent=None,
                        numerical_shapes=None,
                        num_stds_thresh_lst=[0.5, 0.75, 1.0], mu_lst=[0.9, 0.95], num_iters_lst=[10000],
                        path_len_lst=[3], ensembl_folder=ensembl_folder)

    print(color.BLUE + "SPIRAL results are available in this folder: {}".format(data_path) + color.END)
    exit()
