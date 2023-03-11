from SPIRAL_pipeline_funcs import *


def dataset_number(path):
    existing_folders = [d for d in os.listdir(path) if os.path.isdir(os.path.join(path, d))]
    existing_folders = [d for d in existing_folders if d[:4] == 'data' and len(d) > 4]
    if existing_folders:
        new_n = max([int(d[4:]) for d in existing_folders]) + 1
    else:
        new_n = 1
    return new_n


ANALYSIS_FOLDER = '.\\static\\analysis'
STATIC_FOLDER = '.\\static'

if __name__ == '__main__':
    print("#############################################################################")
    print("########                                                             ########")
    print("########                       SPIRAL Local run                      ########")
    print("########                                                             ########")
    print("#############################################################################")

    # create a folder for the new dataset
    data_n = dataset_number(ANALYSIS_FOLDER)
    data_path = os.path.join(ANALYSIS_FOLDER, 'data' + str(data_n))
    os.mkdir(data_path)

    # save dataset name
    print("What is you dataset name?")

    dataset_name = input()
    with open(os.path.join(data_path, 'dataset_name.txt'), 'w') as text_file:
        text_file.write(dataset_name)

    # save species
    print("Select a species:")

    species = {
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

    chosen_species = None
    while chosen_species is None:
        print("choose a species for the following options:")
        for index, val in enumerate(list(species.values())):
            print("{} for {}".format(index + 1, val))
        chosen_index = int(input())
        if 1 <= chosen_index <= len(species):
            chosen_species = list(species)[chosen_index]

    with open(os.path.join(data_path, 'species.txt'), 'w') as text_file:
        text_file.write(chosen_species)

    # save sample name
    samples_names = ["samples", "cells", "spots"]
    chosen_samples_names = None
    while chosen_samples_names is None:
        print("How are your samples called?")
        for index, val in enumerate(samples_names):
            print("{} for {}".format(index + 1, val))
        chosen_index = int(input())
        if 1 <= chosen_index <= len(samples_names):
            chosen_samples_names = samples_names[chosen_index]

    with open(os.path.join(data_path, 'samp_name.txt'), 'w') as text_file:
        text_file.write(chosen_samples_names)

    # save counts
    print("Please upload count matrix file to this folder:{}, and name the file 'counts'".format(data_path))
    print("The file type should be either csv or txt.")
    print("Gene names should appear in the first column. To be compatible with GOrilla (for the enrichment step), "
          "the supported formats are: gene symbol (preferred), gene and protein RefSeq, Uniprot, Unigene and Ensembl.")
    print("Sample names should appear in the first row. If your dataset has labels (for example: 'day0','day1'...) "
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
            "Folder name sould be: {}. The file name should be 'counts'. The file type should be either csv or txt.".
            format(data_path))
        input('when done press enter')

    # labels checkbox
    print("Does your data set has labels?")
    labels_checkbox = input('press 1 if it has labels, anything else otherwise.') == 1
    with open(os.path.join(data_path, 'labels_checkbox.txt'), 'w') as text_file:
        text_file.write(str(labels_checkbox))

    # spatial coordinates
    print("If you have spatial coordinates file, please upload it to this folder:{}, and name the file 'spatial_coors'"
          .format(data_path))
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

    load_data_first_time(analysis_folder=ANALYSIS_FOLDER, data_n=data_n, median_count_normalization_flag=True,
                         with_log=False)
    compute_violin_plots(analysis_folder=ANALYSIS_FOLDER, data_n=data_n, static_path=STATIC_FOLDER, species=species,
                         local_run=True)
    run_SPIRAL_pipeline(analysis_folder=ANALYSIS_FOLDER, data_n=data_n, species=species,
                        min_nFeatures=None, max_nFeatures=None, max_mtpercent=None,
                        numerical_shapes=False,
                        num_stds_thresh_lst=[0.5, 0.75, 1.0], mu_lst=[0.95], num_iters_lst=[10000],
                        path_len_lst=[3])

    print("Results are available in this folder:{}".format(data_path))
    exit()
