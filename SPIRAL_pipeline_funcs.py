#!/home/yaellab/SPIRAL/bin/python3.8

import sys
import os.path
import seaborn as sns
from scipy.stats import pearsonr, norm
from scipy import sparse
from time import time
import pickle5 as pickle
import json
import codecs
import multiprocessing as mp
import SPIRAL_mp_funcs
from shutil import copyfile
from requests.exceptions import ConnectionError
import csv
from zipfile import ZipFile

from SPIRAL_basic_funcs import *
from SPIRAL_visualization_funcs import *
from SPIRAL_enrichment_funcs import *
from SPIRAL_design_excel import *


####################################################################################################################
def load_data_first_time(analysis_folder, data_n, median_count_normalization_flag=False, with_log=False):
    print('load_data_first_time!')

    ############## Load count matrix #############################
    # Load data (and filter out genes and cells\spots with zero reads)
    data_path = data_path_name(analysis_folder, data_n)
    data_loc = glob.glob(os.path.join(data_path, 'counts.*').replace('\\', '/'))[0].replace('\\', '/')

    # Search for the right delimiter
    data = pd.read_csv(data_loc, sep=None, index_col=0, header='infer', engine='python')
    print(data.shape)

    # Check size of count matrix
    num_cells_orig, num_genes_orig = data.shape[1], data.shape[0]

    # if num_cells_orig > 200 and num_cells_orig < 1000:
    #    raise ValueError('Number of spots\cells\samples in count matrix is in the range 200-1000.')

    ############## Load spatial coordinates (if exist) #############################
    spatial = False
    spatial_file_list = glob.glob(os.path.join(data_path, 'spatial_coors.*').replace('\\', '/'))
    if len(spatial_file_list) > 1:
        print('ERROR: too many spatial coordinates files.')
    if spatial_file_list:  # if spatial coordinates file exists
        spatial_loc = spatial_file_list[0].replace('\\', '/')
        spatial = True

    if spatial:
        index_col = 0
        spatial_coors = pd.read_csv(spatial_loc, sep=None, index_col=index_col, header='infer', engine='python')
        print(spatial_coors.shape)

        reload = False
        if spatial_coors.shape[1] == 1:
            index_col = False
            spatial_coors = pd.read_csv(spatial_loc, sep=None, index_col=index_col, header='infer', engine='python')
            print(spatial_coors.shape)

        # Correct for wrong guess of header
        if spatial_coors.shape[0] == data.shape[1] - 1:
            spatial_coors = pd.read_csv(spatial_loc, sep=None, index_col=index_col, header=None)
            print(spatial_coors.shape)
        elif spatial_coors.shape[0] == data.shape[1] + 1:
            spatial_coors = pd.read_csv(spatial_loc, sep=None, index_col=index_col, header=0)
            print(spatial_coors.shape)

        # Check that the sample names in spatial_coors fit the sample names in data
        if spatial_coors.shape[0] != data.shape[1]:
            print(set(list(data)) - set(spatial_coors.index))
            print(set(spatial_coors.index) - set(list(data)))
            raise ValueError(
                'Number of rows in the spatial coordinates file do not fit the number of columns in the count matrix')
        else:
            # If not according to name, then according to order:
            if set(spatial_coors.index) != set(list(data)):
                print('This should be empty:', set(spatial_coors.index).intersection(set(list(data))))
                spatial_coors.index = list(data)

    ############## delete genes and cells with zero reads ##########################
    data = delete_genes_and_cells_with_zero_reads(data)

    if spatial:
        spatial_coors = spatial_coors.loc[list(data), :]

    ############# Check for duplicated genes names #################################
    print(data.shape)
    data = data.groupby(level=0).sum()
    data = data.T.groupby(level=0).sum().T
    print(data.shape)

    ############# median count normalization #######################################
    if median_count_normalization_flag:
        if with_log:
            print('Normalizing to median count and then: log10(x+1)')
        else:
            print('Normalizing to median count (without log normalization)')
        data = median_count_normalization(data, with_log)
    data.index.name = None

    ############# save to file #####################################################
    new_data_loc = data_norm_loc_name(data_path)
    data.to_csv(new_data_loc, sep='\t')

    if spatial:
        new_spatial_loc = spatial_norm_loc_name(data_path)
        spatial_coors.to_csv(new_spatial_loc, sep='\t')

    # Check size of count matrix
    num_cells, num_genes = data.shape[1], data.shape[0]

    ############# check labels #####################################################
    label_set = set(compute_orig_deffs(data_path, data))
    label_set = sort_labels(label_set)
    nlabels = len(label_set)
    if nlabels == num_cells:
        nlabels = 1

    return spatial, num_cells_orig, num_genes_orig, num_cells, num_genes, label_set, nlabels


####################################################################################################################
def compute_violin_plots(analysis_folder, data_n, static_path, species):
    print('compute_violin_plots!')
    # Load data
    data_path = data_path_name(analysis_folder, data_n)
    data_norm_loc = data_norm_loc_name(data_path)
    data = pd.read_csv(data_norm_loc, index_col=0, sep='\t')

    nFeatures = data.astype('bool').sum(axis='index')
    nFeatures.name = 'number of expressed genes'

    # save violin plot of number of expressed genes
    summary = pd.concat([nFeatures], axis=1)
    fig, axs = plt.subplots(1, 1, figsize=(5, 5))
    axs.set_title("number of expressed genes")
    sns.violinplot(y="number of expressed genes", data=summary, ax=axs)

    # save to file
    vln_plot_filename = vln_plot_filename_(static_path=static_path, data_n=data_n)
    plt.tight_layout()
    plt.savefig(vln_plot_filename)
    plt.close()

    # save violin plot of percent of mitochondrial genes
    MTgenes, error = get_MTgenes(data_path, list(data.index), species)
    with_mt = False
    if MTgenes:
        with_mt = True
        mtpercent = np.divide(data.loc[MTgenes, :].sum(), data.sum()) * 100
        mtpercent.name = 'mitochondrial percent of reads'

        summary = pd.concat([mtpercent], axis=1)
        fig, axs = plt.subplots(1, 1, figsize=(5, 5))
        axs.set_title("mitochondrial percent of reads")
        sns.violinplot(y="mitochondrial percent of reads", data=summary, ax=axs)

        # save to file
        vln_plot_mt_filename = vln_plot_mt_filename_(static_path=static_path, data_n=data_n)
        plt.tight_layout()
        plt.savefig(vln_plot_mt_filename)
        plt.close()

    return with_mt, error


####################################################################################################################
def run_SPIRAL_pipeline(analysis_folder, data_n, species=None,
                        min_nFeatures=None, max_nFeatures=None, max_mtpercent=None, numerical_shapes=None,
                        num_stds_thresh_lst=[0.5, 1], mu_lst=[0.95], num_iters_lst=[10000], path_len_lst=[3]):
    # impute_method options: 'agg_wald', 'IaconoClus_dim50', 'IaconoClus', 'agg_wald_opt'

    print('run_SPIRAL_pipeline!')

    print('Algorithm parameters:')
    print('num_stds_thresh_lst:', num_stds_thresh_lst)
    print('mu_lst:', mu_lst)
    print('path_len_lst:', path_len_lst)
    print('num_iters_lst:', num_iters_lst)

    data_path = data_path_name(analysis_folder, data_n)
    data_norm_loc = data_norm_loc_name(data_path)

    ###################### check if spatial coordinates were loaded #############
    spatial = False
    spatial_norm_loc = spatial_norm_loc_name(data_path)
    if os.path.exists(spatial_norm_loc):  # if spatial coordinates file exists
        spatial = True

    ###################### load sample name #####################################
    real_samp_name = open(os.path.join(data_path, 'samp_name.txt'), "r").read()[:-1]

    ######################  set algorithm parameters    #########################
    # no_iterations = False
    use_mp_to_find_structs = True
    errors = []

    seed_mode = 'random_path'  # 'gene' (original method) or 'random_cellpair_set' or 'random_path'

    if seed_mode == 'gene':
        random_cellpair_set, random_path, path_len, num_iters = False, False, None, None
    elif seed_mode == 'random_cellpair_set':
        random_cellpair_set, random_path, path_len, num_iters = True, False, None, 10000
    elif seed_mode == 'random_path':
        random_cellpair_set, random_path, path_len, num_iters = False, True, 3, 10000
    else:
        print('Seed error!')

    min_cells_per_cluster = 10

    use_5000 = False

    data = None
    structs = None
    sigtable_merged = None
    sigtable_filt = None
    sigtable_GO = None

    if species is None:
        species = open(os.path.join(data_path, 'species.txt'), "r").read()
    if species in ['synthetic', 'other']:
        skip_GO = True
    else:
        skip_GO = False

    ######################  run pipeline    #####################################
    # Check if data filtering was already performed, if not- perform data filtering
    data_norm_filt_loc = data_norm_filt_loc_name(data_path)
    spatial_norm_filt_loc = spatial_norm_filt_loc_name(data_path)
    if not os.path.exists(data_norm_filt_loc) or (spatial and (not os.path.exists(spatial_norm_filt_loc))):
        data = filter_data(data_path=data_path,
                           min_nFeatures=min_nFeatures, max_nFeatures=max_nFeatures,
                           max_mtpercent=max_mtpercent, spatial=spatial,
                           species=species,
                           data_norm_loc=data_norm_loc,
                           spatial_norm_loc=spatial_norm_loc,
                           data_norm_filt_loc=data_norm_filt_loc,
                           spatial_norm_filt_loc=spatial_norm_filt_loc)
    else:
        data = pd.read_csv(data_norm_filt_loc, index_col=0, sep='\t')

    # decide on an imputation method
    if data.shape[1] < 200:
        impute_method = 'no_imputation'
    else:
        impute_method = 'agg_wald'

    with open(os.path.join(data_path, 'imputation_method.txt'), 'w') as text_file:
        text_file.write(impute_method)

    # check if labels were already saved to file, if not- save labels to file
    origdeffsfile = os.path.join(data_path, 'orig_deffs.txt')
    if not os.path.exists(origdeffsfile):
        save_orig_deffs_to_file(data=data, data_path=data_path, data_norm_filt_loc=data_norm_filt_loc,
                                origdeffsfile=origdeffsfile)

    # Decide whether to use regular markers or numerical markers in the structures' layouts
    if numerical_shapes is None:
        numerical_shapes = decide_on_numerical_shapes(origdeffsfile)

    # check if layouts of data already exist (if not, compute them)
    norm_filt_pca_coor_file = os.path.join(data_path, 'norm_filt_PCA_coor' + ('5000' * use_5000) + '.csv')
    norm_filt_umap_coor_file = os.path.join(data_path, 'norm_filt_UMAP_coor' + ('5000' * use_5000) + '.csv')
    norm_filt_picfile_pca = os.path.join(data_path, 'norm_filt_PCA' + ('5000' * use_5000) + '.jpg')
    norm_filt_picfile_umap = os.path.join(data_path, 'norm_filt_UMAP' + ('5000' * use_5000) + '.jpg')
    norm_filt_picfile_spatial = os.path.join(data_path, 'norm_filt_spatial.jpg')
    if (not os.path.exists(norm_filt_pca_coor_file)) or (not os.path.exists(norm_filt_umap_coor_file)):
        save_layout_of_orig(norm_filt_pca_coor_file=norm_filt_pca_coor_file,
                            norm_filt_umap_coor_file=norm_filt_umap_coor_file,
                            norm_filt_picfile_pca=norm_filt_picfile_pca,
                            norm_filt_picfile_umap=norm_filt_picfile_umap,
                            norm_filt_picfile_spatial=norm_filt_picfile_spatial,
                            origdeffsfile=origdeffsfile, spatial_norm_filt_loc=spatial_norm_filt_loc,
                            data=data, data_norm_filt_loc=data_norm_filt_loc, data_path=data_path,
                            save_umap=True, save_pca=True, save_spatial=spatial, use_5000=use_5000,
                            numerical_shapes=numerical_shapes)

    # check if repcell computation was already performed, if not- compute repcells
    clustering_file_initial = os.path.join(data_path, impute_method + '_clustering.txt')
    clustering_file_final = os.path.join(data_path, impute_method + '_clustering_final.txt')
    if impute_method != 'no_imputation' and not os.path.exists(clustering_file_final):
        compute_repcells(data=data, impute_method=impute_method,
                         min_cells_per_cluster=min_cells_per_cluster,
                         data_norm_filt_loc=data_norm_filt_loc,
                         clustering_file_initial=clustering_file_initial,
                         clustering_file_final=clustering_file_final)

    # check if an imputed file already exists (if not, compute it)
    repcells_data_loc = os.path.join(data_path, impute_method + '_repcells_counts.txt')
    if impute_method != 'no_imputation':
        if not os.path.exists(repcells_data_loc):
            repcells_data, num_repcells, num_genes = create_imputed_file(data=data,
                                                                         impute_method=impute_method,
                                                                         data_norm_filt_loc=data_norm_filt_loc,
                                                                         clustering_file_final=clustering_file_final,
                                                                         repcells_data_loc=repcells_data_loc)
        else:
            repcells_data = pd.read_csv(repcells_data_loc, index_col=0, sep='\t')
            num_repcells = repcells_data.shape[1]  # number of cells
            num_genes = repcells_data.shape[0]  # number of genes

    # check if layouts of inmputed file already exist (if not, compute them)
    layout_of_imputed_picfile_PCA = os.path.join(data_path,
                                                 impute_method + '_repcells_PCA' + ('5000' * use_5000) + '.jpg')
    layout_of_imputed_picfile_UMAP = os.path.join(data_path,
                                                  impute_method + '_repcells_UMAP' + ('5000' * use_5000) + '.jpg')
    repcelldeffsfile = os.path.join(data_path, impute_method + '_repcells_deffs.txt')
    stepsfile = os.path.join(data_path, impute_method + '_repcells_steps.txt')
    repcellsumapcoorfile = os.path.join(data_path, impute_method + '_repcells_UMAP' + ('5000' * use_5000) + '.csv')
    repcellspcacoorfile = os.path.join(data_path, impute_method + '_repcells_PCA' + ('5000' * use_5000) + '.csv')

    if impute_method != 'no_imputation':
        if not os.path.exists(layout_of_imputed_picfile_PCA) or not os.path.exists(layout_of_imputed_picfile_UMAP):
            save_layout_of_imputed(repcells_data=repcells_data, impute_method=impute_method,
                                   with_legend=True, with_title=False,
                                   repcells_data_loc=repcells_data_loc,
                                   layout_of_imputed_picfile_PCA=layout_of_imputed_picfile_PCA,
                                   layout_of_imputed_picfile_UMAP=layout_of_imputed_picfile_UMAP,
                                   use_5000=use_5000,
                                   data_path=data_path,
                                   repcellsumapcoorfile=repcellsumapcoorfile,
                                   repcellspcacoorfile=repcellspcacoorfile,
                                   repcelldeffsfile=repcelldeffsfile,
                                   origdeffsfile=origdeffsfile,
                                   stepsfile=stepsfile,
                                   with_steps=False,
                                   clustering_file_final=clustering_file_final,
                                   data_norm_filt_loc=data_norm_filt_loc,
                                   numerical_shapes=numerical_shapes)

    # check if visualizations of the repcell partition already exist (on a UMAP, PCA, spatial layouts).
    # If not, create them.
    # repcell partition folder
    repcell_part_folder = os.path.join(data_path, 'repcell_partition')
    if not os.path.exists(repcell_part_folder):
        os.makedirs(repcell_part_folder)
    repcell_partition_UMAP = os.path.join(repcell_part_folder, impute_method + '_repcell_part_UMAP.png')
    repcell_partition_PCA = os.path.join(repcell_part_folder, impute_method + '_repcell_part_PCA.png')
    repcell_partition_spatial = os.path.join(repcell_part_folder, impute_method + '_repcell_part_spatial.png')
    repcell_partition_zipfile = os.path.join(data_path, 'repcell_partition.zip')
    if impute_method != 'no_imputation':
        if (not os.path.exists(repcell_partition_UMAP)) or (not os.path.exists(repcell_partition_PCA)) or (
                spatial and (not os.path.exists(repcell_partition_spatial))
        ):
            visualize_repcell_partition(clustering_file_final=clustering_file_final,
                                        norm_filt_umap_coor_file=norm_filt_umap_coor_file,
                                        norm_filt_pca_coor_file=norm_filt_pca_coor_file,
                                        spatial_norm_filt_loc=spatial_norm_filt_loc,
                                        use_5000=use_5000,
                                        origdeffsfile=origdeffsfile,
                                        repcell_partition_UMAP=repcell_partition_UMAP,
                                        repcell_partition_PCA=repcell_partition_PCA,
                                        repcell_partition_spatial=repcell_partition_spatial,
                                        data=data, data_norm_filt_loc=data_norm_filt_loc,
                                        save_umap=True, save_pca=True, save_spatial=spatial,
                                        with_legend=False, with_title=False,
                                        load_orig_umap_from_file=True, load_orig_pca_from_file=True,
                                        with_deffs=False, numerical_shapes=numerical_shapes)

            # create a zip file of all repcell partition layouts
            zip_all_files_in_folder(zipfile=repcell_partition_zipfile, folder=repcell_part_folder)

    # check if a gene table was already saved, if not- save it
    genetable_file = os.path.join(data_path, impute_method + '_genetable.p')
    if not os.path.exists(genetable_file):
        if impute_method != 'no_imputation':
            save_genetable_to_file(data=repcells_data, genetable_file=genetable_file)
        else:
            save_genetable_to_file(data=data, genetable_file=genetable_file)

    for num_stds_thresh in num_stds_thresh_lst:
        genes_in_sets_coo_file = os.path.join(data_path,
                                              impute_method + '_std_' + str(num_stds_thresh) + '_genes_in_sets_COO.npz')
        genes_in_sets_csc_file = os.path.join(data_path,
                                              impute_method + '_std_' + str(num_stds_thresh) + '_genes_in_sets_CSC.npz')
        genes_in_sets_csr_file = os.path.join(data_path,
                                              impute_method + '_std_' + str(num_stds_thresh) + '_genes_in_sets_CSR.npz')
        genes_in_sets_row_file = os.path.join(data_path,
                                              impute_method + '_std_' + str(num_stds_thresh) + '_genes_in_sets_row.npz')
        genes_in_sets_col_file = os.path.join(data_path,
                                              impute_method + '_std_' + str(num_stds_thresh) + '_genes_in_sets_col.npz')

        if not os.path.exists(genes_in_sets_coo_file):
            if impute_method != 'no_imputation':
                compute_sets_on_standardized_genes(genetable_file=genetable_file,
                                                   genes_in_sets_coo_file=genes_in_sets_coo_file,
                                                   genes_in_sets_csc_file=genes_in_sets_csc_file,
                                                   genes_in_sets_csr_file=genes_in_sets_csr_file,
                                                   genes_in_sets_row_file=genes_in_sets_row_file,
                                                   genes_in_sets_col_file=genes_in_sets_col_file,
                                                   counts=repcells_data, counts_file=repcells_data_loc,
                                                   num_stds_thresh=num_stds_thresh)
            else:
                compute_sets_on_standardized_genes(genetable_file=genetable_file,
                                                   genes_in_sets_coo_file=genes_in_sets_coo_file,
                                                   genes_in_sets_csc_file=genes_in_sets_csc_file,
                                                   genes_in_sets_csr_file=genes_in_sets_csr_file,
                                                   genes_in_sets_row_file=genes_in_sets_row_file,
                                                   genes_in_sets_col_file=genes_in_sets_col_file,
                                                   counts=data, counts_file=data_norm_filt_loc,
                                                   num_stds_thresh=num_stds_thresh)

        for mu in mu_lst:
            for path_len in path_len_lst:
                for num_iters in num_iters_lst:
                    # check if a structures file already exists (if not, compute it)
                    structs_file = structs_filename(data_path=data_path, impute_method=impute_method,
                                                    num_stds_thresh=num_stds_thresh, mu=mu,
                                                    path_len=path_len, num_iters=num_iters)

                    if not os.path.exists(structs_file):
                        if use_mp_to_find_structs:
                            structs = find_structures(genes_in_sets_npz_file=genes_in_sets_csr_file, mu=mu,
                                                      structs_file=structs_file, num_iters=num_iters, path_len=path_len,
                                                      use_mp_to_find_structs=use_mp_to_find_structs)

                    # check if a significance table file already exists (if not, compute it)
                    sigfile = sig_filename(data_path=data_path, impute_method=impute_method,
                                           num_stds_thresh=num_stds_thresh, mu=mu,
                                           path_len=path_len, num_iters=num_iters)

                    if not os.path.exists(sigfile):
                        if impute_method != 'no_imputation':
                            evaluate_significance_of_structs(sigfile=sigfile, genetable_file=genetable_file,
                                                             structs=structs, structs_file=structs_file,
                                                             counts=repcells_data, counts_file=repcells_data_loc,
                                                             impute_method=impute_method, real_samp_name=real_samp_name,
                                                             clustering_file_final=clustering_file_final)
                        else:
                            evaluate_significance_of_structs(sigfile=sigfile, genetable_file=genetable_file,
                                                             structs=structs, structs_file=structs_file,
                                                             counts=data, counts_file=data_norm_filt_loc,
                                                             impute_method=impute_method, real_samp_name=real_samp_name)

    # merge all sigtables
    sigfile_merged = os.path.join(data_path, impute_method + '_sigtable_merged.xlsx')
    if not os.path.exists(sigfile_merged):
        sigtable_merged = merge_sigtables(sigfile_merged=sigfile_merged, data_path=data_path,
                                          impute_method=impute_method, num_stds_thresh_lst=num_stds_thresh_lst,
                                          mu_lst=mu_lst, path_len_lst=path_len_lst, num_iters_lst=num_iters_lst)

    if sigtable_merged.empty:
        return 'No_structures'

    # check if a filtered significance table file already exists (if not, compute it)
    sigfile_filt = os.path.join(data_path, impute_method + '_sigtable_filt.xlsx')

    if not os.path.exists(sigfile_filt):
        sigtable_filt = filter_similar_structures(sigfile=sigfile_merged, significance_table=sigtable_merged,
                                                  sigfile_filt=sigfile_filt, genetable_file=genetable_file,
                                                  data_path=data_path, impute_method=impute_method,
                                                  real_samp_name=real_samp_name, struct_thr=0.8,
                                                  min_nstructs=3, max_nstructs=50)

    if sigtable_filt is None:
        return 'No_structures'

    # check if a filtered significance table file with GO terms and visualizations already exists (if not, compute it)
    sigfile_GO = os.path.join(data_path, impute_method + '_sigtable_filt_GO.xlsx')

    sigfile_GO_temp = sigfile_GO[:-5] + '_temp.xlsx'
    stopflag = False
    start_from_scratch = False

    errors = []
    GO_flag = True
    while not os.path.exists(sigfile_GO):
        try:
            if not skip_GO:
                sigtable_GO = add_GO_terms(sigtable=sigtable_filt, sigtable_file=sigfile_filt,
                                           sigfile_GO=sigfile_GO, sigfile_GO_temp=sigfile_GO_temp,
                                           genetable_file=genetable_file, data_path=data_path,
                                           species=species,
                                           impute_method=impute_method, start_from_scratch=False,
                                           pvals=[0.000001], real_samp_name=real_samp_name)
            else:
                copyfile(sigfile_filt, sigfile_GO)
                sigtable_GO = sigtable_filt
                GO_flag = False
        except ConnectionError as e:
            print(e)
            start_from_scratch = False
        except ConnectionResetError as e:
            print(e)
            start_from_scratch = False
        except:
            print(sys.exc_info()[0])
            errors.append('data ' + str(data_n) + ': ' + str(sys.exc_info()[0]))
            copyfile(sigfile_filt, sigfile_GO)
            sigtable_GO = sigtable_filt
            GO_flag = False
    print(errors)

    with open(os.path.join(data_path, 'GO_flag.txt'), 'w') as text_file:
        text_file.write(str(GO_flag))

    # visualize structures
    sigfile_vis = final_sig_filename(data_path, impute_method)
    picfolder = pic_folder(data_path)

    if os.path.exists(sigfile_GO) and (not os.path.exists(sigfile_vis)):
        sigtable_vis = visualize_structures(sigfile_vis=sigfile_vis, genetable_file=genetable_file,
                                            repcellsumapcoorfile=repcellsumapcoorfile,
                                            repcellspcacoorfile=repcellspcacoorfile,
                                            repcells_data=(
                                                repcells_data if impute_method != 'no_imputation' else None),
                                            repcells_data_loc=(
                                                repcells_data_loc if impute_method != 'no_imputation' else None),
                                            norm_filt_pca_coor_file=norm_filt_pca_coor_file,
                                            repcelldeffsfile=repcelldeffsfile,
                                            origdeffsfile=origdeffsfile,
                                            picfolder=picfolder,
                                            norm_filt_umap_coor_file=norm_filt_umap_coor_file,
                                            data=data, data_path=data_path,
                                            data_norm_filt_loc=data_norm_filt_loc,
                                            clustering_file_final=clustering_file_final,
                                            spatial_norm_filt_loc=spatial_norm_filt_loc,
                                            sigtable=sigtable_GO, sigtable_file=sigfile_GO,
                                            impute_method=impute_method,
                                            save_table=True,
                                            save_network_layouts=True, save_Gnetwork_layouts=True,
                                            save_UMAPs_orig=True, save_PCAs_orig=True,
                                            save_GUMAPs_orig=True, save_GPCAs_orig=True,
                                            save_UMAPs_repcells=(impute_method != 'no_imputation'),
                                            save_PCAs_repcells=(impute_method != 'no_imputation'),
                                            save_GUMAPs_repcells=(impute_method != 'no_imputation'),
                                            save_GPCAs_repcells=(impute_method != 'no_imputation'),
                                            save_spatial=spatial, save_Gspatial=spatial,
                                            with_legend=True, real_samp_name=real_samp_name,
                                            numerical_shapes=numerical_shapes)

    # design excel file
    if os.path.exists(sigfile_vis):
        design_excel(sigfile_vis)

    # create a zip file of all structure layouts
    layouts_zipfile = os.path.join(data_path, 'structure_layouts.zip')
    zip_all_files_in_folder(zipfile=layouts_zipfile, folder=picfolder)
    return 'Success'


####################################################################################################################
def filter_data(data_path, min_nFeatures, max_nFeatures, max_mtpercent, spatial, species,
                data_norm_loc, spatial_norm_loc, data_norm_filt_loc, spatial_norm_filt_loc):
    print('filter_data!')

    ##### Load filtering parameters if they were not received as arguments  ######
    if min_nFeatures is None:
        min_nFeatures = int(open(os.path.join(data_path, 'min_nFeatures.txt'), "r").read())
    if max_nFeatures is None:
        max_nFeatures = int(open(os.path.join(data_path, 'max_nFeatures.txt'), "r").read())

    with_mt = (open(with_mt_filename_(data_path), "r").read().lower() == 'true')

    ##### Filter the data set ######################################################
    data = pd.read_csv(data_norm_loc, index_col=0, sep='\t')
    data = data.loc[:, (data.astype(bool).sum() >= min_nFeatures) & (data.astype(bool).sum() <= max_nFeatures)]
    if with_mt:
        max_mtpercent = int(open(os.path.join(data_path, 'max_mtpercent.txt'), "r").read())
        MTgenes, mt_error = get_MTgenes(data_path, list(data.index), species)
        mtpercent = np.divide(data.loc[MTgenes, :].sum(), data.sum()) * 100
        data = data.loc[:, mtpercent <= max_mtpercent]

    if spatial:
        spatial_coors = pd.read_csv(spatial_norm_loc, index_col=0, sep='\t')
        spatial_coors = spatial_coors.loc[list(data), :]

    ##### save to file ############################################################
    data = median_count_normalization(table=data, with_log=False)
    data.to_csv(data_norm_filt_loc, sep='\t')

    if spatial:
        savetxt(spatial_norm_filt_loc, np.array(spatial_coors), delimiter=',')

    return data


####################################################################################################################
def save_genetable_to_file(data, genetable_file):
    # table with genes' names and indices
    genetable = pd.Series(index=data.index, data=np.array(range(data.shape[0])))
    with open(genetable_file, 'wb') as fp:
        pickle.dump(genetable, fp, protocol=pickle.HIGHEST_PROTOCOL)


####################################################################################################################
def renumber_clusters(labels, min_cells_per_cluster):
    all_clusters = set(labels)
    # remove clusters that represent less than min_cells_per_cluster cells
    small_clusters = [c for c in all_clusters if labels.count(c) < min_cells_per_cluster]
    final_cluster_list = list(all_clusters - set(small_clusters))
    # renumber the clusters such that the small ones will be at the end
    d = dict(zip(final_cluster_list, range(len(final_cluster_list))))
    d.update(dict(zip(small_clusters, np.arange(len(final_cluster_list), len(all_clusters)))))
    labels = [d[l] for l in labels]
    return labels


####################################################################################################################
def compute_lost_percent(orig_scRNA_data, min_cells_per_cluster, n_clusters):
    num_cells = orig_scRNA_data.shape[1]  # number of cells
    clustering = AgglomerativeClustering(linkage='ward', n_clusters=n_clusters).fit(orig_scRNA_data.transpose())
    labels = renumber_clusters(list(clustering.labels_), min_cells_per_cluster)
    all_clusters = set(labels)
    small_clusters = [c for c in all_clusters if labels.count(c) < min_cells_per_cluster]
    lost_percent = len([i for i in range(num_cells) if labels[i] in small_clusters]) / num_cells * 100
    return lost_percent


####################################################################################################################
def binary_search(orig_scRNA_data, minval, maxval, min_cells_per_cluster, max_lost_percent):
    if minval == maxval:
        return minval
    n_clusters = int(np.ceil((minval + maxval) / 2))
    lost_percent = compute_lost_percent(orig_scRNA_data=orig_scRNA_data, min_cells_per_cluster=min_cells_per_cluster,
                                        n_clusters=n_clusters)
    if lost_percent <= max_lost_percent:
        return binary_search(orig_scRNA_data=orig_scRNA_data, minval=n_clusters, maxval=maxval,
                             min_cells_per_cluster=min_cells_per_cluster, max_lost_percent=max_lost_percent)
    else:
        return binary_search(orig_scRNA_data=orig_scRNA_data, minval=minval, maxval=n_clusters - 1,
                             min_cells_per_cluster=min_cells_per_cluster, max_lost_percent=max_lost_percent)


####################################################################################################################
def find_opt_n_clusters(orig_scRNA_data, min_cells_per_cluster, max_lost_percent=3):
    print('find_opt_n_clusters!')
    num_cells = orig_scRNA_data.shape[1]  # number of cells
    max_n_clusters = min(110, int(np.round(num_cells / 30)))

    # first try with the maximal number of clusters
    n_clusters = max_n_clusters
    lost_percent = compute_lost_percent(orig_scRNA_data=orig_scRNA_data, min_cells_per_cluster=min_cells_per_cluster,
                                        n_clusters=n_clusters)
    if lost_percent <= max_lost_percent:
        return n_clusters

    # find a number t, such that if n_clusters==t then (lost_percent <= max_lost_percent), and if n_clusters==t+10
    # then (lost_percent > max_lost_percent)
    while lost_percent > max_lost_percent:
        n_clusters = max(n_clusters - 10, 2)
        lost_percent = compute_lost_percent(orig_scRNA_data=orig_scRNA_data,
                                            min_cells_per_cluster=min_cells_per_cluster,
                                            n_clusters=n_clusters)
        if n_clusters == 2:
            return 2

    # perform a binary search on the range (t, t+10) to find the largest t such that (lost_percent <= max_lost_percent)
    return binary_search(orig_scRNA_data=orig_scRNA_data, minval=n_clusters, maxval=n_clusters + 9,
                         min_cells_per_cluster=min_cells_per_cluster, max_lost_percent=max_lost_percent)


####################################################################################################################
def compute_repcells(clustering_file_initial, clustering_file_final, data=None, impute_method='agg_wald',
                     min_cells_per_cluster=10, data_norm_filt_loc=None):
    # impute_method: 'agg_wald' OR 'agg_wald_opt' OR 'no_imputation'
    # not maintained: 'IaconoClus' OR 'IaconoClus_median' OR 'IaconoClus_dim50' OR 'IaconoClus_dim50_median'
    print('\ncompute_repcells!\n')

    if impute_method == 'no_imputation':
        return None

    if data is None and data_norm_filt_loc is None:
        print('No data path and no data argument.')

    if data is None:
        print('Loading data')
        data = pd.read_csv(data_norm_filt_loc, index_col=0, sep='\t')

    t0 = time()
    if impute_method == 'agg_wald_opt' or impute_method == 'agg_wald':
        num_cells = data.shape[1]  # number of cells
        num_genes = data.shape[0]  # number of genes

        if impute_method == 'agg_wald_opt':
            n_clusters = find_opt_n_clusters(data, min_cells_per_cluster)
        elif impute_method == 'agg_wald':
            if num_cells >= 1000:
                n_clusters = min(100, int(np.round(num_cells / 30)))
            else:
                n_clusters = int(np.round(num_cells / 20))
                min_cells_per_cluster = 5
        print('aiming for', n_clusters, 'clusters')
        print('averagely', np.round(num_cells / n_clusters, 1), 'cells per cluster')

        clustering = AgglomerativeClustering(linkage='ward', n_clusters=n_clusters).fit(data.transpose())
        labels = renumber_clusters(list(clustering.labels_), min_cells_per_cluster)
        clusters = pd.DataFrame(index=range(data.shape[1]), columns=['x'], data=labels)

        clusters.to_csv(clustering_file_initial, sep=' ')

        # change the values for all cells in these clusters in the table to None
        all_clusters = set(labels)
        small_clusters = [c for c in all_clusters if labels.count(c) < min_cells_per_cluster]
        final_cluster_list = list(all_clusters - set(small_clusters))
        clusters.loc[[i for i in clusters.index if clusters.loc[i, 'x'] in small_clusters], 'x'] = None

        print('got', len(final_cluster_list), 'clusters')
        print('averagely', np.round(
            len([i for i in clusters.index if clusters.loc[i, 'x'] in final_cluster_list]) / len(final_cluster_list),
            2), 'cells per cluster')

        clusters.to_csv(clustering_file_final, sep=' ')

    '''    
    elif impute_method == 'IaconoClus' or impute_method == 'IaconoClus_dim50':        
        num_cells = data.shape[1]   #number of cells
        num_genes = data.shape[0]   #number of genes
                
        clusters = pd.read_csv(clustering_file_initial, index_col=0, sep=' ')
        
        clusters.index = range(len(clusters))
        labels = renumber_clusters(list(clusters['x']), min_cells_per_cluster)
        clusters.loc[:, 'x'] = labels
        clusters.to_csv(clustering_file, sep=' ')

        # change the values for all cells in these clusters in the table to None
        all_clusters = set(labels)
        small_clusters = [c for c in all_clusters if labels.count(c) < min_cells_per_cluster]
        final_cluster_list = list(all_clusters - set(small_clusters))
        clusters.loc[[i for i in clusters.index if clusters.loc[i, 'x'] in small_clusters], 'x'] = None

        print('got', len(final_cluster_list), 'clusters')
        print('averagely', np.round(
            len([i for i in clusters.index if clusters.loc[i, 'x'] in final_cluster_list]) / len(final_cluster_list),
            2), 'cells per cluster')
        
        clusters.to_csv(clustering_file_final, sep=' ')
    '''


####################################################################################################################
def create_imputed_file(clustering_file_final, repcells_data_loc,
                        data=None, impute_method='agg_wald', data_norm_filt_loc=None):
    print('create_imputed_file!')

    if data is None and data_norm_filt_loc is None:
        print('No data path and no data argument.')

    if data is None:
        print('Loading data')
        data = pd.read_csv(data_norm_filt_loc, index_col=0, sep='\t')

    if impute_method == 'no_imputation':
        repcells_data = data
    else:
        clusters = pd.read_csv(clustering_file_final, index_col=0, sep=' ').dropna()
        clusters.x = clusters.x.astype(int)

        final_cluster_list = set(clusters['x'])
        repcells_data = pd.DataFrame(data=np.zeros((len(data), len(final_cluster_list))),
                                     index=data.index, columns=np.array(range(len(final_cluster_list))))
        for c in final_cluster_list:
            print('cluster', c, '-', len([i for i in clusters.index if clusters.loc[i, 'x'] == c]), 'cells')
            repcells_data.loc[:, c] = data.iloc[:, [i for i in clusters.index if clusters.loc[i, 'x'] == c]].mean(
                axis=1)

        repcells_data = delete_genes_with_zero_reads(repcells_data)
        repcells_data = median_count_normalization(repcells_data)

    repcells_data.to_csv(repcells_data_loc, sep='\t')

    num_repcells = repcells_data.shape[1]  # number of cells
    num_genes = repcells_data.shape[0]  # number of genes

    return repcells_data, num_repcells, num_genes


####################################################################################################################
def compute_sets_on_standardized_genes(genetable_file,
                                       genes_in_sets_coo_file, genes_in_sets_csc_file, genes_in_sets_csr_file,
                                       genes_in_sets_row_file, genes_in_sets_col_file,
                                       counts=None, counts_file=None, num_stds_thresh=1):
    print('compute_sets_on_standardized_genes!')

    if counts is None and counts_file is None:
        print('No data path and no data argument.')

    if counts is None:
        print('Loading data')
        counts = pd.read_csv(counts_file, index_col=0, sep='\t')

    if not os.path.exists(genetable_file):
        save_genetable_to_file(data=counts, genetable_file=genetable_file)

    with open(genetable_file, 'rb') as fp:
        genes_table = pickle.load(fp)

    # normalize the original matrix
    counts = counts.subtract(counts.mean(axis=1), axis=0).divide(counts.std(axis=1), axis=0)

    # intialize the row and column arrays of the filled indices
    num_genes, num_samples = counts.shape[0], counts.shape[1]
    estimated_num_genes_per_set = num_genes
    genes_in_sets_row = np.empty((int(num_samples * num_samples / 2 * estimated_num_genes_per_set)))
    genes_in_sets_col = np.empty((int(num_samples * num_samples / 2 * estimated_num_genes_per_set)))

    start = time()
    step = 20
    iter_ind = 0
    # v3 = False
    for c1 in range(num_samples):
        for c2 in range(num_samples):
            if c1 != c2:
                gene_inds = list(genes_table.loc[
                                     counts.index[(counts.iloc[:, c2] - counts.iloc[:, c1]) > num_stds_thresh]
                                 ])
                genes_in_sets_row[iter_ind:(iter_ind + len(gene_inds))] = gene_inds
                genes_in_sets_col[iter_ind:(iter_ind + len(gene_inds))] = c1 * num_samples + c2
                iter_ind += len(gene_inds)

        '''
        if ((c1+1) % step == 0) or (c1 == num_cells-1):
            print('finished cell', c1+1)

            # save the row and col indices
            b = genes_in_sets_row.tolist()
            file_path = ('../objects/genes_in_sets/'+impute_method+'/num_stds_thresh_'+str(num_stds_thresh)+
                         '/genes_in_sets_rows_cols_ordered_pairs/'+cell_selection+'/rows/'+data_name+'_'+
                         str(num_genes)+'X'+str(num_samples)+'_'+filter_word+'_c1_'+str(c1+1)+'.json')
            json.dump(b, codecs.open(file_path, 'w', encoding='utf-8'), separators=(',', ':'), sort_keys=True, indent=4)

            b = genes_in_sets_col.tolist()
            file_path = ('../objects/genes_in_sets/'+impute_method+'/num_stds_thresh_'+str(num_stds_thresh)+
                         '/genes_in_sets_rows_cols_ordered_pairs/'+cell_selection+'/cols/'+data_name+'_'+
                         str(num_genes)+'X'+str(num_samples)+'_'+filter_word+'_c1_'+str(c1+1)+'.json')
            json.dump(b, codecs.open(file_path, 'w', encoding='utf-8'), separators=(',', ':'), sort_keys=True, indent=4)

            if 'last_saved_c1' in locals():
                try:
                    os.remove('../objects/genes_in_sets/'+impute_method+'/num_stds_thresh_'+str(num_stds_thresh)+
                              '/genes_in_sets_rows_cols_ordered_pairs/'+cell_selection+'/rows/'+data_name+'_'+
                              str(num_genes)+'X'+str(num_cells)+'_'+filter_word+'_c1_'+str(last_saved_c1+1)+'.json')
                    os.remove('../objects/genes_in_sets/'+impute_method+'/num_stds_thresh_'+str(num_stds_thresh)+
                              '/genes_in_sets_rows_cols_ordered_pairs/'+cell_selection+'/cols/'+data_name+'_'+
                              str(num_genes)+'X'+str(num_samples)+'_'+filter_word+'_c1_'+str(last_saved_c1+1)+'.json')
                except OSError:
                    pass

            last_saved_c1 = c1
        '''

    genes_in_sets_row = genes_in_sets_row[:iter_ind]
    genes_in_sets_col = genes_in_sets_col[:iter_ind]

    # save the row and col indices
    json.dump(genes_in_sets_row.tolist(), codecs.open(genes_in_sets_row_file, 'w', encoding='utf-8'),
              separators=(',', ':'), sort_keys=True, indent=4)
    json.dump(genes_in_sets_col.tolist(), codecs.open(genes_in_sets_col_file, 'w', encoding='utf-8'),
              separators=(',', ':'), sort_keys=True, indent=4)

    genes_in_sets = sparse.coo_matrix((np.ones_like(genes_in_sets_row), (genes_in_sets_row, genes_in_sets_col)),
                                      shape=(num_genes, num_samples * num_samples))

    # save coo_mat to file
    sparse.save_npz(genes_in_sets_coo_file, genes_in_sets)

    # save csc_mat to file
    sparse.save_npz(genes_in_sets_csc_file, genes_in_sets.tocsc())

    # save csr_mat to file
    sparse.save_npz(genes_in_sets_csr_file, genes_in_sets.tocsr())

    t = time() - start

    num_sets = num_samples * (num_samples - 1)

    avg_num_of_genes_per_set = np.absolute(genes_in_sets).sum() / num_sets


####################################################################################################################
def make_structure_list_unique(struct_list):
    struct_list_unique = []
    for s in struct_list:
        # check if exists in unique_list or not
        if s not in struct_list_unique:
            struct_list_unique.append(s)
    return struct_list_unique


####################################################################################################################
def find_structures(genes_in_sets_npz_file, structs_file, mu, num_iters=10000, path_len=3,
                    min_num_genes_in_struct=2, use_mp_to_find_structs=True,
                    # min_num_sets_in_struct=10,
                    ):
    print('find_structures!')

    genes_in_sets = sparse.load_npz(genes_in_sets_npz_file)

    # start with a random path
    if use_mp_to_find_structs:
        with mp.Pool(mp.cpu_count()) as pool:
            struct_list = pool.starmap(SPIRAL_mp_funcs.find_submat_around_random_path,
                                       [(mu, genes_in_sets, path_len, i) for i in range(num_iters)])
    else:
        struct_list = []
        for _ in range(num_iters):
            struct_list.append(SPIRAL_mp_funcs.find_submat_around_random_path(mu, genes_in_sets, path_len))

    relevant_structs = [s for s in struct_list if len(s[0]) >= min_num_genes_in_struct]
    print(len(relevant_structs), 'relevant structures')

    # Remove repeating structures
    relevant_structs = make_structure_list_unique(relevant_structs)
    print(len(relevant_structs), 'relevant structures without repetitions')

    write_struct_list_to_file(filename=structs_file, structlist=relevant_structs)

    return relevant_structs


####################################################################################################################
def evaluate_significance_of_structs(sigfile, genetable_file, impute_method='agg_wald', real_samp_name='sample',
                                     structs=None, structs_file=None, counts=None, counts_file=None,
                                     clustering_file_final=None):
    print('evaluate_significance_of_structs!')

    if structs is None and structs_file is None:
        print('No structure list and no structure list file argument.')

    if structs is None:
        print('Loading data')
        structs = read_struct_list_from_file(structs_file)

    with open(genetable_file, 'rb') as fp:
        genes_table = pickle.load(fp)

    if counts is None and counts_file is None:
        print('No data path and no data argument.')

    if counts is None:
        print('Loading data')
        counts = pd.read_csv(counts_file, index_col=0, sep='\t')

    num_samples, num_genes = counts.shape[1], counts.shape[0]

    normalized_counts = counts.subtract(counts.mean(axis=1), axis=0).divide(counts.std(axis=1), axis=0)
    normalized_counts.columns = [str(x) for x in normalized_counts.columns]

    if impute_method == 'no_imputation':
        samp_name = real_samp_name
        cols = ['num_genes', 'num_' + samp_name + 's',
                'num_genes_in_struct', 'num_' + samp_name + 's_in_struct',
                'num_' + samp_name + '_pairs', 'num_' + samp_name + '_pairs_in_struct',
                'log10_corrected_pval', 'structure_average_std', 'genes', samp_name + '_pairs']

    else:
        samp_name = 'repcell'
        cols = ['num_genes', 'num_' + real_samp_name + 's', 'num_' + samp_name + 's',
                'num_genes_in_struct', 'num_' + real_samp_name + 's_in_struct', 'num_' + samp_name + 's_in_struct',
                'num_' + samp_name + '_pairs', 'num_' + samp_name + '_pairs_in_struct',
                'log10_corrected_pval', 'structure_average_std', 'genes', samp_name + '_pairs']
        cells_to_repcells = pd.read_csv(clustering_file_final, index_col=0, sep=' ')

    significance_table = pd.DataFrame(np.zeros((len(structs), len(cols))), columns=cols)
    significance_table.loc[:, 'num_genes'] = num_genes
    significance_table.loc[:, 'num_' + samp_name + 's'] = num_samples
    significance_table.loc[:, 'num_' + samp_name + '_pairs'] = num_samples * (num_samples - 1)
    if impute_method != 'no_imputation':
        significance_table.loc[:, 'num_' + real_samp_name + 's'] = len(cells_to_repcells)

    for i, s in enumerate(structs):
        # print('structure', i, ':\t#rows:', len(s[0]), '\t#cols:', len(s[1]))

        genes = list(genes_table[s[0]].index)
        significance_table.loc[i, 'num_genes_in_struct'] = len(genes)
        significance_table.loc[i, 'genes'] = write_gene_list(genes)

        sample_pairs = [divmod(x, num_samples) for x in s[1]]
        significance_table.loc[i, 'num_' + samp_name + '_pairs_in_struct'] = len(sample_pairs)
        significance_table.loc[i, samp_name + '_pairs'] = str(sample_pairs)

        # find all the cells that are relevant to this structure:
        samples = list(set([item for sublist in sample_pairs for item in sublist]))
        significance_table.loc[i, 'num_' + samp_name + 's_in_struct'] = len(samples)

        if impute_method != 'no_imputation':
            significance_table.loc[i, 'num_' + real_samp_name + 's_in_struct'] = len(
                [i for i in cells_to_repcells.index if cells_to_repcells.loc[i, 'x'] in samples])

        curr_normalized_counts = normalized_counts.loc[genes, :]
        curr_normalized_counts = curr_normalized_counts.iloc[:, [c for c in samples]]

        # construct the transition matrix B that is relevant to this structure:
        curr_B = np.zeros((len(sample_pairs), len(samples)))
        for k, sett in enumerate(sample_pairs):
            curr_B[k, samples.index(sett[0])] = -1
            curr_B[k, samples.index(sett[1])] = 1

        # multiply the normalized matrix by B to get the vector to compare with
        curr_S = np.dot(curr_B, curr_normalized_counts.transpose())

        # covariance matrix of curr_Y (curr_S is a sample of it)
        cov_mat = np.dot(curr_B, curr_B.transpose())

        # compute structure_average_std
        avg_vec = (1 / len(sample_pairs)) * np.ones((len(sample_pairs)))
        average_diff_of_gene_std = np.sqrt(np.dot(np.dot(avg_vec.transpose(), cov_mat), avg_vec))
        significance_table.loc[i, 'structure_average_std'] = average_diff_of_gene_std

        # average the differences for each gene
        average_diff_of_genes = np.dot(avg_vec.transpose(), curr_S)
        if len(average_diff_of_genes) != len(genes):
            print("PROBLEM 1")

        genes_probs_logs = norm.logsf(average_diff_of_genes / average_diff_of_gene_std) / np.log(10)

        log10pval = np.sum(genes_probs_logs)
        # print('initial probability: 1e', log10pval)

        # Multiple testing correction
        multestcorr = correct_for_multiple_testing(num_genes, len(genes), num_samples * (num_samples - 1),
                                                   len(sample_pairs))
        # print('Multiple testing correction: 1e', multestcorr)

        significance_table.loc[i, 'log10_corrected_pval'] = log10pval + multestcorr
        # print('final probability: 1e', log10pval+multestcorr)

    significance_table = significance_table.sort_values(by='log10_corrected_pval')

    significance_table.to_excel(sigfile)
    return significance_table


####################################################################################################################
def merge_sigtables(sigfile_merged, data_path, impute_method, num_stds_thresh_lst, mu_lst, path_len_lst, num_iters_lst):
    print('merge_sigtables!')

    sigtable_merged = pd.DataFrame()
    for num_stds_thresh in num_stds_thresh_lst:
        for mu in mu_lst:
            for path_len in path_len_lst:
                for num_iters in num_iters_lst:
                    sigfile = sig_filename(data_path=data_path, impute_method=impute_method,
                                           num_stds_thresh=num_stds_thresh, mu=mu,
                                           path_len=path_len, num_iters=num_iters)
                    sigtable = pd.read_excel(sigfile, index_col=0, engine='openpyxl')
                    if not sigtable.empty:
                        sigtable.loc[:, 'old_struct_num'] = sigtable.index
                        sigtable.loc[:, 'num_stds_thresh'] = num_stds_thresh
                        sigtable.loc[:, 'mu'] = mu
                        sigtable.loc[:, 'path_len'] = path_len
                        sigtable.loc[:, 'num_iters'] = num_iters
                        sigtable_merged = pd.concat([sigtable_merged, sigtable], ignore_index=True)
    sigtable_merged.to_excel(sigfile_merged)
    return sigtable_merged


####################################################################################################################
def create_dict_of_struct_lsts(data_path, impute_method, sigtable):
    # create a dictionary of structure lists
    num_std_lst = list(set(sigtable.loc[:, 'num_stds_thresh']))
    mu_lst = list(set(sigtable.loc[:, 'mu']))
    path_len_lst = list(set(sigtable.loc[:, 'path_len']))
    num_iters_lst = list(set(sigtable.loc[:, 'num_iters']))
    struct_dict = dict()
    for num_std in num_std_lst:
        for mu in mu_lst:
            for path_len in path_len_lst:
                for num_iters in num_iters_lst:
                    structs_file = structs_filename(data_path=data_path, impute_method=impute_method,
                                                    num_stds_thresh=num_std,
                                                    mu=mu,
                                                    path_len=path_len,
                                                    num_iters=num_iters)
                    struct_dict[str(num_std) + '_' + str(mu) + '_' + str(path_len) + '_' + str(
                        num_iters)] = read_struct_list_from_file(structs_file)
    return struct_dict


####################################################################################################################
def filter_similar_structures(sigfile, significance_table, sigfile_filt, genetable_file, data_path,
                              struct_thr, min_nstructs, max_nstructs,
                              impute_method='agg_wald', real_samp_name='sample',
                              log_corrpval_thr=-15, Jaccard_thr_genes=0.75, Jaccard_thr_sample_pairs=0.5,
                              # lower_thr=0.4, Jaccard_thr_high_low_repcells=0.5
                              ):
    print('\nfilter_similar_structures!\n')

    if significance_table is None and sigfile is None:
        print('No significance_table path and no significance_table argument.')

    if significance_table is None:
        print('Loading data')
        significance_table = pd.read_excel(sigfile, index_col=0, engine='openpyxl')

    if impute_method == 'no_imputation':
        samp_name = real_samp_name
    else:
        samp_name = 'repcell'

    print('Starting with', len(significance_table), 'structures')

    # save Jaccard_thr_genes to file
    with open(os.path.join(data_path, 'Jaccard_thr_genes.txt'), 'w') as text_file:
        text_file.write(str(Jaccard_thr_genes))

    # save Jaccard_thr_sample_pairs to file
    # with open(os.path.join(data_path, 'Jaccard_thr_sample_pairs.txt'), 'w') as text_file:
    #    text_file.write(str(Jaccard_thr_sample_pairs))

    # filter out structures with log_corrpval>log_corrpval_thr
    significance_table_filt = significance_table[significance_table['log10_corrected_pval'] <= log_corrpval_thr]
    print('There are', len(significance_table_filt), 'structures with log10_corrected_pval<', log_corrpval_thr)

    if len(significance_table_filt) == 0:
        print('NO STRUCTURES FOUND!')
        return None

    #### if there are any structures with log10_corrected_pval<log_corrpval_thr: #####
    # Load genes_table
    with open(genetable_file, 'rb') as fp:
        genes_table = pickle.load(fp)

    struct_dict = create_dict_of_struct_lsts(data_path=data_path, impute_method=impute_method,
                                             sigtable=significance_table_filt)

    # An efficient way to find pairs of similar structures and remove the least favorable structure
    # It counts on the fact that the table is sorted by 'structure_average_std'
    significance_table_filt = significance_table_filt.sort_values(by='structure_average_std', ascending=True)
    inds = list(significance_table_filt.index)
    final_struct_lst = []
    while inds:
        i = inds.pop(0)
        final_struct_lst.append(i)
        similar_to_i = []

        # read i'th gene_list
        if len(significance_table_filt.loc[
                   i, 'genes']) == 32767:  # The genelist was too long to be written in one excel cell
            structs = struct_dict[str(significance_table_filt.loc[i, 'num_stds_thresh']) + '_' +
                                  str(significance_table_filt.loc[i, 'mu']) + '_' +
                                  str(significance_table_filt.loc[i, 'path_len']) + '_' +
                                  str(significance_table_filt.loc[i, 'num_iters'])]
            genes_i = list(genes_table[structs[significance_table_filt.loc[i, 'old_struct_num']][0]].index)
        else:  # read the gene list from the excel file
            genes_i = read_gene_list(significance_table_filt.loc[i, 'genes'])

        # read i'th repcell-pair sets
        sets_in_struct_i = significance_table_filt.loc[i, samp_name + '_pairs']
        sets_in_struct_i = [tuple(list(map(int, s.strip(')').strip('(').split(', ')))) for s in
                            sets_in_struct_i.strip(']').strip('[').split('), (')]

        for j in inds:
            # read j'th gene_list
            if len(significance_table_filt.loc[
                       j, 'genes']) == 32767:  # The genelist was too long to be written in one excel cell
                structs = struct_dict[str(significance_table_filt.loc[j, 'num_stds_thresh']) + '_' +
                                      str(significance_table_filt.loc[j, 'mu']) + '_' +
                                      str(significance_table_filt.loc[j, 'path_len']) + '_' +
                                      str(significance_table_filt.loc[j, 'num_iters'])]
                genes_j = list(genes_table[structs[significance_table_filt.loc[j, 'old_struct_num']][0]].index)
            else:  # read the gene list from the excel file
                genes_j = read_gene_list(significance_table_filt.loc[j, 'genes'])

            # consider (i,j) similar if the Jaccard index of gene lists is above overlap_thresh
            jaccard_i_j_genes = len(set(genes_i).intersection(set(genes_j))) / len(set(genes_i).union(set(genes_j)))
            if jaccard_i_j_genes >= Jaccard_thr_genes:
                similar_to_i.append(j)

            else:  # if not similar by gene sets, let's check similarity by repcell-pair sets
                # read j'th repcell-pair sets
                sets_in_struct_j = significance_table_filt.loc[j, samp_name + '_pairs']
                sets_in_struct_j = [tuple(list(map(int, s.strip(')').strip('(').split(', ')))) for s in
                                    sets_in_struct_j.strip(']').strip('[').split('), (')]

                # consider (i,j) similar if the Jaccard index of sets is above overlap_thresh
                jaccard_i_j_sets = len(set(sets_in_struct_i).intersection(set(sets_in_struct_j))) / len(
                    set(sets_in_struct_i).union(set(sets_in_struct_j)))
                if jaccard_i_j_sets >= Jaccard_thr_sample_pairs:
                    similar_to_i.append(j)
                '''
                else: #if not similar by gene sets or sample pairs, let's check similarity by low and high repcells +
                    # gene set similarity above the lower_thr
                    high_repcells_in_struct_i = set([s[1] for s in sets_in_struct_i])
                    high_repcells_in_struct_j = set([s[1] for s in sets_in_struct_j])
                    jaccard_i_j_high_repcells = len(
                        high_repcells_in_struct_i.intersection(high_repcells_in_struct_j)) / len(
                        high_repcells_in_struct_i.union(high_repcells_in_struct_j))

                    low_repcells_in_struct_i = set([s[0] for s in sets_in_struct_i])
                    low_repcells_in_struct_j = set([s[0] for s in sets_in_struct_j])
                    jaccard_i_j_low_repcells = len(
                        low_repcells_in_struct_i.intersection(low_repcells_in_struct_j)) / len(
                        low_repcells_in_struct_i.union(low_repcells_in_struct_j))
                    
                    if (jaccard_i_j_high_repcells >= Jaccard_thr_high_low_repcells) and (
                            jaccard_i_j_low_repcells >= Jaccard_thr_high_low_repcells) and (
                            jaccard_i_j_genes >= lower_thr):
                        similar_to_i.append(j)
                '''

        inds = [j for j in inds if j not in similar_to_i]

    significance_table_filt = significance_table_filt.loc[final_struct_lst, :]
    print('There are', len(significance_table_filt), 'non-overlapping structures')

    significance_table_filt = significance_table_filt.sort_values(by='structure_average_std', ascending=True)

    # decide on the final structure list
    final_struct_lst = compute_final_struct_list(sigtable=significance_table_filt, struct_thr=struct_thr,
                                                 min_nstructs=min_nstructs, max_nstructs=max_nstructs)
    significance_table_filt = significance_table_filt.loc[final_struct_lst, :]

    print('There are', len(significance_table_filt), 'structures in the final structure list')

    # renumber structures
    significance_table_filt.index = np.arange(1, 1 + len(significance_table_filt))

    # save the matrix of all pair-wise Jaccard indices between structures
    save_Jaccard_mat_genes(sigtable=significance_table_filt, struct_dict=struct_dict, genes_table=genes_table,
                           data_path=data_path)

    significance_table_filt.to_excel(sigfile_filt)
    return significance_table_filt


####################################################################################################################
def compute_final_struct_list(sigtable, struct_thr, min_nstructs, max_nstructs):
    # decide on structures to be visualized and GO termed
    # S is the size of s={all structures}.
    # M is the size of m={structures with structure_average_std<struct_thr}
    # if min_nstructs<=M<=max_nstructs, return m.
    # if M>max_nstructs, return the first max_nstructs.
    # if M<min_nstructs, then: if min_nstructs<S, return the first min_nstructs. Else, return s.
    structs_lst = list(sigtable.index[sigtable['structure_average_std'] <= struct_thr])
    if len(structs_lst) > max_nstructs:
        structs_lst = list(sigtable.sort_values(by='structure_average_std', ascending=True).index[:max_nstructs])
    elif len(structs_lst) < min_nstructs:
        structs_lst = list(
            sigtable.sort_values(by='structure_average_std', ascending=True).index[:min(min_nstructs, len(sigtable))])
    return structs_lst


####################################################################################################################
def save_Jaccard_mat_genes(sigtable, struct_dict, genes_table, data_path):
    inds = sigtable.index
    Jaccard_mat_genes = np.zeros((len(inds), len(inds)))
    for ii, i in enumerate(inds):
        # read i'th gene_list
        if len(sigtable.loc[
                   i, 'genes']) == 32767:  # The genelist was too long to be written in one excel cell
            structs = struct_dict[str(sigtable.loc[i, 'num_stds_thresh']) + '_' +
                                  str(sigtable.loc[i, 'mu']) + '_' +
                                  str(sigtable.loc[i, 'path_len']) + '_' +
                                  str(sigtable.loc[i, 'num_iters'])]
            genes_i = list(genes_table[structs[sigtable.loc[i, 'old_struct_num']][0]].index)
        else:  # read the gene list from the excel file
            genes_i = read_gene_list(sigtable.loc[i, 'genes'])
        for j in inds[ii + 1:]:
            # read j'th gene_list
            if len(sigtable.loc[
                       j, 'genes']) == 32767:  # The genelist was too long to be written in one excel cell
                structs = struct_dict[str(sigtable.loc[j, 'num_stds_thresh']) + '_' +
                                      str(sigtable.loc[j, 'mu']) + '_' +
                                      str(sigtable.loc[j, 'path_len']) + '_' +
                                      str(sigtable.loc[j, 'num_iters'])]
                genes_j = list(genes_table[structs[sigtable.loc[j, 'old_struct_num']][0]].index)
            else:  # read the gene list from the excel file
                genes_j = read_gene_list(sigtable.loc[j, 'genes'])
            jaccard_i_j_genes = len(set(genes_i).intersection(set(genes_j))) / len(set(genes_i).union(set(genes_j)))
            Jaccard_mat_genes[i - 1, j - 1] = jaccard_i_j_genes
            Jaccard_mat_genes[j - 1, i - 1] = jaccard_i_j_genes

    np.save(os.path.join(data_path, 'Jaccard_mat_genes.npy'), Jaccard_mat_genes)


####################################################################################################################


def zip_all_files_in_folder(zipfile, folder):
    print('zip_all_files_in_folder!')

    with ZipFile(zipfile, 'w') as zipObj:
        # Iterate over all the files in directory
        for folderName, subfolders, filenames in os.walk(folder):
            for filename in filenames:
                # create complete filepath of file in directory
                filePath = os.path.join(folderName, filename)
                # Add file to zip
                zipObj.write(filePath, os.path.basename(filePath))
