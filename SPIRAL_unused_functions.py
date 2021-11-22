def choose_cells_to_work_on(data_name, scRNA_data, num_cells, num_genes, cell_selection='use_all'):
    # cell_selection: use_non_GCS_gene_counts or use_first_cells_on_the_list or use_all
    print(data_name)
    print('\nchoose_cells_to_work_on!\n')
    if cell_selection == 'use_non_GCS_gene_counts':
        print('use_non_GCS_gene_counts')
        gene_counts = scRNA_data.astype(bool).sum()

        num_top_genes = 200
        corr_with_gene_counts = pd.Series(index=scRNA_data.index, data=np.empty((len(scRNA_data))))

        for g in scRNA_data.index:
            r, _ = pearsonr(scRNA_data.loc[g, :], gene_counts)
            corr_with_gene_counts.loc[g] = r

        corr_with_gene_counts = corr_with_gene_counts.sort_values(ascending=True)
        genes_for_dropout_computation = corr_with_gene_counts[:-num_top_genes].index.to_list()

        num_cells = 1000  # number of cells

        chosen_cells = scRNA_data.loc[genes_for_dropout_computation, :].astype(bool).sum().sort_values(ascending=False)[
                       :num_cells].index.to_list()
        scRNA_data = scRNA_data.loc[:, chosen_cells]
        print(num_genes, 'X', num_cells, 'expression matrix')

    elif cell_selection == 'use_first_cells_on_the_list':
        print('use_first_cells_on_the_list')
        # Try on a smaller matrix first
        num_cells = 20  # number of cells
        scRNA_data = scRNA_data.iloc[:num_genes, :num_cells]
        print(num_genes, 'X', num_cells, 'expression matrix')

    elif cell_selection == 'use_all':
        print('use all cells')
        print(num_genes, 'X', num_cells, 'expression matrix')
    else:
        print('Did not understand cell selection.')

    return scRNA_data, num_cells, num_genes

####################################################################################################################
def compute_sets(data_name, num_cells, num_genes, cell_selection, impute_method,
                 ordered_pairs=True, save_genes_table_to_file=True,
                 num_stds_thresh=1, read_scRNA_data_from_file=False, scRNA_data=np.nan,
                 filter_word='filtered'):
    # if ordered_pairs=True then S(i,j)!=S(j,i)

    if read_scRNA_data_from_file:
        ind = [m.start() for m in re.finditer('/', data_dict[data_name])][-1]
        imputed_file = data_dict[data_name][:ind + 1] + impute_method + data_dict[data_name][
                                                                        ind:-4] + '_' + filter_word + '_imputed.csv'
        scRNA_data = load_data(imputed_file)
    elif scRNA_data == np.nan:
        print('ERROR: no scRNA_data matrix')

    print(data_name)
    print('\ncompute_sets!\n')

    # table with genes' names and indices
    genes_table = pd.Series(index=scRNA_data.index, data=np.array(range(num_genes)))
    if save_genes_table_to_file:
        genes_table_file = '../objects/genes_tables/' + data_name + '_' + str(num_genes) + 'X' + str(
            num_cells) + '_' + filter_word + '.p'
        validate_folder(genes_table_file)
        with open(genes_table_file, 'wb') as fp:
            pickle.dump(genes_table, fp, protocol=pickle.HIGHEST_PROTOCOL)

    # compute std for every gene
    gene_std = scRNA_data.std(axis=1)

    # intialize the row and column arrays of the filled indices
    estimated_num_genes_per_set = num_genes / 2
    genes_in_sets_row = np.empty((int(num_cells * num_cells / 2 * estimated_num_genes_per_set)))
    genes_in_sets_col = np.empty((int(num_cells * num_cells / 2 * estimated_num_genes_per_set)))
    if not ordered_pairs:
        genes_up_or_down = np.empty((int(num_cells * num_cells / 2 * estimated_num_genes_per_set)))

    start = time()
    step = 20
    iter_ind = 0
    v1, v2, v3 = False, False, False
    for c1 in range(num_cells):
        # print(c1)
        if not ordered_pairs:
            for c2 in np.arange(c1 + 1, num_cells):
                gene_inds = list(genes_table.loc[scRNA_data.index[
                    (scRNA_data.iloc[:, c1] - scRNA_data.iloc[:, c2]) > num_stds_thresh * gene_std]])
                genes_in_sets_row[iter_ind:(iter_ind + len(gene_inds))] = gene_inds
                genes_in_sets_col[iter_ind:(iter_ind + len(gene_inds))] = c1 * num_cells + c2
                genes_up_or_down[iter_ind:(iter_ind + len(gene_inds))] = -1
                iter_ind += len(gene_inds)

                gene_inds = list(genes_table.loc[scRNA_data.index[
                    (scRNA_data.iloc[:, c2] - scRNA_data.iloc[:, c1]) > num_stds_thresh * gene_std]])
                genes_in_sets_row[iter_ind:(iter_ind + len(gene_inds))] = gene_inds
                genes_in_sets_col[iter_ind:(iter_ind + len(gene_inds))] = c1 * num_cells + c2
                genes_up_or_down[iter_ind:(iter_ind + len(gene_inds))] = 1
                iter_ind += len(gene_inds)
        else:
            for c2 in range(num_cells):
                if c1 != c2:
                    gene_inds = list(genes_table.loc[scRNA_data.index[
                        (scRNA_data.iloc[:, c2] - scRNA_data.iloc[:, c1]) > num_stds_thresh * gene_std]])
                    genes_in_sets_row[iter_ind:(iter_ind + len(gene_inds))] = gene_inds
                    genes_in_sets_col[iter_ind:(iter_ind + len(gene_inds))] = c1 * num_cells + c2
                    iter_ind += len(gene_inds)

        if ((c1 + 1) % step == 0) or (c1 == num_cells - 1):
            print('finished cell', c1 + 1)

            # save the row and col indices, and the genes_up_or_down vector to file
            if not ordered_pairs:
                b = genes_in_sets_row.tolist()
                file_path = ('../objects/genes_in_sets/' + impute_method + '/num_stds_thresh_' + str(num_stds_thresh) +
                             '/genes_in_sets_rows_cols_unordered_pairs/' + cell_selection + '/rows/' + data_name + '_' +
                             str(num_genes) + 'X' + str(num_cells) + '_' + filter_word + '_c1_' + str(c1 + 1) + '.json')
                if not v1:
                    validate_folder(file_path)
                    v1 = True
                json.dump(b, codecs.open(file_path, 'w', encoding='utf-8'), separators=(',', ':'), sort_keys=True,
                          indent=4)

                b = genes_in_sets_col.tolist()
                file_path = ('../objects/genes_in_sets/' + impute_method + '/num_stds_thresh_' + str(num_stds_thresh) +
                             '/genes_in_sets_rows_cols_unordered_pairs/' + cell_selection + '/cols/' + data_name + '_' +
                             str(num_genes) + 'X' + str(num_cells) + '_' + filter_word + '_c1_' + str(c1 + 1) + '.json')
                if not v2:
                    validate_folder(file_path)
                    v2 = True
                json.dump(b, codecs.open(file_path, 'w', encoding='utf-8'), separators=(',', ':'), sort_keys=True,
                          indent=4)

                b = genes_up_or_down.tolist()
                file_path = ('../objects/genes_in_sets/' + impute_method + '/num_stds_thresh_' + str(num_stds_thresh) +
                             '/genes_in_sets_rows_cols_unordered_pairs/' + cell_selection + '/genes_up_or_down/' + data_name +
                             '_' + str(num_genes) + 'X' + str(num_cells) + '_' + filter_word + '_c1_' + str(
                            c1 + 1) + '.json')
                if not v3:
                    validate_folder(file_path)
                    v3 = True
                json.dump(b, codecs.open(file_path, 'w', encoding='utf-8'), separators=(',', ':'), sort_keys=True,
                          indent=4)

                try:
                    if (c1 == num_cells - 1) and (np.mod(c1 + 1, step) != 0):
                        last_c1 = int(c1 + 1 / step) * step - 1
                    else:
                        last_c1 = c1 - step
                    os.remove('../objects/genes_in_sets/' + impute_method + '/num_stds_thresh_' + str(num_stds_thresh) +
                              '/genes_in_sets_rows_cols_unordered_pairs/' + cell_selection + '/rows/' + data_name + '_' +
                              str(num_genes) + 'X' + str(num_cells) + '_' + filter_word + '_c1_' + str(
                        last_c1 + 1) + '.json')
                    os.remove('../objects/genes_in_sets/' + impute_method + '/num_stds_thresh_' + str(num_stds_thresh) +
                              '/genes_in_sets_rows_cols_unordered_pairs/' + cell_selection + '/cols/' + data_name + '_' +
                              str(num_genes) + 'X' + str(num_cells) + '_' + filter_word + '_c1_' + str(
                        last_c1 + 1) + '.json')
                    os.remove('../objects/genes_in_sets/' + impute_method + '/num_stds_thresh_' + str(num_stds_thresh) +
                              '/genes_in_sets_rows_cols_unordered_pairs/' + cell_selection + '/genes_up_or_down/' + data_name +
                              '_' + str(num_genes) + 'X' + str(num_cells) + '_' + filter_word + '_c1_' + str(
                        last_c1 + 1) + '.json')
                except OSError:
                    pass
            else:
                b = genes_in_sets_row.tolist()
                file_path = ('../objects/genes_in_sets/' + impute_method + '/num_stds_thresh_' + str(num_stds_thresh) +
                             '/genes_in_sets_rows_cols_ordered_pairs/' + cell_selection + '/rows/' + data_name + '_' +
                             str(num_genes) + 'X' + str(num_cells) + '_' + filter_word + '_c1_' + str(c1 + 1) + '.json')
                if not v1:
                    validate_folder(file_path)
                    v1 = True
                json.dump(b, codecs.open(file_path, 'w', encoding='utf-8'), separators=(',', ':'), sort_keys=True,
                          indent=4)

                b = genes_in_sets_col.tolist()
                file_path = ('../objects/genes_in_sets/' + impute_method + '/num_stds_thresh_' + str(num_stds_thresh) +
                             '/genes_in_sets_rows_cols_ordered_pairs/' + cell_selection + '/cols/' + data_name + '_' +
                             str(num_genes) + 'X' + str(num_cells) + '_' + filter_word + '_c1_' + str(c1 + 1) + '.json')
                if not v2:
                    validate_folder(file_path)
                    v2 = True
                json.dump(b, codecs.open(file_path, 'w', encoding='utf-8'), separators=(',', ':'), sort_keys=True,
                          indent=4)

                try:
                    if (c1 == num_cells - 1) and (np.mod(c1 + 1, step) != 0):
                        last_c1 = int(c1 + 1 / step) * step - 1
                    else:
                        last_c1 = c1 - step
                    os.remove('../objects/genes_in_sets/' + impute_method + '/num_stds_thresh_' + str(num_stds_thresh) +
                              '/genes_in_sets_rows_cols_ordered_pairs/' + cell_selection + '/rows/' + data_name + '_' +
                              str(num_genes) + 'X' + str(num_cells) + '_' + filter_word + '_c1_' + str(
                        last_c1 + 1) + '.json')
                    os.remove('../objects/genes_in_sets/' + impute_method + '/num_stds_thresh_' + str(num_stds_thresh) +
                              '/genes_in_sets_rows_cols_ordered_pairs/' + cell_selection + '/cols/' + data_name + '_' +
                              str(num_genes) + 'X' + str(num_cells) + '_' + filter_word + '_c1_' + str(
                        last_c1 + 1) + '.json')
                except OSError:
                    pass

    if not ordered_pairs:
        file_path = ('../objects/genes_in_sets/' + impute_method + '/num_stds_thresh_' + str(num_stds_thresh) +
                     '/genes_in_sets_rows_cols_unordered_pairs/' + cell_selection + '/rows/' + data_name + '_' +
                     str(num_genes) + 'X' + str(num_cells) + '_' + filter_word + '_c1_' + str(num_cells) + '.json')
        new_file_path = ('../objects/genes_in_sets/' + impute_method + '/num_stds_thresh_' + str(num_stds_thresh) +
                         '/genes_in_sets_rows_cols_unordered_pairs/' + cell_selection + '/rows/' + data_name + '_' +
                         str(num_genes) + 'X' + str(num_cells) + '_' + filter_word + '.json')
        os.rename(file_path, new_file_path)

        file_path = ('../objects/genes_in_sets/' + impute_method + '/num_stds_thresh_' + str(num_stds_thresh) +
                     '/genes_in_sets_rows_cols_unordered_pairs/' + cell_selection + '/cols/' + data_name + '_' +
                     str(num_genes) + 'X' + str(num_cells) + '_' + filter_word + '_c1_' + str(num_cells) + '.json')
        new_file_path = ('../objects/genes_in_sets/' + impute_method + '/num_stds_thresh_' + str(num_stds_thresh) +
                         '/genes_in_sets_rows_cols_unordered_pairs/' + cell_selection + '/cols/' + data_name + '_' +
                         str(num_genes) + 'X' + str(num_cells) + '_' + filter_word + '.json')
        os.rename(file_path, new_file_path)

        file_path = ('../objects/genes_in_sets/' + impute_method + '/num_stds_thresh_' + str(num_stds_thresh) +
                     '/genes_in_sets_rows_cols_unordered_pairs/' + cell_selection + '/genes_up_or_down/' + data_name +
                     '_' + str(num_genes) + 'X' + str(num_cells) + '_' + filter_word + '_c1_' + str(
                    num_cells) + '.json')
        new_file_path = ('../objects/genes_in_sets/' + impute_method + '/num_stds_thresh_' + str(num_stds_thresh) +
                         '/genes_in_sets_rows_cols_unordered_pairs/' + cell_selection + '/genes_up_or_down/' + data_name +
                         '_' + str(num_genes) + 'X' + str(num_cells) + '_' + filter_word + '.json')
        os.rename(file_path, new_file_path)

    else:
        file_path = ('../objects/genes_in_sets/' + impute_method + '/num_stds_thresh_' + str(num_stds_thresh) +
                     '/genes_in_sets_rows_cols_ordered_pairs/' + cell_selection + '/rows/' + data_name + '_' +
                     str(num_genes) + 'X' + str(num_cells) + '_' + filter_word + '_c1_' + str(num_cells) + '.json')
        new_file_path = ('../objects/genes_in_sets/' + impute_method + '/num_stds_thresh_' + str(num_stds_thresh) +
                         '/genes_in_sets_rows_cols_ordered_pairs/' + cell_selection + '/rows/' + data_name + '_' +
                         str(num_genes) + 'X' + str(num_cells) + '_' + filter_word + '.json')
        os.rename(file_path, new_file_path)

        file_path = ('../objects/genes_in_sets/' + impute_method + '/num_stds_thresh_' + str(num_stds_thresh) +
                     '/genes_in_sets_rows_cols_ordered_pairs/' + cell_selection + '/cols/' + data_name + '_' +
                     str(num_genes) + 'X' + str(num_cells) + '_' + filter_word + '_c1_' + str(num_cells) + '.json')
        new_file_path = ('../objects/genes_in_sets/' + impute_method + '/num_stds_thresh_' + str(num_stds_thresh) +
                         '/genes_in_sets_rows_cols_ordered_pairs/' + cell_selection + '/cols/' + data_name + '_' +
                         str(num_genes) + 'X' + str(num_cells) + '_' + filter_word + '.json')
        os.rename(file_path, new_file_path)

    genes_in_sets_row = genes_in_sets_row[:iter_ind]
    genes_in_sets_col = genes_in_sets_col[:iter_ind]
    if not ordered_pairs:
        genes_up_or_down = genes_up_or_down[:iter_ind]
        genes_in_sets = sparse.coo_matrix((genes_up_or_down, (genes_in_sets_row, genes_in_sets_col)),
                                          shape=(num_genes, num_cells * num_cells))
    else:
        genes_in_sets = sparse.coo_matrix((np.ones_like(genes_in_sets_row), (genes_in_sets_row, genes_in_sets_col)),
                                          shape=(num_genes, num_cells * num_cells))

    # save coo_mat to file
    if not ordered_pairs:
        npz_file = ('../objects/genes_in_sets/' + impute_method + '/num_stds_thresh_' + str(num_stds_thresh) +
                    '/genes_in_sets_unordered_pairs/' + cell_selection + '/coo_mat/' + data_name + '_' +
                    str(num_genes) + 'X' + str(num_cells) + '_' + filter_word + '.npz')
        validate_folder(npz_file)
        sparse.save_npz(npz_file, genes_in_sets)
    else:
        npz_file = ('../objects/genes_in_sets/' + impute_method + '/num_stds_thresh_' + str(num_stds_thresh) +
                    '/genes_in_sets_ordered_pairs/' + cell_selection + '/coo_mat/' + data_name + '_' +
                    str(num_genes) + 'X' + str(num_cells) + '_' + filter_word + '.npz')
        validate_folder(npz_file)
        sparse.save_npz(npz_file, genes_in_sets)

    t = time() - start

    if not ordered_pairs:
        num_sets = num_cells * (num_cells - 1) / 2
    else:
        num_sets = num_cells * (num_cells - 1)

    print('it takes', np.round(t / num_sets, 2), 'seconds for a set')
    print('it takes', np.round(t, 2), 'seconds for', num_sets, 'sets')
    print('it would take', np.round((t / 60 / num_sets) * (num_cells * (num_cells - 1)), 2),
          'minutes (=', np.round(t / 60 / 60 / num_sets * num_cells * (num_cells - 1), 2),
          'hours =', np.round(t / 60 / 60 / 24 / num_sets * num_cells * (num_cells - 1), 2),
          'days) for an expression matrix with', num_cells,
          'cells (', num_cells * (num_cells - 1), 'sets)')

    avg_num_of_genes_per_set = np.absolute(genes_in_sets).sum() / num_sets
    print('avg_num_of_genes_per_set is around', np.round(avg_num_of_genes_per_set, 2))

    save_as_csc_mat(data_name, num_cells, num_genes, genes_in_sets, cell_selection, ordered_pairs, filter_word,
                    num_stds_thresh)
    save_as_csr_mat(data_name, num_cells, num_genes, genes_in_sets, cell_selection, ordered_pairs, filter_word,
                    num_stds_thresh)

##############################################################################################################
def save_filtered_normed_to_file(data_loc, new_file, median_count_normalization_flag=False, with_log=False,
                                 filter_out_low_counts_genes=False,
                                 filter_out_low_counts_cells=False,
                                 min_counts_for_gene=20,
                                 min_counts_for_cell=20,
                                 filter_out_low_nFeatures_cells=False,
                                 filter_out_high_nFeatures_cells=False,
                                 max_nFeatures_for_cell=10000,
                                 min_nFeatures_for_cell=300,
                                 filter_mt_percent=False,
                                 mt_thresh=5):
    orig_scRNA_data = load_data(data_loc=data_loc, median_count_normalization_flag=median_count_normalization_flag,
                                with_log=with_log,
                                filter_out_low_counts_genes=filter_out_low_counts_genes,
                                filter_out_low_counts_cells=filter_out_low_counts_cells,
                                min_counts_for_gene=min_counts_for_gene,
                                min_counts_for_cell=min_counts_for_cell,
                                filter_out_low_nFeatures_cells=filter_out_low_nFeatures_cells,
                                filter_out_high_nFeatures_cells=filter_out_high_nFeatures_cells,
                                max_nFeatures_for_cell=max_nFeatures_for_cell,
                                min_nFeatures_for_cell=min_nFeatures_for_cell,
                                filter_mt_percent=filter_mt_percent,
                                mt_thresh=mt_thresh)
    print(orig_scRNA_data[:5])

    ind = [m.start() for m in re.finditer('/', new_file)][-1]
    if not os.path.exists(new_file[:ind]):
        os.makedirs(new_file[:ind])

    orig_scRNA_data.to_csv(new_file, sep='\t')

##################   unused enrichment funcs           ####################################################
def add_ranked_GO_terms(data_name, cell_selection, ordered_pairs, impute_method,
                        read_sigtable_from_file=False, sigtable=np.nan,
                        read_scRNA_data_from_file=False, scRNA_data=np.nan,
                        start_from_scratch=False, pvals=[0.000001, 0.0000001],
                        num_stds_thresh=1, mu=0.9, filt_mode='filt_by_3',
                        filter_word='filtered', struct_thr=0.6, log10_pavg_of_genes_thr=-10,
                        no_iterations=False, random_cellpair_set=False,
                        random_path=False, path_len=None, num_iters=10000):
    # filt_mode = 'filt_by_3' OR 'filt_by_structure_average_std' OR 'filt_by_rule'

    print(data_name)
    print('\nadd_GO_terms!\n')

    gorilla_url = 'http://cbl-gorilla.cs.technion.ac.il/'

    sigfile_filt = ('../results'
                    + random_cellpair_set * ('_random_cellpair_set' + '_num_iters_' + str(num_iters))
                    + random_path * ('_random_path_' + str(path_len) + '_num_iters_' + str(num_iters))
                    + no_iterations * '_no_iterations'
                    + '/sigtables/num_stds_thresh_' + str(num_stds_thresh) + '/mu_' + str(mu) + '/' +
                    (not ordered_pairs) * 'un' + 'ordered_pairs/'
                    + cell_selection + '/' + impute_method + '/' + data_name + '_' + filter_word + '_significance_table_filt.xlsx')
    sigfile_ranked_GO = sigfile_filt[:-5] + '_ranked_GO' + sigfile_filt[-5:]
    sigfile_ranked_GO_temp = sigfile_filt[:-5] + '_ranked_GO_temp' + sigfile_filt[-5:]
    print('\n\n\n', sigfile_filt)
    print('\n\n\n', sigfile_ranked_GO)
    print('\n\n\n', sigfile_ranked_GO_temp)

    if read_sigtable_from_file:
        if start_from_scratch or (not os.path.exists(sigfile_ranked_GO_temp)):
            sigtable = pd.read_excel(sigfile_filt, index_col=0)
            sigtable['Gorilla_access_time'] = 99999
        else:
            # sigtable = pd.read_excel(sigfile_GO, index_col=0)
            sigtable = load_excel_with_openpyxl_and_convert_to_pd_DataFrame(sigfile_ranked_GO_temp)
    elif sigtable == np.nan:
        print('ERROR: No sigtable')

    if len(sigtable) > 0:
        num_cells = int(sigtable.iloc[0, :]['num_cells'])
        num_genes = int(sigtable.iloc[0, :]['num_genes'])
        print('num_cells:', num_cells, ' num_genes:', num_genes)

        # Load genes_table
        with open('../objects/genes_tables/' + data_name + '_' + str(num_genes) + 'X' + str(
                num_cells) + '_' + filter_word + '.p', 'rb') as fp:
            genes_table = pickle.load(fp)

    pvals.sort(reverse=True)
    structs_filename = ('../objects/struct_lists'
                        + random_cellpair_set * ('_random_cellpair_set' + '_num_iters_' + str(num_iters))
                        + random_path * ('_random_path_' + str(path_len) + '_num_iters_' + str(num_iters))
                        + no_iterations * '_no_iterations'
                        + '/num_stds_thresh_' + str(num_stds_thresh) + '/mu_' + str(mu) + '/'
                        + (not ordered_pairs) * 'un' + 'ordered_pairs/'
                        + cell_selection + '/' + impute_method + '/' + data_name + '_' + filter_word + '_relevant_structs.txt')

    if ordered_pairs:
        structs_lst = list(set(sigtable.index[sigtable['structure_average_std'] <= struct_thr]).intersection(
            set(sigtable.index[sigtable['Gorilla_access_time'] == 99999])))
    else:
        structs_lst = list(set(sigtable.index[sigtable['log10_pavg_of_genes'] <= log10_pavg_of_genes_thr]).intersection(
            set(sigtable.index[sigtable['Gorilla_access_time'] == 99999])))

        # Read scRNA_data from file
    if read_scRNA_data_from_file:
        ind = [m.start() for m in re.finditer('/', data_dict[data_name])][-1]
        imputed_file = data_dict[data_name][:ind + 1] + impute_method + data_dict[data_name][
                                                                        ind:-4] + '_' + filter_word + '_imputed.txt'
        scRNA_data = load_data(imputed_file)
    elif scRNA_data == np.nan:
        print('ERROR: no scRNA_data matrix')

    normalized_scRNA_data = scRNA_data.subtract(scRNA_data.mean(axis=1), axis=0).divide(scRNA_data.std(axis=1), axis=0)
    normalized_scRNA_data.columns = [str(x) for x in normalized_scRNA_data.columns]

    for i in structs_lst:
        print('\n\n\nstructure', i)

        sets = sigtable.loc[i, 'sets']
        # print(sets_in_struct)
        sets = [tuple(list(map(int, s.strip(')').strip('(').split(', ')))) for s in
                sets.strip(']').strip('[').split('), (')]
        print(len(sets), 'sets')

        # find all the cells that are relevant to this structure:
        cells = list(set([item for sublist in sets for item in sublist]))
        print(len(cells), 'cells')

        print('cells:', cells)
        print('cells:', [str(c) for c in cells])
        curr_normalized_scRNA_data = normalized_scRNA_data.iloc[:, [c for c in cells]]

        # construct the transition matrix B that is relevant to this structure:
        curr_B = np.zeros((len(sets), len(cells)))
        for k, sett in enumerate(sets):
            curr_B[k, cells.index(sett[0])] = -1
            curr_B[k, cells.index(sett[1])] = 1
        print('curr_B.shape:', curr_B.shape)

        # multiply the normalized matrix by B to get the vector to compare with
        curr_S = np.dot(curr_B, curr_normalized_scRNA_data.transpose())
        print('curr_S.shape:', curr_S.shape)

        if not ordered_pairs:
            curr_S = np.abs(curr_S)

        # covariance matrix of curr_Y (curr_S is a sample of it)
        cov_mat = np.dot(curr_B, curr_B.transpose())
        print('cov_mat.shape:', cov_mat.shape)

        avg_vec = (1 / len(sets)) * np.ones((len(sets)))
        print('avg_vec.shape:', avg_vec.shape)
        average_diff_of_gene_std = np.sqrt(np.dot(np.dot(avg_vec.transpose(), cov_mat), avg_vec))
        print('average_diff_of_gene_std =', average_diff_of_gene_std)

        # average the differences for each gene
        average_diff_of_genes = np.dot(avg_vec.transpose(), curr_S)
        print('average_diff_of_genes[:5] =', average_diff_of_genes[:5])

        if ordered_pairs:
            genes_probs_logs = norm.logsf(average_diff_of_genes / average_diff_of_gene_std) / np.log(10)
        else:
            genes_probs_logs = norm.logsf(average_diff_of_genes / average_diff_of_gene_std) / np.log(10) + np.log10(2)

        genes_probs_logs_table = pd.Series(index=scRNA_data.index, data=genes_probs_logs).sort_values(ascending=True)
        print(genes_probs_logs_table[:10])

        target_list = str(list(genes_probs_logs_table.index)).replace("'", "").replace('[', '').replace(']',
                                                                                                        '').replace(
            ', ', '\n')
        # print(target_list)

        # text_file = open(data_name+'_struct'+str(i)+'_target_list.txt', "w")
        # text_file.write(target_list)
        # text_file.close()

        browser = mechanicalsoup.StatefulBrowser()
        browser.open(gorilla_url)

        # print(browser.url)
        # print(browser.page)

        browser.select_form()

        if data_name in ['Sharon2019', 'Treutlein2014']:
            browser["species"] = "MUS_MUSCULUS"
        elif data_name in ['Wagner2018', 'Wagner2018_']:
            browser["species"] = "DANIO_RERIO"
        else:
            browser["species"] = "HOMO_SAPIENS"

        # <option value="ARABIDOPSIS_THALIANA">Arabidopsis thaliana</option>
        # <option value="SACCHAROMYCES_CEREVISIAE">Saccharomyces cerevisiae</option>
        # <option value="CAENORHABDITIS_ELEGANS">Caenorhabditis elegans</option>
        # <option value="DROSOPHILA_MELANOGASTER">Drosophila melanogaster</option>
        # <option value="DANIO_RERIO">Danio rerio (Zebrafish)</option>
        # <option value="HOMO_SAPIENS" selected="selected">Homo sapiens</option>
        # <option value="MUS_MUSCULUS">Mus musculus</option>
        # <option value="RATTUS_NORVEGICUS">Rattus norvegicus</option>

        browser["run_mode"] = "mhg"
        # "hg" for Two unranked lists of genes (target and background lists)
        # "mhg" for Single ranked list of genes

        browser["target_set"] = target_list

        browser["db"] = "all"
        # "all" OR "proc" (process) OR "func" (function) OR "comp" (component)

        print(pvals[0], str('%f' % (pvals[0])))
        browser["pvalue_thresh"] = str('%f' % (pvals[0]))
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
        print('new_link:', new_link)

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

        print('new_link2:', new_link2)

        ind = [m.start() for m in re.finditer('/', new_link2)][-1]
        for ontology_init, ontology_html in zip(['proc', 'func', 'comp'],
                                                ['GOResultsPROCESS.html', 'GOResultsFUNCTION.html',
                                                 'GOResultsCOMPONENT.html']):
            new_link3 = new_link2[:ind + 1] + ontology_html
            print('new_link3', ontology_init, ':', new_link3)
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
                            [(a + ':' + b) for a, b in zip(df['GO term'], df['Description'])])
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

        # sigtable.to_excel(sigfile_GO)
        sigtable.to_excel(sigfile_ranked_GO_temp)
    # if ordered_pairs:
    #    sigtable = sigtable.sort_values(by='structure_average_std')
    # else:
    #    sigtable = sigtable.sort_values(by='log10_pmax_of_genes')
    sigtable.to_excel(sigfile_ranked_GO)
    try:
        os.remove(sigfile_ranked_GO_temp)
    except:
        pass


####################################################################################################################
def add_ranked_enriched_pathways(data_name, cell_selection, ordered_pairs, impute_method,
                                 read_sigtable_from_file=False, sigtable=np.nan,
                                 read_scRNA_data_from_file=False, scRNA_data=np.nan,
                                 start_from_scratch=False, pvals=[0.001, 0.0001, 0.00001],
                                 num_stds_thresh=1, mu=0.9, filt_mode='filt_by_3',
                                 filter_word='filtered', struct_thr=0.6, log10_pavg_of_genes_thr=-10,
                                 no_iterations=False, random_cellpair_set=False,
                                 random_path=False, path_len=None, num_iters=10000,
                                 load_conversion_table_from_file=True, structs_lst=None, permutation_num=50):
    # filt_mode = 'filt_by_3' OR 'filt_by_structure_average_std' OR 'filt_by_rule'

    print(data_name)
    print('\nadd_ranked_enriched_pathways!\n')

    sigfile_filt = ('../results'
                    + random_cellpair_set * ('_random_cellpair_set' + '_num_iters_' + str(num_iters))
                    + random_path * ('_random_path_' + str(path_len) + '_num_iters_' + str(num_iters))
                    + no_iterations * '_no_iterations'
                    + '/sigtables/num_stds_thresh_' + str(num_stds_thresh) + '/mu_' + str(mu) + '/' +
                    (not ordered_pairs) * 'un' + 'ordered_pairs/'
                    + cell_selection + '/' + impute_method + '/' + data_name + '_' + filter_word + '_significance_table_filt.xlsx')
    sigfile_ranked_gsea = sigfile_filt[:-5] + '_gsea_' + str(permutation_num) + sigfile_filt[-5:]
    sigfile_ranked_gsea_temp = sigfile_filt[:-5] + '_gsea_' + str(permutation_num) + '_temp' + sigfile_filt[-5:]
    print('\n\n\n', sigfile_filt)
    print('\n\n\n', sigfile_ranked_gsea)
    print('\n\n\n', sigfile_ranked_gsea_temp)

    if read_sigtable_from_file:
        if start_from_scratch or (not os.path.exists(sigfile_ranked_gsea_temp)):
            sigtable = pd.read_excel(sigfile_filt, index_col=0)
            sigtable['gsea_comp_time'] = 99999
        else:
            # sigtable = pd.read_excel(sigfile_GO, index_col=0)
            sigtable = load_excel_with_openpyxl_and_convert_to_pd_DataFrame(sigfile_ranked_gsea_temp)
    elif sigtable == np.nan:
        print('ERROR: No sigtable')

    if len(sigtable) > 0:
        num_cells = int(sigtable.iloc[0, :]['num_cells'])
        num_genes = int(sigtable.iloc[0, :]['num_genes'])
        print('num_cells:', num_cells, ' num_genes:', num_genes)

        # Load genes_table
        with open('../objects/genes_tables/' + data_name + '_' + str(num_genes) + 'X' + str(
                num_cells) + '_' + filter_word + '.p', 'rb') as fp:
            genes_table = pickle.load(fp)

    pvals.sort(reverse=True)
    structs_filename = ('../objects/struct_lists'
                        + random_cellpair_set * ('_random_cellpair_set' + '_num_iters_' + str(num_iters))
                        + random_path * ('_random_path_' + str(path_len) + '_num_iters_' + str(num_iters))
                        + no_iterations * '_no_iterations'
                        + '/num_stds_thresh_' + str(num_stds_thresh) + '/mu_' + str(mu) + '/'
                        + (not ordered_pairs) * 'un' + 'ordered_pairs/'
                        + cell_selection + '/' + impute_method + '/' + data_name + '_' + filter_word + '_relevant_structs.txt')

    if structs_lst is None:
        if ordered_pairs:
            structs_lst = list(set(sigtable.index[sigtable['structure_average_std'] <= struct_thr]).intersection(
                set(sigtable.index[sigtable['gsea_comp_time'] == 99999])))
        else:
            structs_lst = list(
                set(sigtable.index[sigtable['log10_pavg_of_genes'] <= log10_pavg_of_genes_thr]).intersection(
                    set(sigtable.index[sigtable['gsea_comp_time'] == 99999])))
    else:
        structs_lst = list(set(structs_lst).intersection(set(sigtable.index[sigtable['gsea_comp_time'] == 99999])))

    # Read scRNA_data from file
    if read_scRNA_data_from_file:
        ind = [m.start() for m in re.finditer('/', data_dict[data_name])][-1]
        imputed_file = data_dict[data_name][:ind + 1] + impute_method + data_dict[data_name][
                                                                        ind:-4] + '_' + filter_word + '_imputed.txt'
        scRNA_data = load_data(imputed_file)
    elif scRNA_data == np.nan:
        print('ERROR: no scRNA_data matrix')

    normalized_scRNA_data = scRNA_data.subtract(scRNA_data.mean(axis=1), axis=0).divide(scRNA_data.std(axis=1), axis=0)
    normalized_scRNA_data.columns = [str(x) for x in normalized_scRNA_data.columns]

    # Get a conversion table for the genes (original gene type -> gene symbol)
    ind = [m.start() for m in re.finditer('/', data_dict[data_name])][-1]
    conversion_table_file = data_dict[data_name][:ind + 1] + impute_method + data_dict[data_name][
                                                                             ind:-4] + '_' + filter_word + '_imputed_conv_table.xlsx'

    if load_conversion_table_from_file and os.path.exists(conversion_table_file):
        print('Loading conversion table from file...')
        conv_table = pd.read_excel(conversion_table_file, index_col=0)
        print('len(conv_table):', len(conv_table))
        print('conv_table[:5]:', conv_table[:5])
    else:
        print('Computing conversion table...')
        # prepare for conversion to gene symbols
        bm = Biomart()
        # view validated marts
        marts = bm.get_marts()
        # view validated dataset
        datasets = bm.get_datasets(mart='ENSEMBL_MART_ENSEMBL')
        # view validated attributes
        attrs = bm.get_attributes(dataset='hsapiens_gene_ensembl')
        # view validated filters
        filters = bm.get_filters(dataset='hsapiens_gene_ensembl')
        # convert to gene symbols
        # query results
        conv_table = bm.query(dataset='hsapiens_gene_ensembl',
                              attributes=['ensembl_gene_id', 'external_gene_name'],
                              filters={'ensembl_gene_id': list(scRNA_data.index)})
        print('len(conv_table):', len(conv_table))
        print('conv_table[:5]:', conv_table[:5])
        conv_table.to_excel(conversion_table_file)

        # results folder
    gsea_folder = ('../results'
                   + random_cellpair_set * '_random_cellpair_set'
                   + random_path * ('_random_path_' + str(path_len))
                   + no_iterations * '_no_iterations'
                   + '/gsea/num_stds_thresh_' + str(num_stds_thresh) + '/mu_' + str(mu) + '/'
                   + (not ordered_pairs) * 'un'
                   + 'ordered_pairs/' + cell_selection + '/' + impute_method + '/' + data_name + '_' + filter_word + '/' + str(
                permutation_num))
    if not os.path.exists(gsea_folder):
        os.makedirs(gsea_folder)

    print(len(structs_lst), 'structures')
    for i in structs_lst:
        print('\n\n\nstructure', i)

        ###################################################################################################
        # compute ranked list of genes for this structure (based on their p-values in relation to the structure)
        sets = sigtable.loc[i, 'sets']
        # print(sets_in_struct)
        sets = [tuple(list(map(int, s.strip(')').strip('(').split(', ')))) for s in
                sets.strip(']').strip('[').split('), (')]
        print(len(sets), 'sets')

        # find all the cells that are relevant to this structure:
        cells = list(set([item for sublist in sets for item in sublist]))
        print(len(cells), 'cells')

        print('cells:', cells)
        print('cells:', [str(c) for c in cells])
        curr_normalized_scRNA_data = normalized_scRNA_data.iloc[:, [c for c in cells]]

        # construct the transition matrix B that is relevant to this structure:
        curr_B = np.zeros((len(sets), len(cells)))
        for k, sett in enumerate(sets):
            curr_B[k, cells.index(sett[0])] = -1
            curr_B[k, cells.index(sett[1])] = 1
        print('curr_B.shape:', curr_B.shape)

        # multiply the normalized matrix by B to get the vector to compare with
        curr_S = np.dot(curr_B, curr_normalized_scRNA_data.transpose())
        print('curr_S.shape:', curr_S.shape)

        if not ordered_pairs:
            curr_S = np.abs(curr_S)

        # covariance matrix of curr_Y (curr_S is a sample of it)
        cov_mat = np.dot(curr_B, curr_B.transpose())
        print('cov_mat.shape:', cov_mat.shape)

        avg_vec = (1 / len(sets)) * np.ones((len(sets)))
        print('avg_vec.shape:', avg_vec.shape)
        average_diff_of_gene_std = np.sqrt(np.dot(np.dot(avg_vec.transpose(), cov_mat), avg_vec))
        print('average_diff_of_gene_std =', average_diff_of_gene_std)

        # average the differences for each gene
        average_diff_of_genes = np.dot(avg_vec.transpose(), curr_S)
        print('average_diff_of_genes[:5] =', average_diff_of_genes[:5])

        if ordered_pairs:
            genes_probs_logs = norm.logsf(average_diff_of_genes / average_diff_of_gene_std) / np.log(10)
        else:
            genes_probs_logs = norm.logsf(average_diff_of_genes / average_diff_of_gene_std) / np.log(10) + np.log10(2)

        genes_probs_logs_table = pd.Series(index=scRNA_data.index, data=genes_probs_logs).sort_values(ascending=True)
        print(genes_probs_logs_table[:10])

        target_list = list(genes_probs_logs_table.index)
        print('target_list[:5]:\n', target_list[:5])

        print(len(target_list), 'genes')
        print(len(set(target_list)), 'genes')

        # text_file = open(data_name+'_struct'+str(i)+'_target_list.txt', "w")
        # text_file.write(target_list)
        # text_file.close()
        ###################################################################################################
        def order(col, target_list=target_list):
            # input:
            # col- pd.Series with ENSEMBL gene names
            # target_list- list of ENSEMBL gene names

            # output: pd.Series in which each of the gene names was replaced by its position in target_list
            key = dict(zip(target_list, range(len(target_list))))
            return (pd.Series([key[element] for element in col]))

        ###################################################################################################
        conv_table_for_struct = conv_table.sort_values(by='ensembl_gene_id', key=order)
        conv_table_for_struct = conv_table_for_struct[~conv_table_for_struct.duplicated(['external_gene_name'])]
        print('len(conv_table_for_struct):', len(conv_table_for_struct))
        print('conv_table_for_struct[:5]:\n', conv_table_for_struct[:5])
        print('conv_table_for_struct[-5:]:\n', conv_table_for_struct[-5:])
        target_list_df = pd.DataFrame(
            {'0': conv_table_for_struct['external_gene_name'], '1': range(len(conv_table_for_struct), 0, -1)})
        target_list_df.index = list(range(len(target_list_df)))
        print('target_list_df[:5]:\n', target_list_df[:5])
        print('target_list_df[-5:]:\n', target_list_df[-5:])

        print(len(target_list_df['0']), 'genes')
        print(len(set(list(target_list_df['0']))), 'genes')

        ###################################################################################################
        # run GSEA
        start = time()
        if not os.path.exists(gsea_folder + '/struct_' + str(i)):
            os.makedirs(gsea_folder + '/struct_' + str(i))
        pre_res = gseapy.prerank(rnk=target_list_df, gene_sets='./gene_sets/c2.cp.v7.2.symbols.gmt',
                                 processes=4,
                                 permutation_num=permutation_num,  # reduce number to speed up testing
                                 outdir=gsea_folder + '/struct_' + str(i), format='png', seed=6)
        print('GSEA computation (with', permutation_num, 'permutations) took', (time() - start) / 60, 'minutes\n')

        # access results through obj.res2d attribute or obj.results
        # print(pre_res.res2d.sort_index().head())
        df = pre_res.res2d
        print(df.head())

        ###################################################################################################
        # save significant pathways
        for pval in pvals:
            df = df[df['fdr'].astype(float) <= pval]
            if len(df) > 0:
                sigtable.loc[i, 'pathways_below_' + str(pval)] = str(
                    [(a + ':' + b) for a, b in zip(df.index, df['fdr'].astype(str))])
            else:
                sigtable.loc[i, 'pathways_below_' + str(pval)] = 'NO PATWAYS'

        sigtable.loc[i, 'gsea_comp_time'] = datetime.now().strftime("%d/%m/%Y, %H:%M:%S")

        # sigtable.to_excel(sigfile_GO)
        sigtable.to_excel(sigfile_ranked_gsea_temp)
    # if ordered_pairs:
    #    sigtable = sigtable.sort_values(by='structure_average_std')
    # else:
    #    sigtable = sigtable.sort_values(by='log10_pmax_of_genes')

    sigtable.to_excel(sigfile_ranked_gsea)
    try:
        os.remove(sigfile_ranked_gsea_temp)
    except:
        pass

    #### list of data sets #####################################################################################################

class Data_dict:
    def __init__(self):
        self.dict = {'Ramos2019_GSM3453214': '../data/Ramos2019/GSE122031_RAW/GSM3453214_A549_4h_MOI20_NormCounts.csv',
'Ramos2019_GSM3453215': '../data/Ramos2019/GSE122031_RAW/GSM3453215_A549_4h_moi02_NormCounts.csv',
'Ramos2019_GSM3453216': '../data/Ramos2019/GSE122031_RAW/GSM3453216_A549_4h_Mock_NormCounts.csv',
'Ramos2019_GSM3453217': '../data/Ramos2019/GSE122031_RAW/GSM3453217_A549_12h_moi20_NormCounts.csv',
'Ramos2019_GSM3453218': '../data/Ramos2019/GSE122031_RAW/GSM3453218_A549_12h_moi02_NormCounts.csv',
'Ramos2019_GSM3453219': '../data/Ramos2019/GSE122031_RAW/GSM3453219_A549_12h_Mock_NormCounts.csv',
'Ramos2019_mock_raw': '../data/Ramos2019/GSE122031_RAW/A549_Mock_RawData_human.txt',
'Ramos2019_MOI02_raw': '../data/Ramos2019/GSE122031_RAW/A549_MOI02_RawData_human.txt',
'Ramos2019_MOI20_raw': '../data/Ramos2019/GSE122031_RAW/A549_MOI20_RawData_human.txt',
'Ramos2019_mock_4h_raw': '../data/Ramos2019/GSE122031_RAW/A549_4h_Mock_RawData_human.txt',
'Ramos2019_mock_12h_raw': '../data/Ramos2019/GSE122031_RAW/A549_12h_Mock_RawData_human.txt',
'Ramos2019_MOI02_4h_raw': '../data/Ramos2019/GSE122031_RAW/A549_4h_MOI02_RawData_human.txt',
'Ramos2019_MOI02_12h_raw': '../data/Ramos2019/GSE122031_RAW/A549_12h_MOI02_RawData_human.txt',
'Ramos2019_MOI20_4h_raw': '../data/Ramos2019/GSE122031_RAW/A549_4h_MOI20_RawData_human.txt',
'Ramos2019_MOI20_12h_raw': '../data/Ramos2019/GSE122031_RAW/A549_12h_MOI20_RawData_human.txt',
'Ramos_all': '../data/Ramos2019/GSE122031_RAW/A549_RawData_human_2000_from_each.txt',
'BenMoshe2019_GSE122084_RAW': '../data/scRNA_Weizmann_data/BenMoshe2019_data/GSE122084_RAW.txt',
'BenMoshe2019_GSE122084_RAW_naive': '../data/BenMoshe2019/GSE122084_RAW_GSM3454528_naive_cells.txt',
'BenMoshe2019_GSE122084_RAW_exposed1': '../data/BenMoshe2019/GSE122084_RAW_GSM3454529_Salmonella_exposed_cells.txt',
'BenMoshe2019_GSE122084_RAW_exposed2': '../data/BenMoshe2019/GSE122084_RAW_GSM3855868_Salmonella_exposed_cells.txt',
'BenMoshe2019_GSE122084_RAW_exposed2_': '../data/BenMoshe2019/GSE122084_RAW_GSM3855868_Salmonella_exposed_cells.txt',
'Zhang2019_inDrop1': '../data/Zhang2019/GSE111912_RAW/GSM3044886_GeneExp.UMIs.inDrop1.txt',
'Zhang2019_inDrop2': '../data/Zhang2019/GSE111912_RAW/GSM3044887_GeneExp.UMIs.inDrop2.txt',
'Zhang2019_DropSeq1': '../data/Zhang2019/GSE111912_RAW/GSM3044888_GeneExp.UMIs.DropSeq1.txt',
'Zhang2019_DropSeq2': '../data/Zhang2019/GSE111912_RAW/GSM3044889_GeneExp.UMIs.DropSeq2.txt',
'Zhang2019_DropSeq3': '../data/Zhang2019/GSE111912_RAW/GSM3044890_GeneExp.UMIs.DropSeq3.txt',
'Zhang2019_10x1': '../data/Zhang2019/GSE111912_RAW/GSM3044891_GeneExp.UMIs.10X1.txt',
'Zhang2019_10x2': '../data/Zhang2019/GSE111912_RAW/GSM3044892_GeneExp.UMIs.10X2.txt',
'Steuerman2018_GSM2884119': '../data/Steuerman2018/GSE107947_RAW/GSM2884119_PBS_CD45p_treated_umis.txt',
'Steuerman2018_GSM2884120': '../data/Steuerman2018/GSE107947_RAW/GSM2884120_PBS_CD45n_treated_umis.txt',
'Steuerman2018_GSM2884121': '../data/Steuerman2018/GSE107947_RAW/GSM2884121_Flu_treated_CD45p_48h_rep1_umis.txt',
'Steuerman2018_GSM2884122': '../data/Steuerman2018/GSE107947_RAW/GSM2884122_Flu_treated_CD45n_48h_rep1_umis.txt',
'Steuerman2018_GSM2884123': '../data/Steuerman2018/GSE107947_RAW/GSM2884123_Flu_treated_CD45p_48h_rep2_umis.txt',
'Steuerman2018_GSM2884124': '../data/Steuerman2018/GSE107947_RAW/GSM2884124_Flu_treated_CD45n_48h_rep2_umis.txt',
'Steuerman2018_GSM2884125': '../data/Steuerman2018/GSE107947_RAW/GSM2884125_Flu_treated_CD45p_72h_umis.txt',
'Steuerman2018_GSM2884126': '../data/Steuerman2018/GSE107947_RAW/GSM2884126_Flu_treated_CD45n_72h_umis.txt',
'Steuerman2018_GSM2884127': '../data/Steuerman2018/GSE107947_RAW/GSM2884127_Flu_treated_CD45p_48h_Irf7KO_umis.txt',
'Steuerman2018_GSM2884128': '../data/Steuerman2018/GSE107947_RAW/GSM2884128_Flu_treated_CD45n_48h_Irf7KO_umis.txt',
'10x_PBMCs': '../data/10x_PBMCs/connect_5k_pbmc_NGSC3_ch1_filtered_feature_bc_matrix.txt',
'Liao2020_C141': '../data/Liao2020/GSE145926_RAW/GSM4339769_C141_filtered_feature_bc_matrix.txt',
'Liao2020_C142': '../data/Liao2020/GSE145926_RAW/GSM4339770_C142_filtered_feature_bc_matrix.txt',
'Liao2020_C143': '../data/Liao2020/GSE145926_RAW/GSM4339771_C143_filtered_feature_bc_matrix.txt',
'Liao2020_C144': '../data/Liao2020/GSE145926_RAW/GSM4339772_C144_filtered_feature_bc_matrix.txt',
'Liao2020_C145': '../data/Liao2020/GSE145926_RAW/GSM4339773_C145_filtered_feature_bc_matrix.txt',
'Liao2020_C146': '../data/Liao2020/GSE145926_RAW/GSM4339774_C146_filtered_feature_bc_matrix.txt',
'Liao2020_C51': '../data/Liao2020/GSE145926_RAW/GSM4475048_C51_filtered_feature_bc_matrix.txt',
'Liao2020_C52': '../data/Liao2020/GSE145926_RAW/GSM4475049_C52_filtered_feature_bc_matrix.txt',
'Liao2020_C100': '../data/Liao2020/GSE145926_RAW/GSM4475050_C100_filtered_feature_bc_matrix.txt',
'Liao2020_C148': '../data/Liao2020/GSE145926_RAW/GSM4475051_C148_filtered_feature_bc_matrix.txt',
'Liao2020_C149': '../data/Liao2020/GSE145926_RAW/GSM4475052_C149_filtered_feature_bc_matrix.txt',
'Liao2020_C152': '../data/Liao2020/GSE145926_RAW/GSM4475053_C152_filtered_feature_bc_matrix.txt',
'Wyler2020_Calu3_mock': '../data/Wyler2020/GSE148729_Calu3_scRNAseq_morethan1000genes_rawcounts_tr/Calu3_scRNAseq_morethan1000genes_rawcounts_tr_mock_without_MIRs.txt',
'Wyler2020_Calu3_s1': '../data/Wyler2020/GSE148729_Calu3_scRNAseq_morethan1000genes_rawcounts_tr/Calu3_scRNAseq_morethan1000genes_rawcounts_tr_s1_without_MIRs.txt',
'Wyler2020_Calu3_s2': '../data/Wyler2020/GSE148729_Calu3_scRNAseq_morethan1000genes_rawcounts_tr/Calu3_scRNAseq_morethan1000genes_rawcounts_tr_s2_without_MIRs.txt',
'Wyler2020_H1299_mock': '../data/Wyler2020/GSE148729_H1299_scRNAseq_morethan2000genes_rawcounts_tr/H1299_scRNAseq_morethan2000genes_rawcounts_tr_mock_without_MIRs.txt',
'Wyler2020_H1299_s1': '../data/Wyler2020/GSE148729_H1299_scRNAseq_morethan2000genes_rawcounts_tr/H1299_scRNAseq_morethan2000genes_rawcounts_tr_s1_without_MIRs.txt',
'Wyler2020_H1299_s2': '../data/Wyler2020/GSE148729_H1299_scRNAseq_morethan2000genes_rawcounts_tr/H1299_scRNAseq_morethan2000genes_rawcounts_tr_s2_without_MIRs.txt',
'Wyler2020_Calu3_s1_4h-A':                      '../data/Wyler2020/GSE148729_Calu3_scRNAseq_morethan1000genes_rawcounts_tr/Calu3_scRNAseq_morethan1000genes_rawcounts_tr_s1_without_MIRs_4h-A.txt',
'Wyler2020_Calu3_s1_8h-B': '../data/Wyler2020/GSE148729_Calu3_scRNAseq_morethan1000genes_rawcounts_tr/Calu3_scRNAseq_morethan1000genes_rawcounts_tr_s1_without_MIRs_8h-B.txt',
 'Wyler2020_Calu3_s1_8h-A': '../data/Wyler2020/GSE148729_Calu3_scRNAseq_morethan1000genes_rawcounts_tr/Calu3_scRNAseq_morethan1000genes_rawcounts_tr_s1_without_MIRs_8h-A.txt',
 'Wyler2020_Calu3_s1_12h-B': '../data/Wyler2020/GSE148729_Calu3_scRNAseq_morethan1000genes_rawcounts_tr/Calu3_scRNAseq_morethan1000genes_rawcounts_tr_s1_without_MIRs_12h-B.txt',
 'Wyler2020_Calu3_s1_4h-B': '../data/Wyler2020/GSE148729_Calu3_scRNAseq_morethan1000genes_rawcounts_tr/Calu3_scRNAseq_morethan1000genes_rawcounts_tr_s1_without_MIRs_4h-B.txt',
 'Wyler2020_Calu3_s1_12h-A': '../data/Wyler2020/GSE148729_Calu3_scRNAseq_morethan1000genes_rawcounts_tr/Calu3_scRNAseq_morethan1000genes_rawcounts_tr_s1_without_MIRs_12h-A.txt',
'Wyler2020_Calu3_s2_4h-B': '../data/Wyler2020/GSE148729_Calu3_scRNAseq_morethan1000genes_rawcounts_tr/Calu3_scRNAseq_morethan1000genes_rawcounts_tr_s2_without_MIRs_4h-B.txt',
 'Wyler2020_Calu3_s2_8h-A': '../data/Wyler2020/GSE148729_Calu3_scRNAseq_morethan1000genes_rawcounts_tr/Calu3_scRNAseq_morethan1000genes_rawcounts_tr_s2_without_MIRs_8h-A.txt',
 'Wyler2020_Calu3_s2_12h-A': '../data/Wyler2020/GSE148729_Calu3_scRNAseq_morethan1000genes_rawcounts_tr/Calu3_scRNAseq_morethan1000genes_rawcounts_tr_s2_without_MIRs_12h-A.txt',
 'Wyler2020_Calu3_s2_4h-A': '../data/Wyler2020/GSE148729_Calu3_scRNAseq_morethan1000genes_rawcounts_tr/Calu3_scRNAseq_morethan1000genes_rawcounts_tr_s2_without_MIRs_4h-A.txt',
 'Wyler2020_Calu3_s2_12h-B': '../data/Wyler2020/GSE148729_Calu3_scRNAseq_morethan1000genes_rawcounts_tr/Calu3_scRNAseq_morethan1000genes_rawcounts_tr_s2_without_MIRs_12h-B.txt',
 'Wyler2020_Calu3_s2_8h-B': '../data/Wyler2020/GSE148729_Calu3_scRNAseq_morethan1000genes_rawcounts_tr/Calu3_scRNAseq_morethan1000genes_rawcounts_tr_s2_without_MIRs_8h-B.txt',
'Zavidij2020_MGUS-1': '../data/Zavidij2020/GSE124310_RAW/GSM3528753_MGUS-1.138N.filtered_gene_bc_matrices.txt',
'Zavidij2020_MGUS-2': '../data/Zavidij2020/GSE124310_RAW/GSM3528755_MGUS-2.138N45P.filtered_gene_bc_matrices.txt',
'Zavidij2020_MGUS-3': '../data/Zavidij2020/GSE124310_RAW/GSM3528757_MGUS-3.138N45P.filtered_gene_bc_matrices.txt',
'Zavidij2020_MGUS-4': '../data/Zavidij2020/GSE124310_RAW/GSM3528759_MGUS-4.138N45P.filtered_gene_bc_matrices.txt',
'Zavidij2020_MGUS-6': '../data/Zavidij2020/GSE124310_RAW/GSM3528762_MGUS-6.138N45P.filtered_gene_bc_matrices.txt',
'Zavidij2020_MM-1': '../data/Zavidij2020/GSE124310_RAW/GSM3528764_MM-1.138N.filtered_gene_bc_matrices.txt',
'Zavidij2020_MM-3': '../data/Zavidij2020/GSE124310_RAW/GSM3528767_MM-3.138N45P.filtered_gene_bc_matrices.txt',
'Zavidij2020_MM-4': '../data/Zavidij2020/GSE124310_RAW/GSM3528769_MM-4.138N45P.filtered_gene_bc_matrices.txt',
'Zavidij2020_MM-5': '../data/Zavidij2020/GSE124310_RAW/GSM3528771_MM-5.138N45P.filtered_gene_bc_matrices.txt',
'Zavidij2020_MM-6': '../data/Zavidij2020/GSE124310_RAW/GSM3528773_MM-6.138N45P.filtered_gene_bc_matrices.txt',
'Zavidij2020_MM-7': '../data/Zavidij2020/GSE124310_RAW/GSM3528775_MM-7.138N45P.filtered_gene_bc_matrices.txt',
'Zavidij2020_MM-8': '../data/Zavidij2020/GSE124310_RAW/GSM3528777_MM-8.138N45P.filtered_gene_bc_matrices.txt',
'Zavidij2020_NBM-1': '../data/Zavidij2020/GSE124310_RAW/GSM3528779_NBM-1.138N45P.filtered_gene_bc_matrices.txt',
'Zavidij2020_NBM-10': '../data/Zavidij2020/GSE124310_RAW/GSM3528781_NBM-10.138N.filtered_gene_bc_matrices.txt',
'Zavidij2020_NBM-11': '../data/Zavidij2020/GSE124310_RAW/GSM3528783_NBM-11.138N.filtered_gene_bc_matrices.txt',
'Zavidij2020_NBM-2': '../data/Zavidij2020/GSE124310_RAW/GSM3528785_NBM-2.138N45P.filtered_gene_bc_matrices.txt',
'Zavidij2020_NBM-3': '../data/Zavidij2020/GSE124310_RAW/GSM3528787_NBM-3.138N45P.filtered_gene_bc_matrices.txt',
'Zavidij2020_NBM-4': '../data/Zavidij2020/GSE124310_RAW/GSM3528789_NBM-4.138N.filtered_gene_bc_matrices.txt',
'Zavidij2020_NBM-5': '../data/Zavidij2020/GSE124310_RAW/GSM3528791_NBM-5.138N.filtered_gene_bc_matrices.txt',
'Zavidij2020_NBM-8': '../data/Zavidij2020/GSE124310_RAW/GSM3528794_NBM-8.138N.filtered_gene_bc_matrices.txt',
'Zavidij2020_NBM-9': '../data/Zavidij2020/GSE124310_RAW/GSM3528796_NBM-9.138N.filtered_gene_bc_matrices.txt',
'Zavidij2020_SMMh-10': '../data/Zavidij2020/GSE124310_RAW/GSM3528798_SMMh-10.138N45P.filtered_gene_bc_matrices.txt',
'Zavidij2020_SMMh-2': '../data/Zavidij2020/GSE124310_RAW/GSM3528800_SMMh-2.138N45P.filtered_gene_bc_matrices.txt',
'Zavidij2020_SMMh-3': '../data/Zavidij2020/GSE124310_RAW/GSM3528802_SMMh-3.138N.filtered_gene_bc_matrices.txt',
'Zavidij2020_SMMh-4': '../data/Zavidij2020/GSE124310_RAW/GSM3528804_SMMh-4.138N.filtered_gene_bc_matrices.txt',
'Zavidij2020_SMMh-6': '../data/Zavidij2020/GSE124310_RAW/GSM3528807_SMMh-6.138N.filtered_gene_bc_matrices.txt',
'Zavidij2020_SMMh-7': '../data/Zavidij2020/GSE124310_RAW/GSM3528809_SMMh-7.138N.filtered_gene_bc_matrices.txt',
'Zavidij2020_SMMh-8': '../data/Zavidij2020/GSE124310_RAW/GSM3528810_SMMh-8.138N45P.filtered_gene_bc_matrices.txt',
'Zavidij2020_SMMh-9': '../data/Zavidij2020/GSE124310_RAW/GSM3528812_SMMh-9.138N45P.filtered_gene_bc_matrices.txt',
'Zavidij2020_SMMl-1': '../data/Zavidij2020/GSE124310_RAW/GSM3528814_SMMl-1.138N.filtered_gene_bc_matrices.txt',
'Zavidij2020_SMMl-2': '../data/Zavidij2020/GSE124310_RAW/GSM3528816_SMMl-2.138N.filtered_gene_bc_matrices.txt',
'Zavidij2020_SMMl-3': '../data/Zavidij2020/GSE124310_RAW/GSM3528818_SMMl-3.138N45P.filtered_gene_bc_matrices.txt',
'Random': '../data/Random/random.txt',
'RinaBulkRNAData': '../data/RinaBulkRNAData/NonAvgCountsAllDataBioInfoCourse.txt',
'yaellab_differentiation_bulk': '../data/yaellab_differentiation/bulk_differentiation_ver1.txt',
'yaeldiffbulk2': '../data/yaellab_differentiation/bulk_differentiation_ver2.csv',
'Soumillon2014_D1': '../data/Soumillon2014/GSE53638_D1_UMI.txt',
'Han2018_H9': '../data/Han2018/H9_all.txt',
'yaellab_covid19': '../data/yaelab_covid19/counts_GRch38_Poria_and_first3_without_cl.txt',
'Sharon2019': '../data/Sharon2019/GSE121416_genes.fpkm_table.txt/GSE121416_genes_fpkm_table_ENSMUSG_only.txt',
'Treutlein2014': '../data/Treutlein2014/GSE52583_RAW/GSE52583_RAW_all.txt',
'Wagner2018': '../data/Wagner2018/GSE112294_RAW/all_1500_from_each_time.txt',
'Wagner2018_': '../data/Wagner2018/GSE112294_RAW/all_1000_from_each_time.txt',
'splatter': '../data/synthetic_splatter/splatter.txt',
'splatter2': '../data/synthetic_splatter/splatter2.txt',
'splatter3': '../data/synthetic_splatter/splatter3.txt',
'splatter5': '../data/synthetic_splatter/splatter5.txt',
'splatter6': '../data/synthetic_splatter/splatter6.txt',
'splatter7': '../data/synthetic_splatter/splatter7.txt',
'synthetic1': '../data/mysynthetic/synthetic1.txt',
'synthetic2': '../data/mysynthetic/synthetic2.txt',
'synthetic3': '../data/mysynthetic/synthetic3.txt',
'synthetic4': '../data/mysynthetic/synthetic4.txt',
'synthetic5': '../data/mysynthetic/synthetic5.txt',
'synthetic6': '../data/mysynthetic/synthetic6.txt',
'synthetic8': '../data/mysynthetic/synthetic8.txt',
'synthetic9': '../data/mysynthetic/synthetic9.txt',
'synthetic10': '../data/mysynthetic/synthetic10.txt',
'synthetic_bulk_1': '../data/mysynthetic/synthetic_bulk_1.txt',
'synthetic_bulk_2': '../data/mysynthetic/synthetic_bulk_2.txt',
'2Dheart': '../data/visium_spatial/V1_Human_Heart/genecounts_raw.txt',
'2Dbreast1': '../data/visium_spatial/V1_Breast_Cancer_Block_A_Section_1/genecounts_raw.txt',
'2Dbreast2': '../data/visium_spatial/V1_Breast_Cancer_Block_A_Section_2/genecounts_filtered.txt',
'2DMouseKidney': '../data/visium_spatial/V1_Mouse_Kidney/genecounts_filtered.txt',
'2Dlymph': '../data/visium_spatial/V1_Human_Lymph_Node/genecounts_filtered.txt',
'2DMouseBrainC': '../data/visium_spatial/V1_Adult_Mouse_Brain/genecounts_filtered.txt',
'2DMouseBrainSA1': '../data/visium_spatial/V1_Mouse_Brain_Sagittal_Anterior/genecounts_filtered.txt',
'2DMouseBrainSA2': '../data/visium_spatial/V1_Mouse_Brain_Sagittal_Anterior_Section_2/genecounts_filtered.txt',
'2DMouseBrainSP1': '../data/visium_spatial/V1_Mouse_Brain_Sagittal_Posterior/genecounts_filtered.txt',
'2DMouseBrainSP2': '../data/visium_spatial/V1_Mouse_Brain_Sagittal_Posterior_Section_2/genecounts_filtered.txt',
'2DFFPE_BreastCancer': '../data/visium_spatial/V1_FFPE_human_breast_cancer/genecounts_filtered.txt',
'2DFFPE_ProstateCancer': '../data/visium_spatial/V1_FFPE_human_prostate_cancer/genecounts_filtered.txt',
'2DFFPE_NormalProstate': '../data/visium_spatial/V1_FFPE_normal_human_prostate/genecounts_filtered.txt',
'2DBerglund2018_P11': '../data/Berglund2018/prostate-twelve/P1.1.tsv',
'2DBerglund2018_P12': '../data/Berglund2018/prostate-twelve/P1.2.tsv',
'2DBerglund2018_P13': '../data/Berglund2018/prostate-twelve/P1.3.tsv',
'2DBerglund2018_P21': '../data/Berglund2018/prostate-twelve/P2.1.tsv',
'2DBerglund2018_P23': '../data/Berglund2018/prostate-twelve/P2.3.tsv',
'2DBerglund2018_P24': '../data/Berglund2018/prostate-twelve/P2.4.tsv',
'2DBerglund2018_P31': '../data/Berglund2018/prostate-twelve/P3.1.tsv',
'2DBerglund2018_P32': '../data/Berglund2018/prostate-twelve/P3.2.tsv',
'2DBerglund2018_P33': '../data/Berglund2018/prostate-twelve/P3.3.tsv',
'2DBerglund2018_P41': '../data/Berglund2018/prostate-twelve/P4.1.tsv',
'2DBerglund2018_P42': '../data/Berglund2018/prostate-twelve/P4.2.tsv',
'2DBerglund2018_P43': '../data/Berglund2018/prostate-twelve/P4.3.tsv',
'2DBerglund2018_Patient2_E2': '../data/Berglund2018/prostate-twelve/Patient2_E2.tsv',
'2DBerglund2018_Patient3_D1': '../data/Berglund2018/prostate-twelve/Patient3_D1.tsv',
'yaellab_diff_sc': '../../Yael_lab_differentiation/data/dataset1/single_cells_all_times.csv',
'yaellab_diff_scF': '../../Yael_lab_differentiation/data/dataset1/single_cells_all_timesF.csv',
'yaellab_diff_sc21F': '../../Yael_lab_differentiation/data/dataset1/single_cells21F.txt',
'yaellab_diff_sc21': '../../Yael_lab_differentiation/data/dataset1/single_cells21.txt',
'yaellab_diff_sc21_96': '../../Yael_lab_differentiation/data/dataset2/scDiffDay21_96cells.csv',
'Iris_Hypo': '../data/Iris/0.1.read_counts_Hypo_0.txt',
'Iris_Norm': '../data/Iris/0.1.read_counts_Norm_V1_1.txt',}
        self.unfit_datasets = ['BenMoshe2019_GSE122084_RAW']

'''        
class Set_num_cells_bidict(object):
    def __init__(self,num_cells):
        cells_set_num_table = pd.DataFrame(data=np.zeros((np.power(num_cells,2),3)), columns=['c1','c2','c1_c2'])
        cells_set_num_table.loc[:,'c1'] = np.repeat(np.array(range(num_cells)),num_cells)
        cells_set_num_table.loc[:,'c2'] = np.stack([np.array(range(num_cells)) for _ in range(num_cells)], axis=0).flatten() 
        cells_set_num_table = cells_set_num_table[cells_set_num_table['c1']!=cells_set_num_table['c2']]
        cells_set_num_table['c1_c2'] = cells_set_num_table['c1'].astype(str) + '_' + cells_set_num_table['c2'].astype(str)
        cells_set_num_table.index = np.array(range(len(cells_set_num_table)))

        set_nums = cells_set_num_table.index.to_list()
        cells = cells_set_num_table['c1_c2'].to_list()

        self.set_num_by_cells = bidict(zip(cells,set_nums))
'''