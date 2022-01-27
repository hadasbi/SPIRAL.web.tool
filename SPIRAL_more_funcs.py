####################################################################################################################
def visualize_multiple_structures(data_name, cell_selection, ordered_pairs, impute_method, structs_to_vis,
                                  read_sigtable_from_file=False, sigtable=np.nan,
                                  read_scRNA_data_from_file=False, scRNA_data=np.nan,
                                  start_from_scratch=False, num_stds_thresh=1, mu=0.9,
                                  filt_mode='filt_by_3', filter_word='filtered',
                                  no_iterations=False):
    # filt_mode = 'filt_by_3' OR 'filt_by_structure_average_std' OR 'filt_by_rule'

    print(data_name)
    print('\nvisualize_multiple_structures!\n')

    if no_iterations:
        if ordered_pairs:
            sigfile_filt = ('../results_no_iterations/sigtables/num_stds_thresh_' + str(num_stds_thresh) + '/mu_' + str(
                mu) + '/ordered_pairs/' +
                            cell_selection + '/' + impute_method + '/' + data_name + '_' + filter_word + '_significance_table_' + filt_mode +
                            '.xlsx')
        else:
            sigfile_filt = ('../results_no_iterations/sigtables/num_stds_thresh_' + str(num_stds_thresh) + '/mu_' + str(
                mu) + '/unordered_pairs/' +
                            cell_selection + '/' + impute_method + '/' + data_name + '_' + filter_word + '_significance_table_' + filt_mode +
                            '.xlsx')
    else:
        if ordered_pairs:
            sigfile_filt = ('../results/sigtables/num_stds_thresh_' + str(num_stds_thresh) + '/mu_' + str(
                mu) + '/ordered_pairs/' +
                            cell_selection + '/' + impute_method + '/' + data_name + '_' + filter_word + '_significance_table_' + filt_mode +
                            '.xlsx')
        else:
            sigfile_filt = ('../results/sigtables/num_stds_thresh_' + str(num_stds_thresh) + '/mu_' + str(
                mu) + '/unordered_pairs/' +
                            cell_selection + '/' + impute_method + '/' + data_name + '_' + filter_word + '_significance_table_' + filt_mode +
                            '.xlsx')
    sigfile_GO = sigfile_filt[:-5] + '_GO' + sigfile_filt[-5:]
    print('\n\n\n', sigfile_filt)
    print('\n\n\n', sigfile_GO)

    if read_sigtable_from_file:
        sigtable = load_excel_with_openpyxl_and_convert_to_pd_DataFrame(sigfile_GO)
    elif sigtable == np.nan:
        print('ERROR: No sigtable')

    if len(sigtable) > 0:
        num_cells = int(sigtable.iloc[0, :]['num_cells'])
        num_genes = int(sigtable.iloc[0, :]['num_genes'])
        print('num_cells:', num_cells, ' num_genes:', num_genes)

        # Load genes_table
        with open('./objects/genes_tables/' + data_name + '_' + str(num_genes) + 'X' + str(
                num_cells) + '_' + filter_word + '.p', 'rb') as fp:
            genes_table = pickle.load(fp)

        # Read scRNA_data from file
        if read_scRNA_data_from_file:
            ind = [m.start() for m in re.finditer('/', data_dict[data_name])][-1]
            imputed_file = data_dict[data_name][:ind + 1] + impute_method + data_dict[data_name][
                                                                            ind:-4] + '_' + filter_word + '_imputed.txt'
            scRNA_data = load_data(imputed_file)
        elif scRNA_data == np.nan:
            print('ERROR: no scRNA_data matrix')

        # cell labels
        general_graphs = False
        if data_name == 'RinaBulkRNAData':
            cell_definitions = list(scRNA_data)
            print(cell_definitions)
            cell_days = np.array([int(d[3:d.find('Rep')]) for d in cell_definitions])
            cell_labels = np.zeros((len(list(scRNA_data)))).astype(int)
            for day in set(cell_days):
                cell_labels[list(np.where(cell_days == day)[0])] = day

        elif data_name == 'yaellab_covid19':
            cell_definitions = list(scRNA_data)
            print(cell_definitions)
            cell_labels = np.zeros((len(list(scRNA_data)))).astype(int)
            for p, patient in enumerate(
                    ['hc', 'cl', 'j', 'n984', 'n4', 'n995', 'n962', 'neg1', 'neg2', 'pos_l', 'pos_h']):
                cell_labels[[i for i, deff in enumerate(cell_definitions) if patient in deff]] = p

        elif data_name in ['Zhang2019_10x2', 'synthetic_splatter_branching_path', 'synthetic_splatter_branching_path2']:
            cell_definitions = [''] * num_cells

        elif data_name in ['Soumillon2014_D1',
                           'Wyler2020_Calu3_mock', 'Wyler2020_Calu3_s1', 'Wyler2020_Calu3_s2',
                           'Wyler2020_H1299_mock', 'Wyler2020_H1299_s1', 'Wyler2020_H1299_s2',
                           'Wagner2018', 'Wagner2018_',
                           'Sharon2019',
                           'Ramos2019_mock_raw', 'Ramos2019_MOI02_raw', 'Ramos2019_MOI20_raw']:
            if (impute_method == 'Iacono_clustering_dimcutoff_50') or (
                    impute_method == 'Iacono_clustering_dimcutoff_50_median'):
                clusters = pd.read_csv(
                    '../representing_clusters/Iacono_clusters_dimcutoff_50/' + data_name + '_' + filter_word + '_clusters.txt',
                    sep=' ', index_col=0)
            else:
                clusters = pd.read_csv(
                    '../representing_clusters/Iacono_clusters/' + data_name + '_' + filter_word + '_clusters.txt',
                    sep=' ', index_col=0)
            clusters.index = np.array(range(len(clusters)))
            clusters.columns = ['Iacono_cluster']
            # print('clusters:', clusters)
            # print('Averagely', np.round(filtered_num_cells/np.max(np.array(clusters)),2), 'cells per cluster')

            # load filtered scRNA matrix
            filtered_scRNA_data = load_data(data_dict[data_name][:-4] + '_' + filter_word + '.txt')
            filtered_num_cells = filtered_scRNA_data.shape[1]  # number of cells
            filtered_num_genes = filtered_scRNA_data.shape[0]  # number of genes
            # print('filtered:', filtered_num_genes, filtered_num_cells)

            # original labels
            orig_labels = list(filtered_scRNA_data)
            print(orig_labels[:5])
            if data_name == 'Soumillon2014_D1':
                labels = ['T0', 'T1', 'T2', 'T3', 'T5', 'T7', 'T9', 'T14']
                for l in labels:
                    cells = [i for i, s in enumerate(orig_labels) if l in s]
                    clusters.loc[cells, 'time_labels'] = l
            elif data_name in ['Wyler2020_Calu3_mock', 'Wyler2020_Calu3_s1', 'Wyler2020_Calu3_s2',
                               'Wyler2020_H1299_mock', 'Wyler2020_H1299_s1', 'Wyler2020_H1299_s2']:
                inds = [[m.start() for m in re.finditer('-', label)] for label in orig_labels]
                labels = list(set([label[ind[1] + 1:ind[2] + 2] for ind, label in zip(inds, orig_labels)]))
                for l in labels:
                    cells = [i for i, s in enumerate(orig_labels) if l in s]
                    clusters.loc[cells, 'time_labels'] = l
            elif data_name in ['Wagner2018', 'Wagner2018_']:
                clusters.loc[:, 'time_labels'] = [s[:s.find('_')] for s in orig_labels]
                labels = list(set(list(clusters.loc[:, 'time_labels'])))
            elif data_name == 'Sharon2019':
                clusters.loc[:, 'time_labels'] = [s[s.find('/') + 1:[i for i, ch in enumerate(s) if ch == '_'][2]].
                                                      replace('_', '.').replace('e', 'E') for s in
                                                  [s.replace('Sample_', '') for s in list(filtered_scRNA_data)]]
                labels = list(set(list(clusters.loc[:, 'time_labels'])))
            elif data_name in ['Ramos2019_mock_raw', 'Ramos2019_MOI02_raw', 'Ramos2019_MOI20_raw']:
                clusters.loc[:, 'time_labels'] = [s[:s.find('_')] for s in list(filtered_scRNA_data)]
                labels = list(set(list(clusters.loc[:, 'time_labels'])))

            clusters_table = pd.DataFrame(columns=labels, data=np.zeros((max(clusters['Iacono_cluster']), len(labels))),
                                          index=(np.array(range(max(clusters['Iacono_cluster']))) + 1))
            day_for_cluster_majority_vote = np.zeros((max(clusters['Iacono_cluster']))).astype(int)
            for c in (np.array(range(max(clusters['Iacono_cluster']))) + 1):
                # print(c)
                best_label_count = 0
                for i, l in enumerate(labels):
                    len_c_l = len(clusters[clusters['Iacono_cluster'] == c][clusters['time_labels'] == l])
                    clusters_table.loc[c, l] = len_c_l
                    if len_c_l > best_label_count:
                        day_for_cluster_majority_vote[c - 1] = i
                        best_label_count = len_c_l
            print('clusters_table:\n', clusters_table)
            cell_labels = day_for_cluster_majority_vote
            cell_definitions = [labels[i] for i in cell_labels]

        elif data_name in ['synthetic_splatter_branching_path3',
                           'synthetic_splatter_branching_path5',
                           'synthetic_splatter_branching_path6']:
            if (impute_method == 'Iacono_clustering_dimcutoff_50') or (
                    impute_method == 'Iacono_clustering_dimcutoff_50_median'):
                clusters = pd.read_csv(
                    '../representing_clusters/Iacono_clusters_dimcutoff_50/' + data_name + '_' + filter_word + '_clusters.txt',
                    sep=' ', index_col=0)
            else:
                clusters = pd.read_csv(
                    '../representing_clusters/Iacono_clusters/' + data_name + '_' + filter_word + '_clusters.txt',
                    sep=' ', index_col=0)
            clusters.index = np.array(range(len(clusters)))
            clusters.columns = ['Iacono_cluster']
            # print('clusters:', clusters)
            # print('Averagely', np.round(filtered_num_cells/np.max(np.array(clusters)),2), 'cells per cluster')

            # load filtered scRNA matrix
            filtered_scRNA_data = load_data(data_dict[data_name][:-4] + '_' + filter_word + '.txt')
            filtered_num_cells = filtered_scRNA_data.shape[1]  # number of cells
            filtered_num_genes = filtered_scRNA_data.shape[0]  # number of genes
            # print('filtered:', filtered_num_genes, filtered_num_cells)

            # load colData
            coldata = pd.read_table(data_dict[data_name][:-4] + '_colData.txt')
            coldata = coldata.loc[list(filtered_scRNA_data), :]

            # label based on group
            orig_labels = list(coldata['Group'])

            # original labels
            labels = list(set(orig_labels))
            for l in labels:
                cells = [i for i, s in enumerate(orig_labels) if l in s]
                clusters.loc[cells, 'time_labels'] = l

            clusters_table = pd.DataFrame(columns=labels, data=np.zeros((max(clusters['Iacono_cluster']), len(labels))),
                                          index=(np.array(range(max(clusters['Iacono_cluster']))) + 1))
            day_for_cluster_majority_vote = np.zeros((max(clusters['Iacono_cluster']))).astype(int)
            for c in (np.array(range(max(clusters['Iacono_cluster']))) + 1):
                # print(c)
                best_label_count = 0
                for i, l in enumerate(labels):
                    len_c_l = len(clusters[clusters['Iacono_cluster'] == c][clusters['time_labels'] == l])
                    clusters_table.loc[c, l] = len_c_l
                    if len_c_l > best_label_count:
                        day_for_cluster_majority_vote[c - 1] = i
                        best_label_count = len_c_l
            print('clusters_table:\n', clusters_table)
            cell_labels = day_for_cluster_majority_vote
            cell_definitions = [labels[i] for i in cell_labels]
            print('cell_labels:', cell_labels)
            print('cell_definitions:', cell_definitions)

        elif data_name == 'Han2018_H9':
            clusters = pd.read_csv(
                '../representing_clusters/Iacono_clusters/' + data_name + '_' + filter_word + '_clusters.txt', sep=' ',
                index_col=0)
            clusters.index = np.array(range(len(clusters)))
            clusters.columns = ['Iacono_cluster']
            # print('clusters:', clusters)
            # print('Averagely', np.round(filtered_num_cells/np.max(np.array(clusters)),2), 'cells per cluster')

            # load filtered scRNA matrix
            filtered_scRNA_data = load_data(data_dict[data_name][:-4] + '_' + filter_word + '.txt')
            filtered_num_cells = filtered_scRNA_data.shape[1]  # number of cells
            filtered_num_genes = filtered_scRNA_data.shape[0]  # number of genes
            # print('filtered:', filtered_num_genes, filtered_num_cells)

            # original labels
            orig_labels = list(filtered_scRNA_data)
            labels = ['day0', 'day10', 'day20']
            for l in labels:
                # print(l)
                cells = [i for i, s in enumerate(orig_labels) if l in s]
                clusters.loc[cells, 'day_cluster'] = l

            clusters_table = pd.DataFrame(columns=labels, data=np.zeros((max(clusters['Iacono_cluster']), len(labels))),
                                          index=(np.array(range(max(clusters['Iacono_cluster']))) + 1))
            day_for_cluster_majority_vote = []
            for c in (np.array(range(max(clusters['Iacono_cluster']))) + 1):
                # print(c)
                best_label_count = 0
                for l in labels:
                    len_c_l = len(clusters[clusters['Iacono_cluster'] == c][clusters['day_cluster'] == l])
                    clusters_table.loc[c, l] = len_c_l
                    if len_c_l > best_label_count:
                        day_for_cluster_majority_vote.append(l)
                        best_label_count = len_c_l
            print('clusters_table:\n', clusters_table)
            cell_labels = [int(l[3:]) for l in day_for_cluster_majority_vote]
            cell_definitions = day_for_cluster_majority_vote

        elif (data_name == 'BenMoshe2019_GSE122084_RAW_naive' and impute_method == 'Iacono_clustering'
              and filter_word == 'filtered'):
            cell_labels = [4, 4, 4, 4, 4, 5, 5, 5, 8, 8, 8, 8, 8, 0, 0, 0, 0, 0, 4, 0, 0, 0,
                           7, 1, 1, 1, 1, 2, 3, 3, 6, 3, 3, 3]
            cell_definitions = [{0: 'CD8 T cell', 1: 'monocyte', 2: 'dendritic cell', 3: 'B cell', 4: 'CD4 T cell',
                                 5: 'unknown', 6: 'dendritic cell', 7: 'monocyte', 8: 'NK cell'}[x] for x in
                                cell_labels]


        elif (data_name == 'BenMoshe2019_GSE122084_RAW_exposed1' and impute_method == 'Iacono_clustering'
              and filter_word == 'filtered'):
            cell_labels = [5, 3, 3, 7, 7, 7, 1, 7, 7, 3, 3, 1, 1, 1, 1, 1, 1, 1, 1, 2, 3, 3,
                           3, 3, 3, 6, 0, 7, 7, 4, 8, 7]
            cell_definitions = [{0: 'monocyte', 1: 'CD8 T cell', 2: 'NK or T cell', 3: 'NK cell', 4: 'B cell',
                                 5: 'CD4 T cell', 6: 'monocyte', 7: 'B cell', 8: 'monocyte or dendritic cell'
                                 }[x] for x in cell_labels]

        elif (data_name == 'BenMoshe2019_GSE122084_RAW_exposed2' and impute_method == 'Iacono_clustering'
              and filter_word == 'filtered'):
            cell_labels = [1, 1, 1, 1, 10, 5, 10, 11, 3, 3, 9, 9, 9, 9, 9, 9, 3,
                           3, 12, 6, 0, 4, 7, 8, 12, 7, 3, 7, 7, 2, 7]
            cell_definitions = [{0: 'B cell', 1: 'monocyte', 2: 'unknown', 3: 'CCR7+ T cell', 4: 'CD8 T cell',
                                 5: 'monocyte', 6: 'dendritic cell', 7: 'CD8 T cell', 8: 'unknown', 9: 'CCR7+ T cell',
                                 10: 'monocyte', 11: 'monocyte', 12: 'NK cell'}[x] for x in cell_labels]

        elif (data_name == 'BenMoshe2019_GSE122084_RAW_exposed2' and impute_method == 'Iacono_clustering_dimcutoff_50'
              and filter_word == 'filtered'):
            cell_labels = [4, 4, 4, 4, 14, 7, 7, 1, 14, 1, 9, 9, 9, 9, 11, 11, 11,
                           9, 9, 11, 11, 9, 9, 11, 9, 9, 9, 9, 9, 11, 11, 9, 10, 5,
                           13, 11, 11, 11, 3, 6, 12, 3, 8, 8, 8, 8, 0, 0, 8, 9, 8,
                           8, 2, 8]
            cell_definitions = [{0: 'NK cell', 1: 'monocyte', 2: 'unknown', 3: 'CD4 T cell', 4: 'monocyte', 5: 'B cell',
                                 6: 'NK cell', 7: 'monocyte', 8: 'NK cell', 9: 'T cell', 10: 'dendritic cell',
                                 11: 'B cell',
                                 12: 'CD8 T cell', 13: 'unknown', 14: 'monocyte'}[x] for x in cell_labels]
        else:
            cell_labels = np.array(range(len(list(scRNA_data)))).astype(np.int)
            general_graphs = True
        print('cell_labels:', cell_labels)

    print('structures to visualize:', structs_to_vis)
    for i in structs_to_vis:
        # For every structure, save a figure and a link to the figure
        print('\nstructure index:', i)
        sets_in_struct = sigtable.loc[i, 'sets']
        # print(sets_in_struct)
        sets_in_struct = [tuple(list(map(int, s.strip(')').strip('(').split(', ')))) for s in
                          sets_in_struct.strip(']').strip('[').split('), (')]
        # print('sets_in_struct:', sets_in_struct)

        # create a graph of the structure
        if ordered_pairs:
            G = nx.DiGraph()
        else:
            G = nx.Graph()

        # sets_in_struct = [(1,2),(1,4),(1,5),(1,8),(3,2),(3,4),(3,5),(3,8)]
        G.add_edges_from(sets_in_struct)
        print(G.nodes)

        plt.figure(figsize=(12, 12))

        if bipartite.is_bipartite(G):
            pos = nx.bipartite_layout(G, list(set([s[0] for s in sets_in_struct])))

        elif ordered_pairs:  # multipartite layout

            # trial 1
            '''
            # Assuming the graph is tripartite, find nodes_0,nodes_1,nodes_2 without knowing them (based on the edges)
            nodes_with_out = set([s[0] for s in G.edges])
            nodes_with_in = set([s[1] for s in G.edges])
            nodes_1 = nodes_with_out.intersection(nodes_with_in)
            nodes_0 = nodes_with_out-nodes_1
            nodes_2 = nodes_with_in-nodes_1
            print(nodes_0)
            print(nodes_1)
            print(nodes_2)
            
            # set the location of the nodes for each set
            pos = dict()
            pos.update( (n, (1, 2*i)) for i, n in enumerate(nodes_0) ) # put nodes from X at x=1
            pos.update( (n, (2, 2*i+1)) for i, n in enumerate(nodes_1) ) # put nodes from Y at x=2
            pos.update( (n, (3, 2*i)) for i, n in enumerate(nodes_2) ) # put nodes from Z at x=3
            '''

            # trial 2
            '''
            nodes_with_out = set([s[0] for s in G.edges])
            nodes_with_in = set([s[1] for s in G.edges])
            nodes_with_in_and_out = nodes_with_out.intersection(nodes_with_in)
            nodes_curr_level = nodes_with_out-nodes_with_in_and_out
            level = 1
            posed_nodes = []
            pos = dict()
            while nodes_curr_level:
                print('level:', level)
                print('nodes_curr_level:', nodes_curr_level)                
                pos.update( (n, (level, 2*j+np.mod(level+1,2))) for j,n in enumerate(nodes_curr_level) ) 
                posed_nodes += list(nodes_curr_level)
                nodes_curr_level = []
                for node in (nodes_with_in-set(posed_nodes)):
                    in_nodes = [s[0] for s in G.edges if s[1]==node]
                    if set(in_nodes).issubset(set(posed_nodes)): 
                        nodes_curr_level.append(node)
                level += 1
            '''

            # trial 3
            nodes_with_out = set([s[0] for s in G.edges])
            nodes_with_in = set([s[1] for s in G.edges])
            nodes_with_in_and_out = nodes_with_out.intersection(nodes_with_in)
            nodes_with_out_only = nodes_with_out - nodes_with_in_and_out
            list_of_nodes_in_levels = [set()] * 100
            posed_nodes = set()
            level = 1
            while (set(G.nodes) - set(posed_nodes)):
                list_of_nodes_in_levels[level] = set(
                    [s[1] for s in G.edges if s[0] in nodes_with_out_only.union(posed_nodes)]
                ) - set(
                    [s[1] for s in G.edges if s[0] not in nodes_with_out_only.union(posed_nodes)]
                ) - posed_nodes

                list_of_nodes_in_levels[level - 1] = list_of_nodes_in_levels[level - 1].union(
                    set([s[0] for s in G.edges if s[1] in list_of_nodes_in_levels[level]]) - posed_nodes)

                posed_nodes = posed_nodes.union(list_of_nodes_in_levels[level])
                posed_nodes = posed_nodes.union(list_of_nodes_in_levels[level - 1])
                level += 1

            pos = dict()
            level = 0
            s = list_of_nodes_in_levels[0]
            while s:
                print('level:', level)
                print('nodes in the current level:', s)
                pos.update((n, (level + 1, 2 * j + np.mod(level, 2))) for j, n in enumerate(s))
                level += 1
                s = list_of_nodes_in_levels[level]

                # pos = nx.layout.spring_layout(G, iterations=200)
            # pos = nx.spectral_layout(G)
            # pos = nx.spiral_layout(G)
            # pos = nx.random_layout(G)
            # pos = nx.kamada_kawai_layout(G)
            # pos = nx.circular_layout(G)
            # pos = nx.planar_layout(G)
            # pos = nx.fruchterman_reingold_layout(G)

        else:
            pos = nx.layout.spring_layout(G, iterations=200)

        if not general_graphs:
            # nx.draw(G, pos, node_color=cell_labels[G.nodes], node_size=800, cmap=plt.cm.Blues, with_labels=True, arrows=True,
            #       labels={list(G.nodes)[i]:'day'+str(cell_labels[G.nodes][i]) for i in range(len(G.nodes))})
            nx.draw(G, pos, node_color=np.array(cell_labels)[G.nodes], node_size=800, cmap=plt.cm.Blues,
                    with_labels=True,
                    arrows=True, labels={
                    list(G.nodes)[i]: str(np.array(cell_definitions)[G.nodes][i]) for i in range(len(G.nodes))})
        # elif data_name=='a': # no labels on nodes
        #    nx.draw(G, pos, node_size=800, cmap=plt.cm.Blues, with_labels=False, arrows=True)
        else:
            nx.draw(G, pos, node_color=np.ones((len(G.nodes))), node_size=800, cmap=plt.cm.Blues, with_labels=True,
                    arrows=True,
                    labels={list(G.nodes)[i]: str(cell_labels[G.nodes][i]) for i in range(len(G.nodes))})

            # save the figure to file
        if ordered_pairs:
            picfolder = ('../results/structure_visualizations/num_stds_thresh_' + str(num_stds_thresh) + '/mu_' + str(
                mu) +
                         '/ordered_pairs/' + cell_selection + '/' + impute_method + '/' + data_name + '_' + filter_word)
        else:
            picfolder = ('../results/structure_visualizations/num_stds_thresh_' + str(num_stds_thresh) + '/mu_' + str(
                mu) +
                         '/unordered_pairs/' + cell_selection + '/' + impute_method + '/' + data_name + '_' + filter_word)
        if not os.path.exists(picfolder):
            os.makedirs(picfolder)
        picfile = picfolder + '/struct_' + str(i) + '.jpg'
        plt.savefig(picfile)
        # plt.clf()
        plt.close()

        # add a hyperlink to the pic file in sigtable
        if ordered_pairs:
            picfolder = ('D:\\Nextcloud\\results\\structure_visualizations\\num_stds_thresh_' + str(
                num_stds_thresh) + '\\mu_' +
                         str(mu) + '\\ordered_pairs\\' + cell_selection + '\\' + impute_method + '\\' + data_name + '_' + filter_word)
        else:
            picfolder = ('D:\\Nextcloud\\results\\structure_visualizations\\num_stds_thresh_' + str(
                num_stds_thresh) + '\\mu_' +
                         str(mu) + '\\unordered_pairs\\' + cell_selection + '\\' + impute_method + '\\' + data_name + '_' + filter_word)
        picfile = os.path.join(picfolder, 'struct_' + str(i) + '.jpg')
        comp_path = picfile.replace('\\', '/')
        hyperlink_txt = '=hyperlink("' + comp_path + '","graph")'
        sigtable.loc[i, 'graphs'] = hyperlink_txt

        # sigtable.to_excel(sigfile_GO)
    sigtable.to_excel(sigfile_GO)


####################################################################################################################
def visualize_multiple_structures_in_one_figure(data_name, cell_selection, ordered_pairs, impute_method,
                                                num_stds_thresh_lst, mu_lst, structs_lst_lst,
                                                read_sigtable_from_file=False, sigtable=np.nan,
                                                read_scRNA_data_from_file=False, scRNA_data=np.nan,
                                                start_from_scratch=False,
                                                num_stds_thresh=1, mu=0.9, filt_mode='filt_by_3',
                                                filter_word='filtered',
                                                no_iterations=False,
                                                struct_thr=0.8, log10_pavg_of_genes_thr=-3, load_umap_from_file=True,
                                                with_legend=True):
    # filt_mode = 'filt_by_3' OR 'filt_by_structure_average_std' OR 'filt_by_rule'

    print(data_name)
    print('\nvisualize_multiple_structures_in_one_figure!\n')

    # umap coordinates
    if data_name != 'yaeldiffbulk2':
        shapes = ['8', 's', "^", '*', 'X', 'D', 'P', '1', '2', '3', '4']
    else:
        shapes = ["$" + str(i) + "$" for i in [0, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21]]
    ind = [m.start() for m in re.finditer('/', data_dict[data_name])][-1]
    umapcoorfile = data_dict[data_name][:ind + 1] + impute_method + data_dict[data_name][
                                                                    ind:-4] + '_' + filter_word + '_imputed_umap.csv'
    if load_umap_from_file and os.path.exists(umapcoorfile):
        print('\nloading UMAP...')
        umap_coor = loadtxt(umapcoorfile, delimiter=',')
    else:
        print('\ncomputing UMAP...')
        umap_coor = umap.UMAP().fit_transform(scRNA_data.transpose())
        savetxt(umapcoorfile, umap_coor, delimiter=',')

    # cell labels
    ind = [m.start() for m in re.finditer('/', data_dict[data_name])][-1]
    celllabelsfile = data_dict[data_name][:ind + 1] + impute_method + data_dict[data_name][
                                                                      ind:-4] + '_' + filter_word + '_imputed_cell_labels.txt'
    celldeffsfile = data_dict[data_name][:ind + 1] + impute_method + data_dict[data_name][
                                                                     ind:-4] + '_' + filter_word + '_imputed_cell_deffs.txt'

    if os.path.exists(celllabelsfile) and os.path.exists(celldeffsfile):
        print('\nloading cell labels and cell definitions from file...')
        cell_labels = []
        cell_definitions = []
        with open(celllabelsfile, "r") as f:
            for line in f:
                cell_labels.append(int(line.strip()))
        with open(celldeffsfile, "r") as f:
            for line in f:
                cell_definitions.append(str(line.strip()))
        general_graphs = False
    else:
        print('\ncomputing cell labels and cell definitions...')
        general_graphs = False
        if data_name == 'RinaBulkRNAData':
            cell_definitions = list(scRNA_data)
            print(cell_definitions)
            cell_days = np.array([int(d[3:d.find('Rep')]) for d in cell_definitions])
            cell_labels = np.zeros((len(list(scRNA_data)))).astype(int)
            for day in set(cell_days):
                cell_labels[list(np.where(cell_days == day)[0])] = day

        elif data_name == 'yaellab_covid19':
            cell_definitions = list(scRNA_data)
            print(cell_definitions)
            cell_labels = np.zeros((len(list(scRNA_data)))).astype(int)
            for p, patient in enumerate(
                    ['hc', 'cl', 'j', 'n984', 'n4', 'n995', 'n962', 'neg1', 'neg2', 'pos_l', 'pos_h']):
                cell_labels[[i for i, deff in enumerate(cell_definitions) if patient in deff]] = p

        elif data_name == 'yaeldiffbulk2':
            cell_definitions = [r[r.find('_') + 1:] for r in list(scRNA_data)]
            b = list(set([r for r in cell_definitions if len(r) == 4]))
            b.sort()
            c = list(set([r for r in cell_definitions if len(r) == 5]))
            c.sort()
            d = b + c
            e = dict(zip(d, list(range(len(d)))))
            cell_labels = [e[r] for r in cell_definitions]

        elif data_name in ['Zhang2019_10x2', 'synthetic_splatter_branching_path', 'synthetic_splatter_branching_path2']:
            cell_definitions = [''] * num_cells
            cell_labels = [0] * num_cells

        elif data_name in ['Soumillon2014_D1',
                           'Wyler2020_Calu3_mock', 'Wyler2020_Calu3_s1', 'Wyler2020_Calu3_s2',
                           'Wyler2020_H1299_mock', 'Wyler2020_H1299_s1', 'Wyler2020_H1299_s2',
                           'Wagner2018', 'Wagner2018_',
                           'Sharon2019',
                           'Ramos2019_mock_raw', 'Ramos2019_MOI02_raw', 'Ramos2019_MOI20_raw']:
            if (impute_method == 'Iacono_clustering_dimcutoff_50') or (
                    impute_method == 'Iacono_clustering_dimcutoff_50_median'):
                clusters = pd.read_csv(
                    '../representing_clusters/Iacono_clusters_dimcutoff_50/' + data_name + '_' + filter_word + '_clusters.txt',
                    sep=' ', index_col=0)
            else:
                clusters = pd.read_csv(
                    '../representing_clusters/Iacono_clusters/' + data_name + '_' + filter_word + '_clusters.txt',
                    sep=' ', index_col=0)
            clusters.index = np.array(range(len(clusters)))
            clusters.columns = ['Iacono_cluster']
            # print('clusters:', clusters)
            # print('Averagely', np.round(filtered_num_cells/np.max(np.array(clusters)),2), 'cells per cluster')

            # load filtered scRNA matrix
            filtered_scRNA_data = load_data(data_dict[data_name][:-4] + '_' + filter_word + '.txt')
            filtered_num_cells = filtered_scRNA_data.shape[1]  # number of cells
            filtered_num_genes = filtered_scRNA_data.shape[0]  # number of genes
            # print('filtered:', filtered_num_genes, filtered_num_cells)

            # clusters is a table of original cells vs. their Iacono cluster (1st column) and their time label (2nd column)

            # original labels
            orig_labels = list(filtered_scRNA_data)
            print(orig_labels[:5])
            if data_name == 'Soumillon2014_D1':
                labels = ['T0', 'T1', 'T2', 'T3', 'T5', 'T7', 'T9', 'T14']
                for l in labels:
                    cells = [i for i, s in enumerate(orig_labels) if l in s]
                    clusters.loc[cells, 'time_labels'] = l
            elif data_name in ['Wyler2020_Calu3_mock', 'Wyler2020_Calu3_s1', 'Wyler2020_Calu3_s2',
                               'Wyler2020_H1299_mock', 'Wyler2020_H1299_s1', 'Wyler2020_H1299_s2']:
                inds = [[m.start() for m in re.finditer('-', label)] for label in orig_labels]
                labels = list(set([label[ind[1] + 1:ind[2] + 2] for ind, label in zip(inds, orig_labels)]))
                for l in labels:
                    cells = [i for i, s in enumerate(orig_labels) if l in s]
                    clusters.loc[cells, 'time_labels'] = l
            elif data_name in ['Wagner2018', 'Wagner2018_']:
                clusters.loc[:, 'time_labels'] = [s[:s.find('_')] for s in orig_labels]
                labels = list(set(list(clusters.loc[:, 'time_labels'])))
            elif data_name == 'Sharon2019':
                clusters.loc[:, 'time_labels'] = [s[s.find('/') + 1:[i for i, ch in enumerate(s) if ch == '_'][2]].
                                                      replace('_', '.').replace('e', 'E') for s in
                                                  [s.replace('Sample_', '') for s in list(filtered_scRNA_data)]]
                labels = list(set(list(clusters.loc[:, 'time_labels'])))
            elif data_name in ['Ramos2019_mock_raw', 'Ramos2019_MOI02_raw', 'Ramos2019_MOI20_raw']:
                clusters.loc[:, 'time_labels'] = [s[:s.find('_')] for s in list(filtered_scRNA_data)]
                labels = list(set(list(clusters.loc[:, 'time_labels'])))

            # clusters_table is a table of Iacono clusters vs. how many of their original cells belong to each of
            # the time labels (columns)
            clusters_table = pd.DataFrame(columns=labels, data=np.zeros((max(clusters['Iacono_cluster']), len(labels))),
                                          index=(np.array(range(max(clusters['Iacono_cluster']))) + 1))
            day_for_cluster_majority_vote = np.zeros((max(clusters['Iacono_cluster']))).astype(int)
            for c in (np.array(range(max(clusters['Iacono_cluster']))) + 1):
                # print(c)
                best_label_count = 0
                for i, l in enumerate(labels):
                    len_c_l = len(clusters[clusters['Iacono_cluster'] == c][clusters['time_labels'] == l])
                    clusters_table.loc[c, l] = len_c_l
                    if len_c_l > best_label_count:
                        day_for_cluster_majority_vote[c - 1] = i
                        best_label_count = len_c_l
            print('clusters_table:\n', clusters_table)
            cell_labels = day_for_cluster_majority_vote
            cell_definitions = [labels[i] for i in cell_labels]

        elif data_name == 'Han2018_H9':
            clusters = pd.read_csv(
                '../representing_clusters/Iacono_clusters/' + data_name + '_' + filter_word + '_clusters.txt', sep=' ',
                index_col=0)
            clusters.index = np.array(range(len(clusters)))
            clusters.columns = ['Iacono_cluster']
            # print('clusters:', clusters)
            # print('Averagely', np.round(filtered_num_cells/np.max(np.array(clusters)),2), 'cells per cluster')

            # load filtered scRNA matrix
            filtered_scRNA_data = load_data(data_dict[data_name][:-4] + '_' + filter_word + '.txt')
            filtered_num_cells = filtered_scRNA_data.shape[1]  # number of cells
            filtered_num_genes = filtered_scRNA_data.shape[0]  # number of genes
            # print('filtered:', filtered_num_genes, filtered_num_cells)

            # original labels
            orig_labels = list(filtered_scRNA_data)
            labels = ['day0', 'day10', 'day20']
            for l in labels:
                # print(l)
                cells = [i for i, s in enumerate(orig_labels) if l in s]
                clusters.loc[cells, 'day_cluster'] = l

            clusters_table = pd.DataFrame(columns=labels, data=np.zeros((max(clusters['Iacono_cluster']), len(labels))),
                                          index=(np.array(range(max(clusters['Iacono_cluster']))) + 1))
            day_for_cluster_majority_vote = []
            for c in (np.array(range(max(clusters['Iacono_cluster']))) + 1):
                # print(c)
                best_label_count = 0
                for l in labels:
                    len_c_l = len(clusters[clusters['Iacono_cluster'] == c][clusters['day_cluster'] == l])
                    clusters_table.loc[c, l] = len_c_l
                    if len_c_l > best_label_count:
                        day_for_cluster_majority_vote.append(l)
                        best_label_count = len_c_l
            print('clusters_table:\n', clusters_table)
            cell_labels = [int(l[3:]) for l in day_for_cluster_majority_vote]
            cell_definitions = day_for_cluster_majority_vote

        elif (data_name == 'BenMoshe2019_GSE122084_RAW_naive' and impute_method == 'Iacono_clustering'
              and filter_word == 'filtered'):
            cell_labels = [4, 4, 4, 4, 4, 5, 5, 5, 8, 8, 8, 8, 8, 0, 0, 0, 0, 0, 4, 0, 0, 0,
                           7, 1, 1, 1, 1, 2, 3, 3, 6, 3, 3, 3]
            cell_definitions = [{0: 'CD8 T cell', 1: 'monocyte', 2: 'dendritic cell', 3: 'B cell', 4: 'CD4 T cell',
                                 5: 'unknown', 6: 'dendritic cell', 7: 'monocyte', 8: 'NK cell'}[x] for x in
                                cell_labels]


        elif (data_name == 'BenMoshe2019_GSE122084_RAW_exposed1' and impute_method == 'Iacono_clustering'
              and filter_word == 'filtered'):
            cell_labels = [5, 3, 3, 7, 7, 7, 1, 7, 7, 3, 3, 1, 1, 1, 1, 1, 1, 1, 1, 2, 3, 3,
                           3, 3, 3, 6, 0, 7, 7, 4, 8, 7]
            cell_definitions = [{0: 'monocyte', 1: 'CD8 T cell', 2: 'NK or T cell', 3: 'NK cell', 4: 'B cell',
                                 5: 'CD4 T cell', 6: 'monocyte', 7: 'B cell', 8: 'monocyte or dendritic cell'
                                 }[x] for x in cell_labels]

        elif (data_name == 'BenMoshe2019_GSE122084_RAW_exposed2' and impute_method == 'Iacono_clustering'
              and filter_word == 'filtered'):
            cell_labels = [1, 1, 1, 1, 10, 5, 10, 11, 3, 3, 9, 9, 9, 9, 9, 9, 3,
                           3, 12, 6, 0, 4, 7, 8, 12, 7, 3, 7, 7, 2, 7]
            cell_definitions = [{0: 'B cell', 1: 'monocyte', 2: 'unknown', 3: 'CCR7+ T cell', 4: 'CD8 T cell',
                                 5: 'monocyte', 6: 'dendritic cell', 7: 'CD8 T cell', 8: 'unknown', 9: 'CCR7+ T cell',
                                 10: 'monocyte', 11: 'monocyte', 12: 'NK cell'}[x] for x in cell_labels]

        elif (data_name == 'BenMoshe2019_GSE122084_RAW_exposed2' and impute_method == 'Iacono_clustering_dimcutoff_50'
              and filter_word == 'filtered'):
            cell_labels = [4, 4, 4, 4, 14, 7, 7, 1, 14, 1, 9, 9, 9, 9, 11, 11, 11,
                           9, 9, 11, 11, 9, 9, 11, 9, 9, 9, 9, 9, 11, 11, 9, 10, 5,
                           13, 11, 11, 11, 3, 6, 12, 3, 8, 8, 8, 8, 0, 0, 8, 9, 8,
                           8, 2, 8]
            cell_definitions = [{0: 'NK cell', 1: 'monocyte', 2: 'unknown', 3: 'CD4 T cell', 4: 'monocyte', 5: 'B cell',
                                 6: 'NK cell', 7: 'monocyte', 8: 'NK cell', 9: 'T cell', 10: 'dendritic cell',
                                 11: 'B cell',
                                 12: 'CD8 T cell', 13: 'unknown', 14: 'monocyte'}[x] for x in cell_labels]
        else:
            cell_labels = np.array(range(len(list(scRNA_data)))).astype(np.int)
            cell_definitions = [' '] * num_cells
            general_graphs = True

        print('cell_definitions:', cell_definitions)
        print('cell_labels:', cell_labels)

        # save cell_labels and cell_definitions to file
        with open(celllabelsfile, "w") as f:
            for s in cell_labels:
                f.write(str(s) + "\n")
        with open(celldeffsfile, "w") as f:
            for s in cell_definitions:
                f.write(str(s) + "\n")

    cell_definitions_sorted_set = []
    lens = list(set([len(r) for r in list(set(cell_definitions))]))
    for l in lens:
        b = list(set([r for r in cell_definitions if len(r) == l]))
        b.sort()
        cell_definitions_sorted_set += b
        ############################################################################################

        # pics folder
        if no_iterations:
            if ordered_pairs:
                picfolder = '../results_no_iterations/structure_visualizations/multiple_structures_in_one_figure/ordered_pairs/'
            else:
                picfolder = '../results_no_iterations/structure_visualizations/multiple_structures_in_one_figure/unordered_pairs/'
        else:
            if ordered_pairs:
                picfolder = '../results/structure_visualizations/multiple_structures_in_one_figure/ordered_pairs/'
            else:
                picfolder = '../results/structure_visualizations/multiple_structures_in_one_figure/unordered_pairs/'
        if not os.path.exists(picfolder):
            os.makedirs(picfolder)

            # picture file name
    picname = (data_name + '_' + filter_word + '_' + impute_method + '_' + str(num_stds_thresh_lst) + '_' + str(
        mu_lst) + '_structs_' +
               str(structs_lst_lst) + '_umap.jpg')

    ########################################################################################################
    fig, ax = plt.subplots(figsize=(15, 15))
    labels = []
    plt.rc('xtick', labelsize=14)
    plt.rc('ytick', labelsize=14)

    # redcolors = ['#FFA07A','#FA8072','#E9967A','#F08080','#CD5C5C','#DC143C','#B22222','#FF0000',
    #             '#8B0000','#800000','#FF6347','#FF4500','#DB7093']
    # redcolors = ['lightcoral','brown','red','tomato','lightsalmon','chocolate','rosybrown','sienna']
    # redcolors = ['c','b','m','r','k']
    redcolors = ['tab:' + c for c in ['blue', 'orange', 'red', 'purple', 'brown', 'pink', 'gray', 'cyan', 'green']] + [
        'magenta', 'teal', 'lightcoral']
    colorind = 0

    all_inds = set()

    for num_stds_thresh, mu, structs_lst in zip(num_stds_thresh_lst, mu_lst, structs_lst_lst):
        if no_iterations:
            if ordered_pairs:
                sigfile_filt = (
                        '../results_no_iterations/sigtables/num_stds_thresh_' + str(num_stds_thresh) + '/mu_' + str(
                    mu) + '/ordered_pairs/' +
                        cell_selection + '/' + impute_method + '/' + data_name + '_' + filter_word + '_significance_table_' + filt_mode +
                        '.xlsx')
            else:
                sigfile_filt = (
                        '../results_no_iterations/sigtables/num_stds_thresh_' + str(num_stds_thresh) + '/mu_' + str(
                    mu) + '/unordered_pairs/' +
                        cell_selection + '/' + impute_method + '/' + data_name + '_' + filter_word + '_significance_table_' + filt_mode +
                        '.xlsx')
        else:
            if ordered_pairs:
                sigfile_filt = ('../results/sigtables/num_stds_thresh_' + str(num_stds_thresh) + '/mu_' + str(
                    mu) + '/ordered_pairs/' +
                                cell_selection + '/' + impute_method + '/' + data_name + '_' + filter_word + '_significance_table_' + filt_mode +
                                '.xlsx')
            else:
                sigfile_filt = ('../results/sigtables/num_stds_thresh_' + str(num_stds_thresh) + '/mu_' + str(
                    mu) + '/unordered_pairs/' +
                                cell_selection + '/' + impute_method + '/' + data_name + '_' + filter_word + '_significance_table_' + filt_mode +
                                '.xlsx')
        sigfile_GO = sigfile_filt[:-5] + '_GO' + sigfile_filt[-5:]
        print('\n\n\n', sigfile_filt)
        print('\n\n\n', sigfile_GO)

        if read_sigtable_from_file:
            sigtable = load_excel_with_openpyxl_and_convert_to_pd_DataFrame(sigfile_GO)
        elif sigtable == np.nan:
            print('ERROR: No sigtable')

        if len(structs_lst) > 0:
            num_cells = int(sigtable.iloc[0, :]['num_cells'])
            num_genes = int(sigtable.iloc[0, :]['num_genes'])
            print('num_cells:', num_cells, ' num_genes:', num_genes)

            # Load genes_table
            with open('./objects/genes_tables/' + data_name + '_' + str(num_genes) + 'X' + str(
                    num_cells) + '_' + filter_word + '.p', 'rb') as fp:
                genes_table = pickle.load(fp)

            # Read scRNA_data from file
            if read_scRNA_data_from_file:
                ind = [m.start() for m in re.finditer('/', data_dict[data_name])][-1]
                imputed_file = data_dict[data_name][:ind + 1] + impute_method + data_dict[data_name][
                                                                                ind:-4] + '_' + filter_word + '_imputed.txt'
                scRNA_data = load_data(imputed_file)
            elif scRNA_data == np.nan:
                print('ERROR: no scRNA_data matrix')

        #########################################################################################

        for i in structs_lst:
            # For every structure, save a figure and a link to the figure
            print('\nstructure index:', i)
            sets_in_struct = sigtable.loc[i, 'sets']
            # print(sets_in_struct)
            sets_in_struct = [tuple(list(map(int, s.strip(')').strip('(').split(', ')))) for s in
                              sets_in_struct.strip(']').strip('[').split('), (')]
            # print('sets_in_struct:', sets_in_struct)

            #############################################################################
            # save merged UMAP of structures
            print('Computing and saving UMAP layout')

            nodes_with_out = set([s[0] for s in sets_in_struct])
            nodes_with_in = set([s[1] for s in sets_in_struct])
            nodes_with_in_and_out = nodes_with_out.intersection(nodes_with_in)
            nodes_curr_level = nodes_with_out - nodes_with_in_and_out
            level = 1
            node_color_dict = dict()
            posed_nodes = []
            while nodes_curr_level:
                # print('level:', level)
                # print('nodes_curr_level:', nodes_curr_level)
                for n in nodes_curr_level:
                    node_color_dict[n] = 'level ' + str(level)
                posed_nodes += list(nodes_curr_level)
                nodes_curr_level = []
                for node in (nodes_with_in - set(posed_nodes)):
                    in_nodes = [s[0] for s in sets_in_struct if s[1] == node]
                    if set(in_nodes).issubset(set(posed_nodes)):
                        nodes_curr_level.append(node)
                level += 1

            for n in range(scRNA_data.shape[1]):
                if n not in nodes_with_out.union(nodes_with_in):
                    node_color_dict[n] = 'not in structure'

            # print('node_color_dict:', node_color_dict)
            cell_levels = [node_color_dict[n] for n in range(scRNA_data.shape[1])]
            # print('cell_levels:', cell_levels)

            # umap.plot.points(mapper, labels=np.array(cell_levels), color_key_cmap='Paired')
            if 'not in structure' in set(cell_levels):
                nLevels = len(set(cell_levels)) - 1
            else:
                nLevels = len(set(cell_levels))

            if nLevels == 2:
                highlevels = ['level 2']
            elif nLevels == 3:
                highlevels = ['level 2', 'level 3']
            elif nLevels == 4:
                highlevels = ['level 3', 'level 4']
            elif nLevels == 5:
                highlevels = ['level 3', 'level 4', 'level 5']
            elif nLevels == 6:
                highlevels = ['level 4', 'level 5', 'level 6']

            inds_highlevels = [n for n, t in enumerate(cell_levels) if t in highlevels]
            for deff, shape in zip(cell_definitions_sorted_set, shapes[:len(set(cell_definitions))]):
                inds_deff = np.where(np.array(cell_definitions) == deff)[0]
                # print('inds_highlevels:', inds_highlevels)
                # print('inds_deff:', inds_deff)
                inds = list(set(inds_highlevels).intersection(set(inds_deff)))
                all_inds = all_inds.union(set(inds))
                if inds:
                    try:
                        l = deff + ' - ' + sigtable.loc[i, 'summary']
                    except:
                        l = deff
                    # print('umap_coor:', umap_coor)
                    # print('inds:', inds)
                    # print('color:', color)
                    # print('shape:', shape)
                    # print(umap_coor[inds, :])
                    ax.scatter(umap_coor[inds, 0], umap_coor[inds, 1], c=redcolors[colorind],
                               label=l, s=130, marker=shape)
                    labels.append(l)
            colorind += 1
    inds_not_in_struct = set(range(len(umap_coor))) - all_inds
    for deff, shape in zip(cell_definitions_sorted_set, shapes[:len(set(cell_definitions))]):
        inds_deff = np.where(np.array(cell_definitions) == deff)[0]
        inds = list(inds_not_in_struct.intersection(set(inds_deff)))
        if inds:
            l = deff + 'not in structure'
            ax.scatter(umap_coor[inds, 0], umap_coor[inds, 1], c='y',
                       label=l, s=130, marker=shape)
            labels.append(l)

    # title = ()
    # ax.set_title(title,fontdict={'fontsize': 22, 'fontweight': 'medium'})
    if with_legend:
        ax.legend(labels, prop={'size': 14})

    ax.set_xlabel('UMAP 1', fontdict={'fontsize': 16, 'fontweight': 'medium'})
    ax.set_ylabel('UMAP 2', fontdict={'fontsize': 16, 'fontweight': 'medium'})

    # save the figure to file
    picfile = picfolder + '/' + picname
    plt.savefig(picfile)
    # plt.clf()
    plt.close()
    # plt.show()


####################################################################################################################
def save_ranked_gene_lists(data_name, cell_selection, ordered_pairs, impute_method,
                           read_sigtable_from_file=False, sigtable=np.nan,
                           read_scRNA_data_from_file=False, scRNA_data=np.nan,
                           start_from_scratch=False,
                           num_stds_thresh=1, mu=0.9, filt_mode='filt_by_3',
                           filter_word='filtered', struct_thr=0.6, log10_pavg_of_genes_thr=-10,
                           no_iterations=False):
    # filt_mode = 'filt_by_3' OR 'filt_by_structure_average_std' OR 'filt_by_rule'

    print(data_name)
    print('\nsave_ranked_gene_lists!\n')

    if no_iterations:
        if ordered_pairs:
            sigfile_filt = ('../results_no_iterations/sigtables/num_stds_thresh_' + str(num_stds_thresh) + '/mu_' + str(
                mu) + '/ordered_pairs/' +
                            cell_selection + '/' + impute_method + '/' + data_name + '_' + filter_word + '_significance_table_' + filt_mode +
                            '.xlsx')
        else:
            sigfile_filt = ('../results_no_iterations/sigtables/num_stds_thresh_' + str(num_stds_thresh) + '/mu_' + str(
                mu) + '/unordered_pairs/' +
                            cell_selection + '/' + impute_method + '/' + data_name + '_' + filter_word + '_significance_table_' + filt_mode +
                            '.xlsx')
    else:
        if ordered_pairs:
            sigfile_filt = ('../results/sigtables/num_stds_thresh_' + str(num_stds_thresh) + '/mu_' + str(
                mu) + '/ordered_pairs/' +
                            cell_selection + '/' + impute_method + '/' + data_name + '_' + filter_word + '_significance_table_' + filt_mode +
                            '.xlsx')
        else:
            sigfile_filt = ('../results/sigtables/num_stds_thresh_' + str(num_stds_thresh) + '/mu_' + str(
                mu) + '/unordered_pairs/' +
                            cell_selection + '/' + impute_method + '/' + data_name + '_' + filter_word + '_significance_table_' + filt_mode +
                            '.xlsx')
    print('\n\n\n', sigfile_filt)

    sigtable = pd.read_excel(sigfile_filt, index_col=0)

    if len(sigtable) > 0:
        num_cells = int(sigtable.iloc[0, :]['num_cells'])
        num_genes = int(sigtable.iloc[0, :]['num_genes'])
        print('num_cells:', num_cells, ' num_genes:', num_genes)

        # Load genes_table
        with open('./objects/genes_tables/' + data_name + '_' + str(num_genes) + 'X' + str(
                num_cells) + '_' + filter_word + '.p', 'rb') as fp:
            genes_table = pickle.load(fp)

    if no_iterations:
        if ordered_pairs:
            structs_filename = (
                    './objects_no_iterations/struct_lists/num_stds_thresh_' + str(num_stds_thresh) + '/mu_' + str(
                mu) +
                    '/ordered_pairs/' + cell_selection + '/' + impute_method + '/' + data_name + '_' + filter_word +
                    '_relevant_structs.txt')
        else:
            structs_filename = (
                    './objects_no_iterations/struct_lists/num_stds_thresh_' + str(num_stds_thresh) + '/mu_' + str(
                mu) +
                    '/unordered_pairs/' + cell_selection + '/' + impute_method + '/' + data_name + '_' + filter_word +
                    '_relevant_structs.txt')
    else:
        if ordered_pairs:
            structs_filename = ('./objects/struct_lists/num_stds_thresh_' + str(num_stds_thresh) + '/mu_' + str(
                mu) + '/ordered_pairs/' +
                                cell_selection + '/' + impute_method + '/' + data_name + '_' + filter_word + '_relevant_structs.txt')
        else:
            structs_filename = ('./objects/struct_lists/num_stds_thresh_' + str(num_stds_thresh) + '/mu_' + str(
                mu) + '/unordered_pairs/' +
                                cell_selection + '/' + impute_method + '/' + data_name + '_' + filter_word + '_relevant_structs.txt')

    if ordered_pairs:
        structs_lst = list(sigtable.index[sigtable['structure_average_std'] <= struct_thr])
    else:
        structs_lst = list(sigtable.index[sigtable['log10_pavg_of_genes'] <= log10_pavg_of_genes_thr])

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

    if no_iterations:
        if ordered_pairs:
            ranked_gene_list_folder = (
                    '../results_no_iterations/ranked_gene_lists/num_stds_thresh_' + str(num_stds_thresh) +
                    '/mu_' + str(mu) + '/ordered_pairs/' + cell_selection + '/' + impute_method + '/' + data_name)
        else:
            ranked_gene_list_folder = (
                    '../results_no_iterations/ranked_gene_lists/num_stds_thresh_' + str(num_stds_thresh) +
                    '/mu_' + str(mu) + '/unordered_pairs/' + cell_selection + '/' + impute_method + '/' + data_name)
    else:
        if ordered_pairs:
            ranked_gene_list_folder = (
                    '../results/ranked_gene_lists/num_stds_thresh_' + str(num_stds_thresh) + '/mu_' + str(mu) +
                    '/ordered_pairs/' + cell_selection + '/' + impute_method + '/' + data_name)
        else:
            ranked_gene_list_folder = (
                    '../results/ranked_gene_lists/num_stds_thresh_' + str(num_stds_thresh) + '/mu_' + str(mu) +
                    '/unordered_pairs/' + cell_selection + '/' + impute_method + '/' + data_name)
    print('\n\n\n', ranked_gene_list_folder)
    if not os.path.exists(ranked_gene_list_folder):
        os.makedirs(ranked_gene_list_folder)

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

        genes_probs_logs_table1 = pd.DataFrame(data={'ID': genes_probs_logs_table.index,
                                                     't': np.flip(np.array(range(len(genes_probs_logs_table))))})
        print(genes_probs_logs_table1[:10])

        genes_probs_logs_table1.to_csv(ranked_gene_list_folder + '/struct_' + str(i) + '_ranked_genes.txt', sep='\t',
                                       index=False)


####################################################################################################################
def compute_structure_average_std(data_name, cell_selection, ordered_pairs, impute_method,
                                  read_sigtable_from_file=False, sigtable=np.nan,
                                  read_scRNA_data_from_file=False, scRNA_data=np.nan, filter_word='filtered'):
    print(data_name)
    print('\ncompute_structure_average_std!\n')

    # Read scRNA_data from file
    if read_scRNA_data_from_file:
        imputed_file = data_dict[data_name][:ind + 1] + impute_method + data_dict[data_name][
                                                                        ind:-4] + '_' + filter_word + '_imputed.txt'
        scRNA_data = load_data(imputed_file)
    elif scRNA_data == np.nan:
        print('ERROR: no scRNA_data matrix')

    normalized_scRNA_data = scRNA_data.subtract(scRNA_data.mean(axis=1), axis=0).divide(scRNA_data.std(axis=1), axis=0)
    normalized_scRNA_data.columns = [str(x) for x in normalized_scRNA_data.columns]
    print(normalized_scRNA_data.shape)
    # print('index:', list(normalized_scRNA_data.index))
    # print(normalized_scRNA_data.index[0], type(normalized_scRNA_data.index[0]))

    if ordered_pairs:
        sigfile = ('../results/sigtables/num_stds_thresh_' + str(num_stds_thresh) + '/mu_' + str(
            mu) + '/ordered_pairs/' +
                   cell_selection + '/' + impute_method + '/' + data_name + '_' + filter_word + '_significance_table.xlsx')
    else:
        sigfile = ('../results/sigtables/num_stds_thresh_' + str(num_stds_thresh) + '/mu_' + str(
            mu) + '/unordered_pairs/' +
                   cell_selection + '/' + impute_method + '/' + data_name + '_' + filter_word + '_significance_table.xlsx')

    if read_sigtable_from_file:
        sigtable = pd.read_excel(sigfile, index_col=0)
    elif sigtable == np.nan:
        print('ERROR: no sigtable')

    num_cells = int(sigtable.iloc[0, :]['num_cells'])
    num_genes = int(sigtable.iloc[0, :]['num_genes'])
    print('num_cells:', num_cells, ' num_genes:', num_genes)

    structs_loaded = False
    for i in list(sigtable.index):
        print('\n\n\nstructure', i)

        if len(sigtable.loc[i, 'genes']) == 32767:  # The genelist was too long to be written in one excel cell
            print('reading gene list from structures file')
            if structs_loaded:
                genes = list(genes_table[relevant_structs[i][0]].index)
            else:
                structs_filename = './objects/struct_lists/ordered_pairs/' + cell_selection + '/' + impute_method + '/' + data_name + '_relevant_structs.txt'
                relevant_structs = read_struct_list_from_file(structs_filename)
                genes = list(genes_table[relevant_structs[i][0]].index)
                structs_loaded = True
        else:  # read the gene list from the excel file
            print('reading gene list from excel')
            genes = read_gene_list(sigtable.loc[i, 'genes'])

        # print('\ngenes:', genes, '\n')
        # print('genes:', genes)
        # print(genes[0], type(genes[0]))
        print('number of genes in struct is', len(genes))
        non_id_genes = [g for g in genes if g not in list(normalized_scRNA_data.index)]
        print('number of un-identified genes is', len(non_id_genes))
        print(non_id_genes)

        sets_in_struct = sigtable.loc[i, 'sets']
        sets = [tuple(list(map(int, s.strip(')').strip('(').split(', ')))) for s in
                sets_in_struct.strip(']').strip('[').split('), (')]
        print('sets:', sets)

        # find all the cells that are relevant to this structure:
        cells = list(set([item for sublist in sets for item in sublist]))
        print(len(cells), 'cells')
        print('cells:', cells)
        print('cells:', [str(c) for c in cells])

        # construct the transition matrix B that is relevant to this structure:
        curr_B = np.zeros((len(sets), len(cells)))
        for k, sett in enumerate(sets):
            curr_B[k, cells.index(sett[0])] = -1
            curr_B[k, cells.index(sett[1])] = 1
        print('curr_B.shape:', curr_B.shape)

        # covariance matrix of curr_Y (curr_S is a sample of it)
        cov_mat = np.dot(curr_B, curr_B.transpose())
        print('cov_mat.shape:', cov_mat.shape)

        avg_vec = (1 / len(sets)) * np.ones((len(sets)))
        print('avg_vec.shape:', avg_vec.shape)
        average_diff_of_gene_std = np.sqrt(np.dot(np.dot(avg_vec.transpose(), cov_mat), avg_vec))
        print('average_diff_of_gene_std =', average_diff_of_gene_std)
        sigtable.loc[i, 'structure_average_std'] = average_diff_of_gene_std

        # pick the relevant slices of the normalized data matrix
        curr_normalized_scRNA_data = normalized_scRNA_data.loc[genes, :]
        # curr_normalized_scRNA_data = curr_normalized_scRNA_data.loc[:,[str(c) for c in cells]]
        curr_normalized_scRNA_data = curr_normalized_scRNA_data.iloc[:, [c for c in cells]]

        # multiply the normalized matrix by B to get the vector to compare with
        curr_S = np.dot(curr_B, curr_normalized_scRNA_data.transpose())
        print('curr_S.shape:', curr_S.shape)

        # average the differences for each gene
        average_diff_of_genes = np.dot(avg_vec.transpose(), curr_S)
        if len(average_diff_of_genes) != len(genes):
            print("PROBLEM 1")
        print('average_diff_of_genes[:5] =', average_diff_of_genes[:5])
        genes_probs_logs = norm.logsf(average_diff_of_genes / average_diff_of_gene_std) / np.log(10)
        if len(genes_probs_logs) != len(genes):
            print("PROBLEM 2")
        print('genes_probs_logs[:5] =', genes_probs_logs[:5])
        print('np.max(genes_probs_logs) =', np.max(genes_probs_logs))

        sigtable.loc[i, 'log10_pmax_of_genes'] = np.max(genes_probs_logs)
        sigtable.loc[i, 'log10_pavg_of_genes'] = np.log10(np.mean(np.power(10, genes_probs_logs)))
        print(np.log10(np.mean(np.power(10, genes_probs_logs))))
    sigtable.to_excel(sigfile)
    return sigtable, num_genes, num_cells


####################################################################################################################
def find_opposite_structures(data_name, cell_selection, ordered_pairs, impute_method, read_sigtable_from_file=False,
                             sigtable=np.nan, overlap_thresh=0.9, filt_mode='filt_by_3', filter_word='filtered'):
    print(data_name)
    print('\nfilter_similar_structures!\n')

    if ordered_pairs:
        sigfile_filt = ('../results/sigtables/num_stds_thresh_' + str(num_stds_thresh) + '/mu_' + str(
            mu) + '/ordered_pairs/' +
                        cell_selection + '/' + impute_method + '/' + data_name + '_' + filter_word + '_significance_table_' + filt_mode +
                        '.xlsx')
    else:
        sigfile_filt = ('../results/sigtables/num_stds_thresh_' + str(num_stds_thresh) + '/mu_' + str(
            mu) + '/unordered_pairs/' +
                        cell_selection + '/' + impute_method + '/' + data_name + '_' + filter_word + '_significance_table_' + filt_mode +
                        '.xlsx')
    sigfile_GO = sigfile_filt[:-5] + '_GO' + sigfile_filt[-5:]
    sigtable = load_excel_with_openpyxl_and_convert_to_pd_DataFrame(sigfile_GO)

    # Load genes_table
    with open('./objects/genes_tables/' + data_name + '_' + str(num_genes) + 'X' + str(
            num_cells) + '_' + filter_word + '.p', 'rb') as fp:
        genes_table = pickle.load(fp)

    num_cells = int(sigtable.iloc[0, :]['num_cells'])
    num_genes = int(sigtable.iloc[0, :]['num_genes'])
    print('num_cells:', num_cells, ' num_genes:', num_genes)

    opporiste_structs_table = pd.DataFrame(np.zeros((100, 24)),
                                           columns=[
                                               ['num_genes', 'num_cells', 'num_sets', 'struct1_ind', 'struct2_ind'] +
                                               ['struct1_' + s for s in ['num_rows_in_struct', 'num_cols_in_struct',
                                                                         'num_cells_in_struct',
                                                                         '#1a_log10_multiplication_of_Chebyshev_bound_of_all_genes',
                                                                         '#1b_log10_factorial_of_num_rows_in_struct',
                                                                         '#2_log10_init_prob',
                                                                         '#3_log10_init_prob',
                                                                         '#4_log10_init_prob',
                                                                         '#5_log10_init_prob',
                                                                         'log10_mul_test_corr',
                                                                         '#1_log10_final_prob',
                                                                         '#2_log10_final_prob',
                                                                         '#3_log10_final_prob',
                                                                         '#4_log10_final_prob',
                                                                         '#5_log10_final_prob',
                                                                         'log10_pmax_of_genes', 'log10_pavg_of_genes',
                                                                         'structure_average_std',
                                                                         'genes', 'sets', 'graphs']] +
                                               ['struct2_' + s for s in ['num_rows_in_struct', 'num_cols_in_struct',
                                                                         'num_cells_in_struct',
                                                                         '#1a_log10_multiplication_of_Chebyshev_bound_of_all_genes',
                                                                         '#1b_log10_factorial_of_num_rows_in_struct',
                                                                         '#2_log10_init_prob',
                                                                         '#3_log10_init_prob',
                                                                         '#4_log10_init_prob',
                                                                         '#5_log10_init_prob',
                                                                         'log10_mul_test_corr',
                                                                         '#1_log10_final_prob',
                                                                         '#2_log10_final_prob',
                                                                         '#3_log10_final_prob',
                                                                         '#4_log10_final_prob',
                                                                         '#5_log10_final_prob',
                                                                         'log10_pmax_of_genes', 'log10_pavg_of_genes',
                                                                         'structure_average_std',
                                                                         'genes', 'sets', 'graphs']]])

    opp_ind = 0
    structs_loaded = False
    struct_inds = list(sigtable.index)
    for k, i in enumerate(struct_inds):
        # print('\nk=',k,' i=',i)   #k is the serial number of i in inds
        for t, j in enumerate(struct_inds[k + 1:]):
            # print('t=',t,' t+k+1=',t+k+1,' j=', j)  #t+k+1 is the serial number of j in inds

            sets_i = sigtable.loc[i, 'sets']
            sets_i = [tuple(list(map(int, s.strip(')').strip('(').split(', ')))) for s in
                      sets_i.strip(']').strip('[').split('), (')]
            print('sets_i:', sets_i)

            sets_j = sigtable.loc[j, 'sets']
            sets_j = [tuple(list(map(int, s.strip(')').strip('(').split(', ')))) for s in
                      sets_j.strip(']').strip('[').split('), (')]
            print('sets_j:', sets_j)

            sets_j_opposite = [(s[1], s[0]) for s in sets_j]
            print('sets_j_opposite:', sets_j_opposite)

            # Jaccard
            if len(set(sets_i).intersection(set(sets_j_opposite))) / len(set(sets_i).union(set(sets_j_opposite))) > 0.9:
                opporiste_structs_table.loc[opp_ind, ['num_genes', 'num_cells', 'num_sets']] = sigtable.loc[
                    i, ['num_genes', 'num_cells', 'num_sets']]
                opporiste_structs_table.loc[opp_ind, 'struct1_ind'] = i
                opporiste_structs_table.loc[opp_ind, 'struct2_ind'] = j
                opporiste_structs_table.loc[
                    opp_ind, ['struct1_' + s for s in ['num_rows_in_struct', 'num_cols_in_struct',
                                                       'num_cells_in_struct',
                                                       '#1a_log10_multiplication_of_Chebyshev_bound_of_all_genes',
                                                       '#1b_log10_factorial_of_num_rows_in_struct',
                                                       '#2_log10_init_prob',
                                                       '#3_log10_init_prob',
                                                       '#4_log10_init_prob',
                                                       '#5_log10_init_prob',
                                                       'log10_mul_test_corr',
                                                       '#1_log10_final_prob',
                                                       '#2_log10_final_prob',
                                                       '#3_log10_final_prob',
                                                       '#4_log10_final_prob',
                                                       '#5_log10_final_prob',
                                                       'log10_pmax_of_genes', 'log10_pavg_of_genes',
                                                       'structure_average_std',
                                                       'genes', 'sets', 'graphs']]] = sigtable.loc[
                    i, ['num_rows_in_struct', 'num_cols_in_struct',
                        'num_cells_in_struct',
                        '#1a_log10_multiplication_of_Chebyshev_bound_of_all_genes',
                        '#1b_log10_factorial_of_num_rows_in_struct',
                        '#2_log10_init_prob',
                        '#3_log10_init_prob',
                        '#4_log10_init_prob',
                        '#5_log10_init_prob',
                        'log10_mul_test_corr',
                        '#1_log10_final_prob',
                        '#2_log10_final_prob',
                        '#3_log10_final_prob',
                        '#4_log10_final_prob',
                        '#5_log10_final_prob',
                        'log10_pmax_of_genes', 'log10_pavg_of_genes',
                        'structure_average_std',
                        'genes', 'sets', 'graphs']]
                opporiste_structs_table.loc[
                    opp_ind, ['struct2_' + s for s in ['num_rows_in_struct', 'num_cols_in_struct',
                                                       'num_cells_in_struct',
                                                       '#1a_log10_multiplication_of_Chebyshev_bound_of_all_genes',
                                                       '#1b_log10_factorial_of_num_rows_in_struct',
                                                       '#2_log10_init_prob',
                                                       '#3_log10_init_prob',
                                                       '#4_log10_init_prob',
                                                       '#5_log10_init_prob',
                                                       'log10_mul_test_corr',
                                                       '#1_log10_final_prob',
                                                       '#2_log10_final_prob',
                                                       '#3_log10_final_prob',
                                                       '#4_log10_final_prob',
                                                       '#5_log10_final_prob',
                                                       'log10_pmax_of_genes', 'log10_pavg_of_genes',
                                                       'structure_average_std',
                                                       'genes', 'sets', 'graphs']]] = sigtable.loc[
                    j, ['num_rows_in_struct', 'num_cols_in_struct',
                        'num_cells_in_struct',
                        '#1a_log10_multiplication_of_Chebyshev_bound_of_all_genes',
                        '#1b_log10_factorial_of_num_rows_in_struct',
                        '#2_log10_init_prob',
                        '#3_log10_init_prob',
                        '#4_log10_init_prob',
                        '#5_log10_init_prob',
                        'log10_mul_test_corr',
                        '#1_log10_final_prob',
                        '#2_log10_final_prob',
                        '#3_log10_final_prob',
                        '#4_log10_final_prob',
                        '#5_log10_final_prob',
                        'log10_pmax_of_genes', 'log10_pavg_of_genes',
                        'structure_average_std',
                        'genes', 'sets', 'graphs']]
                opp_ind += 1

    oppfile = sigfile[:-5] + '_opposite_structures' + sigfile[-5:]
    opporiste_structs_table.to_excel(oppfile)
    return opporiste_structs_table


####################################################################################################################
def visualize_multiple_structures_with_arrows(data_name, cell_selection, ordered_pairs, impute_method, structs_lst,
                                              read_sigtable_from_file=False, sigtable=np.nan,
                                              read_scRNA_data_from_file=False, scRNA_data=np.nan,
                                              num_stds_thresh=1, mu=0.9, filt_mode='filt_by_3', filter_word='filtered',
                                              no_iterations=False, mode='PCA', load_coor_from_file=True,
                                              with_legend=True, red_only=False):
    # filt_mode = 'filt_by_3' OR 'filt_by_structure_average_std' OR 'filt_by_rule'

    print(data_name)
    print('\nvisualize_multiple_structures_with_arrows!\n')

    if no_iterations:
        if ordered_pairs:
            sigfile_filt = ('../results_no_iterations/sigtables/num_stds_thresh_' + str(num_stds_thresh) + '/mu_' + str(
                mu) + '/ordered_pairs/' +
                            cell_selection + '/' + impute_method + '/' + data_name + '_' + filter_word + '_significance_table_' + filt_mode +
                            '.xlsx')
        else:
            sigfile_filt = ('../results_no_iterations/sigtables/num_stds_thresh_' + str(num_stds_thresh) + '/mu_' + str(
                mu) + '/unordered_pairs/' +
                            cell_selection + '/' + impute_method + '/' + data_name + '_' + filter_word + '_significance_table_' + filt_mode +
                            '.xlsx')
    else:
        if ordered_pairs:
            sigfile_filt = ('../results/sigtables/num_stds_thresh_' + str(num_stds_thresh) + '/mu_' + str(
                mu) + '/ordered_pairs/' +
                            cell_selection + '/' + impute_method + '/' + data_name + '_' + filter_word + '_significance_table_' + filt_mode +
                            '.xlsx')
        else:
            sigfile_filt = ('../results/sigtables/num_stds_thresh_' + str(num_stds_thresh) + '/mu_' + str(
                mu) + '/unordered_pairs/' +
                            cell_selection + '/' + impute_method + '/' + data_name + '_' + filter_word + '_significance_table_' + filt_mode +
                            '.xlsx')
    sigfile_GO = sigfile_filt[:-5] + '_GO' + sigfile_filt[-5:]
    print('\n\n\n', sigfile_filt)
    print('\n\n\n', sigfile_GO)

    if read_sigtable_from_file:
        sigtable = load_excel_with_openpyxl_and_convert_to_pd_DataFrame(sigfile_GO)
    elif sigtable == np.nan:
        print('ERROR: No sigtable')

    num_cells = int(sigtable.iloc[0, :]['num_cells'])
    num_genes = int(sigtable.iloc[0, :]['num_genes'])
    print('num_cells:', num_cells, ' num_genes:', num_genes)

    # Load genes_table
    with open('./objects/genes_tables/' + data_name + '_' + str(num_genes) + 'X' + str(
            num_cells) + '_' + filter_word + '.p', 'rb') as fp:
        genes_table = pickle.load(fp)

    # Read scRNA_data from file
    if read_scRNA_data_from_file:
        ind = [m.start() for m in re.finditer('/', data_dict[data_name])][-1]
        imputed_file = data_dict[data_name][:ind + 1] + impute_method + data_dict[data_name][
                                                                        ind:-4] + '_' + filter_word + '_imputed.txt'
        scRNA_data = load_data(imputed_file)
    elif scRNA_data == np.nan:
        print('ERROR: no scRNA_data matrix')

    # umap and PCA shapes and sample names
    if data_name != 'yaeldiffbulk2':
        shapes = ['8', 's', "^", '*', 'X', 'D', 'P', '1', '2', '3', '4']
        samps = 'repcells'
    else:
        shapes = ["$" + str(i) + "$" for i in [0, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21]]
        samps = 'samples'

    # umap and PCA coordinates
    if mode == 'umap':
        ind = [m.start() for m in re.finditer('/', data_dict[data_name])][-1]
        umapcoorfile = data_dict[data_name][:ind + 1] + impute_method + data_dict[data_name][
                                                                        ind:-4] + '_' + filter_word + '_imputed_umap.csv'
        if load_coor_from_file and os.path.exists(umapcoorfile):
            print('\nloading UMAP...')
            coors = loadtxt(umapcoorfile, delimiter=',')
        else:
            print('\ncomputing UMAP...')
            coors = umap.UMAP().fit_transform(scRNA_data.transpose())
            savetxt(umapcoorfile, coors, delimiter=',')
    elif mode == 'PCA':
        ind = [m.start() for m in re.finditer('/', data_dict[data_name])][-1]
        pcacoorfile = data_dict[data_name][:ind + 1] + impute_method + data_dict[data_name][
                                                                       ind:-4] + '_' + filter_word + '_imputed_pca.csv'
        if load_coor_from_file and os.path.exists(pcacoorfile):
            print('\nloading PCA...')
            coors = loadtxt(pcacoorfile, delimiter=',')
        else:
            print('\ncomputing PCA...')
            pca = PCA(n_components=2)
            datamatrix = normalize(scRNA_data.transpose())
            coors = pca.fit_transform(datamatrix)
            savetxt(pcacoorfile, coors, delimiter=',')

            # cell labels
    ind = [m.start() for m in re.finditer('/', data_dict[data_name])][-1]
    celllabelsfile = data_dict[data_name][:ind + 1] + impute_method + data_dict[data_name][
                                                                      ind:-4] + '_' + filter_word + '_imputed_cell_labels.txt'
    celldeffsfile = data_dict[data_name][:ind + 1] + impute_method + data_dict[data_name][
                                                                     ind:-4] + '_' + filter_word + '_imputed_cell_deffs.txt'

    if os.path.exists(celllabelsfile) and os.path.exists(celldeffsfile):
        print('\nloading cell labels and cell definitions from file...')
        cell_labels = []
        cell_definitions = []
        with open(celllabelsfile, "r") as f:
            for line in f:
                cell_labels.append(int(line.strip()))
        with open(celldeffsfile, "r") as f:
            for line in f:
                cell_definitions.append(str(line.strip()))
        general_graphs = False
    else:
        print('\ncomputing cell labels and cell definitions...')
        general_graphs = False
        if data_name == 'RinaBulkRNAData':
            cell_definitions = list(scRNA_data)
            print(cell_definitions)
            cell_days = np.array([int(d[3:d.find('Rep')]) for d in cell_definitions])
            cell_labels = np.zeros((len(list(scRNA_data)))).astype(int)
            for day in set(cell_days):
                cell_labels[list(np.where(cell_days == day)[0])] = day

        elif data_name == 'yaellab_covid19':
            cell_definitions = list(scRNA_data)
            print(cell_definitions)
            cell_labels = np.zeros((len(list(scRNA_data)))).astype(int)
            for p, patient in enumerate(
                    ['hc', 'cl', 'j', 'n984', 'n4', 'n995', 'n962', 'neg1', 'neg2', 'pos_l', 'pos_h']):
                cell_labels[[i for i, deff in enumerate(cell_definitions) if patient in deff]] = p

        elif data_name == 'yaeldiffbulk2':
            cell_definitions = [r[r.find('_') + 1:] for r in list(scRNA_data)]
            b = list(set([r for r in cell_definitions if len(r) == 4]))
            b.sort()
            c = list(set([r for r in cell_definitions if len(r) == 5]))
            c.sort()
            d = b + c
            e = dict(zip(d, list(range(len(d)))))
            cell_labels = [e[r] for r in cell_definitions]

        elif data_name in ['Zhang2019_10x2', 'synthetic_splatter_branching_path', 'synthetic_splatter_branching_path2']:
            cell_definitions = [''] * num_cells
            cell_labels = [0] * num_cells

        elif data_name in ['Soumillon2014_D1',
                           'Wyler2020_Calu3_mock', 'Wyler2020_Calu3_s1', 'Wyler2020_Calu3_s2',
                           'Wyler2020_H1299_mock', 'Wyler2020_H1299_s1', 'Wyler2020_H1299_s2',
                           'Wagner2018', 'Wagner2018_',
                           'Sharon2019',
                           'Ramos2019_mock_raw', 'Ramos2019_MOI02_raw', 'Ramos2019_MOI20_raw']:
            if (impute_method == 'Iacono_clustering_dimcutoff_50') or (
                    impute_method == 'Iacono_clustering_dimcutoff_50_median'):
                clusters = pd.read_csv(
                    '../representing_clusters/Iacono_clusters_dimcutoff_50/' + data_name + '_' + filter_word + '_clusters.txt',
                    sep=' ', index_col=0)
            else:
                clusters = pd.read_csv(
                    '../representing_clusters/Iacono_clusters/' + data_name + '_' + filter_word + '_clusters.txt',
                    sep=' ', index_col=0)
            clusters.index = np.array(range(len(clusters)))
            clusters.columns = ['Iacono_cluster']
            # print('clusters:', clusters)
            # print('Averagely', np.round(filtered_num_cells/np.max(np.array(clusters)),2), 'cells per cluster')

            # load filtered scRNA matrix
            filtered_scRNA_data = load_data(data_dict[data_name][:-4] + '_' + filter_word + '.txt')
            filtered_num_cells = filtered_scRNA_data.shape[1]  # number of cells
            filtered_num_genes = filtered_scRNA_data.shape[0]  # number of genes
            # print('filtered:', filtered_num_genes, filtered_num_cells)

            # clusters is a table of original cells vs. their Iacono cluster (1st column) and their time label (2nd column)

            # original labels
            orig_labels = list(filtered_scRNA_data)
            print(orig_labels[:5])
            if data_name == 'Soumillon2014_D1':
                labels = ['T0', 'T1', 'T2', 'T3', 'T5', 'T7', 'T9', 'T14']
                for l in labels:
                    cells = [i for i, s in enumerate(orig_labels) if l in s]
                    clusters.loc[cells, 'time_labels'] = l
            elif data_name in ['Wyler2020_Calu3_mock', 'Wyler2020_Calu3_s1', 'Wyler2020_Calu3_s2',
                               'Wyler2020_H1299_mock', 'Wyler2020_H1299_s1', 'Wyler2020_H1299_s2']:
                inds = [[m.start() for m in re.finditer('-', label)] for label in orig_labels]
                labels = list(set([label[ind[1] + 1:ind[2] + 2] for ind, label in zip(inds, orig_labels)]))
                for l in labels:
                    cells = [i for i, s in enumerate(orig_labels) if l in s]
                    clusters.loc[cells, 'time_labels'] = l
            elif data_name in ['Wagner2018', 'Wagner2018_']:
                clusters.loc[:, 'time_labels'] = [s[:s.find('_')] for s in orig_labels]
                labels = list(set(list(clusters.loc[:, 'time_labels'])))
            elif data_name == 'Sharon2019':
                clusters.loc[:, 'time_labels'] = [s[s.find('/') + 1:[i for i, ch in enumerate(s) if ch == '_'][2]].
                                                      replace('_', '.').replace('e', 'E') for s in
                                                  [s.replace('Sample_', '') for s in list(filtered_scRNA_data)]]
                labels = list(set(list(clusters.loc[:, 'time_labels'])))
            elif data_name in ['Ramos2019_mock_raw', 'Ramos2019_MOI02_raw', 'Ramos2019_MOI20_raw']:
                clusters.loc[:, 'time_labels'] = [s[:s.find('_')] for s in list(filtered_scRNA_data)]
                labels = list(set(list(clusters.loc[:, 'time_labels'])))

            # clusters_table is a table of Iacono clusters vs. how many of their original cells belong to each of
            # the time labels (columns)
            clusters_table = pd.DataFrame(columns=labels, data=np.zeros((max(clusters['Iacono_cluster']), len(labels))),
                                          index=(np.array(range(max(clusters['Iacono_cluster']))) + 1))
            day_for_cluster_majority_vote = np.zeros((max(clusters['Iacono_cluster']))).astype(int)
            for c in (np.array(range(max(clusters['Iacono_cluster']))) + 1):
                # print(c)
                best_label_count = 0
                for i, l in enumerate(labels):
                    len_c_l = len(clusters[clusters['Iacono_cluster'] == c][clusters['time_labels'] == l])
                    clusters_table.loc[c, l] = len_c_l
                    if len_c_l > best_label_count:
                        day_for_cluster_majority_vote[c - 1] = i
                        best_label_count = len_c_l
            print('clusters_table:\n', clusters_table)
            cell_labels = day_for_cluster_majority_vote
            cell_definitions = [labels[i] for i in cell_labels]

        elif data_name == 'Han2018_H9':
            clusters = pd.read_csv(
                '../representing_clusters/Iacono_clusters/' + data_name + '_' + filter_word + '_clusters.txt', sep=' ',
                index_col=0)
            clusters.index = np.array(range(len(clusters)))
            clusters.columns = ['Iacono_cluster']
            # print('clusters:', clusters)
            # print('Averagely', np.round(filtered_num_cells/np.max(np.array(clusters)),2), 'cells per cluster')

            # load filtered scRNA matrix
            filtered_scRNA_data = load_data(data_dict[data_name][:-4] + '_' + filter_word + '.txt')
            filtered_num_cells = filtered_scRNA_data.shape[1]  # number of cells
            filtered_num_genes = filtered_scRNA_data.shape[0]  # number of genes
            # print('filtered:', filtered_num_genes, filtered_num_cells)

            # original labels
            orig_labels = list(filtered_scRNA_data)
            labels = ['day0', 'day10', 'day20']
            for l in labels:
                # print(l)
                cells = [i for i, s in enumerate(orig_labels) if l in s]
                clusters.loc[cells, 'day_cluster'] = l

            clusters_table = pd.DataFrame(columns=labels, data=np.zeros((max(clusters['Iacono_cluster']), len(labels))),
                                          index=(np.array(range(max(clusters['Iacono_cluster']))) + 1))
            day_for_cluster_majority_vote = []
            for c in (np.array(range(max(clusters['Iacono_cluster']))) + 1):
                # print(c)
                best_label_count = 0
                for l in labels:
                    len_c_l = len(clusters[clusters['Iacono_cluster'] == c][clusters['day_cluster'] == l])
                    clusters_table.loc[c, l] = len_c_l
                    if len_c_l > best_label_count:
                        day_for_cluster_majority_vote.append(l)
                        best_label_count = len_c_l
            print('clusters_table:\n', clusters_table)
            cell_labels = [int(l[3:]) for l in day_for_cluster_majority_vote]
            cell_definitions = day_for_cluster_majority_vote

        elif (data_name == 'BenMoshe2019_GSE122084_RAW_naive' and impute_method == 'Iacono_clustering'
              and filter_word == 'filtered'):
            cell_labels = [4, 4, 4, 4, 4, 5, 5, 5, 8, 8, 8, 8, 8, 0, 0, 0, 0, 0, 4, 0, 0, 0,
                           7, 1, 1, 1, 1, 2, 3, 3, 6, 3, 3, 3]
            cell_definitions = [{0: 'CD8 T cell', 1: 'monocyte', 2: 'dendritic cell', 3: 'B cell', 4: 'CD4 T cell',
                                 5: 'unknown', 6: 'dendritic cell', 7: 'monocyte', 8: 'NK cell'}[x] for x in
                                cell_labels]


        elif (data_name == 'BenMoshe2019_GSE122084_RAW_exposed1' and impute_method == 'Iacono_clustering'
              and filter_word == 'filtered'):
            cell_labels = [5, 3, 3, 7, 7, 7, 1, 7, 7, 3, 3, 1, 1, 1, 1, 1, 1, 1, 1, 2, 3, 3,
                           3, 3, 3, 6, 0, 7, 7, 4, 8, 7]
            cell_definitions = [{0: 'monocyte', 1: 'CD8 T cell', 2: 'NK or T cell', 3: 'NK cell', 4: 'B cell',
                                 5: 'CD4 T cell', 6: 'monocyte', 7: 'B cell', 8: 'monocyte or dendritic cell'
                                 }[x] for x in cell_labels]

        elif (data_name == 'BenMoshe2019_GSE122084_RAW_exposed2' and impute_method == 'Iacono_clustering'
              and filter_word == 'filtered'):
            cell_labels = [1, 1, 1, 1, 10, 5, 10, 11, 3, 3, 9, 9, 9, 9, 9, 9, 3,
                           3, 12, 6, 0, 4, 7, 8, 12, 7, 3, 7, 7, 2, 7]
            cell_definitions = [{0: 'B cell', 1: 'monocyte', 2: 'unknown', 3: 'CCR7+ T cell', 4: 'CD8 T cell',
                                 5: 'monocyte', 6: 'dendritic cell', 7: 'CD8 T cell', 8: 'unknown', 9: 'CCR7+ T cell',
                                 10: 'monocyte', 11: 'monocyte', 12: 'NK cell'}[x] for x in cell_labels]

        elif (data_name == 'BenMoshe2019_GSE122084_RAW_exposed2' and impute_method == 'Iacono_clustering_dimcutoff_50'
              and filter_word == 'filtered'):
            cell_labels = [4, 4, 4, 4, 14, 7, 7, 1, 14, 1, 9, 9, 9, 9, 11, 11, 11,
                           9, 9, 11, 11, 9, 9, 11, 9, 9, 9, 9, 9, 11, 11, 9, 10, 5,
                           13, 11, 11, 11, 3, 6, 12, 3, 8, 8, 8, 8, 0, 0, 8, 9, 8,
                           8, 2, 8]
            cell_definitions = [{0: 'NK cell', 1: 'monocyte', 2: 'unknown', 3: 'CD4 T cell', 4: 'monocyte', 5: 'B cell',
                                 6: 'NK cell', 7: 'monocyte', 8: 'NK cell', 9: 'T cell', 10: 'dendritic cell',
                                 11: 'B cell',
                                 12: 'CD8 T cell', 13: 'unknown', 14: 'monocyte'}[x] for x in cell_labels]
        else:
            cell_labels = np.array(range(len(list(scRNA_data)))).astype(np.int)
            cell_definitions = [' '] * num_cells
            general_graphs = True

        print('cell_definitions:', cell_definitions)
        print('cell_labels:', cell_labels)

        # save cell_labels and cell_definitions to file
        with open(celllabelsfile, "w") as f:
            for s in cell_labels:
                f.write(str(s) + "\n")
        with open(celldeffsfile, "w") as f:
            for s in cell_definitions:
                f.write(str(s) + "\n")

                # sort cell_definitions
    cell_definitions_sorted_set = []
    lens = list(set([len(r) for r in list(set(cell_definitions))]))
    for l in lens:
        b = list(set([r for r in cell_definitions if len(r) == l]))
        b.sort()
        cell_definitions_sorted_set += b
    print('cell_definitions_sorted_set:', cell_definitions_sorted_set)
    ############################################################################################
    # pics folder
    if no_iterations:
        if ordered_pairs:
            picfolder = ('../results_no_iterations/structure_visualizations/num_stds_thresh_' + str(
                num_stds_thresh) + '/mu_' + str(mu) +
                         '/ordered_pairs/' + cell_selection + '/' + impute_method + '/' + data_name + '_' + filter_word)
            picfolderlink = ('D:\\Nextcloud\\results_no_iterations\\structure_visualizations\\num_stds_thresh_' + str(
                num_stds_thresh) + '\\mu_' +
                             str(mu) + '\\ordered_pairs\\' + cell_selection + '\\' + impute_method + '\\' + data_name + '_' + filter_word)
        else:
            picfolder = ('../results_no_iterations/structure_visualizations/num_stds_thresh_' + str(
                num_stds_thresh) + '/mu_' + str(mu) +
                         '/unordered_pairs\\' + cell_selection + '/' + impute_method + '/' + data_name + '_' + filter_word)
            picfolderlink = ('D:\\Nextcloud\\results_no_iterations\\structure_visualizations\\num_stds_thresh_' + str(
                num_stds_thresh) + '\\mu_' +
                             str(mu) + '\\unordered_pairs\\' + cell_selection + '\\' + impute_method + '\\' + data_name + '_' + filter_word)
    else:
        if ordered_pairs:
            picfolder = ('../results/structure_visualizations/num_stds_thresh_' + str(num_stds_thresh) + '/mu_' + str(
                mu) +
                         '/ordered_pairs/' + cell_selection + '/' + impute_method + '/' + data_name + '_' + filter_word)
            picfolderlink = ('D:\\Nextcloud\\results\\structure_visualizations\\num_stds_thresh_' + str(
                num_stds_thresh) + '\\mu_' +
                             str(mu) + '\\ordered_pairs\\' + cell_selection + '\\' + impute_method + '\\' + data_name + '_' + filter_word)
        else:
            picfolder = ('../results/structure_visualizations/num_stds_thresh_' + str(num_stds_thresh) + '/mu_' + str(
                mu) +
                         '/unordered_pairs\\' + cell_selection + '/' + impute_method + '/' + data_name + '_' + filter_word)
            picfolderlink = ('D:\\Nextcloud\\results\\structure_visualizations\\num_stds_thresh_' + str(
                num_stds_thresh) + '\\mu_' +
                             str(mu) + '\\unordered_pairs\\' + cell_selection + '\\' + impute_method + '\\' + data_name + '_' + filter_word)
    if not os.path.exists(picfolder):
        os.makedirs(picfolder)

        ############################################################################################
    class_array = np.ones((num_cells, len(structs_lst))) * 1000  # 1000 means that the cell is not in the structure

    ############################################################################################
    for j, i in enumerate(structs_lst):
        # For every structure, save a figure and a link to the figure
        print('\nstructure index:', i)
        sets_in_struct = sigtable.loc[i, 'sets']
        # print(sets_in_struct)
        sets_in_struct = [tuple(list(map(int, s.strip(')').strip('(').split(', ')))) for s in
                          sets_in_struct.strip(']').strip('[').split('), (')]
        # print('sets_in_struct:', sets_in_struct)

        #############################################################################
        nodes_with_out = list(set([s[0] for s in sets_in_struct]))
        nodes_with_in = list(set([s[1] for s in sets_in_struct]))
        class_array[nodes_with_out, j] = 0
        class_array[nodes_with_in, j] = 1

        #############################################################################

    # classification
    unique = np.unique(class_array, axis=0)
    print(len(unique), 'clusters')
    class_vector = np.zeros((num_cells))
    for r, row in enumerate(unique):
        for c in range(num_cells):
            if np.array_equal(class_array[c, :], row):
                class_vector[c] = r

    # colorss = ['g','c','b','m','r','k','y', ]
    colorss = ['tab:' + c for c in ['blue', 'orange', 'red', 'purple', 'brown', 'pink', 'gray', 'cyan', 'green']] + [
        'magenta', 'teal', 'lightcoral']

    # figure
    plt.rc('xtick', labelsize=20)
    plt.rc('ytick', labelsize=20)

    fig, ax = plt.subplots(figsize=(20, 20))
    labels = []

    ########################################
    # Add GO terms
    lowx = np.min(coors[:, 0])
    highx = np.max(coors[:, 0])
    lowy = np.min(coors[:, 1])
    highy = np.max(coors[:, 1])

    ########################################
    # x,y positions of GO terms
    # create an ellipse of GO terms
    npos = len(structs_lst) + (len(structs_lst) % 2)
    mid = [(lowx + highx) / 2, (lowy + highy) / 2]
    radius = [(highx - lowx) * 4 / 5, (highy - lowy) * 4 / 5]
    tempx = np.linspace(mid[0] - radius[0], mid[0] + radius[0], int(npos / 2 + 1))
    x = np.concatenate((tempx, tempx[1:-1]))
    temp = np.sqrt(
        np.power(radius[1], 2) * (np.power(radius[0], 2) - np.power(tempx - mid[0], 2)) / (np.power(radius[0], 2)))
    y = np.concatenate((temp + mid[1], -temp[1:-1] + mid[1]))
    print('x:', x)
    print('y:', y)

    # Reorganize the text positions such that they be close to the upregulating cells of these terms
    newx = np.zeros((len(structs_lst)))
    newy = np.zeros((len(structs_lst)))
    for j, i in enumerate(structs_lst):
        upcells = np.where(class_array[:, j] == 1)[0]
        midupcells = [np.mean(coors[upcells, 0]), np.mean(coors[upcells, 1])]
        dists = np.power(midupcells[0] - x, 2) + np.power(midupcells[1] - y, 2)
        indmin = np.argmin(dists)
        newx[j], newy[j] = x[indmin], y[indmin]
        x = np.concatenate((x[:indmin], x[indmin + 1:]))
        y = np.concatenate((y[:indmin], y[indmin + 1:]))
    x = newx
    y = newy
    print('x:', x)
    print('y:', y)

    # set axis limits
    stdx = np.std(coors[:, 0])
    stdy = np.std(coors[:, 1])
    ax.set_xlim(lowx - stdx / 4, highx + stdx / 4)
    ax.set_ylim(lowy - stdy / 4, highy + stdy / 4)

    for i, xx, yy in zip(structs_lst, x, y):
        goterm = sigtable.loc[i, 'good_summary']
        # split to 2 rows if necessary
        spaces = [i for i, letter in enumerate(goterm) if letter == ' ']
        if len(spaces) > 3:
            goterm = goterm[:spaces[2]] + '\n' + goterm[spaces[2] + 1:]
        # move the left terms a bit more to the left
        if xx < (lowx - stdx / 4):
            ax.text(xx - 0.5, yy, goterm, fontsize=35)
        else:
            ax.text(xx, yy, goterm, fontsize=35)

    for deff, shape in zip(cell_definitions_sorted_set, shapes[:len(set(cell_definitions))]):
        inds_deff = np.where(np.array(cell_definitions) == deff)[0]
        # for r,color in zip(range(len(unique)),colorss[:len(unique)]):
        for r in range(len(unique)):
            inds_cluster = np.where(class_vector == r)[0]
            inds = list(set(inds_cluster).intersection(set(inds_deff)))
            if data_name == 'yaeldiffbulk2':
                if len(deff) == 4:
                    size = 130
                else:
                    size = 250
            else:
                size = 270
            ax.scatter(coors[inds, 0], coors[inds, 1],
                       # c=color,
                       label=deff, s=size, marker=shape)
        labels.append(deff)
    if with_legend:
        ax.legend(labels, prop={'size': 18})

    # Add arrows
    for r, row in enumerate(unique):
        inds_cluster = np.where(class_vector == r)[0]
        for j, val in enumerate(row):
            if val == 1:
                cluster_center = [np.mean(coors[inds_cluster, 0]), np.mean(coors[inds_cluster, 1])]
                arrow1 = mpatches.Arrow(x=cluster_center[0], y=cluster_center[1],
                                        dx=x[j] - cluster_center[0], dy=y[j] - cluster_center[1],
                                        width=0.03, capstyle='projecting', color='red')
                arrow1.set_clip_on(False)
                ax.add_patch(arrow1)
            if val == 0 and (not red_only):
                cluster_center = [np.mean(coors[inds_cluster, 0]), np.mean(coors[inds_cluster, 1])]
                arrow1 = mpatches.Arrow(x=cluster_center[0], y=cluster_center[1],
                                        dx=x[j] - cluster_center[0], dy=y[j] - cluster_center[1],
                                        width=0.1, capstyle='projecting', color='green')
                arrow1.set_clip_on(False)
                ax.add_patch(arrow1)

    # ax.set_title(title,fontdict={'fontsize': 22, 'fontweight': 'medium'})
    ax.set_xlabel(mode + ' 1', fontdict={'fontsize': 22, 'fontweight': 'medium'})
    ax.set_ylabel(mode + ' 2', fontdict={'fontsize': 22, 'fontweight': 'medium'})

    # picture file name
    if red_only:
        picname = data_name + '_' + filter_word + '_' + impute_method + '_' + str(num_stds_thresh) + '_' + str(
            mu) + str(structs_lst) + '_' + mode + '_arrows_red_only.jpg'
    else:
        picname = data_name + '_' + filter_word + '_' + impute_method + '_' + str(num_stds_thresh) + '_' + str(
            mu) + str(structs_lst) + '_' + mode + '_arrows.jpg'

    # save the figure to file
    picfile = picfolder + '/' + picname
    plt.savefig(picfile, bbox_inches='tight')
    # plt.clf()
    plt.close()
