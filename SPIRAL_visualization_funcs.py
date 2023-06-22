#!/var/www/html/SPIRAL.web.tool/spiral_venv/bin/python3.9

import os.path
import pickle
import networkx as nx
import umap.umap_ as umap
from numpy import savetxt, loadtxt
from SPIRAL_basic_funcs import *
import itertools
import matplotlib.pyplot as plt
import matplotlib
import importlib

####################################################################################################################
def visualize_repcell_partition(clustering_file_final,
                                norm_filt_umap_coor_file, norm_filt_pca_coor_file,
                                use_5000,
                                origdeffsfile,
                                repcell_partition_UMAP, repcell_partition_PCA, repcell_partition_spatial,
                                spatial_norm_filt_loc=None,
                                data=None, data_norm_filt_loc=None,
                                save_umap=True, save_pca=True, save_spatial=True,
                                with_legend=False, with_title=False,
                                load_orig_umap_from_file=True, load_orig_pca_from_file=True,
                                with_deffs=False, numerical_shapes=False,
                                ):
    print('visualize_repcell_partition!')

    if data is None and data_norm_filt_loc is None:
        print('No data path and no data argument.')

    if data is None:
        print('Loading data')
        data = pd.read_csv(data_norm_filt_loc, index_col=0, sep='\t')

    # use 5000 most variable genes for coors computation
    if use_5000:
        data = data.loc[data.var(axis=1).sort_values(ascending=False)[:5000].index, :]

    cells_to_repcells = pd.read_csv(clustering_file_final, index_col=0, sep=' ').dropna()

    cell_definitions_sorted_set = None
    shapes = None
    if with_deffs:
        print('loading cell labels and cell definitions from file...')
        cell_definitions_orig = []
        with open(origdeffsfile, "r") as f:
            for line in f:
                cell_definitions_orig.append(str(line.strip()))

        cells_to_repcells.loc[:, 'deffs'] = cell_definitions_orig
        cell_definitions_sorted_set = sort_labels(set(cell_definitions_orig))

        # shapes
        if numerical_shapes:
            shapes = ["$\mathsf{" + str(i) + "}$" for i in cell_definitions_sorted_set]
            with_legend = False
        else:
            shapes = ['8', 's', "^", '*', 'X', 'D', 'P', '1', '2', '3', '4']

    def save_layout(axislabel, cells_to_repcells, coor, picfile, numerical_shapes,
                    with_legend=False, with_title=False,
                    with_deffs=False, cell_definitions_sorted_set=None, shapes=None):
        plt.rc('xtick', labelsize=18)
        plt.rc('ytick', labelsize=18)

        cm = plt.get_cmap('gist_rainbow')

        fig, ax = plt.subplots(figsize=(15, 15))
        num_repcells = len(set(cells_to_repcells['x']))
        ax.set_prop_cycle('color', np.random.permutation([cm(1. * i / num_repcells) for i in range(num_repcells)]))

        marker_size = get_marker_size(num_points=num_repcells, numerical_shapes=numerical_shapes)

        if with_deffs:
            for deff, shape in zip(cell_definitions_sorted_set, shapes[:len(cell_definitions_sorted_set)]):
                inds_deff = [i for i in cells_to_repcells.index if cells_to_repcells.loc[i, 'deffs'] == deff]
                for repcell in set(cells_to_repcells['x']):
                    inds_repcell = [i for i in cells_to_repcells.index if cells_to_repcells.loc[i, 'x'] == repcell]
                    inds = set(inds_repcell).intersection(set(inds_deff))
                    ax.scatter(coor[inds, 0], coor[inds, 1],
                               label=str(deff) + ', repcell ' + str(int(repcell)), s=marker_size,
                               marker=shape)
        else:
            for repcell in set(cells_to_repcells['x']):
                inds = [i for i in cells_to_repcells.index if cells_to_repcells.loc[i, 'x'] == repcell]
                ax.scatter(coor[inds, 0], coor[inds, 1], label='repcell ' + str(int(repcell)), s=marker_size,
                           marker="$" + str(int(repcell)) + "$")

        if with_title:
            ax.set_title('Repcell partition', fontdict={'fontsize': 30, 'fontweight': 'medium'})
        ax.set_xlabel(axislabel + ' 1', fontdict={'fontsize': 22, 'fontweight': 'medium'})
        ax.set_ylabel(axislabel + ' 2', fontdict={'fontsize': 22, 'fontweight': 'medium'})

        if with_legend:
            ax.legend()
        plt.tight_layout()
        plt.savefig(picfile)
        plt.close()

        # umap, PCA, spatial coordinates

    if save_umap:
        # get coordinates
        if load_orig_umap_from_file and os.path.exists(norm_filt_umap_coor_file):
            print('loading UMAP of original data set...')
            umap_coor_orig = loadtxt(norm_filt_umap_coor_file, delimiter=',')
        else:
            print('computing UMAP of original data set...')
            umap_coor_orig = umap.UMAP(random_state=42).fit_transform(data.transpose())
            savetxt(norm_filt_umap_coor_file, umap_coor_orig, delimiter=',')

        # save repcell partition
        save_layout(axislabel='UMAP', cells_to_repcells=cells_to_repcells, coor=umap_coor_orig,
                    picfile=repcell_partition_UMAP,
                    with_legend=with_legend, with_title=with_title, numerical_shapes=numerical_shapes,
                    with_deffs=with_deffs, cell_definitions_sorted_set=cell_definitions_sorted_set, shapes=shapes)

    if save_pca:
        # get coordinates
        if load_orig_pca_from_file and os.path.exists(norm_filt_pca_coor_file):
            print('loading PCA of original data set...')
            pca_coor_orig = loadtxt(norm_filt_pca_coor_file, delimiter=',')
        else:
            print('computing PCA of original data set...')
            pca = PCA(n_components=2)
            datamatrix = normalize(data.transpose())
            pca_coor_orig = pca.fit_transform(datamatrix)
            savetxt(norm_filt_pca_coor_file, pca_coor_orig, delimiter=',')

        # save repcell partition
        save_layout(axislabel='PC', cells_to_repcells=cells_to_repcells, coor=pca_coor_orig,
                    picfile=repcell_partition_PCA,
                    with_legend=with_legend, with_title=with_title, numerical_shapes=numerical_shapes,
                    with_deffs=with_deffs, cell_definitions_sorted_set=cell_definitions_sorted_set, shapes=shapes)

    if save_spatial:
        # get coordinates
        print('loading spatial coordinates of data set...')
        spatial_coor = loadtxt(spatial_norm_filt_loc, delimiter=',')

        # save repcell partition
        save_layout(axislabel='spatial axis', cells_to_repcells=cells_to_repcells, coor=spatial_coor,
                    picfile=repcell_partition_spatial,
                    with_legend=with_legend, with_title=with_title, numerical_shapes=numerical_shapes,
                    with_deffs=with_deffs, cell_definitions_sorted_set=cell_definitions_sorted_set, shapes=shapes)


##################################################################################################################
def compute_orig_deffs(data_path, data):
    labels_checkbox = open(os.path.join(data_path, 'labels_checkbox.txt'), "r").read()
    if (labels_checkbox == 'True') or (not all('_' in s for s in list(data))):
        return ['0'] * data.shape[1]
    return [s[:s.find('_')] for s in list(data)]


##################################################################################################################
def save_orig_deffs_to_file(data, data_path, data_norm_filt_loc, origdeffsfile):
    # load count matrix
    if data is None and data_norm_filt_loc is None:
        print('No data path and no data argument.')

    if data is None:
        print('Loading data')
        data = pd.read_csv(data_norm_filt_loc, index_col=0, sep='\t')

    orig_deffs = compute_orig_deffs(data_path, data)
    with open(origdeffsfile, "w") as f:
        for s in orig_deffs:
            f.write(str(s) + "\n")


##################################################################################################################
def decide_on_numerical_shapes(origdeffsfile):
    # Decide whether to use regular markers or numerical markers in the structures' layouts
    print('loading cell definitions from file...')
    cell_definitions_orig = []
    with open(origdeffsfile, "r") as f:
        for line in f:
            cell_definitions_orig.append(str(line.strip()))

    # if there are at least 5 different labels, and each label is no longer than 3, then use numerical_shapes.
    if len(set(cell_definitions_orig)) >= 5 and max([len(d) for d in cell_definitions_orig]) <= 3:
        return True
    return False


##################################################################################################################
def get_marker_size(num_points, numerical_shapes):
    if not numerical_shapes:
        return 130
    else:
        if num_points <= 150:
            return 500
        else:
            return 300


##################################################################################################################
def compute_deffs_for_imputed(impute_method, repcelldeffsfile, origdeffsfile, clustering_file_final,
                              data_path,
                              repcells_data=None, repcells_data_loc=None, stepsfile=None, with_steps=False,
                              data=None, data_norm_filt_loc=None):
    print('compute_deffs_for_imputed!')

    if stepsfile == None:
        with_steps = False

    if repcells_data is None and repcells_data_loc is None:
        print('No repcells_data path and no repcells_data argument.')

    if repcells_data is None:
        print('Loading repcells_data')
        repcells_data = pd.read_csv(repcells_data_loc, index_col=0, sep='\t')

    if impute_method == 'no_imputation':
        repcell_deffs = [s[:s.find('_')] for s in list(repcells_data)]
    else:
        # Load clustering file
        clusters = pd.read_csv(clustering_file_final, index_col=0, sep=' ')
        clusters.index = range(len(clusters))
        clusters.columns = ['repcell']
        clusters = clusters.dropna()

        if not os.path.exists(origdeffsfile):
            save_orig_deffs_to_file(data=data, data_path=data_path, data_norm_filt_loc=data_norm_filt_loc,
                                    origdeffsfile=origdeffsfile)

        # Load orig_deffs
        orig_deffs = []
        with open(origdeffsfile, "r") as f:
            for line in f:
                orig_deffs.append(str(line.strip()))

        clusters.loc[:, 'deffs'] = [orig_deffs[i] for i in clusters.index]
        deffs_list = list(set(clusters['deffs']))
        final_cluster_list = list(set(clusters['repcell']))
        print('len(final_cluster_list):', len(final_cluster_list))

        # clusters_table = pd.DataFrame(columns=labels, data=np.zeros((len(final_cluster_list),len(labels))),
        #                             index=final_cluster_list)
        deff_for_cluster_majority_vote = np.zeros((len(final_cluster_list))).astype(int)
        if with_steps:
            step_for_cluster = np.zeros((len(final_cluster_list))).astype(int)
            # Load filtered scRNA_data_orig
            if data is None and data_norm_filt_loc is None:
                print('No data path and no data argument.')

            if data is None:
                print('Loading data')
                data = pd.read_csv(data_norm_filt_loc, index_col=0, sep='\t')

            # load colData
            coldata = pd.read_table(stepsfile, index_col=0, sep=',')
            coldata = coldata.loc[list(data), :]
            clusters['steps'] = np.array(coldata.iloc[clusters.index, :].loc[:, 'Step'].astype(int))

        for ci, c in enumerate(final_cluster_list):
            # print(c)
            best_label_count = 0
            for li, l in enumerate(deffs_list):
                len_c_l = len([ind for ind in clusters.index if
                               clusters.loc[ind, 'repcell'] == c and clusters.loc[ind, 'deffs'] == l])
                # clusters_table.loc[c,l] = len_c_l
                if len_c_l > best_label_count:
                    deff_for_cluster_majority_vote[ci] = li
                    best_label_count = len_c_l
            if with_steps:
                step_for_cluster[ci] = np.int(np.mean(clusters.loc[[i for i in clusters.index if
                                                                    clusters.loc[i, 'repcell'] == c and clusters.loc[
                                                                        i, 'deffs'] == deffs_list[
                                                                        deff_for_cluster_majority_vote[ci]]], 'steps']))

        repcell_deffs = [deffs_list[li] for li in deff_for_cluster_majority_vote]

        if with_steps:
            with open(stepsfile, "w") as f:
                for s in step_for_cluster:
                    f.write(str(s) + "\n")

    with open(repcelldeffsfile, "w") as f:
        for s in repcell_deffs:
            f.write(str(s) + "\n")


########################################################################################################
def save_simple_layout(axislabel, coors, cell_definitions, cell_definitions_sorted_set, deff_to_shape_dict,
                       picfile, numerical_shapes,
                       steps=None, with_steps=False, with_legend=True, with_title=False, title=''):
    plt.rc('xtick', labelsize=18)
    plt.rc('ytick', labelsize=18)

    cm = plt.get_cmap('gist_rainbow')

    marker_size = get_marker_size(num_points=len(cell_definitions), numerical_shapes=numerical_shapes)

    fig, ax = plt.subplots(figsize=(15, 15))
    for deff in cell_definitions_sorted_set:
        inds = [i for i, d in enumerate(cell_definitions) if d == deff]
        if inds:
            if with_steps and (steps is not None):
                im = ax.scatter(coors[inds, 0], coors[inds, 1], c=steps[inds], cmap=plt.cm.jet,
                                marker=deff_to_shape_dict[deff], label=deff, s=marker_size)
                im.set_clim(0, 1000 * int(np.ceil(np.max(steps) / 1000)))
            else:
                im = ax.scatter(coors[inds, 0], coors[inds, 1],
                                marker=deff_to_shape_dict[deff], label=deff, s=marker_size)

    if with_legend and len(cell_definitions_sorted_set) > 1:
        # Add a legend
        ax.legend(cell_definitions_sorted_set, prop={'size': 18})
    if with_steps and (steps is not None):
        # Add a colorbar
        cbar = fig.colorbar(im, ax=ax)
        cbar.set_label('steps along a path', size=18)

    if with_title:
        ax.set_title(title, fontdict={'fontsize': 30, 'fontweight': 'medium'})
    ax.set_xlabel(axislabel + ' 1', fontdict={'fontsize': 22, 'fontweight': 'medium'})
    ax.set_ylabel(axislabel + ' 2', fontdict={'fontsize': 22, 'fontweight': 'medium'})

    plt.tight_layout()
    plt.savefig(picfile)
    plt.close()


########################################################################################################
def sort_labels(label_set):
    label_set_sorted = []
    lens = list(set([len(r) for r in list(set(label_set))]))
    for l in lens:
        b = list(set([r for r in set(label_set) if len(r) == l]))
        b.sort()
        label_set_sorted += b
    return label_set_sorted


####################################################################################################################
def save_layout_of_imputed(repcells_data, impute_method, repcells_data_loc,
                           layout_of_imputed_picfile_PCA, layout_of_imputed_picfile_UMAP,
                           repcelldeffsfile, origdeffsfile, stepsfile,
                           repcellsumapcoorfile, repcellspcacoorfile,
                           clustering_file_final, data_path,
                           data_norm_filt_loc,
                           with_steps=False,
                           save_UMAP=True, save_PCA=True, use_5000=False,
                           load_coor_from_file=True, with_legend=True, with_title=True,
                           numerical_shapes=False):
    print('save_layout_of_imputed!')

    if repcells_data is None and repcells_data_loc is None:
        print('No repcells_data path and no repcells_data argument.')

    if repcells_data is None:
        print('Loading repcells_data')
        repcells_data = pd.read_csv(repcells_data_loc, index_col=0, sep='\t')

    if use_5000:
        repcells_data = repcells_data.loc[repcells_data.var(axis=1).sort_values(ascending=False)[:5000].index, :]

    # cell labels
    if with_steps:
        if not os.path.exists(repcelldeffsfile) or not os.path.exists(stepsfile):
            compute_deffs_for_imputed(impute_method=impute_method, repcelldeffsfile=repcelldeffsfile,
                                      origdeffsfile=origdeffsfile, clustering_file_final=clustering_file_final,
                                      repcells_data=repcells_data, stepsfile=stepsfile, with_steps=True,
                                      data_norm_filt_loc=data_norm_filt_loc, data_path=data_path)
    elif not os.path.exists(repcelldeffsfile):
        compute_deffs_for_imputed(impute_method=impute_method, repcelldeffsfile=repcelldeffsfile,
                                  origdeffsfile=origdeffsfile, clustering_file_final=clustering_file_final,
                                  repcells_data=repcells_data, with_steps=False,
                                  data_norm_filt_loc=data_norm_filt_loc, data_path=data_path)

    repcell_definitions = []
    with open(repcelldeffsfile, "r") as f:
        for line in f:
            repcell_definitions.append(str(line.strip()))
    if with_steps:
        step_for_cluster = []
        with open(stepsfile, "r") as f:
            for line in f:
                step_for_cluster.append(int(line.strip()))
        step_for_cluster = np.array(step_for_cluster)
    else:
        step_for_cluster = None

    ############################################# cell_definitions_sorted #################################
    repcell_definitions_sorted_set = sort_labels(repcell_definitions)

    # print('cell_definitions:', cell_definitions)
    print('cell_definitions_sorted_set:', repcell_definitions_sorted_set)

    # shapes
    if numerical_shapes:
        shapes = ["$\mathsf{" + str(i) + "}$" for i in repcell_definitions_sorted_set]
        with_legend = False
    else:
        shapes = ['8', 's', "^", '*', 'X', 'D', 'P', '1', '2', '3', '4']

    deff_to_shape_dict = dict(zip(repcell_definitions_sorted_set, shapes[:len(set(repcell_definitions))]))
    # print('deff_to_shape_dict:', deff_to_shape_dict)

    ###############################################################################################################
    # umap and PCA
    if save_UMAP:
        if load_coor_from_file and os.path.exists(repcellsumapcoorfile):
            print('loading UMAP...')
            coors = loadtxt(repcellsumapcoorfile, delimiter=',')
        else:
            print('computing UMAP...')
            coors = umap.UMAP(random_state=42).fit_transform(repcells_data.transpose())
            savetxt(repcellsumapcoorfile, coors, delimiter=',')

        print('Saving UMAP layout')
        save_simple_layout(axislabel='UMAP', coors=coors, cell_definitions=repcell_definitions,
                           cell_definitions_sorted_set=repcell_definitions_sorted_set,
                           deff_to_shape_dict=deff_to_shape_dict, picfile=layout_of_imputed_picfile_UMAP,
                           steps=step_for_cluster, with_steps=with_steps,
                           with_legend=with_legend, with_title=with_title, numerical_shapes=numerical_shapes)
    if save_PCA:
        if load_coor_from_file and os.path.exists(repcellspcacoorfile):
            print('loading PCA...')
            coors = loadtxt(repcellspcacoorfile, delimiter=',')
        else:
            print('computing PCA...')
            pca = PCA(n_components=2)
            datamatrix = normalize(repcells_data.transpose())
            coors = pca.fit_transform(datamatrix)
            savetxt(repcellspcacoorfile, coors, delimiter=',')

        print('Saving PCA layout')
        save_simple_layout(axislabel='PC', coors=coors, cell_definitions=repcell_definitions,
                           cell_definitions_sorted_set=repcell_definitions_sorted_set,
                           deff_to_shape_dict=deff_to_shape_dict, picfile=layout_of_imputed_picfile_PCA,
                           steps=step_for_cluster, with_steps=with_steps,
                           with_legend=with_legend, with_title=with_title, numerical_shapes=numerical_shapes)


####################################################################################################################
def save_layout_of_orig(norm_filt_pca_coor_file, norm_filt_umap_coor_file,
                        norm_filt_picfile_pca, norm_filt_picfile_umap, norm_filt_picfile_spatial,
                        origdeffsfile, data_path,
                        spatial_norm_filt_loc=None,
                        data=None, data_norm_filt_loc=None,
                        save_umap=True, save_pca=True, save_spatial=False,
                        use_5000=False, numerical_shapes=False,
                        with_steps=False, stepsfile=None,
                        load_coor_from_file=True, with_legend=True, with_title=True):
    print('save_layout_of_orig!')
    if data is None and data_norm_filt_loc is None:
        print('No data path and no data argument.')

    if data is None:
        print('Loading data')
        data = pd.read_csv(data_norm_filt_loc, index_col=0, sep='\t')

    # use 5000 most variable genes for coors computation
    if use_5000:
        data = data.loc[data.var(axis=1).sort_values(ascending=False)[:5000].index, :]

    # cell labels
    if not os.path.exists(origdeffsfile):
        save_orig_deffs_to_file(data=data, data_path=data_path, data_norm_filt_loc=data_norm_filt_loc,
                                origdeffsfile=origdeffsfile)

    print('loading cell labels and cell definitions from file...')
    cell_definitions_orig = []
    with open(origdeffsfile, "r") as f:
        for line in f:
            cell_definitions_orig.append(str(line.strip()))

    ############################################# cell_definitions_sorted #################################
    cell_definitions_orig_sorted_set = sort_labels(set(cell_definitions_orig))
    print('cell_definitions_orig_sorted_set:', cell_definitions_orig_sorted_set)

    # shapes
    if numerical_shapes:
        shapes = ["$\mathsf{" + str(i) + "}$" for i in cell_definitions_orig_sorted_set]
        with_legend = False
    else:
        shapes = ['8', 's', "^", '*', 'X', 'D', 'P', '1', '2', '3', '4']

    deff_to_shape_dict = dict(zip(cell_definitions_orig_sorted_set, shapes[:len(set(cell_definitions_orig))]))
    print('deff_to_shape_dict:', deff_to_shape_dict)

    ########################################################################################################
    if with_steps:
        # load colData
        coldata = pd.read_table(stepsfile, index_col=0, sep=',')
        coldata = coldata.loc[list(data), :]
        steps = np.array(coldata.loc[:, 'Step'].astype(int))
    else:
        steps = None

    ###############################################################################################################

    if save_umap:
        if load_coor_from_file and os.path.exists(norm_filt_umap_coor_file):
            print('loading UMAP of original data set...')
            coors = loadtxt(norm_filt_umap_coor_file, delimiter=',')
        else:
            print('computing UMAP of original data set...')
            coors = umap.UMAP(random_state=42).fit_transform(data.transpose())
            savetxt(norm_filt_umap_coor_file, coors, delimiter=',')

        print('Saving UMAP layout')
        save_simple_layout(axislabel='UMAP', coors=coors, cell_definitions=cell_definitions_orig,
                           cell_definitions_sorted_set=cell_definitions_orig_sorted_set,
                           deff_to_shape_dict=deff_to_shape_dict, picfile=norm_filt_picfile_umap,
                           steps=steps, with_steps=with_steps,
                           with_legend=with_legend, with_title=with_title, numerical_shapes=numerical_shapes)

    if save_pca:
        if load_coor_from_file and os.path.exists(norm_filt_pca_coor_file):
            print('loading PCA of original data set...')
            coors = loadtxt(norm_filt_pca_coor_file, delimiter=',')
        else:
            print('computing PCA of original data set...')
            pca = PCA(n_components=2)
            datamatrix = normalize(data.transpose())
            coors = pca.fit_transform(datamatrix)
            savetxt(norm_filt_pca_coor_file, coors, delimiter=',')

        print('Saving PCA layout')
        save_simple_layout(axislabel='PC', coors=coors, cell_definitions=cell_definitions_orig,
                           cell_definitions_sorted_set=cell_definitions_orig_sorted_set,
                           deff_to_shape_dict=deff_to_shape_dict, picfile=norm_filt_picfile_pca,
                           steps=steps, with_steps=with_steps,
                           with_legend=with_legend, with_title=with_title, numerical_shapes=numerical_shapes)

    if save_spatial:
        print('loading spatial coordinates of original data set...')
        coors = loadtxt(spatial_norm_filt_loc, delimiter=',')

        print('Saving PCA layout')
        save_simple_layout(axislabel='spatial axis', coors=coors, cell_definitions=cell_definitions_orig,
                           cell_definitions_sorted_set=cell_definitions_orig_sorted_set,
                           deff_to_shape_dict=deff_to_shape_dict, picfile=norm_filt_picfile_spatial,
                           steps=steps, with_steps=with_steps,
                           with_legend=with_legend, with_title=with_title, numerical_shapes=numerical_shapes)


####################################################################################################################
def cluster_to_levels_based_on_sets_in_struct(all_nodes, sets_in_struct):
    # clustering- right to left (arrows are as short as possible)
    nodes_in_struct = set(itertools.chain.from_iterable(sets_in_struct))
    nodes_with_out = set([s[0] for s in sets_in_struct])
    nodes_with_in = set([s[1] for s in sets_in_struct])
    nodes_with_in_and_out = nodes_with_out.intersection(nodes_with_in)
    nodes_with_out_only = nodes_with_out - nodes_with_in_and_out
    list_of_nodes_in_levels = [set()] * 100
    clustered_nodes = set()
    level = 1
    while (nodes_in_struct - clustered_nodes):
        list_of_nodes_in_levels[level] = set(
            [s[1] for s in sets_in_struct if s[0] in nodes_with_out_only.union(clustered_nodes)]
        ) - set(
            [s[1] for s in sets_in_struct if s[0] not in nodes_with_out_only.union(clustered_nodes)]
        ) - clustered_nodes

        list_of_nodes_in_levels[level - 1] = list_of_nodes_in_levels[level - 1].union(
            set([s[0] for s in sets_in_struct if s[1] in list_of_nodes_in_levels[level]]) - clustered_nodes)

        clustered_nodes = clustered_nodes.union(list_of_nodes_in_levels[level])
        clustered_nodes = clustered_nodes.union(list_of_nodes_in_levels[level - 1])
        level += 1

    list_of_nodes_in_levels_dict = {i: list(l) for i, l in enumerate(list_of_nodes_in_levels)}
    node_level_dict = {v: k for k, l in list_of_nodes_in_levels_dict.items() for v in l}
    node_level_name_dict = {v: ('level ' + str(k + 1)) for k, l in list_of_nodes_in_levels_dict.items() for v in l}
    # print('node_level_dict:', node_level_dict)

    for n in set(all_nodes) - nodes_in_struct:
        node_level_name_dict[n] = 'not in structure'

    '''
    # Another method of clustering- left to right
    nodes_with_out_only = nodes_with_out-nodes_with_in_and_out
    nodes_curr_level = nodes_with_out_only
    level = 1
    node_color_dict = dict()
    clustered_nodes = []
    while nodes_curr_level:
        #print('level:', level)
        #print('nodes_curr_level:', nodes_curr_level)
        for n in nodes_curr_level:
            node_color_dict[n] = 'level '+str(level)
        clustered_nodes += list(nodes_curr_level)
        nodes_curr_level = []
        for node in (nodes_with_in-set(clustered_nodes)):
            in_nodes = [s[0] for s in sets_in_struct if s[1]==node]
            if set(in_nodes).issubset(set(clustered_nodes)):
                nodes_curr_level.append(node)
        level += 1
    '''
    return list_of_nodes_in_levels, node_level_dict, node_level_name_dict


####################################################################################################################
def save_layout_gradual(cell_table, picfile, data, genelist, cell_definitions_sorted_set, colormap, axislabel,
                        coor, shapes, title, color_by_log_expression, numerical_shapes):
    plt.rc('xtick', labelsize=18)
    plt.rc('ytick', labelsize=18)

    fig, ax = plt.subplots(figsize=(15, 15))
    genelist_ = list(set(genelist).intersection(set(data.index)))
    if len(genelist_) < len(genelist):
        print('missing genes:', set(genelist) - set(genelist_))

    marker_size = get_marker_size(num_points=coor.shape[0], numerical_shapes=numerical_shapes)

    for deff, shape in zip(cell_definitions_sorted_set, shapes[:len(cell_definitions_sorted_set)]):
        inds = [i for i in cell_table.index if cell_table.loc[i, 'deffs'] == deff]
        if inds:
            size = 130
            if len(cell_definitions_sorted_set) == 1:
                if color_by_log_expression:
                    h = ax.scatter(coor[inds, 0], coor[inds, 1], cmap=colormap,
                                   vmin=np.log10(data.loc[genelist_, :].mean().min()),
                                   vmax=np.log10(data.loc[genelist_, :].mean().max()),
                                   s=marker_size,
                                   marker=shape, c=np.log10(data.iloc[:, inds].loc[genelist_, :].mean()))
                else:
                    h = ax.scatter(coor[inds, 0], coor[inds, 1], cmap=colormap,
                                   vmin=data.loc[genelist_, :].mean().min(),
                                   vmax=data.loc[genelist_, :].mean().max(),
                                   s=marker_size,
                                   marker=shape, c=data.iloc[:, inds].loc[genelist_, :].mean())
            else:
                if color_by_log_expression:
                    h = ax.scatter(coor[inds, 0], coor[inds, 1], cmap=colormap,
                                   vmin=np.log10(data.loc[genelist_, :].mean().min()),
                                   vmax=np.log10(data.loc[genelist_, :].mean().max()),
                                   s=marker_size, label=deff,
                                   marker=shape, c=np.log10(data.iloc[:, inds].loc[genelist_, :].mean()))
                else:
                    h = ax.scatter(coor[inds, 0], coor[inds, 1], cmap=colormap,
                                   vmin=data.loc[genelist_, :].mean().min(),
                                   vmax=data.loc[genelist_, :].mean().max(),
                                   s=marker_size, label=deff,
                                   marker=shape, c=data.iloc[:, inds].loc[genelist_, :].mean())

    cbar = plt.colorbar(h)
    # cbar.ax.tick_params(labelsize=30)
    cbar.set_label((color_by_log_expression * 'log10 ') + 'average expression level', size=20)

    # ax.legend()
    ax.set_title(title, fontdict={'fontsize': 30, 'fontweight': 'medium'})
    ax.set_xlabel(axislabel + ' 1', fontdict={'fontsize': 22, 'fontweight': 'medium'})
    ax.set_ylabel(axislabel + ' 2', fontdict={'fontsize': 22, 'fontweight': 'medium'})

    plt.tight_layout()
    plt.savefig(picfile)
    plt.close()


####################################################################################################################
def save_struct_layout(cell_table, picfile, cell_definitions_sorted_set, curr_colors,
                       axislabel, coor, shapes, title, level_desc_dict, numerical_shapes, size_by_degree=False,
                       num_arrows_per_cell=None, with_legend=False):
    # cell_table has the columns "levels" and "cell_definitions".
    # Its indices are the indices of the rows in coors that will be visualized.

    matplotlib.pyplot.rc('xtick', labelsize=18)
    matplotlib.pyplot.rc('ytick', labelsize=18)

    marker_size = get_marker_size(num_points=coor.shape[0], numerical_shapes=numerical_shapes)

    fig, ax = matplotlib.pyplot.subplots(figsize=(15, 15))
    hs = []
    labels = []
    for level, color in zip(sorted(set(cell_table['levels'])), curr_colors):
        inds_level = [i for i in cell_table.index if cell_table.loc[i, 'levels'] == level]
        for deff, shape in zip(cell_definitions_sorted_set, shapes[:len(cell_definitions_sorted_set)]):
            inds_deff = [i for i in cell_table.index if cell_table.loc[i, 'deffs'] == deff]
            inds = list(set(inds_level).intersection(set(inds_deff)))
            if inds:
                if level == 'not in structure':
                    h = ax.scatter(coor[inds, 0], coor[inds, 1], color=color,
                                   label=(len(cell_definitions_sorted_set) != 1) * (deff + ' - ') + level_desc_dict[
                                       level],
                                   s=marker_size, marker=shape)
                    hs.append(h)
                    labels.append((len(cell_definitions_sorted_set) != 1) * (deff + ' - ') + level_desc_dict[level])
                else:
                    if size_by_degree:
                        poss_vals = [num_arrows_per_cell[ind] for ind in inds]
                        for num_arrows in np.arange(1, 1 + max(poss_vals)):
                            curr_inds = [ind for ind in inds if num_arrows_per_cell[ind] == num_arrows]
                            h = ax.scatter(coor[curr_inds, 0], coor[curr_inds, 1], color=color,
                                           label=(len(cell_definitions_sorted_set) != 1) * (deff + ' - ') +
                                                 level_desc_dict[
                                                     level], s=130 + 15 * (num_arrows - 1), marker=shape)
                            if num_arrows == 1:
                                hs.append(h)
                                labels.append(
                                    (len(cell_definitions_sorted_set) != 1) * (deff + ' - ') + level_desc_dict[level])
                    else:
                        h = ax.scatter(coor[inds, 0], coor[inds, 1], color=color,
                                       label=(len(cell_definitions_sorted_set) != 1) * (deff + ' - ') + level_desc_dict[
                                           level], s=marker_size, marker=shape)
                        hs.append(h)
                        labels.append((len(cell_definitions_sorted_set) != 1) * (deff + ' - ') + level_desc_dict[level])

    ax.set_title(title, fontdict={'fontsize': 30, 'fontweight': 'medium'})
    ax.set_xlabel(axislabel + ' 1', fontdict={'fontsize': 22, 'fontweight': 'medium'})
    ax.set_ylabel(axislabel + ' 2', fontdict={'fontsize': 22, 'fontweight': 'medium'})

    if with_legend:
        ax.legend(hs, labels, prop={'size': 22})

    matplotlib.pyplot.tight_layout()
    matplotlib.pyplot.savefig(picfile)
    matplotlib.pyplot.close()


####################################################################################################################
def visualize_structures(sigfile_vis, genetable_file, repcellsumapcoorfile, repcellspcacoorfile,
                         norm_filt_pca_coor_file, repcelldeffsfile, origdeffsfile, picfolder,
                         norm_filt_umap_coor_file, data_norm_filt_loc, clustering_file_final, spatial_norm_filt_loc,
                         data_path,
                         numerical_shapes=False, data=None,
                         repcells_data=None, repcells_data_loc=None,
                         sigtable=None, sigtable_file=None,
                         impute_method='agg_wald',
                         save_table=True, round_pvals=True,
                         save_network_layouts=True, save_Gnetwork_layouts=True,
                         save_UMAPs_repcells=True, save_PCAs_repcells=True,
                         save_UMAPs_orig=True, save_PCAs_orig=True,
                         save_GUMAPs_repcells=True, save_GUMAPs_orig=True,
                         save_GPCAs_repcells=True, save_GPCAs_orig=True,
                         save_spatial=True, save_Gspatial=True,
                         with_legend=True,
                         color_by_log_expression=False,
                         with_links_in_sigtable=False,
                         real_samp_name='sample',
                         save_layers_to_excel=False):
    print('visualize_structures!')

    if sigtable is None and sigtable_file is None:
        print('No significance_table path and no significance_table argument.')

    if sigtable is None:
        print('Loading data')
        sigtable = load_excel_with_openpyxl_and_convert_to_pd_DataFrame(sigtable_file)

    if data is None and data_norm_filt_loc is None:
        print('No data path and no data_norm_filt_loc argument.')

    if data is None:
        print('Loading data')
        data = pd.read_csv(data_norm_filt_loc, sep='\t', index_col=0)

    if impute_method != 'no_imputation':
        if repcells_data is None and repcells_data_loc is None:
            print('No repcells_data path and no repcells_data argument.')

        if repcells_data is None:
            print('Loading data')
            repcells_data = pd.read_csv(repcells_data_loc, sep='\t', index_col=0)

    if impute_method == 'no_imputation':
        samp_name = real_samp_name
    else:
        samp_name = 'repcell'

    num_samples = int(sigtable.iloc[0, :]['num_' + samp_name + 's'])
    num_genes = int(sigtable.iloc[0, :]['num_genes'])
    print('num_samples:', num_samples, ' num_genes:', num_genes)

    # Load genes_table
    with open(genetable_file, 'rb') as fp:
        genes_table = pickle.load(fp)

    '''
    def read_imputed():
        # Read imputed scRNA_data from file
        if read_scRNA_data_from_file:
            print('Loading imputed scRNA data...')
            ind = [m.start() for m in re.finditer('/', data_dict[data_name])][-1]
            repcells_data_loc = data_dict[data_name][:ind + 1] + impute_method + data_dict[data_name][
                                                                            ind:-4] + '_' + filter_word + '_imputed.csv'
            scRNA_data = load_data(repcells_data_loc)
        elif scRNA_data==np.nan:
            print('ERROR: no scRNA_data matrix')
        # Take top 5000 variable genes
        if use_5000:
            print('Filtering for top 5000 variable genes...')
            scRNA_data = scRNA_data.loc[scRNA_data.var(axis=1).sort_values(ascending=False)[:5000].index, :]
        return scRNA_data
    
    def read_orig():     
        #Read original data set from file
        print('Loading original scRNA data...')
        scRNA_data_orig = load_data(data_dict[data_name][:-4]+'_'+filter_word+'.csv')
        if use_5000:
            print('Filtering for top 5000 variable genes...')
            scRNA_data_orig = scRNA_data_orig.loc[
                              scRNA_data_orig.var(axis=1).sort_values(ascending=False)[:5000].index, :]
        return scRNA_data_orig
    
    
    scRNA_data = read_imputed()
    if save_UMAPs_orig or save_PCAs_orig or save_spatial or save_Gspatial or save_GUMAPs_orig or save_GPCAs_orig:
        scRNA_data_orig = read_orig()
    '''

    # umap and PCA coordinates
    if save_UMAPs_repcells:
        if os.path.exists(repcellsumapcoorfile):
            print('loading UMAP...')
            umap_coor = loadtxt(repcellsumapcoorfile, delimiter=',')
        else:
            print('computing UMAP...')
            umap_coor = umap.UMAP(random_state=42).fit_transform(repcells_data.transpose())
            savetxt(repcellsumapcoorfile, umap_coor, delimiter=',')
    if save_PCAs_repcells:
        if os.path.exists(repcellspcacoorfile):
            print('loading PCA...')
            pca_coor = loadtxt(repcellspcacoorfile, delimiter=',')
        else:
            print('computing PCA...')
            pca = PCA(n_components=2)
            datamatrix = normalize(repcells_data.transpose())
            pca_coor = pca.fit_transform(datamatrix)
            savetxt(repcellspcacoorfile, pca_coor, delimiter=',')
    if save_UMAPs_orig:
        if os.path.exists(norm_filt_umap_coor_file):
            print('loading UMAP of original data set...')
            umap_coor_orig = loadtxt(norm_filt_umap_coor_file, delimiter=',')
        else:
            print('computing UMAP of original data set...')
            umap_coor_orig = umap.UMAP(random_state=42).fit_transform(data.transpose())
            savetxt(norm_filt_umap_coor_file, umap_coor_orig, delimiter=',')
    if save_PCAs_orig:
        if os.path.exists(norm_filt_pca_coor_file):
            print('loading PCA of original data set...')
            pca_coor_orig = loadtxt(norm_filt_pca_coor_file, delimiter=',')
        else:
            print('computing PCA of original data set...')
            pca = PCA(n_components=2)
            datamatrix = normalize(data.transpose())
            pca_coor_orig = pca.fit_transform(datamatrix)
            savetxt(norm_filt_pca_coor_file, pca_coor_orig, delimiter=',')

    '''
    if save_GUMAPs_repcells or save_GPCAs_repcells:
        repcells_data_norm = norm_0_max(repcells_data,1000)
    if save_GUMAPs_orig or save_GPCAs_orig:
        data_norm = norm_0_max(data,1000)
    '''

    # load mapping of original cells to repcells
    if impute_method != 'no_imputation':
        if save_UMAPs_orig or save_PCAs_orig or save_spatial or save_Gspatial or save_GUMAPs_orig or save_GPCAs_orig:
            cells_to_repcells = pd.read_csv(clustering_file_final, index_col=0, sep=' ').dropna()

    # Loading spatial coors
    if save_spatial or save_Gspatial:
        spatial_coor = loadtxt(spatial_norm_filt_loc, delimiter=',')

    # definitions
    if not os.path.exists(origdeffsfile):
        save_orig_deffs_to_file(data=data, data_path=data_path, data_norm_filt_loc=data_norm_filt_loc,
                                origdeffsfile=origdeffsfile)

    if impute_method != 'no_imputation' and not os.path.exists(repcelldeffsfile):
        compute_deffs_for_imputed(impute_method=impute_method, repcelldeffsfile=repcelldeffsfile,
                                  origdeffsfile=origdeffsfile, clustering_file_final=clustering_file_final,
                                  repcells_data_loc=repcells_data_loc, data_norm_filt_loc=data_norm_filt_loc,
                                  data_path=data_path)

    print('loading cell definitions from file...')

    if save_UMAPs_orig or save_PCAs_orig or save_spatial or save_Gspatial or save_GUMAPs_orig or save_GPCAs_orig:
        cell_definitions_orig = []
        with open(origdeffsfile, "r") as f:
            for line in f:
                cell_definitions_orig.append(str(line.strip()))
        if impute_method != 'no_imputation':
            cells_to_repcells.loc[:, 'deffs'] = [cell_definitions_orig[i] for i in cells_to_repcells.index]

    if save_UMAPs_repcells or save_GUMAPs_repcells or save_PCAs_repcells or save_GPCAs_repcells:
        repcell_definitions = []
        with open(repcelldeffsfile, "r") as f:
            for line in f:
                repcell_definitions.append(str(line.strip()))

    ############################################################################################
    # pics folder
    if not os.path.exists(picfolder):
        os.makedirs(picfolder)

    # colormap
    colormap = plt.get_cmap('cool')

    ############################################# cell_definitions_sorted #################################
    cell_definitions_sorted_set = sort_labels(set(cell_definitions_orig))
    # print('cell_definitions:', cell_definitions)
    print('cell_definitions_sorted_set:', cell_definitions_sorted_set)

    # umap and PCA shapes and sample names
    if save_UMAPs_repcells or save_PCAs_repcells or save_UMAPs_orig or save_PCAs_orig:
        if numerical_shapes:
            shapes = ["$\mathsf{" + str(i) + "}$" for i in cell_definitions_sorted_set]
            with_legend = False
        else:
            shapes = ['8', 's', "^", '*', 'X', 'D', 'P', '1', '2', '3', '4']

    deff_to_shape_dict = dict(zip(cell_definitions_sorted_set, shapes[:len(set(cell_definitions_orig))]))
    # print('deff_to_shape_dict:', deff_to_shape_dict)

    #######################################################################################################

    if impute_method == 'no_imputation':
        deffs = cell_definitions_orig
        curr_data = data
    else:
        deffs = repcell_definitions
        curr_data = repcells_data

    for i in sigtable.index:
        # reloading matplotlib once in a while solves the following issue in Windows:
        # after several iterations the error "failed to allocate bitmap" appears
        importlib.reload(matplotlib)
        matplotlib.use('Agg')

        # For every structure, save a figure and a link to the figure
        print('structure ' + str(i))
        sets_in_struct = sigtable.loc[i, samp_name + '_pairs']
        # print(sets_in_struct)
        sets_in_struct = [tuple(list(map(int, s.strip(')').strip('(').split(', ')))) for s in
                          sets_in_struct.strip(']').strip('[').split('), (')]
        # print('sets_in_struct:', sets_in_struct)

        ################################## count number of pairs (=columns in S, =arrows) each cell has ###############
        l = list(itertools.chain.from_iterable(sets_in_struct))
        num_arrows_per_cell = dict([[x, l.count(x)] for x in set(l)])
        # print('num_arrows_per_cell:', num_arrows_per_cell)

        ################################## Decide on clustering to levels of expression #####################
        list_of_nodes_in_levels, node_level_dict, node_level_name_dict = cluster_to_levels_based_on_sets_in_struct(
            range(curr_data.shape[1]), sets_in_struct)

        if save_GUMAPs_repcells or save_GPCAs_repcells or save_GUMAPs_orig or save_GPCAs_orig or save_Gspatial or save_Gnetwork_layouts:
            if len(sigtable.loc[i, 'genes']) == 32767:  # The genelist was too long to be written in one excel cell
                structs_file = structs_filename(data_path=data_path, impute_method=impute_method,
                                                num_stds_thresh=sigtable.loc[i, 'num_stds_thresh'],
                                                mu=sigtable.loc[i, 'mu'],
                                                path_len=sigtable.loc[i, 'path_len'],
                                                num_iters=sigtable.loc[i, 'num_iters'])
                structs = read_struct_list_from_file(structs_file)
                genelist = list(genes_table[structs[sigtable.loc[i, 'old_struct_num']][0]].index)
            else:  # read the gene list from the excel file
                genelist = read_gene_list(sigtable.loc[i, 'genes'])

        ################################## round p-value for title ############################
        if round_pvals:
            origp = sigtable.loc[i, 'log10_corrected_pval']
            if origp > -100:
                p = int(np.around(origp, 0))
            elif origp > -1000:
                p = int(np.around(np.ceil(origp / 10) * 10, -1))
            else:
                p = int(np.around(np.ceil(origp / 100) * 100, -2))
        else:
            p = np.around(sigtable.loc[i, 'log10_corrected_pval'], 1)

        ################################## title ##############################################
        if impute_method == 'no_imputation':
            title = ('Algorithm parameters: ' + r'$\alpha$' + '=' + str(sigtable.loc[i, 'num_stds_thresh']) + ', ' +
                     r'$\mu$' + '=' + str(sigtable.loc[i, 'mu']) + ', ' +
                     'L=' + str(sigtable.loc[i, 'path_len']) + ', ' +
                     'T=' + str(sigtable.loc[i, 'num_iters']) +
                     '\nStructure scores: p-value=1E' + str(p) +
                     ', ' + r'$\widetilde{\sigma}$' + '=' + str(
                        np.around(sigtable.loc[i, 'structure_average_std'], 2)) +
                     '\n' + str(int(sigtable.loc[i, 'num_genes_in_struct'])) + '/' + str(
                        sigtable.loc[i, 'num_genes']) +
                     ' genes, ' + str(int(sigtable.loc[i, 'num_' + samp_name + 's_in_struct'])) + '/' +
                     str(sigtable.loc[i, 'num_' + samp_name + 's']) + ' ' + samp_name + 's')
        else:
            title = ('Algorithm parameters: ' + r'$\alpha$' + '=' + str(sigtable.loc[i, 'num_stds_thresh']) + ', ' +
                     r'$\mu$' + '=' + str(sigtable.loc[i, 'mu']) + ', ' +
                     'L=' + str(sigtable.loc[i, 'path_len']) + ', ' +
                     'T=' + str(sigtable.loc[i, 'num_iters']) +
                     '\nStructure scores: p-value=1E' + str(p) +
                     ', ' + r'$\widetilde{\sigma}$' + '=' + str(
                        np.around(sigtable.loc[i, 'structure_average_std'], 2)) +
                     '\n' + str(int(sigtable.loc[i, 'num_genes_in_struct'])) + '/' + str(
                        sigtable.loc[i, 'num_genes']) + ' genes, ' +
                     str(int(sigtable.loc[i, 'num_' + real_samp_name + 's_in_struct'])) + '/' +
                     str(int(sigtable.loc[i, 'num_' + real_samp_name + 's'])) + ' ' + real_samp_name + 's, ' +
                     str(int(sigtable.loc[i, 'num_' + samp_name + 's_in_struct'])) + '/' +
                     str(sigtable.loc[i, 'num_' + samp_name + 's']) + ' ' + samp_name + 's')

        ############################### save network layout (regular and GRADUAL) ########################################
        def compute_pos_for_network_layout():
            # multipartite layout
            pos = dict()
            level = 0
            s = list_of_nodes_in_levels[0]
            while s:
                # print('level:', level)
                # print('nodes in the current level:', s)
                # print('definitions of nodes in the current level:', np.array(cell_definitions)[list(s)])

                s_sorted = []
                all_defs_in_s = list(set(list(np.array(deffs)[list(s)])))
                # print('all_defs_in_s:', all_defs_in_s)
                lens = list(set([len(r) for r in all_defs_in_s]))
                # print('lens:', lens)
                for l in lens:
                    b = list(set([r for r in all_defs_in_s if len(r) == l]))
                    b.sort()
                    for defb in b:
                        s_sorted += [r for r in s if np.array(deffs)[r] == defb]

                pos.update((n, (level + 1, 2 * j + np.mod(level, 2))) for j, n in enumerate(s_sorted))
                level += 1
                s = list_of_nodes_in_levels[level]

                # print('list_of_nodes_in_levels:', list_of_nodes_in_levels)
                # print('pos:', pos)

                # pos = nx.layout.spring_layout(G, iterations=200)
                # pos = nx.spectral_layout(G)
                # pos = nx.spiral_layout(G)
                # pos = nx.random_layout(G)
                # pos = nx.kamada_kawai_layout(G)
                # pos = nx.circular_layout(G)
                # pos = nx.planar_layout(G)
                # pos = nx.fruchterman_reingold_layout(G)

            return pos

        if save_network_layouts or save_Gnetwork_layouts:
            print('Computing and saving network layout')

            # create a graph of the structure
            G = nx.DiGraph()
            G.add_edges_from(sets_in_struct)
            # print('G.nodes:', G.nodes)

            pos = compute_pos_for_network_layout()

            if save_layers_to_excel:
                sigtable.loc[i, 'layers'] = str([s for s in list_of_nodes_in_levels if len(s) != 0])

            if save_network_layouts:
                fig, ax = plt.subplots(figsize=(15, 15))
                ax.set_title(title, fontdict={'fontsize': 30, 'fontweight': 'medium'})

                if len(set(deffs)) == 1:
                    # print('labels:', {list(G.nodes)[i]:str(cell_labels[G.nodes][i]) for i in range(len(G.nodes))})
                    nx.draw(G, pos, node_color=[node_level_dict[node] for node in G.nodes],
                            node_size=700 + np.array([50 * (num_arrows_per_cell[node] - 1) for node in G.nodes]),
                            cmap=plt.cm.cool, with_labels=True, arrows=True,
                            labels={node: node for node in G.nodes})
                else:
                    for deff, shape in zip(cell_definitions_sorted_set, shapes[:len(set(deffs))]):
                        inds_deff = np.where(np.array(deffs) == deff)[0]
                        nodelist = list(set(inds_deff).intersection(set(G.nodes)))
                        # print('node_shape', shape)
                        # print('node_color:',[node_level_dict[node] for node in nodelist])
                        nx.draw_networkx_nodes(G, pos, nodelist=nodelist,
                                               node_size=700 + np.array(
                                                   [50 * (num_arrows_per_cell[node] - 1) for node in nodelist]),
                                               node_color=[node_level_dict[node] for node in nodelist],
                                               cmap=plt.cm.cool, vmin=0, vmax=max(node_level_dict.values()),
                                               node_shape=shape)
                    nx.draw_networkx_edges(G, pos, arrows=True)
                    nx.draw_networkx_labels(G, pos, labels={node: deffs[node] for node in list(G.nodes)})

                    # if data_name=='yaeldiffbulk2'
                    # nx.draw(G, pos, node_color=[node_level_dict[node] for node in G.nodes],
                    #        node_size=700+np.array([50*(num_arrows_per_cell[node]-1) for node in G.nodes]),
                    #        cmap=plt.cm.cool, with_labels=True, arrows=True,
                    #        #labels={list(G.nodes)[i]:str(np.array(cell_definitions)[G.nodes][i]) for i in range(len(G.nodes))})
                    #        labels={node:repcell_definitions[node] for node in G.nodes})

                # picture file name
                picname = 'struct_' + str(i) + '_net.jpg'

                # save the figure to file
                picfile = os.path.join(picfolder, picname)
                plt.tight_layout()
                plt.savefig(picfile)
                # plt.clf()
                plt.close()

                if with_links_in_sigtable:
                    # add a hyperlink to the pic file in sigtable
                    comp_path = picfile.replace('\\', '/')
                    hyperlink_txt = '=hyperlink("' + comp_path + '","network")'
                    sigtable.loc[i, 'networks'] = hyperlink_txt

            if save_Gnetwork_layouts:

                genelist_ = list(set(genelist).intersection(set(curr_data.index)))
                if len(genelist_) < len(genelist):
                    print('missing genes:', set(genelist) - set(genelist_))

                vmin = curr_data.loc[genelist_, :].mean().min()
                vmax = curr_data.loc[genelist_, :].mean().max()

                fig, ax = plt.subplots(figsize=(15, 15))
                ax.set_title(title, fontdict={'fontsize': 30, 'fontweight': 'medium'})

                cdict = dict(
                    zip(list(set(G.nodes)), list(curr_data.iloc[:, list(set(G.nodes))].loc[genelist_, :].mean())))
                if len(set(deffs)) == 1:
                    # print('labels:', {list(G.nodes)[i]:str(cell_labels[G.nodes][i]) for i in range(len(G.nodes))})
                    nx.draw(G, pos, node_color=[cdict[node] for node in G.nodes],
                            node_size=700 + np.array([50 * (num_arrows_per_cell[node] - 1) for node in G.nodes]),
                            cmap=plt.cm.cool, with_labels=True, arrows=True,
                            labels={node: node for node in G.nodes}, vmin=vmin, vmax=vmax)
                else:
                    for deff, shape in zip(cell_definitions_sorted_set, shapes[:len(set(deffs))]):
                        inds_deff = np.where(np.array(deffs) == deff)[0]
                        nodelist = list(set(inds_deff).intersection(set(G.nodes)))
                        # print('node_shape', shape)
                        # print('node_color:',[node_level_dict[node] for node in nodelist])
                        nx.draw_networkx_nodes(G, pos, nodelist=nodelist,
                                               node_size=700 + np.array(
                                                   [50 * (num_arrows_per_cell[node] - 1) for node in nodelist]),
                                               node_color=[cdict[node] for node in nodelist],
                                               cmap=plt.cm.cool,
                                               node_shape=shape, vmin=vmin, vmax=vmax)
                    nx.draw_networkx_edges(G, pos, arrows=True)
                    nx.draw_networkx_labels(G, pos, labels={node: deffs[node] for node in list(G.nodes)})

                    # if data_name=='yaeldiffbulk2'
                    # nx.draw(G, pos, node_color=[cdict[node] for node in G.nodes],
                    #        node_size=700+np.array([50*(num_arrows_per_cell[node]-1) for node in G.nodes]),
                    #        cmap=plt.cm.cool, with_labels=True, arrows=ordered_pairs,
                    #        #labels={list(G.nodes)[i]:str(np.array(cell_definitions)[G.nodes][i]) for i in range(len(G.nodes))})
                    #        labels={node: cell_definitions[node] for node in G.nodes}, vmin=vmin, vmax=vmax)

                sm = plt.cm.ScalarMappable(cmap=colormap, norm=plt.Normalize(vmin=vmin, vmax=vmax))
                sm.set_array([])
                cbar = plt.colorbar(sm)
                cbar.set_label('average expression level', size=20)

                # picture file name
                picname = 'struct_' + str(i) + '_Gnet.jpg'

                # save the figure to file
                picfile = os.path.join(picfolder, picname)
                plt.tight_layout()
                plt.savefig(picfile)
                # plt.clf()
                plt.close()

                if with_links_in_sigtable:
                    # add a hyperlink to the pic file in sigtable
                    comp_path = picfile.replace('\\', '/')
                    hyperlink_txt = '=hyperlink("' + comp_path + '","Gnetwork")'
                    sigtable.loc[i, 'Gnetworks'] = hyperlink_txt

        ############################################# UMAPs\PCAs\spatial ########################################
        if save_UMAPs_repcells or save_PCAs_repcells or save_UMAPs_orig or save_PCAs_orig or save_spatial:
            print('Computing and saving UMAP\PCA\spatial layout')

            levels = np.array([node_level_name_dict[n] for n in range(curr_data.shape[1])])

            if (save_UMAPs_orig or save_PCAs_orig or save_spatial) and impute_method != 'no_imputation':
                # Compute levels_orig
                nLevels = np.max(list(node_level_dict.values())) + 1
                for node in set(cells_to_repcells['x']) - set(node_level_dict.keys()):
                    node_level_dict[node] = 'not in structure'
                    node_level_name_dict[node] = 'not in structure'
                cell_to_level = {i: node_level_dict[cells_to_repcells.loc[i, 'x']] for i in cells_to_repcells.index}
                levels_orig = [node_level_name_dict[int(cells_to_repcells.loc[i, 'x'])] for i in
                               cells_to_repcells.index]
                cells_to_repcells['levels'] = levels_orig

            # umap.plot.points(mapper, labels=np.array(levels), color_key_cmap='Paired')
            level_desc_dict = dict()
            if 'not in structure' in set(levels):
                # print('"not in structure" in set(levels)')
                level_desc_dict['not in structure'] = 'not in structure'
                nLevels = len(set(levels)) - 1
                curr_colors = list(colormap(np.linspace(0, 1.0, nLevels))) + ['tab:gray']
                # print('curr_colors:', curr_colors)
            else:
                # print('"not in structure" NOT in set(levels)')
                nLevels = len(set(levels))
                curr_colors = colormap(np.linspace(0, 1.0, nLevels))
                # print('curr_colors:', curr_colors)

            if nLevels == 2:
                level_desc_dict['level 1'] = 'low expression'
                level_desc_dict['level 2'] = 'high expression'
                # curr_colors = ['g','r']
            elif nLevels == 3:
                level_desc_dict['level 1'] = 'low expression'
                level_desc_dict['level 2'] = 'medium expression'
                level_desc_dict['level 3'] = 'high expression'
                # curr_colors = ['g','b','r']
            elif nLevels == 4:
                level_desc_dict['level 1'] = 'low expression'
                level_desc_dict['level 2'] = 'medium expression'
                level_desc_dict['level 3'] = 'high expression'
                level_desc_dict['level 4'] = 'very high expression'
                # curr_colors = ['g','c','b','r']
            elif nLevels == 5:
                level_desc_dict['level 1'] = 'very low expression'
                level_desc_dict['level 2'] = 'low expression'
                level_desc_dict['level 3'] = 'medium expression'
                level_desc_dict['level 4'] = 'high expression'
                level_desc_dict['level 5'] = 'very high expression'
                # curr_colors = ['g','c','b','m','r']
            else:  # more than 5 levels
                level_desc_dict['level 1'] = 'lowest expression'
                for l in np.arange(1, nLevels - 1):
                    level_desc_dict['level ' + str(l + 1)] = 'level ' + str(l + 1)
                level_desc_dict['level ' + str(nLevels)] = 'highest expression'

            '''
            elif nLevels==6:
                level_desc_dict['level 1'] = 'very low expression'
                level_desc_dict['level 2'] = 'low expression'
                level_desc_dict['level 3'] = 'medium expression'
                level_desc_dict['level 4'] = 'high expression'
                level_desc_dict['level 5'] = 'very high expression'
                level_desc_dict['level 6'] = 'highest expression'
                curr_colors = ['g','c','b','m','r','k']
            elif nLevels==7:
                level_desc_dict['level 1'] = 'lowest expression'
                level_desc_dict['level 2'] = 'very low expression'
                level_desc_dict['level 3'] = 'low expression'
                level_desc_dict['level 4'] = 'medium expression'
                level_desc_dict['level 5'] = 'high expression'
                level_desc_dict['level 6'] = 'very high expression'
                level_desc_dict['level 7'] = 'highest expression'
                curr_colors = ['g','c','b','m','r','#A52A2A','k']
            elif nLevels==8:
                level_desc_dict['level 1'] = 'lowest expression'
                level_desc_dict['level 2'] = 'very low expression'
                level_desc_dict['level 3'] = 'low expression'
                level_desc_dict['level 4'] = 'medium expression'
                level_desc_dict['level 5'] = 'high expression'
                level_desc_dict['level 6'] = 'very high expression'
                level_desc_dict['level 7'] = 'highest expression'
                curr_colors = ['g','c','b','#00008B','m','r','#A52A2A','k']
            '''

            if impute_method != 'no_imputation':
                repcell_table = pd.DataFrame(columns=['levels', 'deffs'])
                repcell_table.loc[:, 'levels'] = levels
                repcell_table.loc[:, 'deffs'] = repcell_definitions
            else:
                orig_table = pd.DataFrame(columns=['levels', 'deffs'])
                orig_table.loc[:, 'levels'] = levels
                orig_table.loc[:, 'deffs'] = cell_definitions_orig

            if save_UMAPs_repcells:
                # picture file name
                picname = 'struct_' + str(i) + '_repcells_UMAP.jpg'

                # save the figure to file
                picfile = os.path.join(picfolder, picname)
                save_struct_layout(cell_table=repcell_table, picfile=picfile,
                                   cell_definitions_sorted_set=cell_definitions_sorted_set,
                                   curr_colors=curr_colors, axislabel='UMAP', coor=umap_coor, shapes=shapes,
                                   title=title, level_desc_dict=level_desc_dict,
                                   size_by_degree=False, num_arrows_per_cell=None, with_legend=with_legend,
                                   numerical_shapes=numerical_shapes)

                if with_links_in_sigtable:
                    # add a hyperlink to the pic file in sigtable
                    comp_path = picfile.replace('\\', '/')
                    hyperlink_txt = '=hyperlink("' + comp_path + '","repcells_UMAP")'
                    sigtable.loc[i, 'repcells_UMAPs'] = hyperlink_txt

            if save_PCAs_repcells:
                # picture file name
                picname = 'struct_' + str(i) + '_repcells_PCA.jpg'

                # save the figure to file
                picfile = picfolder + '/' + picname
                save_struct_layout(cell_table=repcell_table, picfile=picfile,
                                   cell_definitions_sorted_set=cell_definitions_sorted_set,
                                   curr_colors=curr_colors, axislabel='PC', coor=pca_coor, shapes=shapes,
                                   title=title, level_desc_dict=level_desc_dict,
                                   size_by_degree=False, num_arrows_per_cell=None, with_legend=with_legend,
                                   numerical_shapes=numerical_shapes)

                if with_links_in_sigtable:
                    # add a hyperlink to the pic file in sigtable
                    comp_path = picfile.replace('\\', '/')
                    hyperlink_txt = '=hyperlink("' + comp_path + '","repcells_PCA")'
                    sigtable.loc[i, 'repcells_PCAs'] = hyperlink_txt

            if save_UMAPs_orig:
                # picture file name
                picname = 'struct_' + str(i) + '_orig_UMAP.jpg'

                # save the figure to file
                picfile = os.path.join(picfolder, picname)
                if impute_method != 'no_imputation':
                    table = cells_to_repcells
                else:
                    table = orig_table

                save_struct_layout(cell_table=table, picfile=picfile,
                                   cell_definitions_sorted_set=cell_definitions_sorted_set,
                                   curr_colors=curr_colors, axislabel='UMAP', coor=umap_coor_orig, shapes=shapes,
                                   title=title, level_desc_dict=level_desc_dict,
                                   size_by_degree=False, num_arrows_per_cell=None, with_legend=with_legend,
                                   numerical_shapes=numerical_shapes)

                if with_links_in_sigtable:
                    # add a hyperlink to the pic file in sigtable
                    comp_path = picfile.replace('\\', '/')
                    hyperlink_txt = '=hyperlink("' + comp_path + '","orig_UMAP")'
                    sigtable.loc[i, 'orig_UMAPs'] = hyperlink_txt

            if save_PCAs_orig:
                # picture file name
                picname = 'struct_' + str(i) + '_orig_PCA.jpg'

                # save the figure to file
                picfile = picfolder + '/' + picname
                if impute_method != 'no_imputation':
                    table = cells_to_repcells
                else:
                    table = orig_table

                save_struct_layout(cell_table=table, picfile=picfile,
                                   cell_definitions_sorted_set=cell_definitions_sorted_set,
                                   curr_colors=curr_colors, axislabel='PC', coor=pca_coor_orig, shapes=shapes,
                                   title=title, level_desc_dict=level_desc_dict,
                                   size_by_degree=False, num_arrows_per_cell=None, with_legend=with_legend,
                                   numerical_shapes=numerical_shapes)

                if with_links_in_sigtable:
                    # add a hyperlink to the pic file in sigtable
                    comp_path = picfile.replace('\\', '/')
                    hyperlink_txt = '=hyperlink("' + comp_path + '","orig_PCA")'
                    sigtable.loc[i, 'orig_PCAs'] = hyperlink_txt

            if save_spatial:
                # picture file name
                picname = 'struct_' + str(i) + '_spatial.jpg'

                # save the figure to file
                picfile = os.path.join(picfolder, picname)
                if impute_method != 'no_imputation':
                    table = cells_to_repcells
                else:
                    table = orig_table

                save_struct_layout(cell_table=table, picfile=picfile,
                                   cell_definitions_sorted_set=cell_definitions_sorted_set,
                                   curr_colors=curr_colors, axislabel='spatial axis', coor=spatial_coor, shapes=shapes,
                                   title=title, level_desc_dict=level_desc_dict,
                                   size_by_degree=False, num_arrows_per_cell=None, with_legend=with_legend,
                                   numerical_shapes=numerical_shapes)

                if with_links_in_sigtable:
                    # add a hyperlink to the pic file in sigtable
                    comp_path = picfile.replace('\\', '/')
                    hyperlink_txt = '=hyperlink("' + comp_path + '","spatial")'
                    sigtable.loc[i, 'spatial'] = hyperlink_txt

            ############################################# GRADUAL UMAPs\PCAs\spatial ########################################
            print('Computing and saving gradual UMAP\PCA\spatial layout')

            if save_GUMAPs_repcells:
                # picture file name
                picname = 'struct_' + str(i) + '_repcells_GUMAP.jpg'

                # save the figure to file
                picfile = os.path.join(picfolder, picname)
                save_layout_gradual(cell_table=repcell_table, picfile=picfile, data=repcells_data,
                                    genelist=genelist, cell_definitions_sorted_set=cell_definitions_sorted_set,
                                    colormap=colormap, axislabel='UMAP', coor=umap_coor, shapes=shapes, title=title,
                                    color_by_log_expression=color_by_log_expression,
                                    numerical_shapes=numerical_shapes)

                if with_links_in_sigtable:
                    # add a hyperlink to the pic file in sigtable
                    comp_path = picfile.replace('\\', '/')
                    hyperlink_txt = '=hyperlink("' + comp_path + '","repcells_Gumap")'
                    sigtable.loc[i, 'repcells_Gumaps'] = hyperlink_txt

            if save_GPCAs_repcells:
                # picture file name
                picname = 'struct_' + str(i) + '_repcells_GPCA.jpg'

                # save the figure to file
                picfile = os.path.join(picfolder, picname)
                save_layout_gradual(cell_table=repcell_table, picfile=picfile, data=repcells_data,
                                    genelist=genelist, cell_definitions_sorted_set=cell_definitions_sorted_set,
                                    colormap=colormap, axislabel='PC', coor=pca_coor, shapes=shapes, title=title,
                                    color_by_log_expression=color_by_log_expression,
                                    numerical_shapes=numerical_shapes)

                if with_links_in_sigtable:
                    # add a hyperlink to the pic file in sigtable
                    comp_path = picfile.replace('\\', '/')
                    hyperlink_txt = '=hyperlink("' + comp_path + '","repcells_GPCA")'
                    sigtable.loc[i, 'repcells_GPCAs'] = hyperlink_txt

            if save_GUMAPs_orig:
                # picture file name
                picname = 'struct_' + str(i) + '_orig_GUMAP.jpg'
                if impute_method != 'no_imputation':
                    table = cells_to_repcells
                else:
                    table = orig_table

                # save the figure to file
                picfile = picfolder + '/' + picname
                save_layout_gradual(cell_table=table, picfile=picfile,
                                    data=data, genelist=genelist,
                                    cell_definitions_sorted_set=cell_definitions_sorted_set, colormap=colormap,
                                    axislabel='UMAP', coor=umap_coor_orig, shapes=shapes, title=title,
                                    color_by_log_expression=color_by_log_expression,
                                    numerical_shapes=numerical_shapes)

                if with_links_in_sigtable:
                    # add a hyperlink to the pic file in sigtable
                    comp_path = picfile.replace('\\', '/')
                    hyperlink_txt = '=hyperlink("' + comp_path + '","orig_GUMAP")'
                    sigtable.loc[i, 'orig_GUMAPs'] = hyperlink_txt

            if save_GPCAs_orig:
                # picture file name
                picname = 'struct_' + str(i) + '_orig_GPCA.jpg'
                if impute_method != 'no_imputation':
                    table = cells_to_repcells
                else:
                    table = orig_table

                # save the figure to file
                picfile = picfolder + '/' + picname
                # save_layout_gradual(scRNA_data_orig_norm,axislabel='PC',coor=pca_coor_orig,cell_table=cells_to_repcells,picfile=picfile)
                save_layout_gradual(cell_table=table, picfile=picfile, data=data,
                                    genelist=genelist, cell_definitions_sorted_set=cell_definitions_sorted_set,
                                    colormap=colormap, axislabel='PC', coor=pca_coor_orig, shapes=shapes, title=title,
                                    color_by_log_expression=color_by_log_expression,
                                    numerical_shapes=numerical_shapes)

                if with_links_in_sigtable:
                    # add a hyperlink to the pic file in sigtable
                    comp_path = picfile.replace('\\', '/')
                    hyperlink_txt = '=hyperlink("' + comp_path + '","orig_GPCA")'
                    sigtable.loc[i, 'orig_GPCAs'] = hyperlink_txt

            if save_Gspatial:
                # picture file name
                picname = 'struct_' + str(i) + '_Gspatial.jpg'
                if impute_method != 'no_imputation':
                    table = cells_to_repcells
                else:
                    table = orig_table

                # save the figure to file
                picfile = os.path.join(picfolder, picname)
                save_layout_gradual(cell_table=table, picfile=picfile,
                                    data=data, genelist=genelist,
                                    cell_definitions_sorted_set=cell_definitions_sorted_set, colormap=colormap,
                                    axislabel='spatial axis', coor=spatial_coor, shapes=shapes, title=title,
                                    color_by_log_expression=color_by_log_expression,
                                    numerical_shapes=numerical_shapes)

                if with_links_in_sigtable:
                    # add a hyperlink to the pic file in sigtable
                    comp_path = picfile.replace('\\', '/')
                    hyperlink_txt = '=hyperlink("' + comp_path + '","Gspatial")'
                    sigtable.loc[i, 'Gspatial'] = hyperlink_txt

    if save_table:
        sigtable.to_excel(sigfile_vis)

    return sigtable
