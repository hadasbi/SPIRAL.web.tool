#!/var/www/html/SPIRAL.web.tool/spiral_venv/bin/python3.9

import numpy as np


def find_submat_around_g(g, mu, abs_genes_in_sets):
    rows_g = [g]
    cols_g = list(abs_genes_in_sets.nonzero()[1][np.where(abs_genes_in_sets.nonzero()[0] == g)[0]])
    converged = False
    while not converged:
        new_rows_g = list(np.where(np.ravel(np.sum(abs_genes_in_sets[:, cols_g], axis=1)) >= mu * len(cols_g))[0])
        new_cols_g = list(
            np.where(np.ravel(np.sum(abs_genes_in_sets[new_rows_g, :], axis=0)) >= mu * len(new_rows_g))[0])
        if new_rows_g != rows_g or new_cols_g != cols_g:
            rows_g = new_rows_g
            cols_g = new_cols_g
        else:
            converged = True
    return (new_rows_g, new_cols_g)


def find_submat_around_g_no_iterations(g, mu, abs_genes_in_sets):
    rows_g = [g]
    cols_g = list(abs_genes_in_sets.nonzero()[1][np.where(abs_genes_in_sets.nonzero()[0] == g)[0]])
    new_rows_g = list(np.where(np.ravel(np.sum(abs_genes_in_sets[:, cols_g], axis=1)) >= mu * len(cols_g))[0])
    return (new_rows_g, cols_g)


def find_submat_around_random_cellpair_set(mu, abs_genes_in_sets):
    # remove the possibility that both the pair (i,j) and the pair (j,i) would be picked for the seed (for any i,j)
    num_cells = int(np.sqrt(abs_genes_in_sets.shape[1]))
    cell_pairs = [divmod(x, num_cells) for x in range(abs_genes_in_sets.shape[1])]
    cell_pairs = [x for x in cell_pairs if x[0] != x[1]]
    rand_order = np.random.choice(np.array(range(num_cells)), size=num_cells)
    cell_pairs = [x for x in cell_pairs if rand_order[x[0]] < rand_order[x[1]]]
    columns = [x[0] * num_cells + x[1] for x in cell_pairs]

    # randomly pick a seed of cell-pairs
    cols_g = list(np.random.choice(columns, size=max(int(np.round(abs_genes_in_sets.shape[1] / 100)), 10)))
    rows_g = list(np.where(np.ravel(np.sum(abs_genes_in_sets[:, cols_g], axis=1)) >= mu * len(cols_g))[0])

    converged = False
    while not converged:
        new_cols_g = list(np.where(np.ravel(np.sum(abs_genes_in_sets[rows_g, :], axis=0)) >= mu * len(rows_g))[0])
        new_rows_g = list(
            np.where(np.ravel(np.sum(abs_genes_in_sets[:, new_cols_g], axis=1)) >= mu * len(new_cols_g))[0])
        if new_rows_g != rows_g or new_cols_g != cols_g:
            cols_g = new_cols_g
            rows_g = new_rows_g
        else:
            converged = True
    return (new_rows_g, new_cols_g)


def find_submat_around_random_path(mu, abs_genes_in_sets, path_len, seed):
    num_cells = int(np.sqrt(abs_genes_in_sets.shape[1]))

    # randomly pick (path_len+1) cells
    rng = np.random.default_rng(seed)
    cells_in_path = rng.choice(range(num_cells), size=path_len + 1, replace=False)
    cols_g = []
    for i, c1 in enumerate(cells_in_path):
        for c2 in cells_in_path[(i + 1):]:
            cols_g.append(c1 * num_cells + c2)

    rows_g = list(np.where(np.ravel(np.sum(abs_genes_in_sets[:, cols_g], axis=1)) >= mu * len(cols_g))[0])

    converged = False
    while not converged:
        new_cols_g = list(np.where(np.ravel(np.sum(abs_genes_in_sets[rows_g, :], axis=0)) >= mu * len(rows_g))[0])
        new_rows_g = list(
            np.where(np.ravel(np.sum(abs_genes_in_sets[:, new_cols_g], axis=1)) >= mu * len(new_cols_g))[0])
        if new_rows_g != rows_g or new_cols_g != cols_g:
            cols_g = new_cols_g
            rows_g = new_rows_g
        else:
            converged = True
    return (new_rows_g, new_cols_g)


def compute_sets(genes_table, scRNA_data, gene_std, cell_pair, num_cells, ordered_pairs):
    c1 = cell_pair[0]
    c2 = cell_pair[1]
    if not ordered_pairs:
        gene_inds1 = np.array(
            list(genes_table.loc[scRNA_data.index[(scRNA_data.iloc[:, c1] - scRNA_data.iloc[:, c2]) > gene_std]]))
        gene_inds2 = np.array(
            list(genes_table.loc[scRNA_data.index[(scRNA_data.iloc[:, c2] - scRNA_data.iloc[:, c1]) > gene_std]]))
        genes_in_sets_row = np.concatenate((gene_inds1, gene_inds2))
        genes_in_sets_col = (np.ones(len(gene_inds1 + len(gene_inds1))) * (c1 * num_cells + c2)).astype(int)
        genes_up_or_down = (np.concatenate((np.ones(len(gene_inds1)) * (-1), np.ones(len(gene_inds2))))).astype(int)
        return (genes_in_sets_row, genes_in_sets_col, genes_up_or_down)
    else:
        gene_inds = np.array(
            list(genes_table.loc[scRNA_data.index[(scRNA_data.iloc[:, c2] - scRNA_data.iloc[:, c1]) > gene_std]]))
        genes_in_sets_col = (np.ones(len(gene_inds)) * (c1 * num_cells + c2)).astype(int)
        return (gene_inds, genes_in_sets_col)


def compute_sets_efficient(genes_table, scRNA_data, gene_std, cell_pair, num_cells):
    c1 = cell_pair[0]
    c2 = cell_pair[1]
    gene_inds = np.array(
        list(genes_table.loc[scRNA_data.index[(scRNA_data.iloc[:, c2] - scRNA_data.iloc[:, c1]) > gene_std]]))
    return (np.concatenate((np.array([c1 * num_cells + c2]), gene_inds)))
