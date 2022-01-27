import pandas as pd
from SPIRAL_basic_funcs import *

dataloc = './static/analysis/data150/'
sigtable = pd.read_excel(dataloc + 'agg_wald_sigtable_filt_GO_vis.xlsx', index_col=0)
samp_name = 'repcell'

for i in sigtable.index[:15]:
    print('\n')
    for j in sigtable.index[:15]:
        genes_i = read_gene_list(sigtable.loc[i, 'genes'])
        genes_j = read_gene_list(sigtable.loc[j, 'genes'])
        jaccard_genes_i_j = len(set(genes_i).intersection(set(genes_j))) / len(set(genes_i).union(set(genes_j)))

        ###############
        sets_in_struct_i = sigtable.loc[i, samp_name + '_pairs']
        sets_in_struct_i = [tuple(list(map(int, s.strip(')').strip('(').split(', ')))) for s in
                            sets_in_struct_i.strip(']').strip('[').split('), (')]

        sets_in_struct_j = sigtable.loc[j, samp_name + '_pairs']
        sets_in_struct_j = [tuple(list(map(int, s.strip(')').strip('(').split(', ')))) for s in
                            sets_in_struct_j.strip(']').strip('[').split('), (')]

        jaccard_sets_i_j = len(set(sets_in_struct_i).intersection(set(sets_in_struct_j))) / len(
            set(sets_in_struct_i).union(set(sets_in_struct_j)))

        ###############
        low_repcells_in_struct_i = set([s[0] for s in sets_in_struct_i])
        low_repcells_in_struct_j = set([s[0] for s in sets_in_struct_j])
        jaccard_low_repcells_i_j = len(low_repcells_in_struct_i.intersection(low_repcells_in_struct_j)) / len(
            low_repcells_in_struct_i.union(low_repcells_in_struct_j))

        ###############
        high_repcells_in_struct_i = set([s[1] for s in sets_in_struct_i])
        high_repcells_in_struct_j = set([s[1] for s in sets_in_struct_j])
        jaccard_high_repcells_i_j = len(high_repcells_in_struct_i.intersection(high_repcells_in_struct_j)) / len(
            high_repcells_in_struct_i.union(high_repcells_in_struct_j))

        if jaccard_high_repcells_i_j > 0.3:
            print(i, j, jaccard_genes_i_j, jaccard_sets_i_j,
                  (jaccard_genes_i_j + jaccard_sets_i_j) / 2, jaccard_low_repcells_i_j, jaccard_high_repcells_i_j)
