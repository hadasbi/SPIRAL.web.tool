#!/var/www/html/SPIRAL.web.tool/spiral_venv/bin/python3.10

from SPIRAL_pipeline_funcs import *

ANALYSIS_FOLDER = './static/analysis'
static_path = './static'
species = "MUS_MUSCULUS"
# species = "HOMO_SAPIENS"
# species = 'synthetic'
# species = 'DANIO_RERIO'
ensembl_folder = ensembl_loc(production=False, hostname=None)
for data_n in [6]:
    print('\n\n\n\n\n\n\n', data_n)
    load_data_first_time(analysis_folder=ANALYSIS_FOLDER, data_n=data_n, median_count_normalization_flag=True,
                         with_log=False)
    compute_violin_plots(analysis_folder=ANALYSIS_FOLDER, data_n=data_n, static_path=static_path, species=species,
                         ensembl_folder=ensembl_folder)
