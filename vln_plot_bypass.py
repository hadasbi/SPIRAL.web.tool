from SPIRAL_pipeline_funcs import *

ANALYSIS_FOLDER = './static/analysis'
static_path = './static'
species = "MUS_MUSCULUS"
for data_n in range(126, 148):
    print('\n\n\n\n\n\n\n', data_n)
    try:
        load_data_first_time(analysis_folder=ANALYSIS_FOLDER, data_n=data_n, median_count_normalization_flag=True,
                             with_log=False)
        compute_violin_plots(analysis_folder=ANALYSIS_FOLDER, data_n=data_n, static_path=static_path, species=species)
    except:
        pass


