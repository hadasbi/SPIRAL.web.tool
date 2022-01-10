from SPIRAL_pipeline_funcs import *

ANALYSIS_FOLDER = './static/analysis'
# for data_n in range(120, 121):
# for data_n in range(121, 148):
for data_n in [122]:
    if os.path.isdir(os.path.join(ANALYSIS_FOLDER, 'data' + str(data_n))):
        try:
            print('\n\n\n\n\n\n\n', data_n)
            # load_data_first_time(analysis_folder=ANALYSIS_FOLDER, data_n=data_n, median_count_normalization_flag=True, with_log=False)
            run_SPIRAL_pipeline(analysis_folder=ANALYSIS_FOLDER, data_n=data_n, species=None,
                                min_nFeatures=None, max_nFeatures=None, max_mtpercent=None,
                                num_stds_thresh_lst=[1.0], mu_lst=[0.95], num_iters_lst=[10000], path_len_lst=[3])
        except:
            pass
