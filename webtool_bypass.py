#!/var/www/html/SPIRAL.web.tool/spiral_venv/bin/python3.9

from SPIRAL_pipeline_funcs import *

ANALYSIS_FOLDER = './static/analysis'
if __name__ == '__main__':
    for data_n in [5000]:
        if os.path.isdir(os.path.join(ANALYSIS_FOLDER, 'data' + str(data_n))):
            try:
                print('\n\n\n\n\n\n\n', data_n)
                if data_n in [3, 4, 6]:
                    numerical_shapes = True
                else:
                    numerical_shapes = False
                analysis_folder = ANALYSIS_FOLDER
                data_path = data_path_name(analysis_folder, data_n)
                species = open(os.path.join(data_path, 'species.txt'), "r").read()
                #load_data_first_time(analysis_folder=ANALYSIS_FOLDER, data_n=data_n, median_count_normalization_flag=True,
                #                     with_log=False)
                # compute_violin_plots(analysis_folder=ANALYSIS_FOLDER, data_n=data_n, static_path='./static',
                #                     species=species)
                outcome = run_SPIRAL_pipeline(analysis_folder=ANALYSIS_FOLDER, data_n=data_n, species=None,
                                              min_nFeatures=None, max_nFeatures=None, max_mtpercent=None,
                                              numerical_shapes=numerical_shapes,
                                              num_stds_thresh_lst=[0.5, 0.75, 1.0], mu_lst=[0.9, 0.95],
                                              num_iters_lst=[10000],
                                              path_len_lst=[3],
                                              save_layers_to_excel=True)
                with open(outcome_path(data_path), 'w') as text_file:
                    text_file.write(str(outcome))
            except Exception as error:
                print(error)
                outcome = 'Error'
                with open(outcome_path(data_path), 'w') as text_file:
                    text_file.write(str(outcome))
