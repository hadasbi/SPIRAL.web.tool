import os, io
from time import sleep
from flask import Flask, render_template, request, redirect, url_for, flash, send_from_directory, jsonify
from flask_bootstrap import Bootstrap
from flask_wtf import FlaskForm
from wtforms import StringField, IntegerField, SubmitField, SelectField, SelectMultipleField
from wtforms.validators import DataRequired, Email
from flask_wtf.file import FileField, FileRequired, FileAllowed
from flask_uploads import UploadSet, configure_uploads, patch_request_class, IMAGES
from flask_mail import Mail, Message
import json
import pandas as pd
import plotly
import plotly.express as px
import base64

# from SPIRAL_basic_funcs import *
from SPIRAL_pipeline_funcs import *
# from SPIRAL_visualization_funcs import *
# from SPIRAL_enrichment_funcs import *
# from SPIRAL_design_excel import *

from threading import Thread

# UPLOADS_FOLDER = 'D:/Nextcloud/SPIRAL/uploads/'

app = Flask(__name__)

app.config['MAIL_SERVER'] = 'smtp.gmail.com'
app.config['MAIL_PORT'] = 465
app.config['MAIL_USERNAME'] = 'SPIRAL.web.tool@gmail.com'
app.config['MAIL_PASSWORD'] = '123SPIRAL'
app.config['MAIL_USE_TLS'] = False
app.config['MAIL_USE_SSL'] = True
mail = Mail(app)

app.config['MAX_CONTENT_LENGTH'] = 1024 * 1024 * 1024

ANALYSIS_FOLDER = os.path.join(app._static_folder, 'analysis')

# Flask-WTF requires an encryption key - the string can be anything
app.config['SECRET_KEY'] = 'kncwhgJAVKBJAHFvlv,Klz'
app.config['UPLOADED_TABLES_DEST'] = ANALYSIS_FOLDER

app.config['SERVER_NAME'] = '132.69.181.166:5000'

# SHREK:
# app.config['SERVER_NAME'] = '132.68.107.4:5000'


# app.config['UPLOADED_PHOTOS_DEST'] = ANALYSIS_FOLDER

# app._static_folder = 'D:/Nextcloud/SPIRAL/code/static'
# app.config['UPLOAD_FOLDER'] = 'D:/Nextcloud/SPIRAL/uploads/'
# app.debug = True

# Flask-Bootstrap requires this line
Bootstrap(app)

tables = UploadSet(name='tables', extensions=('txt', 'csv'))
configure_uploads(app, (tables))

# photos = UploadSet('photos', IMAGES)
# configure_uploads(app, (photos))

patch_request_class(app, 1024 * 1024 * 1024)

'''
##########################  test   ###########################################
class NameForm(FlaskForm):
    name = StringField('Which actor is your favorite?', validators=[DataRequired()])
    pic = FileField('Table of actor:',
                    validators=[FileRequired(), FileAllowed(tables,  'csv or txt files only!')]
                    )
#    pic = FileField('Picture of actor:',
#                    validators=[FileRequired(), FileAllowed(photos, 'images only!')]
#                    )
    submit = SubmitField('Submit')

@app.route('/test', methods=['GET', 'POST'])
def test():
    # you must tell the variable 'form' what you named the class, above
    # 'form' is the variable name used in this template: index.html
    form = NameForm(request.form)
    print('form.errors:', form.errors)
    print('form.pic.errors:', form.pic.errors)
    print(form.validate_on_submit())

    #if form.validate_on_submit():
    if request.method == 'POST' and 'pic' in request.files:
        print('Hello')
        # create a folder for the new dataset
        #new_n = dataset_number(UPLOADS_FOLDER)
        #data_prefix = os.path.join(UPLOADS_FOLDER, 'data' + str(new_n) + '_')
        #print(data_prefix)

        #analysis_path = os.path.join(ANALYSIS_FOLDER, 'data' + str(new_n))
        #os.mkdir(analysis_path)

        new_n = dataset_number(ANALYSIS_FOLDER)
        data_path = os.path.join(ANALYSIS_FOLDER, 'data' + str(new_n))
        print(data_path)
        os.mkdir(data_path)

        # save e-mail
        #with open(data_prefix + 'name_of_actor.txt', 'w') as text_file:
        #    text_file.write(form.name.data)
        with open(os.path.join(data_path, 'name_of_actor.txt'), 'w') as text_file:
            text_file.write(form.name.data)


        #uploaded_file = request.files['pic']
        #print(uploaded_file)
        #if uploaded_file:
        #    print(uploaded_file.filename)
        #    uploaded_file.save(os.path.join(data_path, 'file.png'),)


        filename = tables.save(storage=request.files['pic'], folder='data' + str(new_n), name='count_matrix.')
        #filename = photos.save(storage=request.files['pic'], folder='data' + str(new_n), name='pic.')
        print(filename)

        return redirect(url_for('home'))
    return render_template('test.html', form=form)
'''


#############################   forms   ################################################
def dataset_number(path):
    existing_folders = [d for d in os.listdir(path) if os.path.isdir(os.path.join(path, d))]
    existing_folders = [d for d in existing_folders if d[:4] == 'data' and len(d) > 4]
    print(existing_folders)
    if existing_folders:
        new_n = max([int(d[4:]) for d in existing_folders]) + 1
    else:
        new_n = 1
    '''
    onlyfiles = [f for f in os.listdir(path) if os.path.isfile(os.path.join(path, f))]
    prefixes = set([s[:s.find('_')] for s in onlyfiles])
    if prefixes:
        new_n = max([int(p[4:]) for p in prefixes]) + 1
    else:
        new_n = 1
    '''
    print(new_n)
    return new_n


class HomePage(FlaskForm):
    submit = SubmitField('Run SPIRAL')


class LoadData(FlaskForm):
    dataset_name = StringField('Name of dataset:')
    email = StringField('* E-mail address:', validators=[DataRequired(),
                                                         Email("This field requires a valid email address")])
    count_matrix = FileField('* Gene expression matrix:',
                             validators=[FileRequired(), FileAllowed(tables, 'csv or txt files only!')]
                             )
    spatial_coors = FileField('Spatial coordinates:',
                              validators=[FileAllowed(tables, 'csv or txt files only!')]
                              )
    species = SelectField('Species:', choices=[("ARABIDOPSIS_THALIANA", 'Arabidopsis thaliana'),
                                               ("SACCHAROMYCES_CEREVISIAE", 'Saccharomyces cerevisiae'),
                                               ("CAENORHABDITIS_ELEGANS", 'Caenorhabditis elegans'),
                                               ("DROSOPHILA_MELANOGASTER", 'Drosophila melanogaster'),
                                               ("DANIO_RERIO", 'Danio rerio (Zebrafish)'),
                                               ("HOMO_SAPIENS", 'Homo sapiens'),
                                               ("MUS_MUSCULUS", 'Mus musculus'),
                                               ("RATTUS_NORVEGICUS", 'Rattus norvegicus')],
                          validators=[DataRequired()])
    samp_name = SelectField('How are your samples called?',
                            choices=[("samples", 'samples'), ("cells", 'cells'), ("spots", 'spots')],
                            validators=[DataRequired()], default="samples")
    submit = SubmitField('Submit')


class CheckData(FlaskForm):
    keep_going = SubmitField("Sounds about right, let's proceed")
    go_back = SubmitField("Go back and reload files")


class CheckFilteringParams(FlaskForm):
    keep_going = SubmitField("Sounds good, next step please")
    go_back = SubmitField("Go back and change the filtering parameters")


class vln_plot_form(FlaskForm):
    min_nFeatures = IntegerField('* Minimal number of expressed genes:', validators=[DataRequired()])
    max_nFeatures = IntegerField('* Maximal number of expressed genes:', validators=[DataRequired()])
    max_mtpercent = IntegerField('Maximal percent of mitochondrial gene expression:')
    submit = SubmitField('Submit')


class alg_params_form(FlaskForm):
    num_std_thr = SelectMultipleField('\u03B1:',
                                      choices=[('0.5', '0.5'), ('0.75', '0.75'), ('1', '1'), ('1.25', '1.25'),
                                               ('1.5', '1.5')])
    mu = SelectMultipleField('\u03bc:', choices=[('0.9', '0.9'), ('0.95', '0.95')])
    path_len = SelectField('Path length:', choices=[('3', '3')])
    num_iter = SelectField('Number of iterations:', choices=[('10000', '10000')])
    submit = SubmitField('Submit')


########################################################################################
# @app.route('/static/vln_plots/<filename>')
# def get_vln_route(filename):
#    return '/static/vln_plots/' + filename

@app.route("/static/vln_plots/<path:filename>")
def get_vln_route(filename):
    # return send_from_directory('/static/vln_plots/', filename, as_attachment=True)
    return '/static/vln_plots/' + filename


@app.route('/', methods=['POST', 'GET'])
def home():
    print('home!!!')
    form = HomePage()
    if form.validate_on_submit():
        return redirect(url_for('load_data_form'))
    return render_template('index.html', form=form)


@app.route('/run_SPIRAL_step1', methods=['POST', 'GET'])
def load_data_form():
    print('load_data_form!!!')
    form = LoadData(request.form)
    print('hi')
    print('form.errors:', form.errors)
    print('form.email.errors:', form.email.errors)
    print('form.validate_on_submit():', form.validate_on_submit())

    if request.method == 'POST' and 'count_matrix' in request.files and not form.email.errors:
        print('Hi')
        print(form.errors)

        # create a folder for the new dataset
        data_n = dataset_number(ANALYSIS_FOLDER)
        data_path = os.path.join(ANALYSIS_FOLDER, 'data' + str(data_n))
        print(data_path)
        os.mkdir(data_path)

        # save dataset_name
        if form.dataset_name.data != '':
            with open(os.path.join(data_path, 'dataset_name.txt'), 'w') as text_file:
                text_file.write(form.dataset_name.data)

        # save e-mail
        with open(os.path.join(data_path, 'email.txt'), 'w') as text_file:
            text_file.write(form.email.data)

        filename = tables.save(storage=request.files['count_matrix'], folder='data' + str(data_n), name='counts.')
        print(filename)

        # if 'spatial_coors' in request.files:
        if request.files['spatial_coors'].filename:
            filename = tables.save(storage=request.files['spatial_coors'], folder='data' + str(data_n),
                                   name='spatial_coors.')
            print(filename)

        # save species
        species = form.species.data
        with open(os.path.join(data_path, 'species.txt'), 'w') as text_file:
            text_file.write(species)

        # save sample name
        samp_name = form.samp_name.data
        with open(os.path.join(data_path, 'samp_name.txt'), 'w') as text_file:
            text_file.write(samp_name)

        return redirect(url_for('check_data_form', data_n=data_n, samp_name=samp_name, species=species))

        # flash("This field requires a valid email address")
    return render_template('load_data_form.html', form=form)


@app.route('/run_SPIRAL_step1.5', methods=['POST', 'GET'])
def check_data_form():
    print('check_data_form!!!')
    data_n = request.args['data_n']
    samp_name = request.args['samp_name']
    species = request.args['species']
    spatial, num_cells_orig, num_genes_orig, num_cells, num_genes, label_set, nlabels = load_data_first_time(
        analysis_folder=ANALYSIS_FOLDER, data_n=data_n, median_count_normalization_flag=True, with_log=False)

    form = CheckData(request.form)
    if form.validate_on_submit():
        if form.keep_going.data:
            return redirect(url_for('violin_plots', data_n=data_n, spatial=spatial, species=species,
                                    num_cells=num_cells, num_genes=num_genes, samp_name=samp_name))
        elif form.go_back.data:
            return redirect(url_for('load_data_form'))

    return render_template('check_data_form.html', form=form, spatial=spatial,
                           num_cells_orig=num_cells_orig, num_genes_orig=num_genes_orig,
                           num_cells=num_cells, num_genes=num_genes, label_set=label_set, nlabels=nlabels,
                           samp_name=samp_name)


@app.route('/run_SPIRAL_step2', methods=['POST', 'GET'])
def violin_plots():
    print('violin_plots!!!')
    data_n = request.args['data_n']
    spatial = request.args.get('spatial', type=lambda v: v.lower() == 'true')
    num_cells = request.args['num_cells']
    num_genes = request.args['num_genes']
    samp_name = request.args['samp_name']
    species = request.args['species']

    data_path = os.path.join(ANALYSIS_FOLDER, 'data' + str(data_n))
    vln_plot_filename = vln_plot_filename_(static_path=app._static_folder, data_n=data_n)
    with_mt_filename = os.path.join(data_path, 'with_mt.txt')
    if not os.path.exists(vln_plot_filename) or not os.path.exists(with_mt_filename):
        with_mt = compute_violin_plots(analysis_folder=ANALYSIS_FOLDER, data_n=data_n, static_path=app._static_folder,
                                       species=species)
        with open(with_mt_filename, 'w') as text_file:
            text_file.write(str(with_mt))
    else:
        with_mt = (open(with_mt_filename, "r").read().lower() == 'true')

    form = vln_plot_form()
    if request.method == 'POST' and not form.min_nFeatures.errors and not form.max_nFeatures.errors:
        ###### WRITE THIS !!!! ###############

        # if not check_nFeatures_filtering_values(data_n=data_n,
        #                                        min_nFeatures=form.min_nFeatures.data,
        #                                        max_nFeatures=form.max_nFeatures.data):
        #    flash('NO NO NO')
        # if not check_mtpercent_filtering_values(data_n=data_n, max_mtpercent=form.max_mtpercent.data):
        #    flash('NO NO NO')

        data_path = os.path.join(ANALYSIS_FOLDER, 'data' + str(data_n))

        # save min_nFeatures
        min_nFeatures = form.min_nFeatures.data
        with open(os.path.join(data_path, 'min_nFeatures.txt'), 'w') as text_file:
            text_file.write(str(min_nFeatures))

        # save max_nFeatures
        max_nFeatures = form.max_nFeatures.data
        with open(os.path.join(data_path, 'max_nFeatures.txt'), 'w') as text_file:
            text_file.write(str(max_nFeatures))

        # save max_mtpercent
        max_mtpercent = form.max_mtpercent.data
        if max_mtpercent is not None:
            with open(os.path.join(data_path, 'max_mtpercent.txt'), 'w') as text_file:
                text_file.write(str(max_mtpercent))
        print(max_mtpercent)
        return redirect(url_for('check_data_filtering_params',
                                data_n=data_n, min_nFeatures=min_nFeatures, species=species,
                                max_nFeatures=max_nFeatures, max_mtpercent=max_mtpercent,
                                spatial=spatial, num_cells=num_cells, num_genes=num_genes,
                                samp_name=samp_name))

    return render_template('violin_plots.html', form=form, vln_plot_file='/static/vln_data' + str(data_n) + '.jpg',
                           with_mt=with_mt)


@app.route('/run_SPIRAL_step2.5', methods=['POST', 'GET'])
def check_data_filtering_params():
    print('check_data_filtering_params!!!')
    data_n = request.args['data_n']
    min_nFeatures = int(request.args['min_nFeatures'])
    max_nFeatures = int(request.args['max_nFeatures'])
    try:
        max_mtpercent = int(request.args['max_mtpercent'])
    except:
        max_mtpercent = None
    spatial = request.args.get('spatial', type=lambda v: v.lower() == 'true')
    num_cells_orig2 = request.args['num_cells']
    num_genes_orig2 = request.args['num_genes']
    samp_name = request.args['samp_name']
    species = request.args['species']

    data_path = os.path.join(ANALYSIS_FOLDER, 'data' + str(data_n))
    data_norm_loc = data_norm_loc_name(data_path)
    spatial_norm_loc = spatial_norm_loc_name(data_path)
    data_norm_filt_loc = data_norm_filt_loc_name(data_path)
    spatial_norm_filt_loc = spatial_norm_filt_loc_name(data_path)

    num_genes_filename = os.path.join(data_path, 'num_genes.txt')
    num_cells_filename = os.path.join(data_path, 'num_cells.txt')
    if not os.path.exists(num_genes_filename) or not os.path.exists(num_cells_filename):
        data = filter_data(data_path=data_path, min_nFeatures=min_nFeatures, max_nFeatures=max_nFeatures,
                           max_mtpercent=max_mtpercent, spatial=spatial, species=species,
                           data_norm_loc=data_norm_loc, spatial_norm_loc=spatial_norm_loc,
                           data_norm_filt_loc=data_norm_filt_loc, spatial_norm_filt_loc=spatial_norm_filt_loc)
        num_genes, num_cells = data.shape[0], data.shape[1]
        with open(num_genes_filename, 'w') as text_file:
            text_file.write(str(num_genes))
        with open(num_cells_filename, 'w') as text_file:
            text_file.write(str(num_cells))
    else:
        num_genes = open(num_genes_filename, "r").read()
        num_cells = open(num_cells_filename, "r").read()

    form = CheckFilteringParams(request.form)
    if form.validate_on_submit():
        if form.keep_going.data:
            return redirect(url_for('alg_params', data_n=data_n))
        elif form.go_back.data:
            return redirect(url_for('violin_plots', data_n=data_n, spatial=spatial, species=species,
                                    num_cells=num_cells_orig2, num_genes=num_genes_orig2, samp_name=samp_name))

    return render_template('check_data_filtering_params.html', form=form,
                           num_cells_orig=num_cells_orig2, num_genes_orig=num_genes_orig2,
                           num_cells=num_cells, num_genes=num_genes, samp_name=samp_name)


@app.route('/run_SPIRAL_step3', methods=['POST', 'GET'])
def alg_params():
    print('alg_params!!!')
    data_n = request.args['data_n']
    form = alg_params_form()
    if form.validate_on_submit():
        num_stds_thresh_lst = [float(s) for s in form.num_std_thr.data]
        mu_lst = [float(s) for s in form.mu.data]
        path_len_lst = [int(form.path_len.data)]
        num_iters_lst = [int(form.num_iter.data)]
        print(num_stds_thresh_lst, mu_lst, path_len_lst, num_iters_lst)
        thread = Thread(target=run_SPIRAL,
                        kwargs={'num_stds_thresh_lst': num_stds_thresh_lst, 'mu_lst': mu_lst,
                                'num_iters_lst': num_iters_lst, 'path_len_lst': path_len_lst,
                                'analysis_folder': ANALYSIS_FOLDER, 'data_n': data_n})
        thread.start()

        return redirect(url_for('all_done', data_n=data_n))
    return render_template('alg_params.html', form=form)


def run_SPIRAL(analysis_folder, data_n, num_stds_thresh_lst, mu_lst, num_iters_lst, path_len_lst):
    run_SPIRAL_pipeline(analysis_folder=analysis_folder, data_n=data_n,
                        num_stds_thresh_lst=num_stds_thresh_lst, mu_lst=mu_lst,
                        num_iters_lst=num_iters_lst, path_len_lst=path_len_lst)
    with app.app_context():
        email(data_n=data_n)


@app.route('/all_done', methods=['POST', 'GET'])
def all_done():
    print('all_done!!!')
    data_n = request.args['data_n']
    return 'All done, wait for an email.'


@app.route("/email")
def email(data_n):
    data_path = os.path.join(ANALYSIS_FOLDER, 'data' + str(data_n))
    recipient_email = open(os.path.join(data_path, 'email.txt'), "r").read()
    msg = Message('SPIRAL results are ready for you', sender='SPIRAL.web.tool@gmail.com', recipients=[recipient_email])
    msg.body = "Hello, \nThe results link is " + url_for('results_panel', url=data_n_to_url(data_n), _external=True)
    mail.send(msg)
    return "Message sent!"


#################### results panel ##########################################################################
'''
@app.route('/callback', methods=['POST', 'GET'])
def cb():
    return gm(request.args.get('data'))

@app.route('/results')
def results():
    return render_template('results_panel.html', img_data=gm())


def gm(country='Israel', data_n=48):
    #df = pd.DataFrame(px.data.gapminder())
    #fig = px.line(df[df['country'] == country], x="year", y="gdpPercap")
    #print(type(fig))

    #graphJSON = json.dumps(fig, cls=plotly.utils.PlotlyJSONEncoder)
    #print(type(graphJSON))
    #print(graphJSON)
    #return graphJSON

    if country == 'Israel':
        picfile = os.path.join(ANALYSIS_FOLDER, 'data' + str(data_n), 'structure_layouts', 'struct_33_Gnet.jpg')
    elif country == 'Germany':
        picfile = os.path.join(ANALYSIS_FOLDER, 'data' + str(data_n), 'structure_layouts', 'struct_33_Gspatial.jpg')

    #data = {}
    #with open(picfile, mode='rb') as file:
    #    img = file.read()
    #data['img'] = base64.encodebytes(img).decode('utf-8')
    #imgJSON = json.dumps(data)
    #print(type(imgJSON))
    #return imgJSON

    im = Image.open(picfile)
    data = io.BytesIO()
    im.save(data, "JPEG")
    encoded_img_data = base64.b64encode(data.getvalue())
    return encoded_img_data.decode('utf-8')

@app.route('/results')
def results():
    return render_template('results_panel.html', data_n='48')
'''


############################################################################################################
@app.route('/_pic_names_and_GO_terms')
def pic_names_and_GO_terms():
    # get file names of the wanted figures
    struct = request.args.get('struct', 'no_struct', type=str)
    vis_type = request.args.get('vis_type', 'no_vis_type', type=str)
    color_type = request.args.get('color_type', 'no_color_type', type=str)
    data_n = request.args.get('data_n', '000', type=str)
    print(struct, vis_type, color_type, data_n)

    data_path = os.path.join(ANALYSIS_FOLDER, 'data' + str(data_n))
    pic_path = pic_folder(data_path)
    print(pic_path)

    spatial = False
    if [file for file in os.listdir(pic_path) if file.endswith("_spatial.jpg")]:
        spatial = True

    if vis_type == 'no_vis_type':  # impute_method = 'no_imputation'
        impute_method = 'no_imputation'
        if color_type == 'avg_expression':
            struct_pics = ['struct_' + struct + '_Gnet.jpg',
                           'struct_' + struct + '_orig_GPCA.jpg',
                           'struct_' + struct + '_orig_GUMAP.jpg']
            if spatial:
                struct_pics.append('struct_' + struct + '_Gspatial.jpg')
        elif color_type == 'layers':
            struct_pics = ['struct_' + struct + '_net.jpg',
                           'struct_' + struct + '_orig_PCA.jpg',
                           'struct_' + struct + '_orig_UMAP.jpg']
            if spatial:
                struct_pics.append('struct_' + struct + '_spatial.jpg')
    else:
        impute_method = 'agg_wald'
        if vis_type == 'repcells':
            if color_type == 'avg_expression':
                struct_pics = ['struct_' + struct + '_Gnet.jpg',
                               'struct_' + struct + '_repcells_GPCA.jpg',
                               'struct_' + struct + '_repcells_GUMAP.jpg']
            elif color_type == 'layers':
                struct_pics = ['struct_' + struct + '_net.jpg',
                               'struct_' + struct + '_repcells_PCA.jpg',
                               'struct_' + struct + '_repcells_UMAP.jpg']
        elif vis_type == 'orig':
            if color_type == 'avg_expression':
                struct_pics = ['struct_' + struct + '_orig_GPCA.jpg',
                               'struct_' + struct + '_orig_GUMAP.jpg']
                if spatial:
                    struct_pics.append('struct_' + struct + '_Gspatial.jpg')
            elif color_type == 'layers':
                struct_pics = ['struct_' + struct + '_orig_PCA.jpg',
                               'struct_' + struct + '_orig_UMAP.jpg']
                if spatial:
                    struct_pics.append('struct_' + struct + '_spatial.jpg')
    # pic_path = '/static/analysis/data' + data_n + '/structure_layouts/'

    struct_pics = ['/' + os.path.join(pic_path, file).replace('\\', '/') for file in struct_pics]
    print(struct_pics)

    # get GO enrichment terms
    # impute_method = request.args.get('impute_method', 'no_impute_method', type=str)
    sigfile = final_sig_filename(data_path, impute_method)
    # sigtable = pd.read_csv(sigfile, index_col=0, sep='\t')
    sigtable = load_excel_with_openpyxl_and_convert_to_pd_DataFrame(sigfile)
    GO_terms = [sigtable.loc[int(struct), 'proc_GOterms_below_1e-06'],
                sigtable.loc[int(struct), 'func_GOterms_below_1e-06'],
                sigtable.loc[int(struct), 'comp_GOterms_below_1e-06']]

    # get gene list
    gene_list = sigtable.loc[int(struct), 'genes'].replace("['", "").replace("']", "").replace("\n", "").split("', '")
    print(gene_list)
    return jsonify(struct_pics=struct_pics, gene_list=gene_list, GO_terms=GO_terms)


@app.route('/results/result_table_<data_n>', methods=['GET', 'POST'])
def download_sigfile(data_n):
    data_path = os.path.join(ANALYSIS_FOLDER, 'data' + str(data_n))
    impute_method = open(os.path.join(data_path, 'imputation_method.txt'), "r").read()
    sigfile = os.path.basename(final_sig_filename(data_path, impute_method))
    return send_from_directory(data_path, sigfile)


@app.route('/results/structure_layout_zipfile_<data_n>', methods=['GET', 'POST'])
def download_structure_layout_zipfile(data_n):
    data_path = os.path.join(ANALYSIS_FOLDER, 'data' + str(data_n))
    structure_layout_zipfile = 'structure_layouts.zip'
    return send_from_directory(data_path, structure_layout_zipfile)


@app.route('/results/download_repcell_partition_<data_n>', methods=['GET', 'POST'])
def download_repcell_partition(data_n):
    data_path = os.path.join(ANALYSIS_FOLDER, 'data' + str(data_n))
    repcell_partition_zipfile = 'repcell_partition.zip'
    return send_from_directory(data_path, repcell_partition_zipfile)


@app.route('/results/<url>')
def results_panel(url):
    data_n = url_to_data_n(url)

    # get structure list
    full_pic_path = os.path.join(ANALYSIS_FOLDER, 'data' + str(data_n), 'structure_layouts')
    pic_suffs = [file[file.find('_') + 1:] for file in os.listdir(full_pic_path) if
                 file.endswith(".jpg") and 'struct' in file]
    struct_lst = list(set([int(file[:file.find('_')]) for file in pic_suffs]))
    struct_lst = [str(s) for s in struct_lst]

    # get imputation method
    data_path = os.path.join(ANALYSIS_FOLDER, 'data' + str(data_n))
    impute_method = open(os.path.join(data_path, 'imputation_method.txt'), "r").read()
    return render_template('results_panel.html', data_n=data_n, struct_lst=struct_lst, impute_method=impute_method)


def data_n_to_url(data_n):
    sample_string = str(1000000 + int(data_n))
    sample_string_bytes = sample_string.encode("ascii")
    base64_bytes = base64.b64encode(sample_string_bytes)
    base64_string = base64_bytes.decode("ascii")
    return base64_string


def url_to_data_n(url):
    base64_string = url
    base64_bytes = base64_string.encode("ascii")
    sample_string_bytes = base64.b64decode(base64_bytes)
    sample_string = int(sample_string_bytes.decode("ascii")) - 1000000
    return sample_string


############################################################################################################
if __name__ == '__main__':
    # app.run(host='0.0.0.0', port=5000, debug=True)
    app.run(host='0.0.0.0', port=5000)
