#!/var/www/html/SPIRAL.web.tool/spiral_venv/bin/python3.10

from numba.core.errors import NumbaDeprecationWarning, NumbaPendingDeprecationWarning
import warnings
warnings.simplefilter('ignore', category=NumbaDeprecationWarning)
warnings.simplefilter('ignore', category=NumbaPendingDeprecationWarning)

import urllib
from flask import Flask, render_template, request, redirect, url_for, flash, send_from_directory, jsonify, Markup
from flask_bootstrap import Bootstrap
from flask_wtf import FlaskForm
from wtforms import StringField, IntegerField, SubmitField, SelectField, SelectMultipleField, IntegerRangeField, \
    ValidationError, BooleanField
from wtforms.validators import DataRequired, Email, NumberRange, EqualTo
from flask_wtf.file import FileField, FileRequired, FileAllowed
from flask_uploads import UploadSet, configure_uploads, IMAGES
from flask_mail import Mail, Message

from SPIRAL_pipeline_funcs import *

from threading import Thread
import gzip

import socket

# get hostname.
hostname = socket.gethostname()
if hostname in ['bi-fiona', 'Spin']:  # production (linux)
    production = True
else:
    production = False  # development (Windows)

if production:
    if hostname == 'bi-fiona':
        app = Flask(__name__, static_folder='/home/yaellab/SPIRAL/static')
    else:
        app = Flask(__name__, static_folder='/var/www/html/SPIRAL.web.tool/static')
    app.config['SERVER_NAME'] = 'spiral.technion.ac.il'
else:
    app = Flask(__name__)
    app.config['SERVER_NAME'] = '132.69.236.238:5000'

app.config['MAIL_SERVER'] = 'smtp.gmail.com'
app.config['MAIL_PORT'] = 465
app.config['MAIL_USERNAME'] = 'SPIRAL.web.tool@gmail.com'
app.config['MAIL_USE_TLS'] = False
app.config['MAIL_USE_SSL'] = True

# This is a special password for gmail from this app, as explained here https://stackoverflow.com/questions/72478573/how-to-send-an-email-using-python-after-googles-policy-update-on-not-allowing-j
if production:
    if hostname == 'bi-fiona':
        app.config['MAIL_PASSWORD'] = open('/home/yaellab/SPIRAL/mail_password.txt', 'r').read()
    else:
        app.config['MAIL_PASSWORD'] = open('/var/www/html/SPIRAL.web.tool/mail_password.txt', 'r').read()
else:
    app.config['MAIL_PASSWORD'] = open('./mail_password.txt', 'r').read()

mail = Mail(app)

app.config['MAX_CONTENT_LENGTH'] = 1024 * 1024 * 1024

ANALYSIS_FOLDER = os.path.join(app._static_folder, 'analysis')

# Flask-WTF requires an encryption key - the string can be anything
if production:
    if hostname == 'bi-fiona':
        app.config['SECRET_KEY'] = open('/home/yaellab/SPIRAL/Flask-WTF_encription_key.txt', 'r').read()
    else:
        app.config['SECRET_KEY'] = open('/var/www/html/SPIRAL.web.tool/Flask-WTF_encription_key.txt', 'r').read()
else:
    app.config['SECRET_KEY'] = open('./Flask-WTF_encription_key.txt', 'r').read()

app.config['UPLOADED_TABLES_DEST'] = ANALYSIS_FOLDER
app.config['UPLOADED_MATRIX_DEST'] = ANALYSIS_FOLDER
app.config['UPLOADED_BARCODES_DEST'] = ANALYSIS_FOLDER

# Flask-Bootstrap requires this line
Bootstrap(app)

tables = UploadSet(name='tables', extensions=('txt', 'csv'))
matrix = UploadSet(name='matrix', extensions='mtx.gz')
barcodes_features = UploadSet(name='barcodes', extensions='tsv.gz')
configure_uploads(app, (tables, matrix, barcodes_features))

#patch_request_class(app, 1024 * 1024 * 1024)
app.config['MAX_CONTENT_LENGTH'] = 1024 * 1024 * 1024

ensembl_folder = ensembl_loc(production, hostname)

#############################   forms   ################################################
def dataset_number(path):
    existing_folders = [d for d in os.listdir(path) if os.path.isdir(os.path.join(path, d))]
    existing_folders = [d for d in existing_folders if d[:4] == 'data' and len(d) > 4]
    if existing_folders:
        new_n = max([int(d[4:]) for d in existing_folders]) + 1
    else:
        new_n = 1
    return new_n


class HomePage(FlaskForm):
    style = {
        'style': 'background-color: #005959; color: white; padding: 15px 32px; text-align: center; font-size: 16px;'}
    submit = SubmitField('Run SPIRAL on my data set', render_kw=style)


class DataError(FlaskForm):
    style = {
        'style': 'background-color: #005959; color: white; padding: 15px 32px; text-align: center; font-size: 16px;'}
    submit = SubmitField('Go back and reload data', render_kw=style)


class LoadData(FlaskForm):
    dataset_name = StringField('Name of dataset: ')
    email = StringField('E-mail address: ', validators=[Email("Please enter a valid email address")])
    count_matrix = FileField('Gene expression matrix: ',
                             validators=[FileRequired(), FileAllowed(tables, 'csv or txt files only!')])
    visium_10x_matrix = FileField('Vision 10x matrix: ', validators=[FileAllowed(matrix, 'mtx.gz files only!')])
    visium_10x_features = FileField('Vision 10x features: ',
                                    validators=[FileAllowed(barcodes_features, 'tsv.gz files only!')])
    visium_10x_barcodes = FileField('Vision 10x barcodes: ',
                                    validators=[FileAllowed(barcodes_features, 'tsv.gz files only!')])
    spatial_coors = FileField('Spatial coordinates: ', validators=[FileAllowed(tables, 'csv or txt files only!')])
    species = SelectField('Species: ', choices=[("ARABIDOPSIS_THALIANA", 'Arabidopsis thaliana'),
                                                ("SACCHAROMYCES_CEREVISIAE", 'Saccharomyces cerevisiae'),
                                                ("CAENORHABDITIS_ELEGANS", 'Caenorhabditis elegans'),
                                                ("DROSOPHILA_MELANOGASTER", 'Drosophila melanogaster'),
                                                ("DANIO_RERIO", 'Danio rerio (Zebrafish)'),
                                                ("HOMO_SAPIENS", 'Homo sapiens'),
                                                ("MUS_MUSCULUS", 'Mus musculus'),
                                                ("RATTUS_NORVEGICUS", 'Rattus norvegicus'),
                                                ('other', 'other (GO enrichment analysis will not be performed)'),
                                                ('synthetic',
                                                 'synthetic (GO enrichment analysis will not be performed)')],
                          default=("HOMO_SAPIENS", 'Homo sapiens'), validators=[DataRequired()])
    samp_name = SelectField('How are your samples called? ',
                            choices=[("samples", 'samples'), ("cells", 'cells'), ("spots", 'spots')],
                            validators=[DataRequired()], default="samples")
    labels_checkbox = BooleanField(Markup('My data set does <strong>not</strong> have labels. '))
    submit = SubmitField('Submit')


class CheckData(FlaskForm):
    keep_going_style = {
        'style': 'background-color: #2a7d4f; color: white; padding: 15px 32px; text-align: center; font-size: 16px;'}
    keep_going = SubmitField("Sounds about right, let's proceed", render_kw=keep_going_style)
    go_back_style = {
        'style': 'background-color: #800000; color: white; padding: 15px 32px; text-align: center; font-size: 16px;'}
    go_back = SubmitField("Go back and reload files", render_kw=go_back_style)


class CheckFilteringParams(FlaskForm):
    keep_going_style = {
        'style': 'background-color: #2a7d4f; color: white; padding: 15px 32px; text-align: center; font-size: 16px;'}
    keep_going = SubmitField("Sounds good, next step please", render_kw=keep_going_style)
    go_back_style = {
        'style': 'background-color: #800000; color: white; padding: 15px 32px; text-align: center; font-size: 16px;'}
    go_back = SubmitField("Go back and change the filtering parameters", render_kw=go_back_style)


class MoreThan(object):

    def __init__(self, fieldname, message=None):
        self.fieldname = fieldname
        self.message = message

    def __call__(self, form, field):
        try:
            other = form[self.fieldname]
        except KeyError:
            raise ValidationError(field.gettext("Invalid field name '%s'.") % self.fieldname)
        if field.data <= other.data:
            d = {
                'other_label': hasattr(other, 'label') and other.label.text or self.fieldname,
                'other_name': self.fieldname
            }
            message = self.message
            if message is None:
                message = field.gettext(
                    'The maximal number of expressed genes has to be larger than the minimal number of expressed genes')

            raise ValidationError(message)


class vln_plot_form(FlaskForm):
    min_nFeatures = IntegerField('Minimal number of expressed genes: ', validators=[DataRequired()])
    max_nFeatures = IntegerField('Maximal number of expressed genes: ',
                                 validators=[DataRequired(), MoreThan('min_nFeatures')])
    max_mtpercent = IntegerField('Maximal percent of mitochondrial gene expression: ', default=100,
                                 validators=[NumberRange(1, 100)])
    submit = SubmitField('Submit')


class alg_params_form(FlaskForm):
    fill_default_vals = BooleanField('Fill recommended values')
    num_std_thr = SelectMultipleField(
        'Number of standard deviations threshold \u03B1: ',
        choices=[('0.5', '0.5'), ('0.75', '0.75'), ('1', '1'), ('1.25', '1.25'),
                 ('1.5', '1.5')])
    mu = SelectMultipleField('Density parameter \u03bc: ', choices=[('0.9', '0.9'), ('0.95', '0.95')])
    path_len = SelectField('Path length L: ', choices=[('3', '3')])
    num_iter = SelectField('Number of iterations T: ', choices=[('10000', '10000')])
    submit = SubmitField('Submit')


def flash_errors(form):
    """Flashes form errors"""
    for field, errors in form.errors.items():
        for error in errors:
            flash(u"Error in the %s field - %s" % (
                getattr(form, field).label.text,
                error
            ), 'error')


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
    print(app.instance_path)
    form = HomePage()
    Zhang2019_link = '/results/' + data_n_to_url(1)
    Wagner2018_link = '/results/' + data_n_to_url(3)
    MouseSP2_link = '/results/' + data_n_to_url(2)
    normal_prostate_link = '/results/' + data_n_to_url(5)
    yaellab_differentiation_link = '/results/' + data_n_to_url(4)
    bulk_mouse_Bcells_link = '/results/' + data_n_to_url(6)
    if form.validate_on_submit():
        return redirect(url_for('load_data_form'))
    return render_template('index.html', form=form, Zhang2019_link=Zhang2019_link, Wagner2018_link=Wagner2018_link,
                           MouseSP2_link=MouseSP2_link, normal_prostate_link=normal_prostate_link,
                           yaellab_differentiation_link=yaellab_differentiation_link,
                           bulk_mouse_Bcells_link=bulk_mouse_Bcells_link)


@app.route('/run_SPIRAL_step1', methods=['POST', 'GET'])
def load_data_form():
    print('load_data_form!!!')
    form = LoadData(request.form)
    if request.method == 'POST' and not form.email.errors and ('count_matrix' in request.files or all(
            x in request.files for x in
            ['visium_10x_matrix', 'visium_10x_features', 'visium_10x_barcodes'])):
        # create a folder for the new dataset
        data_n = dataset_number(ANALYSIS_FOLDER)
        data_path = os.path.join(ANALYSIS_FOLDER, 'data' + str(data_n))
        os.mkdir(data_path)

        # save dataset_name
        if form.dataset_name.data != '':
            with open(os.path.join(data_path, 'dataset_name.txt'), 'w') as text_file:
                text_file.write(form.dataset_name.data)

        # save e-mail
        with open(os.path.join(data_path, 'email.txt'), 'w') as text_file:
            text_file.write(form.email.data)

        # save species
        species = form.species.data
        with open(os.path.join(data_path, 'species.txt'), 'w') as text_file:
            text_file.write(species)

        # save sample name
        samp_name = form.samp_name.data
        with open(os.path.join(data_path, 'samp_name.txt'), 'w') as text_file:
            text_file.write(samp_name)

        # labels_checkbox
        labels_checkbox = form.labels_checkbox.data
        with open(os.path.join(data_path, 'labels_checkbox.txt'), 'w') as text_file:
            text_file.write(str(labels_checkbox))

        if request.files['count_matrix'].filename:
            tables.save(storage=request.files['count_matrix'], folder='data' + str(data_n), name='counts.')
            # if 'spatial_coors' in request.files:
            if request.files['spatial_coors'].filename:
                tables.save(storage=request.files['spatial_coors'], folder='data' + str(data_n), name='spatial_coors.')
        else:
            matrix_dir = os.path.join('data' + str(data_n), 'temp')
            data_path = os.path.join(ANALYSIS_FOLDER, matrix_dir)
            os.mkdir(data_path)
            matrix.save(storage=request.files['visium_10x_matrix'], name=matrix_dir + '/matrix.mtx.gz')
            barcodes_features.save(storage=request.files['visium_10x_features'], name=matrix_dir + '/features.tsv.gz')
            barcodes_features.save(storage=request.files['visium_10x_barcodes'], name=matrix_dir + '/barcodes.tsv.gz')
            # read in MEX format matrix as table
            mat_filtered = scipy.io.mmread(os.path.join(data_path, "matrix.mtx.gz"))
            # list of transcript ids, e.g. 'ENSG00000187634'
            features_path = os.path.join(data_path, "features.tsv.gz")
            # list of gene names, e.g. 'SAMD11'
            gene_names = [row[1] for row in csv.reader(gzip.open(features_path, mode="rt"), delimiter="\t")]
            # list of barcodes, e.g. 'AAACATACAAAACG-1'
            barcodes_path = os.path.join(data_path, "barcodes.tsv.gz")
            barcodes = [row[0] for row in csv.reader(gzip.open(barcodes_path, mode="rt"), delimiter="\t")]
            # transform table to pandas dataframe and label rows and columns
            matrix_df = pd.DataFrame.sparse.from_spmatrix(mat_filtered)
            matrix_df.columns = barcodes
            matrix_df.insert(loc=0, column="gene", value=gene_names)
            # save the table as a CSV (note the CSV will be a very large file)
            matrix_df.to_csv(os.path.join(ANALYSIS_FOLDER, 'data' + str(data_n), "counts.csv"), index=False)
            # if 'spatial_coors' in request.files:
            if request.files['spatial_coors'].filename:
                # tables.save(storage=request.files['spatial_coors'], name=matrix_dir + '/spatial_coors.')
                spatial_coors = pd.read_csv(request.files['spatial_coors'])
                spatial_coors.drop(spatial_coors.iloc[:, 1:4], axis=1, inplace=True)
                spatial_coors.to_csv(os.path.join(ANALYSIS_FOLDER, 'data' + str(data_n), "spatial_coors.csv"),
                                     index=False)
            # os.remove(matrix_dir)

        return redirect(url_for('check_data_form', data_n=data_n, samp_name=samp_name, species=species))

        # flash("This field requires a valid email address")
    return render_template('load_data_form.html', form=form)


@app.route('/results/download_count_matrix_<data_n>', methods=['GET', 'POST'])
def download_count_matrix(data_n):
    data_path = os.path.join(ANALYSIS_FOLDER, 'data' + str(data_n))
    count_matrix_file = 'counts.csv'
    return send_from_directory(data_path, count_matrix_file)


@app.route('/results/download_spatial_coors_<data_n>', methods=['GET', 'POST'])
def download_spatial_coors(data_n):
    data_path = os.path.join(ANALYSIS_FOLDER, 'data' + str(data_n))
    spatial_coors_file = 'spatial_coors.csv'
    return send_from_directory(data_path, spatial_coors_file)


@app.route('/data_error', methods=['POST', 'GET'])
def bad_data_file():
    print('bad_data_file!!!')
    form = DataError()
    return render_template('bad_data_file.html', form=form)


@app.route('/run_SPIRAL_step2', methods=['POST', 'GET'])
def check_data_form():
    print('check_data_form!!!')
    data_n = request.args['data_n']
    samp_name = request.args['samp_name']
    species = request.args['species']
    try:
        spatial, num_cells_orig, num_genes_orig, num_cells, num_genes, label_set, nlabels = load_data_first_time(
            analysis_folder=ANALYSIS_FOLDER, data_n=data_n, median_count_normalization_flag=True, with_log=False)
    except:
        return redirect(url_for('bad_data_file'))

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


@app.route('/run_SPIRAL_step3', methods=['POST', 'GET'])
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
    vln_plot_mt_filename = vln_plot_mt_filename_(static_path=app._static_folder, data_n=data_n)
    with_mt_filename = with_mt_filename_(data_path)
    mt_error_filename = os.path.join(data_path, 'mt_error.txt')
    if not os.path.exists(vln_plot_filename) or not os.path.exists(with_mt_filename):
        with_mt, mt_error = compute_violin_plots(analysis_folder=ANALYSIS_FOLDER, data_n=data_n,
                                                 static_path=app._static_folder, ensembl_folder=ensembl_folder,
                                                 species=species)
        with open(with_mt_filename, 'w') as text_file:
            text_file.write(str(with_mt))
        with open(mt_error_filename, 'w') as text_file:
            text_file.write(mt_error)
    else:
        with_mt = (open(with_mt_filename, "r").read().lower() == 'true')
        mt_error = open(mt_error_filename, "r").read()

    form = vln_plot_form()
    if request.method == 'POST':
        if form.validate_on_submit():
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
            return redirect(url_for('check_data_filtering_params',
                                    data_n=data_n, min_nFeatures=min_nFeatures, species=species,
                                    max_nFeatures=max_nFeatures, max_mtpercent=max_mtpercent,
                                    spatial=spatial, num_cells=num_cells, num_genes=num_genes,
                                    samp_name=samp_name))
        else:
            flash_errors(form)  # doesn't work
    return render_template('violin_plots.html', form=form, vln_plot_file='/static/vln_data' + str(data_n) + '.jpg',
                           vln_plot_mt_file='/static/vln_mt_data' + str(data_n) + '.jpg',
                           with_mt=with_mt, mt_error=mt_error)


@app.route('/run_SPIRAL_step4', methods=['POST', 'GET'])
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
                           data_norm_filt_loc=data_norm_filt_loc, spatial_norm_filt_loc=spatial_norm_filt_loc,
                           ensembl_folder=ensembl_folder)
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
            os.remove(num_genes_filename)
            os.remove(num_cells_filename)
            return redirect(url_for('violin_plots', data_n=data_n, spatial=spatial, species=species,
                                    num_cells=num_cells_orig2, num_genes=num_genes_orig2, samp_name=samp_name))

    return render_template('check_data_filtering_params.html', form=form,
                           num_cells_orig=num_cells_orig2, num_genes_orig=num_genes_orig2,
                           num_cells=num_cells, num_genes=num_genes, samp_name=samp_name)


@app.route('/run_SPIRAL_step5', methods=['POST', 'GET'])
def alg_params():
    print('alg_params!!!')
    data_n = request.args['data_n']
    form = alg_params_form()
    if form.validate_on_submit() and form.submit.data:
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

        return redirect(url_for('results_panel', url=data_n_to_url(data_n)))
    return render_template('alg_params.html', form=form)


def run_SPIRAL(analysis_folder, data_n, num_stds_thresh_lst, mu_lst, num_iters_lst, path_len_lst):
    try:
        outcome = run_SPIRAL_pipeline(analysis_folder=analysis_folder, data_n=data_n,
                                      num_stds_thresh_lst=num_stds_thresh_lst, mu_lst=mu_lst,
                                      num_iters_lst=num_iters_lst, path_len_lst=path_len_lst,
                                      ensembl_folder=ensembl_folder)
        data_path = os.path.join(ANALYSIS_FOLDER, 'data' + str(data_n))
        with open(outcome_path(data_path), 'w') as text_file:
            text_file.write(str(outcome))
        with app.app_context():
            email(data_n=data_n, outcome=outcome)
    except Exception as error:
        print(error)
        outcome = 'Error'
        with open(outcome_path(data_path), 'w') as text_file:
            text_file.write(str(outcome))
        with app.app_context():
            email(data_n=data_n, outcome=outcome, error=error)


@app.route("/email")
def email(data_n, outcome, error=None):
    # if outcome!='Error': send an email to user and a (bcc) copy to self
    # if outcome=='Error': send an email to user and a (bcc) copy to self, also send an error notification to self
    data_path = os.path.join(ANALYSIS_FOLDER, 'data' + str(data_n))
    recipient_email = open(os.path.join(data_path, 'email.txt'), "r").read()
    if recipient_email:
        dataset_name = open(os.path.join(data_path, 'dataset_name.txt'), "r").read()
        if outcome == 'Success':
            msg = Message('SPIRAL results are ready', sender='SPIRAL.web.tool@gmail.com', recipients=[recipient_email],
                          bcc=['SPIRAL.web.tool@gmail.com'])
            msg.body = ("Hello, \nSPIRAL finished processing" + (dataset_name != '') * (' ' + dataset_name) + '.' +
                        "\nThe results link is " + urllib.parse.unquote(
                        url_for('results_panel', url=data_n_to_url(data_n), _external=True, _scheme='https')))
        elif outcome == 'No_structures':
            msg = Message('SPIRAL results', sender='SPIRAL.web.tool@gmail.com', recipients=[recipient_email],
                          bcc=['SPIRAL.web.tool@gmail.com'])
            msg.body = ("Hello, \nUnfortunately SPIRAL was not able to find any structures in" + (
                    dataset_name != '') * (
                                ' ' + dataset_name) +
                        '. You might want to check again the format of the data you uploaded. Also, we encourage you to contact us with any questions.')
        elif outcome == 'Error':
            msg = Message('SPIRAL results', sender='SPIRAL.web.tool@gmail.com', recipients=[recipient_email],
                          bcc=['SPIRAL.web.tool@gmail.com'])
            msg.body = ("Hello, \nUnfortunately SPIRAL encountered an error while computing the structures for" + (
                    dataset_name != '') * (
                                ' ' + dataset_name) +
                        '. We will look into this shortly. Also, you may reply this email with any questions.')
        mail.send(msg)

    if outcome == 'Error':
        email_to_inform = "spiral.web.tool@gmail.com"
        msg = Message('SPIRAL: An error occurred when computing data ' + str(data_n),
                      sender='SPIRAL.web.tool@gmail.com', recipients=[email_to_inform])
        msg.body = 'An error occurred:' + str(error) + '.\ndata_n:' + str(data_n) + '.'
        mail.send(msg)

    if recipient_email:
        if outcome != 'Error':
            return "Message to user sent!"
        else:
            return "Message to user sent! Error notification to SPIRAL sent!"
    else:
        if outcome != 'Error':
            return "No user email address"
        else:
            return "Error notification to SPIRAL sent!"


#################### results panel ##########################################################################
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

    if production:
        struct_pics = ['/static/analysis/data' + str(data_n) + '/structure_layouts/' + file for file in struct_pics]
    else:
        struct_pics = ['/' + os.path.join(pic_path, file).replace('\\', '/') for file in struct_pics]

    # get GO enrichment terms
    sigfile = final_sig_filename(data_path)
    sigtable = load_excel_with_openpyxl_and_convert_to_pd_DataFrame(sigfile)
    proc_cols = [col for col in list(sigtable) if 'proc_GOterms_below' in col]
    func_cols = [col for col in list(sigtable) if 'func_GOterms_below' in col]
    comp_cols = [col for col in list(sigtable) if 'comp_GOterms_below' in col]
    if len(proc_cols) >= 1 and len(func_cols) >= 1 and len(comp_cols) >= 1:
        proc_col = proc_cols[0]
        func_col = func_cols[0]
        comp_col = comp_cols[0]
        GO_terms = []
        for col in [proc_col, func_col, comp_col]:
            if sigtable.loc[int(struct), col] == 'NO TERMS':
                GO_terms.append('NO TERMS')
            else:
                GO_terms.append([m.group(0) for m in
                                 re.finditer(
                                     r"(GO:){1}[0-9]{7}(:){1}[A-Za-z0-9\- ]+( \(){1}(qval){1}[0-9E\-\.]+(\)){1}",
                                     sigtable.loc[int(struct), col])])
    else:
        GO_terms = [['Not available'], ['Not available'], ['Not available']]

    # get gene list
    gene_list = read_gene_list(sigtable.loc[int(struct), 'genes'])
    return jsonify(struct_pics=struct_pics, gene_list=gene_list, GO_terms=GO_terms)


@app.route('/results/result_table_<data_n>', methods=['GET', 'POST'])
def download_sigfile(data_n):
    data_path = os.path.join(ANALYSIS_FOLDER, 'data' + str(data_n))
    sigfile = os.path.basename(final_sig_filename(data_path))
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


@app.route('/_filter_struct_lst_for_result_panel')
def filter_struct_lst_for_result_panel():
    # produces a list of structures, in which there is no pair of structures with Jaccard index larger than
    # Jaccard_thr_genes. The function relies on the fact that the structure numbers are
    # range(1, len(Jaccard_mat_genes) + 1), and that structure i corresponds to row i-1 in Jaccard_mat_genes.
    print('filter_struct_lst_for_result_panel!!!')

    data_n = request.args.get('data_n', None, type=str)
    new_Jaccard_thr_genes = float(request.args.get('new_Jaccard_thr_genes', None, type=str))

    data_path = os.path.join(ANALYSIS_FOLDER, 'data' + str(data_n))
    # get the Jaccard matrix of gene lists
    Jaccard_mat_genes = np.load(os.path.join(data_path, 'Jaccard_mat_genes.npy'))

    if new_Jaccard_thr_genes is None:
        if Jaccard_mat_genes is None:
            return jsonify(thr_struct_lst=[])
        else:
            return jsonify(thr_struct_lst=range(1, len(Jaccard_mat_genes) + 1))

    inds = list(range(1, len(Jaccard_mat_genes) + 1))
    thr_struct_lst = []

    while inds:
        i = inds.pop(0)
        thr_struct_lst.append(i)
        similar_to_i = []

        for j in inds:
            if Jaccard_mat_genes[i - 1][j - 1] >= new_Jaccard_thr_genes:
                similar_to_i.append(j)

        inds = [j for j in inds if j not in similar_to_i]

    return jsonify(thr_struct_lst=thr_struct_lst)


@app.route('/results/<url>/<struct>/<go_file>')
def go_webpage(url, struct, go_file):
    data_n = url_to_data_n(url)
    struct_path = 'go/data' + str(data_n) + '/' + struct

    if go_file == 'GOResultsPROCESS':
        return render_template(struct_path + '/GOResultsPROCESS.html')
    if go_file == 'GOResultsFUNCTION':
        return render_template(struct_path + '/GOResultsFUNCTION.html')
    if go_file == 'GOResultsCOMPONENT':
        return render_template(struct_path + '/GOResultsCOMPONENT.html')

    return render_template('cant_find_data.html')


@app.route('/results/<url>')
def results_panel(url):
    data_n = url_to_data_n(url)
    data_path = os.path.join(ANALYSIS_FOLDER, 'data' + str(data_n))

    if not os.path.exists(data_path):
        return render_template('cant_find_data.html')

    dataset_name = open(os.path.join(data_path, 'dataset_name.txt'), "r").read()
    full_pic_path = os.path.join(ANALYSIS_FOLDER, 'data' + str(data_n), 'structure_layouts')
    GO_flag_path = os.path.join(data_path, 'GO_flag.txt')
    imputation_method_path = os.path.join(data_path, 'imputation_method.txt')
    Jaccard_thr_genes_path = os.path.join(data_path, 'Jaccard_thr_genes.txt')

    if os.path.exists(outcome_path(data_path)):
        outcome = open(outcome_path(data_path), "r").read()
        print('outcome:', outcome)
        if outcome == 'Success':
            # get structure list
            pic_suffs = [file[file.find('_') + 1:] for file in os.listdir(full_pic_path) if
                         file.endswith(".jpg") and 'struct' in file]
            struct_lst = list(set([int(file[:file.find('_')]) for file in pic_suffs]))
            struct_lst = [str(s) for s in struct_lst]

            # get GO_flag
            GO_flag = (open(GO_flag_path, "r").read() == 'True')

            # get imputation method
            impute_method = open(imputation_method_path, "r").read()

            # get the Jaccard_thr_genes and Jaccard_thr_sample_pairs that were used
            Jaccard_thr_genes = open(Jaccard_thr_genes_path, "r").read()
            # Jaccard_thr_sample_pairs = open(os.path.join(data_path, 'Jaccard_thr_sample_pairs.txt'), "r").read()

            return render_template('results_panel.html', data_n=data_n, struct_lst=struct_lst,
                                   impute_method=impute_method,
                                   dataset_name=dataset_name, GO_flag=GO_flag, Jaccard_thr_genes=Jaccard_thr_genes)

        elif outcome == 'No_structures':
            return render_template('no_structures.html')

    return render_template('all_done.html')


############################################################################################################
if __name__ == '__main__':
    app.run(host='0.0.0.0', port=5000)
