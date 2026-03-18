from flask import Flask, render_template, request, Response, send_file, make_response, Blueprint, after_this_request, request, redirect, jsonify
from function_set import *
from io import StringIO
import csv
import uuid
from werkzeug.utils import secure_filename
import re
import logging

main = Blueprint('main', __name__, url_prefix='/HelixHarbor')


logging.basicConfig(stream=sys.stdout, level=logging.DEBUG)
logger = logging.getLogger()


def _has_comparable_values(values):
    if values is None:
        return False
    if isinstance(values, pd.DataFrame):
        return not values.empty
    return len(values.dropna()) > 0



@main.route('/', methods=['GET', 'POST'])
def index():
    
    
    try:
        option_organism = request.form.get("option4")
        option_background = request.form.get("option")
        option_region = request.form.get("option2")
        option_feature = request.form.get("option3")
        
        input_type = request.form.get('input_type')

        orientation = request.form.get("orientation")

        if request.method == 'POST':
            results = {}
            sequence = request.form.get('sequence', '')
            sequence = sequence.replace('\n', ' ')
            list1 = request.form.get('list_1', '')
            list1 = list1.replace(' ', '').replace('\n', '')
            list2 = request.form.get('list_2', '')
            list2 = list2.replace(' ', '').replace('\n', '')
            

            if input_type == 'sequence':

                if is_valid_amino_acid_sequence(sequence):
                    df = run_for_sequence(sequence)

                    # Full path for the file
                    file_path = os.path.join('Gene_list_predictions', 'sequence_output.tsv')

                    # Save the DataFrame
                    df.to_csv(file_path, sep='\t', index=False)

                    results['sequence'] = df.to_html(classes='myDataFrame')

                    
                    return render_template('home.html', results=results)
                else:
                    results['Error'] = f"Please Enter A Valid Amino Acid Sequence"
                    return render_template('home.html', results=results)  

            elif input_type == 'ps_aac':
                
                uniprot_ids = re.split(r'[,\s]+', sequence.strip())
                uniprot_ids = list(set(uniprot_ids))
                for id in uniprot_ids:
                    if is_valid_uniprot_id(id) != True:
                        results['Error'] = "uniprot id " + id + " is not valid. "
                        return render_template('home.html', results=results)

                temp = ps_aac(uniprot_ids, orientation)
                if temp:
                    results['plot_hydrophobicity'] = temp
                    return render_template('home.html', results=results) 
                results['Error'] = f"An error occurred: No proteins were found in database"
                return render_template('home.html', results=results)          

            elif input_type=='uniprot':
                
                df_background = load_tmh_dbs(option_background, option_organism)
                
                uniprot_ids = re.split(r'[,\s]+', sequence.strip())
                uniprot_ids = list(set(uniprot_ids))
                ''''
                for id in uniprot_ids:
                    if is_valid_uniprot_id(id) != True:
                        results['Error'] = "uniprot id " + id + " is not valid. "
                        return render_template('home.html', results=results)
                '''
                df, other_features = run_for_tmh_list(uniprot_ids, df_background)

                if len(other_features) > 0:
                    results['non_transmembrane'] = other_features
                    
                
                df = df.rename(columns={'ID': 'id'})
                df.rename(columns={'Type': 'type'}, inplace=True)
                df_background.rename(columns={'Type': 'type'}, inplace=True)


                
                
                if option_region == 'option1':
                    selected_second_row_values = extract_second_row_values(df, option_feature, 'list')
                    all_second_row_values = extract_second_row_values(df_background, option_feature, 'bk')

                    if not _has_comparable_values(selected_second_row_values) or not _has_comparable_values(all_second_row_values):
                        results['Error'] = "No comparable transmembrane regions were found for your list or the selected background."
                        return render_template('home.html', results=results)
                   
                    if option_feature == 'option6':
                        results ['plot_hydrophobicity']= aac_density_plot(selected_second_row_values, all_second_row_values)
                    else:    
                        save_density_raw_data(df, df_background, "List of Interest", "Background", option_feature, option_region)
                        results['plot_hydrophobicity'] = denisty_plot(selected_second_row_values, all_second_row_values,"List of Interest", "Background")

                elif option_region == 'option2':

                    selected_tmh_values = extract_tmh_values(df, option_feature, 'list')
                    all_tmh_values = extract_tmh_values(df_background, option_feature, 'bk')

                    if not _has_comparable_values(selected_tmh_values) or not _has_comparable_values(all_tmh_values):
                        results['Error'] = "No comparable transmembrane regions were found for your list or the selected background."
                        return render_template('home.html', results=results)

                    if option_feature == 'option6':
                        results['plot_hydrophobicity'] = aac_density_plot(selected_tmh_values, all_tmh_values)    
                    else:
                        save_density_raw_data(df, df_background, "List of Interest", "Background", option_feature, option_region)
                        results['plot_hydrophobicity'] = denisty_plot(selected_tmh_values, all_tmh_values, "List of Interest", "Background")

                return render_template('home.html', results=results)
            
            elif input_type == 'setcmp':
                uniprot_ids_1 = re.split(r'[,\s]+', list1.strip())
                uniprot_ids_2 = re.split(r'[,\s]+', list2.strip())
                uniprot_ids_1 = [x for x in uniprot_ids_1 if x]
                uniprot_ids_2 = [x for x in uniprot_ids_2 if x]
                if not uniprot_ids_1 or not uniprot_ids_2:
                    results['Error'] = "Please provide both List 1 and List 2 with at least one UniProt ID each."
                    return render_template('home.html', results=results)
                uniprot_ids_1 = list(set(uniprot_ids_1))
                uniprot_ids_2 = list(set(uniprot_ids_2))
                # Concatenate lists to iterate over both
                all_uniprot_ids = uniprot_ids_1 + uniprot_ids_2
                ''''
                for id in all_uniprot_ids:
                    if not is_valid_uniprot_id(id):
                        results['Error'] = "uniprot id " + id + " is not valid."
                        return render_template('home.html', results=results)
                '''
                df_background = load_tmh_dbs('option1', 'option1')
                list1_df, list2_df = run_for_list_cmp(uniprot_ids_1, uniprot_ids_2, df_background)
                
                

            

                if option_region == 'option1':
                    selected_second_row_values = extract_second_row_values(list1_df, option_feature, 'list1')
                    all_second_row_values = extract_second_row_values(list2_df, option_feature, 'list2')

                    if not _has_comparable_values(selected_second_row_values) or not _has_comparable_values(all_second_row_values):
                        results['Error'] = "No comparable transmembrane regions were found for one or both lists."
                        return render_template('home.html', results=results)
                    
                    
                    if option_feature == 'option5':
                        results['plot_hydrophobicity'] = heat_plot(selected_second_row_values)

                    elif option_feature == 'option6':
                        results ['plot_hydrophobicity']= aac_density_plot(selected_second_row_values, all_second_row_values)
                    else:    
                        save_density_raw_data(list1_df, list2_df, "List 1", "List 2", option_feature, option_region)
                        results['plot_hydrophobicity'] = denisty_plot(selected_second_row_values, all_second_row_values, "List 1", "List 2")

                elif option_region == 'option2':

                    selected_tmh_values = extract_tmh_values(list1_df, option_feature, 'list1')
                    all_tmh_values = extract_tmh_values(list2_df, option_feature, 'list2')

                    if not _has_comparable_values(selected_tmh_values) or not _has_comparable_values(all_tmh_values):
                        results['Error'] = "No comparable transmembrane regions were found for one or both lists."
                        return render_template('home.html', results=results)
                   
                    if option_feature == 'option5':
                        results['plot_hydrophobicity'] = heat_plot(selected_tmh_values)
                    elif option_feature == 'option6':
                        results['plot_hydrophobicity'] = aac_density_plot(selected_tmh_values, all_tmh_values)    
                    else:
                        save_density_raw_data(list1_df, list2_df, "List 1", "List 2", option_feature, option_region)
                        results['plot_hydrophobicity'] = denisty_plot(selected_tmh_values, all_tmh_values, "List 1", "List 2")

                return render_template('home.html', results=results)
            
        else:
            return render_template('home.html')

    except Exception as e:
        results = {}
        results['Error'] = f"An error occurred: {str(e)}"
        return render_template('home.html', results=results)        

@main.route('/download_csv')
def download_csv():
    file_path = os.path.join(os.getcwd(), 'Gene_list_predictions', 'sequence_output.tsv')
    df = pd.read_csv(file_path, delimiter= '\t')
    output = StringIO()
    df.to_csv(output,  sep='\t', index=False)
    response = make_response(output.getvalue())
    response.headers["Content-Disposition"] = "attachment; filename=helix_harbor_output.tsv"
    response.headers["Content-type"] = "text/csv"
    return response

@main.route('/generate_job_id', methods=['GET'])
def generate_job_id():
    job_id = str(uuid.uuid4())
    return job_id

@main.route('/download_background', methods=['GET'])
def download_background():
    # Get the selected organism option from the request parameters
    option = request.args.get('option')
    
    # Load the background set data using the load_background function
    df = load_background(option)
    
    # Convert the DataFrame to a CSV string
    csv_str = df.to_csv(index=False)
    
    # Create a response with the CSV string
    response = make_response(csv_str)
    response.headers['Content-Disposition'] = 'attachment; filename=background_set.csv'
    response.headers['Content-Type'] = 'text/csv'
    
    return response

@main.route('/download_custom_sheet', methods=['GET'])
def download_custom_sheet():
    file_path = os.path.join(os.getcwd(), 'DCS', 'amino_acid_values.csv')
    df = pd.read_csv(file_path, delimiter= '\t')
    output = StringIO()
    df.to_csv(output,  sep='\t', index=False)
    response = make_response(output.getvalue())
    response.headers["Content-Disposition"] = "attachment; filename=amino_acid_values.csv"
    response.headers["Content-type"] = "text/csv"
    return response

@main.route('/download_aac', methods=['GET'])
def download_aac():
    file_path = os.path.join(os.getcwd(), 'AAC_values', 'AAC_values_helixharbor.tsv')
    df = pd.read_csv(file_path, delimiter= '\t')
    output = StringIO()
    df.to_csv(output,  sep='\t', index=False)
    response = make_response(output.getvalue())
    response.headers["Content-Disposition"] = "attachment; filename=AAC_values_helixharbor.tsv"
    response.headers["Content-type"] = "text/csv"
    return response    

@main.route('/manual_pdf', methods=['GET'])
def manual_pdf():
    file_path = os.path.join(os.getcwd(), 'HelixHarbor_installation_guide.pdf')
    if not os.path.exists(file_path):
        response = make_response("Manual PDF not found.")
        response.status_code = 404
        return response
    return send_file(
        file_path,
        mimetype='application/pdf',
        as_attachment=False,
        download_name='HelixHarbor_installation_guide.pdf'
    )

@main.route('/download_plot_raw_data', methods=['GET'])
def download_plot_raw_data():
    file_path = os.path.join(os.getcwd(), 'list_report', 'plot_raw_data.tsv')
    if not os.path.exists(file_path):
        response = make_response("No plot raw data available. Generate a plot first.")
        response.status_code = 404
        return response

    df = pd.read_csv(file_path, delimiter='\t')
    output = StringIO()
    df.to_csv(output, sep='\t', index=False)
    response = make_response(output.getvalue())
    response.headers["Content-Disposition"] = "attachment; filename=plot_raw_data.tsv"
    response.headers["Content-type"] = "text/csv"
    return response

@main.route('/upload_custom_sheet', methods=['POST'])
def upload_custom_sheet():
    if request.method == 'POST':
        filename = 'amino_acid_values.csv'
        f = request.files['file']
        if f:
            upload_folder = 'UCS'
            
            os.makedirs(upload_folder, exist_ok=True)

            file_path = os.path.join(upload_folder, secure_filename(filename))
            f.save(file_path)
            
            return 'File uploaded successfully'
        else:
            return 'No file uploaded'

@main.route('/lower_upper_range', methods=['POST'])
def lower_upper_range():
    # Extract values from the request form
    upper_range = request.form.get('upper_range')
    lower_range = request.form.get('lower_range')
    selected_feature = request.form.get('selected_feature')
    selected_method = request.form.get('selected_method')
    print(selected_feature, selected_method)

    feature_column_map = {
        'option1': 'volume',
        'option2': 'SA',  
        'option3': 'bulkiness',
        'option4': 'hydrophobicity'
    }

    column_name = feature_column_map.get(selected_feature, 'volume')  # Default to 'volume' if not found


    if selected_method == 'setcmp':
        file_path = os.path.join(os.getcwd(), 'list_report', 'list_report1.tsv')
        df1 = pd.read_csv(file_path, delimiter= '\t')
        filtered_df1_1 = df1[(df1[column_name] >= int(upper_range))]
        filtered_df2_1 = filtered_df1_1[(filtered_df1_1[column_name] <= int(lower_range))]
        unique_ids_1 = filtered_df2_1['id'].unique().tolist()

        file_path = os.path.join(os.getcwd(), 'list_report', 'list_report2.tsv')
        df2 = pd.read_csv(file_path, delimiter= '\t')
        filtered_df1_2 = df2[(df2[column_name] >= int(upper_range))]
        filtered_df2_2 = filtered_df1_2[(filtered_df1_2[column_name] <= int(lower_range))]
        unique_ids_2 = filtered_df2_2['id'].unique().tolist()

        return jsonify(unique_ids_1=unique_ids_1 , unique_ids_2=unique_ids_2)
        
    elif selected_method == 'uniprot':
        
        file_path = os.path.join(os.getcwd(), 'list_report', 'list_report.tsv')
        df = pd.read_csv(file_path, delimiter= '\t')

        filtered_df1 = df[(df[column_name] >= int(upper_range))]
        filtered_df2 = filtered_df1[(filtered_df1[column_name] <= int(lower_range))]

        # Get unique IDs from the filtered DataFrame
        unique_ids = filtered_df2['id'].unique().tolist()
        
        # Return the list of protein IDs as JSON
        return jsonify(unique_ids=unique_ids)

app = Flask(__name__, template_folder='templates')
app.register_blueprint(main)
