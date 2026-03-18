import io
import logging
import base64
import xml.etree.ElementTree as ET
import requests
from collections import Counter
from Bio import AlignIO
from Bio.Align.Applications import ClustalOmegaCommandline
from flask import Flask, render_template, request
import glob
import matplotlib.pyplot as plt
import sys
import json
import os
import subprocess
import pandas as pd
import seaborn as sns
import scipy.stats as stats
from collections import OrderedDict
import numpy as np
import threading
import matplotlib
matplotlib.use('Agg')  # or 'Qt5Agg', 'Qt4Agg', etc. depending on your backend
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry
from scipy import stats
logging.basicConfig(stream=sys.stdout, level=logging.INFO)
logger = logging.getLogger()

WEBSITE_API = "https://rest.uniprot.org/"
PROTEINS_API = "https://www.ebi.ac.uk/proteins/api"

def is_valid_amino_acid_sequence(sequence):
    valid_amino_acids = set("ACDEFGHIKLMNPQRSTVWYUX")
    
    sequence = sequence.upper()
    
    for char in sequence:
        if char not in valid_amino_acids:
            return False
    
    return True

def ps_aac(protein_ids_of_interest, orientation):
    def calculate_aa_frequencies_varying_lengths(sequences):
        if not sequences:
            return {}  # Return an empty dict if there are no sequences

        # Initialize aa_counts dictionary with empty lists
        aa_counts = {aa: [] for aa in 'ACDEFGHIKLMNPQRSTVWY'}
        max_length = 0  # Keep track of the maximum sequence length

        # Update counts and extend lists as needed
        for seq in sequences:
            seq_length = len(seq)
            max_length = max(max_length, seq_length)

            for i, aa in enumerate(seq):
                if aa in aa_counts:
                    # Extend the counts list with zeros if necessary
                    while len(aa_counts[aa]) <= i:
                        aa_counts[aa].append(0)
                    aa_counts[aa][i] += 1

        # Extend all counts to max_length with zeros for uniformity
        for aa in aa_counts:
            aa_counts[aa].extend([0] * (max_length - len(aa_counts[aa])))

        # Calculate frequencies
        total_seqs = len(sequences)
        aa_frequencies = {}
        for aa, counts in aa_counts.items():
            frequencies = []
            for i in range(max_length):
                # Count sequences that include this position
                seqs_including_position = sum(1 for seq in sequences if len(seq) > i)
                if seqs_including_position > 0:
                    freq = counts[i] / seqs_including_position
                else:
                    freq = 0
                frequencies.append(freq)
            aa_frequencies[aa] = frequencies

        return aa_frequencies

    def process_sequences(df):
        df = df.sort_values(by=['id', 'begin'])

        df['reversed_sequence'] = df['sequence']

        df['prev_orientation'] = df['orientation'].shift(1)
        df['next_orientation'] = df['orientation'].shift(-1)


        if orientation == 'outside':
            mask = (df['orientation'] == 'membrane') & (df['prev_orientation'] == 'inside') & (df['next_orientation'] == 'outside')
        elif orientation == 'inside':
            mask = (df['orientation'] == 'membrane') & (df['prev_orientation'] == 'outside') & (df['next_orientation'] == 'inside')



        df.loc[mask, 'reversed_sequence'] = df.loc[mask, 'sequence'].apply(lambda x: x[::-1])

        df.drop(['prev_orientation', 'next_orientation'], axis=1, inplace=True)

        return df

    # Load the dataset
    data = pd.read_csv(os.path.join('backgrounds_tsv', 'TMalpha_homo_sapiens.csv'))
    
    data['sequence'] = data['sequence'].str.replace('\n', '')

    processed = process_sequences(data)
    processed = processed[processed['id'].isin(protein_ids_of_interest)]
    if processed.empty:
        return False
        
    
    membrane_only = processed[processed['orientation'] == 'membrane']
    if membrane_only.empty:
        return False
    aa_frequencies = calculate_aa_frequencies_varying_lengths(membrane_only['reversed_sequence'].tolist())
    if not aa_frequencies:
        return False

    aa_frequencies_df = pd.DataFrame(aa_frequencies)
    aa_frequencies_df = aa_frequencies_df.T
    aa_frequencies_df.index.name = "amino_acid"
    ps_aac_file_path = os.path.join('list_report', 'plot_raw_data.tsv')
    aa_frequencies_df.to_csv(ps_aac_file_path, sep='\t')

    # Save the plot
    plt.figure(figsize=(14, 10))
    ax = sns.heatmap(aa_frequencies_df, cmap='viridis', annot=False)
    ax.set_xticklabels(range(1, aa_frequencies_df.shape[1] + 1))
    plt.title('Heatmap of Amino Acid Frequencies by Position')
    plt.xlabel('Position')
    plt.ylabel('Amino Acid')

    img_bytes = io.BytesIO()
    plt.savefig(img_bytes, format='png')
    img_bytes.seek(0)

    # encode the byte stream as a base64 string
    img_base64 = base64.b64encode(img_bytes.getvalue()).decode('ascii')

    # Clear the current figure to free memory
    plt.clf()
    plt.close()

    return img_base64


def fetch_uniprot_sequence(uniprot_id):
    url = f"https://www.uniprot.org/uniprot/{uniprot_id}.fasta"
    response = requests.get(url)
    sequence = "".join(response.text.split("\n")[1:])
    return sequence


def retrieve_protein_features(accession_number):
    url = f"https://rest.uniprot.org/uniprot/{accession_number}.xml"
    
    # set up retry strategy
    retry_strategy = Retry(
        total=5,  # Total number of retries
        status_forcelist=[429, 500, 502, 503, 504],
        method_whitelist=["HEAD", "GET", "OPTIONS"],
        backoff_factor=1  # Wait time between retries (exponential backoff)
    )

    # create an adapter with the retry strategy
    adapter = HTTPAdapter(max_retries=retry_strategy)

    # create a session and mount the adapter
    http = requests.Session()
    http.mount("https://", adapter)
    
    try:
        response = http.get(url, timeout=10)  # Adjust the timeout value as needed
        response.raise_for_status()  # Check for any HTTP errors

        root = ET.fromstring(response.content)
        namespace = {'uniprot': 'http://uniprot.org/uniprot'}
        protein_features = []
        other_features = set()
        sequence = None

     
        sequence = fetch_uniprot_sequence(accession_number)
        # Look for the 'feature' tags
        for feature in root.findall('.//uniprot:feature', namespace):
            feature_type = feature.get('type')
            if feature_type in ['transmembrane region', 'topological domain', 'intramembrane region']:
                print("hello")
                begin_element = feature.find('uniprot:location/uniprot:begin', namespace)
                end_element = feature.find('uniprot:location/uniprot:end', namespace)
                
                if begin_element is not None and end_element is not None:
                    begin_position = begin_element.get('position')
                    end_position = end_element.get('position')

                    trimmed_sequence = sequence[int(begin_position)-1:int(end_position)]

                    protein_features.append({
                        "id": accession_number,
                        "type": feature_type,
                        "begin": begin_position,
                        "end": end_position,
                        "full_sequence": sequence  # Add sequence to each feature
                    })
            else:
                other_features.add(accession_number)
        
        return protein_features, list(other_features)
    
    except requests.exceptions.RequestException as e:
        print(f"An error occurred: {e}")
        return [], []




def read_aa_values_from_csv(file_path):
    # Read the CSV file into a DataFrame
    df = pd.read_csv(file_path)

    # Ensure the DataFrame has the expected columns
    if 'AA' in df.columns and 'Value' in df.columns:
        # Convert DataFrame to a dictionary
        aa_values_dict = pd.Series(df.Value.values, index=df.AA).to_dict()
        return aa_values_dict
    else:
        raise ValueError("CSV must contain 'AA' and 'Value' columns")


def sliding_window_scores(sequence, index_values, window_size=1, normalize_values=False):
    scores = []
    window_size = len(sequence)
    for i in range(len(sequence) - window_size + 1):
        window = sequence[i:i+window_size]
        score = sum(index_values.get(aa, 0) for aa in window) / window_size
        scores.append(score)

    return normalize(scores) if normalize_values else scores

# Function to normalize bulkiness values


def normalize(values):
    min_val = min(values)
    max_val = max(values)
    range_val = max_val - min_val
    return [(val - min_val) / range_val for val in values]


def extract_second_row_values(df, feature, label):
    if not df.index.is_unique:
        df = df.reset_index(drop=True)
    
    # Filter the DataFrame for rows where 'type' is 'transmembrane regions'
    transmembrane_df = df.query("type == 'transmembrane region' or type == 'intramembrane region'")


    # For each 'id', get the first row where 'type' is 'transmembrane regions'
    first_transmembrane_rows = transmembrane_df.groupby('id').first().reset_index()

    # Save the DataFrame to a file based on the label
    if label == 'list1':
        file_path = os.path.join('list_report', 'list_report1.tsv')
    elif label == 'list2':
        file_path = os.path.join('list_report', 'list_report2.tsv')
    else:
        file_path = os.path.join('list_report', 'list_report.tsv')
    first_transmembrane_rows.to_csv(file_path, sep='\t', index=False)

    if first_transmembrane_rows.empty:
        if feature == 'option6':
            return pd.DataFrame(columns=['id', 'sequence'])
        return pd.Series(dtype=float)

    # Depending on the feature requested, select and return the appropriate column(s)
    if feature == 'option1':
        return first_transmembrane_rows['volume'].dropna()
    elif feature == 'option2':
        return first_transmembrane_rows['SA'].dropna()
    elif feature == 'option3':
        return first_transmembrane_rows['bulkiness'].dropna()
    elif feature == 'option4':
        return first_transmembrane_rows['hydrophobicity'].dropna()
    elif feature == 'option6':
        return first_transmembrane_rows[['id', 'sequence']].dropna()
    else:
        # Apply a custom scale based on 'sequence'
        first_transmembrane_rows['custom_scale'] = first_transmembrane_rows['sequence'].apply(
            lambda x: sliding_window_scores(
                x, 
                read_aa_values_from_csv(os.path.join(os.getcwd(), 'UCS', 'amino_acid_values.csv')), 
                window_size=1, 
                normalize_values=False
            )
        )
        first_transmembrane_rows['custom_scale'] = first_transmembrane_rows['custom_scale'].apply(lambda x: x[0] if x else None)
        return first_transmembrane_rows['custom_scale'].dropna()



def is_valid_uniprot_id(uniprot_id):
    url = f"https://www.uniprot.org/uniprot/{uniprot_id}.fasta"
    response = requests.get(url)
    return response.ok


def run_for_sequence(sequence):
    results = {}
    run_tmbed_sequnce(sequence)
    # Parse the TMbed tool output text into structured data using the function
    parsed_data = []
    file_path = os.path.join(os.getcwd(), 'tmbed', 'sample.pred')
    with open(file_path) as file:
        file_content = file.read()
        entries = file_content.strip().split('>')[1:]
    for entry in entries:
        lines = entry.strip().split('\n')
        if len(lines) < 3:
            print(f"Skipping an entry due to unexpected format: {lines[0]}")
            continue  # Skip this entry
        header = lines[0]
        sequence = lines[1]
        annotation = lines[2]
        if header != 'user_sequence':
            protein_id = header.split('|')[1]
        else:
            protein_id = "user_sequence"

        
        parsed_data.extend(parse_tmbed_annotation(
            protein_id, sequence, annotation))

    # Create a DataFrame from the parsed data
    result_df = pd.DataFrame(parsed_data)

    result_df['SA'] = result_df['sequence'].apply(lambda x: sliding_window_scores(
        x, surface_area_index_theo, window_size=1, normalize_values=False))
    result_df['volume'] = result_df['sequence'].apply(lambda x: sliding_window_scores(
        x, volume_index, window_size=1, normalize_values=False))
    result_df['bulkiness'] = result_df['sequence'].apply(lambda x: sliding_window_scores(
        x, bulkiness_index, window_size=1, normalize_values=False))
    result_df['hydrophobicity'] = result_df['sequence'].apply(kd_scale_avg)
    result_df['volume'] = result_df['volume'].apply(
        lambda x: x[0] if x else None)
    result_df['SA'] = result_df['SA'].apply(lambda x: x[0] if x else None)
    result_df['bulkiness'] = result_df['bulkiness'].apply(
        lambda x: x[0] if x else None)

    # Reset the index of the result DataFrame if needed
    result_df.reset_index(drop=True, inplace=True)

    results = result_df

    return results


def get_annotations():
    annotaions = list()
    with open('seq_pred_web/region_output.gff3') as f:
        for line in f:
            if not (line.startswith('#')):
                parts = line.split('\t')
                start = int(parts[3])
                end = int(parts[4])
                annotaions.append(start)
                annotaions.append(end)
    return annotaions




def denisty_plot(data1, data2, label1, label2):
    if data1.empty:
        print("data1 is empty after dropping NaNs.")
    if data2.empty:
        print("data2 is empty after dropping NaNs.")    




    data1 = data1.dropna()
    data2 = data2.dropna()
    
    # Print out the number of NaN values in each dataset
    print(f"NaNs in {label1}: {data1.isna().sum()}")
    print(f"NaNs in {label2}: {data2.isna().sum()}")
    
    
    # Set the clip parameter
    clip = (-10, None)
    
    # Create the plot
    sns.kdeplot(data1, clip=clip, label=label1)
    sns.kdeplot(data2, clip=clip, label=label2)
    
    # Perform Kolmogorov-Smirnov (KS) test
    ks_statistic, ks_p_value = stats.ks_2samp(data1, data2)
    
    # Perform t-test
    t_statistic, t_p_value = stats.ttest_ind(data1, data2, equal_var=False)
    
    # Perform Wilcoxon rank-sum test (Mann-Whitney U test)
    wilcoxon_statistic, wilcoxon_p_value = stats.mannwhitneyu(data1, data2, alternative='two-sided')
    
    # Annotate the plot with the test results
    plt.text(0.05, 0.95, f'KS Test: p-value = {ks_p_value:.2e}', transform=plt.gca().transAxes, fontsize=9, verticalalignment='top')
    plt.text(0.05, 0.90, f'T-Test: p-value = {t_p_value:.2e}', transform=plt.gca().transAxes, fontsize=9, verticalalignment='top')
    plt.text(0.05, 0.85, f'Wilcoxon Test: p-value = {wilcoxon_p_value:.2e}', transform=plt.gca().transAxes, fontsize=9, verticalalignment='top')
    
    plt.legend()
    
    # Store the plot in a byte stream
    img_bytes = io.BytesIO()
    plt.savefig(img_bytes, format='png')
    img_bytes.seek(0)
    
    # Encode the byte stream as a base64 string
    img_base64 = base64.b64encode(img_bytes.getvalue()).decode('ascii')
    
    plt.clf()
    
    return img_base64


def _prepare_density_raw_data(df, feature, region):
    work_df = df.copy()
    if not work_df.index.is_unique:
        work_df = work_df.reset_index(drop=True)

    tmh_rows = work_df.query("type == 'transmembrane region' or type == 'intramembrane region'").copy()
    if region == 'option1':
        tmh_rows = tmh_rows.groupby('id').first().reset_index()

    feature_map = {
        'option1': 'volume',
        'option2': 'SA',
        'option3': 'bulkiness',
        'option4': 'hydrophobicity'
    }

    if feature in feature_map:
        tmh_rows['value'] = tmh_rows[feature_map[feature]]
    else:
        tmh_rows['value'] = tmh_rows['sequence'].apply(
            lambda x: sliding_window_scores(
                x,
                read_aa_values_from_csv(os.path.join(os.getcwd(), 'UCS', 'amino_acid_values.csv')),
                window_size=1,
                normalize_values=False
            )
        )
        tmh_rows['value'] = tmh_rows['value'].apply(lambda x: x[0] if x else None)

    output_columns = ['id', 'type', 'begin', 'end', 'value']
    for col in output_columns:
        if col not in tmh_rows.columns:
            tmh_rows[col] = None

    return tmh_rows[output_columns].dropna(subset=['value'])


def save_density_raw_data(df1, df2, label1, label2, feature, region):
    mapped_df1 = _prepare_density_raw_data(df1, feature, region)
    mapped_df2 = _prepare_density_raw_data(df2, feature, region)

    mapped_df1.insert(0, 'group', label1)
    mapped_df2.insert(0, 'group', label2)

    raw_df = pd.concat([mapped_df1, mapped_df2], ignore_index=True)
    file_path = os.path.join('list_report', 'plot_raw_data.tsv')
    raw_df.to_csv(file_path, sep='	', index=False)


def heat_plot(selected_second_row_values):

    df_aac = get_aac(selected_second_row_values[['id', 'sequence']])
    file_path = os.path.join('list_report', 'plot_raw_data.tsv')
    df_aac.to_csv(file_path, sep='	', index=False)
    filtered_ids = df_aac[(df_aac['S'] >= 0.5)]['id']
    print(filtered_ids)
    df_heat_plot = df_aac.set_index('id')

    df_heat_plot.index = df_heat_plot.index.map(str)

    # Create the heatmap
    plt.figure(figsize=(15, 10))  # Adjust the figure size as needed
    ax = sns.heatmap(df_heat_plot, annot=False, cmap="viridis", yticklabels=True)
    ax.set_xticklabels(range(1, df_heat_plot.shape[1] + 1))
    plt.xticks(rotation=90)  # Rotate x labels if needed
    plt.yticks(rotation=0)  # Ensure y labels are horizontal for readability

    # Aesthetics
    plt.title('Amino Acid Composition Heatmap')
    plt.xlabel('Amino Acids')
    plt.ylabel('IDs')

    img_bytes = io.BytesIO()
    plt.savefig(img_bytes, format='png')
    img_bytes.seek(0)

    # encode the byte stream as a base64 string
    img_base64 = base64.b64encode(img_bytes.getvalue()).decode('ascii')

    # Clear the current figure to free memory
    plt.clf()
    plt.close()
    return img_base64


def get_sequence(entry):
    with open("TMH_sequence_fasta/sequence.fasta", "w") as file:
        r = get_url(f"https://rest.uniprot.org/uniprotkb/{entry}.fasta")
        file.write(r.text)


def read_fasta(fasta_file_name):
    fasta_data = pd.DataFrame(columns=['discrp', 'seq'])
    current_sequence = ""
    last_header = ""
    i = 0
    with open(fasta_file_name) as fasta_file:
        for line in fasta_file:
            line = line.strip()
            if line.startswith(">"):
                if current_sequence == "":
                    last_header = line
                    continue
                else:
                    fasta_data.loc[i] = [last_header.strip('>').split(" ", 1)[
                        0]] + [current_sequence]
                    i += 1
                    current_sequence = ""
                    last_header = line
            else:
                current_sequence += line
        if current_sequence != "":
            fasta_data.loc[i] = [last_header.strip('>').split(" ", 1)[
                0]] + [current_sequence]
            i += 1
        return fasta_data


def get_url(url, **kwargs):
    response = requests.get(url, **kwargs)

    if not response.ok:
        print(response.text)
        response.raise_for_status()
        sys.exit()

    return response





def kd_scale(seq):
    kd = {'A': 1.8, 'R': -4.5, 'N': -3.5, 'D': -3.5, 'C': 2.5,
          'Q': -3.5, 'E': -3.5, 'G': -0.4, 'H': -3.2, 'I': 4.5,
          'L': 3.8, 'K': -3.9, 'M': 1.9, 'F': 2.8, 'P': -1.6,
          'S': -0.8, 'T': -0.7, 'W': -0.9, 'Y': -1.3, 'V': 4.2}

    kd_sum = 0.0
    for aa in seq:
        if aa in kd:
            kd_sum += kd[aa]

    return round(kd_sum, 2)


def kd_scale_avg(seq):
    kd = {'A': 1.8, 'R': -4.5, 'N': -3.5, 'D': -3.5, 'C': 2.5,
          'Q': -3.5, 'E': -3.5, 'G': -0.4, 'H': -3.2, 'I': 4.5,
          'L': 3.8, 'K': -3.9, 'M': 1.9, 'F': 2.8, 'P': -1.6,
          'S': -0.8, 'T': -0.7, 'W': -0.9, 'Y': -1.3, 'V': 4.2}

    kd_sum = 0.0
    for aa in seq:
        if aa in kd:
            kd_sum += kd[aa]

    return round((kd_sum/len(seq)), 2)


def ges_scale_avg(seq):
    ges = {'A': 1.6, 'R': -12.3, 'N': -4.8, 'D': -9.2, 'C': 2.0,
           'Q': -4.10, 'E': -8.2, 'G': 1.0, 'H': -3.0, 'I': 3.10,
           'L': 2.8, 'K': -8.8, 'M': 3.4, 'F': 3.7, 'P': -2.0,
           'S': -0.6, 'T': 1.2, 'W': 1.9, 'Y': -0.7, 'V': 2.6}

    ges_sum = 0.0
    for aa in seq:
        if aa in ges:
            ges_sum += ges[aa]

    return round((ges_sum/len(seq)), 2)


def net_charge(seq):
    charge = {'A': 0, 'R': 1, 'N': 0, 'D': -1, 'C': 0,
              'Q': 0, 'E': -1, 'G': 0, 'H': 0, 'I': 0,
              'L': 0, 'K': 1, 'M': 0, 'F': 0, 'P': 0,
              'S': 0, 'T': 0, 'W': 0, 'Y': 0, 'V': 0}

    total_charge = 0.0
    for aa in seq:
        if aa in charge:
            total_charge += charge[aa]

    return round(total_charge, 1)


def zm_polarity(seq):
    zm = {'A': 0, 'R': 53.00, 'N': 3.38, 'D': 49.70, 'C': 1.48,
          'Q': 3.35, 'E': 49.90, 'G': 0, 'H': 51.60, 'I': 0.13,
          'L': 0.13, 'K': 49.50, 'M': 1.43, 'F': 0.35, 'P': 1.58,
          'S': 1.67, 'T': 1.66, 'W': 2.10, 'Y': 1.61, 'V': 0.13}

    zm_sum = 0.0
    for aa in seq:
        if aa in zm:
            zm_sum += zm[aa]

    return round((zm_sum/len(seq)), 2)


def zm_bulkiness(seq):
    zm = {
        'A': 11.50, 'R': 14.28, 'D': 11.68, 'N': 12.82, 'C': 13.46,
        'E': 13.57, 'Q': 14.45, 'G': 3.40, 'H': 13.69, 'L': 21.40,
        'I': 21.40, 'K': 15.71, 'M': 16.25, 'F': 19.80, 'P': 17.43,
        'S': 9.47, 'T': 15.77, 'W': 21.67, 'Y': 18.03, 'V': 21.57}

    zm_sum = 0.0
    for aa in seq:
        if aa in zm:
            zm_sum += zm[aa]

    return round((zm_sum/len(seq)), 2)


def calc_aac(seq_id_pair):
    seq_id, seq = seq_id_pair
    aac_dict = OrderedDict({x: 0 for x in list("ACDEFGHIKLMNPQRSTVWY")})
    for c in seq:
        aac_dict[c] += 1
    aac_values = np.divide(np.array(list(aac_dict.values())), len(seq))
    return [seq_id] + aac_values.tolist()


def get_aac(df):
    aac_data = df.apply(calc_aac, axis=1).tolist()
    aac_df = pd.DataFrame(
        data=aac_data,
        index=df.index,
        columns=['id'] + list("ACDEFGHIKLMNPQRSTVWY")
    )
    return aac_df


def normalize_tmh_dataframe(df):
    normalized = df.copy()

    unnamed_columns = [col for col in normalized.columns if str(col).startswith('Unnamed:')]
    if unnamed_columns:
        normalized = normalized.drop(columns=unnamed_columns)

    rename_map = {
        'ID': 'id',
        'Type': 'type',
        'Orientation': 'orientation',
        'Begin': 'begin',
        'End': 'end',
        'Full_sequence': 'full_sequence',
        'Sequence': 'sequence',
    }
    normalized = normalized.rename(columns=rename_map)

    for canonical in set(rename_map.values()):
        duplicate_columns = [col for col in normalized.columns if col == canonical]
        if len(duplicate_columns) > 1:
            merged = normalized.loc[:, duplicate_columns].bfill(axis=1).iloc[:, 0]
            normalized = normalized.drop(columns=duplicate_columns)
            normalized[canonical] = merged

    return normalized



def load_background(option):
    option = str(option)
    if option == 'option1':
        file_path = os.path.join(os.getcwd(), 'backgrounds_tsv', 'PositiveNegative_homo_sapiens.csv')
    elif option == 'option2':
        file_path = os.path.join(os.getcwd(), 'backgrounds_tsv', 'PositiveNegative.csv')
    else:
        file_path = os.path.join(os.getcwd(), 'backgrounds_tsv', 'PositiveNegative_homo_sapiens.csv')
    df = normalize_tmh_dataframe(pd.read_csv(file_path))
    if 'sequence' in df.columns:
        df['sequence'] = df['sequence'].astype(str)
    return df


def load_tmh_dbs(option_background, option_organism):
    option_background = str(option_background)
    option_organism = str(option_organism)
    if option_organism == 'option1':
        if option_background == 'option1':
            file_path = os.path.join(os.getcwd(), 'backgrounds_tsv', 'PositiveNegative_homo_sapiens.csv')
            df_uniprot = normalize_tmh_dataframe(pd.read_csv(file_path))
            df_uniprot['sequence'] = df_uniprot['sequence'].astype(str)
            return df_uniprot
        elif option_background == 'option2':
            file_path = os.path.join(os.getcwd(), 'backgrounds_tsv', 'TopDB.csv')
            df_topdb = normalize_tmh_dataframe(pd.read_csv(file_path))
            df_topdb['sequence'] = df_topdb['sequence'].astype(str)
            return df_topdb
        elif option_background == 'option3':
            file_path = os.path.join(os.getcwd(), 'backgrounds_tsv', 'TMalpha_homo_sapiens.csv')
            df_tmaplha = normalize_tmh_dataframe(pd.read_csv(file_path))
            df_tmaplha['sequence'] = df_tmaplha['sequence'].astype(str)
            return df_tmaplha
    elif option_organism == 'option2':
        if option_background == 'option1':
            file_path = os.path.join(os.getcwd(), 'backgrounds_tsv', 'PositiveNegative.csv')
            df_uniprot = normalize_tmh_dataframe(pd.read_csv(file_path))
            df_uniprot['sequence'] = df_uniprot['sequence'].astype(str)
            return df_uniprot
        elif option_background == 'option2':
            file_path = os.path.join(os.getcwd(), 'backgrounds_tsv', 'TopDB.csv')
            df_topdb = normalize_tmh_dataframe(pd.read_csv(file_path))
            df_topdb['sequence'] = df_topdb['sequence'].astype(str)
            return df_topdb
        elif option_background == 'option3':
            file_path = os.path.join(os.getcwd(), 'backgrounds_tsv', 'TMAlpha.csv')
            df_tmaplha = normalize_tmh_dataframe(pd.read_csv(file_path))
            df_tmaplha['sequence'] = df_tmaplha['sequence'].astype(str)
            return df_tmaplha

def run_for_tmh_list(uniprot_ids, df_background):

    def extract_sequence(row):
        # Convert 'Begin' and 'End' to integers
        begin_index = int(row['begin']) - 1  # Adjust for 1-based indexing
        end_index = int(row['end'])
        return row['full_sequence'][begin_index:end_index]

    result_df = pd.DataFrame(columns=["id", "begin", "end", "full_sequence", "sequence", "SA", "volume", "bulkiness", "hydrophobicity"])
    proteins_count_bk = 0
    proteins_count_up = 0
    other_features = set()
    for uniprot_id in uniprot_ids:
        if uniprot_id in df_background['id'].values:
            rows = normalize_tmh_dataframe(df_background[df_background['id'] == uniprot_id])
            result_df = pd.concat([result_df, rows], ignore_index=True)
            logger.info("found in background")
            logger.info("-------------------")
            proteins_count_bk +=1
            
        else:
            protein_features, other_features = retrieve_protein_features(uniprot_id)
            proteins_count_up += 1
            if protein_features:
                
                new_rows = normalize_tmh_dataframe(pd.DataFrame(protein_features))
                result_df = pd.concat([result_df, new_rows], ignore_index=True)
                logger.info("extracted from UniProt")
                logger.info("-------------------")
                
                
    logger.info('in background: %s, in uniprot: %s', proteins_count_bk, proteins_count_up)

                
    result_df.dropna(subset=['begin'], inplace=True)
    result_df.dropna(subset=['end'], inplace=True)
    
    result_df.dropna(subset=['full_sequence'], inplace=True)
    # Apply the function to each row
    result_df['sequence'] = result_df.apply(extract_sequence, axis=1)
    result_df.dropna(subset=['sequence'], inplace=True)
    result_df = result_df[result_df['sequence'] != '']
    result_df['SA'] = result_df['sequence'].apply(lambda x: sliding_window_scores(x, surface_area_index_theo, window_size=1, normalize_values=False))
    result_df['volume'] = result_df['sequence'].apply(lambda x: sliding_window_scores( x, volume_index, window_size=1, normalize_values=False))
    result_df['bulkiness'] = result_df['sequence'].apply(lambda x: sliding_window_scores(x, bulkiness_index, window_size=1, normalize_values=False))
    result_df['hydrophobicity'] = result_df['sequence'].apply(kd_scale_avg)
    result_df['volume'] = result_df['volume'].apply(lambda x: x[0] if x else None)
    result_df['SA'] = result_df['SA'].apply(lambda x: x[0] if x else None)
    result_df['bulkiness'] = result_df['bulkiness'].apply(lambda x: x[0] if x else None)
    result_df.rename(columns={'Type': 'type'}, inplace=True)
    result_df = result_df.dropna(subset=['sequence'])
    result_df = result_df.dropna(subset=['full_sequence'])


    result_df.reset_index(drop=True, inplace=True)
    return result_df, other_features


def run_for_list_cmp(uniprot_ids_1, uniprot_ids_2, df_background):

    def extract_sequence(row):
        full_sequence = row.get('full_sequence')
        if pd.isna(full_sequence) or not str(full_sequence):
            return None

        try:
            begin_index = int(row['begin']) - 1  # Adjust for 1-based indexing
            end_index = int(row['end'])
        except (TypeError, ValueError):
            return None

        if begin_index < 0 or end_index <= begin_index:
            return None

        return full_sequence[begin_index:end_index]

    def populate_missing_full_sequences(result_df):
        if result_df.empty or 'full_sequence' not in result_df.columns:
            return result_df

        missing_full_sequence = result_df['full_sequence'].isna() | (result_df['full_sequence'].astype(str) == '')
        missing_ids = result_df.loc[missing_full_sequence, 'id'].dropna().unique()
        if len(missing_ids) == 0:
            return result_df

        sequences = {uniprot_id: fetch_uniprot_sequence(uniprot_id) for uniprot_id in missing_ids}
        for uniprot_id, seq in sequences.items():
            mask = (result_df['id'] == uniprot_id) & missing_full_sequence
            result_df.loc[mask, 'full_sequence'] = seq
        return result_df

    result_df_1 = pd.DataFrame(columns=["id", "begin", "end", "full_sequence", "sequence", "SA", "volume", "bulkiness", "hydrophobicity"])
    result_df_2 = pd.DataFrame(columns=["id", "begin", "end", "full_sequence", "sequence", "SA", "volume", "bulkiness", "hydrophobicity"])
    


    for uniprot_id in uniprot_ids_1:
        if uniprot_id in df_background['id'].values:
            rows = normalize_tmh_dataframe(df_background[df_background['id'] == uniprot_id])
            result_df_1 = pd.concat([result_df_1, rows], ignore_index=True)
            print("found in background")
            print('----------------')
            
        else:
            protein_features, other_features = retrieve_protein_features(uniprot_id)
            if protein_features:
                new_rows = normalize_tmh_dataframe(pd.DataFrame(protein_features))
                result_df_1 = pd.concat([result_df_1, new_rows], ignore_index=True)
                print("extracted from UniProt")
                print('----------------')
                
    
    
    for uniprot_id in uniprot_ids_2:
        if uniprot_id in df_background['id'].values:
            rows = normalize_tmh_dataframe(df_background[df_background['id'] == uniprot_id])
            result_df_2 = pd.concat([result_df_2, rows], ignore_index=True)
            print("found in background")
            print('----------------')
           
        else:
            protein_features, other_features = retrieve_protein_features(uniprot_id)
            if protein_features:
                new_rows = normalize_tmh_dataframe(pd.DataFrame(protein_features))
                result_df_2 = pd.concat([result_df_2, new_rows], ignore_index=True)
                print("extracted from UniProt")
                print('----------------')
                
    result_df_1 = populate_missing_full_sequences(result_df_1)
    result_df_2 = populate_missing_full_sequences(result_df_2)

    result_df_1['sequence'] = result_df_1.apply(extract_sequence, axis=1)
    result_df_1 = result_df_1.dropna(subset=['sequence'])
    result_df_1 = result_df_1[result_df_1['sequence'] != '']
    result_df_1['SA'] = result_df_1['sequence'].apply(lambda x: sliding_window_scores(x, surface_area_index_theo, window_size=1, normalize_values=False))
    result_df_1['volume'] = result_df_1['sequence'].apply(lambda x: sliding_window_scores( x, volume_index, window_size=1, normalize_values=False))
    result_df_1['bulkiness'] = result_df_1['sequence'].apply(lambda x: sliding_window_scores(x, bulkiness_index, window_size=1, normalize_values=False))
    result_df_1['hydrophobicity'] = result_df_1['sequence'].apply(kd_scale_avg)
    result_df_1['volume'] = result_df_1['volume'].apply(lambda x: x[0] if x else None)
    result_df_1['SA'] = result_df_1['SA'].apply(lambda x: x[0] if x else None)
    result_df_1['bulkiness'] = result_df_1['bulkiness'].apply(lambda x: x[0] if x else None)
    result_df_1.reset_index(drop=True, inplace=True)

    result_df_2['sequence'] = result_df_2.apply(extract_sequence, axis=1)
    result_df_2 = result_df_2.dropna(subset=['sequence'])
    result_df_2 = result_df_2[result_df_2['sequence'] != '']
    result_df_2['SA'] = result_df_2['sequence'].apply(lambda x: sliding_window_scores(x, surface_area_index_theo, window_size=1, normalize_values=False))
    result_df_2['volume'] = result_df_2['sequence'].apply(lambda x: sliding_window_scores( x, volume_index, window_size=1, normalize_values=False))
    result_df_2['bulkiness'] = result_df_2['sequence'].apply(lambda x: sliding_window_scores(x, bulkiness_index, window_size=1, normalize_values=False))
    result_df_2['hydrophobicity'] = result_df_2['sequence'].apply(kd_scale_avg)
    result_df_2['volume'] = result_df_2['volume'].apply(lambda x: x[0] if x else None)
    result_df_2['SA'] = result_df_2['SA'].apply(lambda x: x[0] if x else None)
    result_df_2['bulkiness'] = result_df_2['bulkiness'].apply(lambda x: x[0] if x else None)
    result_df_2.reset_index(drop=True, inplace=True)
    result_df_1.rename(columns={'Type': 'type'}, inplace=True)
    result_df_2.rename(columns={'Type': 'type'}, inplace=True)
    print(result_df_2.columns, result_df_1.columns)

    return result_df_1, result_df_2    



# Define a function to parse the annotation into structured data
def parse_tmbed_annotation(protein_id, sequence, annotation):

    # Define the mappings for type and orientation based on the annotation characters
    type_mapping = {'B': 'Transmembrane beta strand', 'H': 'Transmembrane alpha helix',
                    'S': 'Signal peptide', 'i': 'Non-Transmembrane', 'o': 'Non-Transmembrane'}
    orientation_mapping = {'i': 'inside', 'o': 'outside',
                           'B': 'membrane', 'H': 'membrane', 'S': 'membrane'}

    # Initialize variables to keep track of the current annotation type and its start position
    current_type = ''
    current_start = 0
    parsed_data = []

    # Iterate over the annotation characters and sequence simultaneously
    for i, (anno_char, seq_char) in enumerate(zip(annotation, sequence)):
        anno_type = type_mapping.get(anno_char, 'Unknown')
        # Check if the current character is different from the last one, indicating a new region
        if anno_char != current_type:
            # If there was a previous region, save it
            if current_type:
                parsed_data.append({
                    'id': protein_id,
                    'sequence': sequence[current_start:i],
                    'type': type_mapping[current_type],
                    'orientation': orientation_mapping[current_type],
                    'begin': current_start + 1,  # 1-based indexing
                    'end': i,
                    'sequence_length': i - current_start
                })
            # Update the current region's start position and type
            current_start = i
            current_type = anno_char

    # Add the last region if the sequence ended
    if current_type and current_start < len(sequence):
        parsed_data.append({
            'id': protein_id,
            'sequence': sequence[current_start:],
            'type': type_mapping[current_type],
            'orientation': orientation_mapping[current_type],
            'begin': current_start + 1,  # 1-based indexing
            'end': len(sequence),
            'sequence_length': len(sequence) - current_start
        })

    return parsed_data


def run_tmbed_sequnce(sequence):

    file_path = os.path.join(os.getcwd(), 'user_seq', 'user_seq.fasta')
    with open(file_path, "w") as file:
        file.write(">user_sequence\n")  # Default sequence name
        file.write(sequence + "\n")

    working_directory = os.path.join(os.getcwd(), 'tmbed')
    cmd = f"python -m tmbed predict -f {file_path} -p sample.pred --out-format=1"

    try:
        # Run the command with a timeout (e.g., 5 minutes)
        output = subprocess.run(
            cmd, shell=True, capture_output=True, text=True, timeout=300, cwd=working_directory)

        # Print standard output and standard error
        print("Standard Output:\n", output.stdout)
        print("Standard Error:\n", output.stderr)
    except subprocess.TimeoutExpired:
        print("The command took too long to run and was terminated.")
    except subprocess.CalledProcessError as e:
        # This will catch other errors related to the subprocess
        print(f"An error occurred: {e}")


bulkiness_index = {
    'A': 11.500, 'R': 14.280, 'N': 12.820,
    'D': 11.680, 'C': 13.460, 'Q': 14.450,
    'E': 13.570, 'G': 3.400,  'H': 13.690,
    'I': 21.400, 'L': 21.400, 'K': 15.710,
    'M': 16.250, 'F': 19.800, 'P': 17.430,
    'S': 9.470,  'T': 15.770, 'W': 21.670,
    'Y': 18.030, 'V': 21.570
}

volume_index = {
    'A': 88.6,   # Alanine
    'R': 173.4,  # Arginine
    'N': 114.1,  # Asparagine
    'D': 111.1,  # Aspartic acid
    'C': 108.5,  # Cysteine
    'Q': 143.8,  # Glutamine
    'E': 138.4,  # Glutamic Acid
    'G': 60.1,   # Glycine
    'H': 153.2,  # Histidine
    'I': 166.7,  # Isoleucine
    'L': 166.7,  # Leucine
    'K': 168.6,  # Lysine
    'M': 162.9,  # Methionine
    'F': 189.9,  # Phenylalanine
    'P': 112.7,  # Proline
    'S': 89.0,   # Serine
    'T': 116.1,  # Threonine
    'W': 227.8,  # Tryptophan
    'Y': 193.6,  # Tyrosine
    'V': 140.0,  # Valine
    'B': 111.1,  # Asparagine or Aspartic acid (ambiguous)
    'Z': 143.8   # Glutamine or Glutamic acid (ambiguous)
}
# Dictionary for Tien et al. 2013 (Theoretical)
surface_area_index_theo = {
    "A": 129.0,  # Alanine
    "R": 274.0,  # Arginine
    "N": 195.0,  # Asparagine
    "D": 193.0,  # Aspartate
    "C": 167.0,  # Cysteine
    "E": 223.0,  # Glutamate
    "Q": 225.0,  # Glutamine
    "G": 104.0,  # Glycine
    "H": 224.0,  # Histidine
    "I": 197.0,  # Isoleucine
    "L": 201.0,  # Leucine
    "K": 236.0,  # Lysine
    "M": 224.0,  # Methionine
    "F": 240.0,  # Phenylalanine
    "P": 159.0,  # Proline
    "S": 155.0,  # Serine
    "T": 172.0,  # Threonine
    "W": 285.0,  # Tryptophan
    "Y": 263.0,  # Tyrosine
    "V": 174.0,  # Valine
}

# Dictionary for Tien et al. 2013 (Empirical)
surface_area_index_emp = {
    "A": 121.0,  # Alanine
    "R": 265.0,  # Arginine
    "N": 187.0,  # Asparagine
    "D": 187.0,  # Aspartate
    "C": 148.0,  # Cysteine
    "E": 214.0,  # Glutamate
    "Q": 214.0,  # Glutamine
    "G": 97.0,   # Glycine
    "H": 216.0,  # Histidine
    "I": 195.0,  # Isoleucine
    "L": 191.0,  # Leucine
    "K": 230.0,  # Lysine
    "M": 203.0,  # Methionine
    "F": 228.0,  # Phenylalanine
    "P": 154.0,  # Proline
    "S": 143.0,  # Serine
    "T": 163.0,  # Threonine
    "W": 264.0,  # Tryptophan
    "Y": 255.0,  # Tyrosine
    "V": 165.0,  # Valine
}


def extract_tmh_values(df, feature, label):
    if not df.index.is_unique:
        df = df.reset_index(drop=True)

    tmh_rows = df.query("type == 'transmembrane region' or type == 'intramembrane region'")


    # Save the DataFrame to a file based on the label
    if label == 'list1':
        file_path = os.path.join('list_report', 'list_report1.tsv')
    elif label == 'list2':
        file_path = os.path.join('list_report', 'list_report2.tsv')
    else:
        file_path = os.path.join('list_report', 'list_report.tsv')
    tmh_rows.to_csv(file_path, sep='\t', index=False)

    if tmh_rows.empty:
        if feature == 'option6':
            return pd.DataFrame(columns=['id', 'sequence'])
        return pd.Series(dtype=float)

    # Extract the specified feature from the TMH rows
    if feature == 'option1':
        return tmh_rows['volume'].dropna()
    elif feature == 'option2':
        return tmh_rows['SA'].dropna()
    elif feature == 'option3':
        return tmh_rows['bulkiness'].dropna()
    elif feature == 'option4':
        return tmh_rows['hydrophobicity'].dropna()
    elif feature == 'option6':
        return tmh_rows[['id', 'sequence']].dropna()
    else:
        # Apply a custom scale based on 'sequence'
        tmh_rows['custom_scale'] = tmh_rows['sequence'].apply(
            lambda x: sliding_window_scores(
                x, 
                read_aa_values_from_csv(os.path.join(os.getcwd(), 'UCS', 'amino_acid_values.csv')), 
                window_size=1, 
                normalize_values=False
            )
        )
        tmh_rows['custom_scale'] = tmh_rows['custom_scale'].apply(lambda x: x[0] if x else None)
        return tmh_rows['custom_scale'].dropna()

def aac_density_plot(selected, backrground):

    
   
    df_aac = get_aac(selected[['id', 'sequence']])
    df_aac_2 = get_aac(backrground[['id', 'sequence']])

    raw_selected = df_aac.copy()
    raw_selected.insert(1, 'group', 'List of interest')
    raw_background = df_aac_2.copy()
    raw_background.insert(1, 'group', 'Background Set')
    raw_combined = pd.concat([raw_selected, raw_background], ignore_index=True)
    raw_file_path = os.path.join('list_report', 'plot_raw_data.tsv')
    raw_combined.to_csv(raw_file_path, sep='	', index=False)


    file_path = os.path.join('AAC_values', 'AAC_values_helixharbor.tsv')
    df_aac.to_csv(file_path, sep='\t', index=False)

    
    df_aac.drop('id', axis=1, inplace=True)
    df_aac_2.drop('id', axis=1, inplace=True)
    
    
    clip = (-10, None)
    bin_size = 0.1
    max_bins = 9

    # Creating a 5x4 grid of plots
    fig, axes = plt.subplots(5, 4, figsize=(20, 15))
    fig.tight_layout(pad=5.0)

    # Loop through each amino acid
    for i, amino_acid in enumerate(df_aac.columns):
        ax = axes[i // 4, i % 4]

        # Extracting data for the current amino acid
        data1 = df_aac[amino_acid]
        data2 = df_aac_2[amino_acid]

         
        # Compute the number of bins
        n_bins = min(max_bins, int((data1.max() - data1.min()) / bin_size))
        n_bins_2 = min(max_bins, int((data2.max() - data2.min()) / bin_size))

        # Compute the bandwidth and adjustment factor
        bw = bin_size / 2
        bw_adjust = bw / (data1.max() - data1.min()) * n_bins
        bw_adjust_2 = bw / (data2.max() - data2.min()) * n_bins_2
        
        

        # Create the density plot for the current amino acid
        sns.kdeplot(data1, clip=clip, label='List of interest', ax=ax)
        sns.kdeplot(data2, clip=clip, label='Background Set', ax=ax)
        ax.set_title(f'Density Plot for {amino_acid}')

        # Perform statistical tests
        ks_statistic, ks_p_value = stats.ks_2samp(data1, data2)
        t_statistic, t_p_value = stats.ttest_ind(data1, data2, equal_var=False)
        wilcoxon_statistic, wilcoxon_p_value = stats.mannwhitneyu(
            data1, data2, alternative='two-sided')

        # Annotate the plot with test results
        ax.text(0.95, 0.95, f'KS: p={ks_p_value:.2e}', transform=ax.transAxes, fontsize=9, 
                verticalalignment='top', horizontalalignment='right')
        ax.text(0.95, 0.90, f'T-Test: p={t_p_value:.2e}', transform=ax.transAxes, fontsize=9, 
                verticalalignment='top', horizontalalignment='right')
        ax.text(0.95, 0.85, f'Wilcoxon: p={wilcoxon_p_value:.2e}', transform=ax.transAxes, fontsize=9, 
                verticalalignment='top', horizontalalignment='right')

       
      
    
    plt.legend()
    plt.show()



    img_bytes = io.BytesIO()
    plt.savefig(img_bytes, format='png')
    img_bytes.seek(0)

    # encode the byte stream as a base64 string
    img_base64 = base64.b64encode(img_bytes.getvalue()).decode('ascii')

   

    return img_base64
