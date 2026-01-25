import os
import shutil
from utils.convert_to_cif import atom2cif
from utils import parse_hhr
from utils import parse_seq_temp
import csv
import ast
from utils.utils import makedir_if_not_exists
import pandas as pd
from common.config import pdb_mmcif_database_dir
from utils.extract_cif_chains import filter_cif_by_chain

def _copy(source_dir_file,destination_dir):
    shutil.copy(source_dir_file,destination_dir)

def copy_templates(source_dir,destination_dir):
    for f in os.listdir(source_dir):
        if f.split(".")[-1] == "atom":
            input_path = os.path.join(source_dir,f)
            _copy(input_path,destination_dir)
            destination_atom = os.path.join(destination_dir,f)
            destination_cif = os.path.join(destination_dir,f.split(".")[0]+".cif")
            atom2cif(destination_atom,destination_cif)

def handle_pdb_templates(predictor_path, predictor, destination_path, dest_template_csv):
    template_path = os.path.join(predictor_path, "msas", "pdb_hits.hhr")
    output_csv_path = prepare_template_csv_for_af3(template_path, dest_template_csv)

    template_cif_destination_path = os.path.join(destination_path, "pdb_templates")
    makedir_if_not_exists(template_cif_destination_path)

    if output_csv_path:
        df = pd.read_csv(output_csv_path)
        for template_name in df["template_name"]:
            input_cif = os.path.join(pdb_mmcif_database_dir, template_name[:4].lower() + ".cif")
            temp_cif_path = os.path.join(template_cif_destination_path, template_name[:4] + ".cif")
            final_cif_path = os.path.join(template_cif_destination_path, template_name + ".cif")

            if  os.path.exists(final_cif_path):
                continue

            if os.path.exists(input_cif):
                shutil.copy(input_cif, temp_cif_path)
                filter_cif_by_chain(temp_cif_path, template_name[4:].upper(), final_cif_path)

    return return_template_info(dest_template_csv, "pdb_templates", template_cif_destination_path)




def handle_inhouse_templates(target_path,chain, predictor_path, predictor, destination_path, dest_template_csv):
    if chain:
        template_path = os.path.join(target_path, "N2_monomer_template_search",chain, "sequence_templates.csv")
    else:
        template_path = os.path.join(target_path, "N2_monomer_template_search", "sequence_templates.csv")
    if not os.path.exists(template_path):
        return []
    output_csv_path = prepare_template_csv_for_af3(template_path, dest_template_csv)

    template_cif_destination_path = os.path.join(destination_path, "inhouse_templates")
    makedir_if_not_exists(template_cif_destination_path)
    inhouse_template_dataset_dir = os.path.join(predictor_path, "templates")

    if output_csv_path:
        df = pd.read_csv(output_csv_path)
        for template_name in df["template_name"]:
            atom_input = os.path.join(inhouse_template_dataset_dir, template_name + ".atom")
            atom_output = os.path.join(template_cif_destination_path, template_name + ".atom")

            cif_output = atom_output[:-4] + "cif"
            if os.path.exists(cif_output):
                continue
            if os.path.exists(atom_input):
                shutil.copy(atom_input, atom_output)
                atom2cif(atom_output, cif_output)

    return return_template_info(dest_template_csv, "inhouse_templates", template_cif_destination_path)



def prepare_template_csv_for_af3(input_file,output_csv):
    if input_file.split(".")[-1] == "hhr":
        parse_hhr.run_hhr_parser(hhr_file=input_file,output=output_csv)
        return output_csv

    elif input_file.split(".")[-1] == "csv":
        parse_seq_temp.run_seq_temp_parser(template_file=input_file,output=output_csv)
        return output_csv
    
    else:
        print("Invalid template format\nRunning Template free")
        return 0


def return_template_info(af3_template_csv, template_base_dir, template_cif_path):
    """
    Returns cleaned template info for AlphaFold3. Removes missing template-residue mappings.

    Args:
        af3_template_csv (str): Path to CSV with template_name, query_indices, template_indices
        template_cif_path (str): Directory with .cif files
        template_base_dir (str): Path prefix to be added to mmcifPath

    Returns:
        List[Dict]: List of cleaned template info dicts
    """
    template_info_list = []

    with open(af3_template_csv, 'r') as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            template_name = row['template_name']
            try:
                query_indices = ast.literal_eval(row['query_indices'])
                template_indices = ast.literal_eval(row['template_indices'])
            except Exception as e:
                print(f"Warning: Skipping {template_name} due to parse error: {e}")
                continue

            cif_file = os.path.join(template_cif_path, f"{template_name}.cif")
            if not os.path.exists(cif_file):
                # print(f"Warning: Skipping {template_name}, CIF not found at {cif_file}")
                continue

            try:
                residue_indices = set()
                atom_lines_found = False

                with open(cif_file, 'r') as cif:
                    for line in cif:
                        if line.startswith("ATOM"):
                            atom_lines_found = True
                            tokens = line.strip().split()
                            try:
                                residue_index = int(tokens[8])  # column with residue number
                                residue_indices.add(residue_index)
                            except (IndexError, ValueError):
                                continue

                if not atom_lines_found:
                    continue

                # Filter out missing template_indices (adjusted +1)
                adjusted_template_indices = [idx + 1 for idx in template_indices]
                cleaned_query_indices = []
                cleaned_template_indices = []

                for qi, ti, ti_adj in zip(query_indices, template_indices, adjusted_template_indices):
                    if ti_adj in residue_indices:
                        cleaned_query_indices.append(qi)
                        cleaned_template_indices.append(ti)

                if not cleaned_template_indices:
                    continue

                template_info = {
                    "mmcifPath": f"{template_base_dir}/{template_name}.cif",
                    "queryIndices": cleaned_query_indices,
                    "templateIndices": cleaned_template_indices
                }
                template_info_list.append(template_info)

            except Exception as e:
                print(f"Error: Skipping {template_name}, failed reading CIF: {e}")

    return template_info_list[:4]