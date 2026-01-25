
import os
import re
import pandas as pd
import string
import shutil
from common.config import pdb_mmcif_database_dir
from templates import return_template_info
import os
from pathlib import Path
from utils.convert_to_cif import atom2cif
from utils import parse_hhr
from utils import parse_seq_temp
from utils import parse_seq_temp_tmsearch
from utils import parse_seq_temp_foldseek
import csv
import ast
from utils.utils import makedir_if_not_exists
import pandas as pd
from utils.extract_cif_chains import filter_cif_by_chain
from utils.generate_jsons import read_fasta
from template_search.pipeline import _process_single_chain
from utils.repair_atom import repair_atom


def transform_csv_all(input_csv: str, output_csv: str):
    # Read input CSV
    df = pd.read_csv(input_csv)

    # Extract template (0th element when splitting "name" on space)
    df["target"] = df["template"]

    # aligned_length
    df["alnlen"] = df["aligned_length"]

    # aln_temp
    df["taln"] = df["aln_temp"]

    # aln_query
    df["qaln"] = df["aln_query"]
    df["tstart"] = df["tstart"]
    df["tend"] = df["tend"]
    df["qstart"] = df["qstart"]
    df["qend"] = df["qend"]
    df["sum_probs"] = df["sum_probs"] if "sum_probs" in df else "not_applicable"
    df["tmscore"]  = df["tmscore"]  if "tmscore"  in df else "not_applicable"

    # Select and reorder the required columns
    out_df = df[[
        "target",
        "alnlen",
        "qaln",
        "qstart",
        "qend",
        "taln",
        "tstart",
        "tend",
        "sum_probs",
        "tmscore"
    ]]

    # Save to new CSV
    out_df.to_csv(output_csv, index=False)


def transform_csv_pdb70_seq(input_csv: str, output_csv: str):
    # Read input CSV
    df = pd.read_csv(input_csv)

    # Extract template (0th element when splitting "name" on space)
    df["template"] = df["name"].str.split(" ").str[0]

    # aligned_length
    df["aligned_length"] = df["aligned_cols"]

    # aln_temp
    df["aln_temp"] = df["hit_sequence"]

    # aln_query
    df["aln_query"] = df["query"]

    # Find tstart and tend from indices_hit
    # shoukd I add 1 to start and end based 
    # on /bmlfast/casp16_ts_qs_qa/TS_run/valid/H1213_human/N5_monomer_templates_concatenation/pdb70_seq/sequence_templates.csv?
    def find_first_and_last_valid(values):
        vals = [int(v) for v in values.split("_") if v != "-1"]
        if not vals:
            return None, None
        return vals[0], vals[-1]
    def find_first_and_last_valid(values):
        vals = [int(v) for v in values.split("_") if v != "-1"]
        if not vals:
            return None, None
        return vals[0]+1, vals[-1]+1

    df[["tstart", "tend"]] = df["indices_hit"].apply(lambda x: pd.Series(find_first_and_last_valid(x)))

    # Find qstart and qend from indices_query
    df[["qstart", "qend"]] = df["indices_query"].apply(lambda x: pd.Series(find_first_and_last_valid(x)))

    # sum_probs stays the same
    df["sum_probs"] = df["sum_probs"]

    # Select and reorder the required columns
    out_df = df[[
        "template",
        "aligned_length",
        "aln_temp",
        "tstart",
        "tend",
        "aln_query",
        "qstart",
        "qend",
        "sum_probs"
    ]]

    # Save to new CSV
    out_df.to_csv(output_csv, index=False)


def split_chains_csv(input_csv, output_dir,template_type_dir,template_type,is_pdb70_seq=False):
    os.makedirs(output_dir, exist_ok=True)

    # Load CSV, drop dummy first index column if present
    df = pd.read_csv(input_csv)
    df = df.loc[:, ~df.columns.str.match("Unnamed")]  # remove unnamed cols like "Unnamed: 0"

    # Find all numbered columns
    chain_groups = {}
    for col in df.columns:
        match = re.match(r"(.+?)(\d+)$", col)  # match base name + number
        if match:
            base, num = match.groups()
            chain_groups.setdefault(num, []).append((col, base))
        # else: skip non-numbered columns like "index", "tpdbcode"
    # Create per-chain CSVs
    templates_atom_dir = os.path.join(output_dir,"templates")
    os.makedirs(templates_atom_dir, exist_ok=True)
    for idx, (chain_num, col_pairs) in enumerate(sorted(chain_groups.items(), key=lambda x: int(x[0]))):
        selected_cols = [c[0] for c in col_pairs]
        base_names = [c[1] for c in col_pairs]

        sub_df = df[selected_cols].copy()
        sub_df.columns = base_names  # strip numbers

        # Use letters A, B, C...
        chain_letter = string.ascii_uppercase[idx]
        # output_dir
        output_path_temp = os.path.join(output_dir, f"{chain_letter}_temp.csv")
        sub_df.to_csv(output_path_temp, index=False)
        output_path_final = os.path.join(output_dir, f"{chain_letter}.csv")
        dest_template_csv = os.path.join(output_dir, f"{chain_letter}_templates.csv")
        if template_type == "pdb70_seq":
            output_path_temp_pdb70_seq = os.path.join(output_dir, f"{chain_letter}_temp_pdb70_seq.csv")
            transform_csv_pdb70_seq(output_path_temp,output_path_temp_pdb70_seq)
            output_path_temp = output_path_temp_pdb70_seq
            transform_csv_all(output_path_temp,output_path_final)
            handle_pdb70_templates(output_path_final,"",output_dir,dest_template_csv)
            # shutil.copytree
        else:
            for template_atom in os.listdir(os.path.join(template_type_dir,"templates")):
                if template_atom.endswith(".atom") and not os.path.exists(os.path.join(templates_atom_dir,template_atom)):
                    shutil.copy(os.path.join(template_type_dir,"templates",template_atom),os.path.join(templates_atom_dir,template_atom))
         
            transform_csv_all(output_path_temp,output_path_final)
        
            # dest_template_csv = os.path.join(output_dir, f"{chain_letter}_templates.csv")
            if template_type=="tmsearch":
                handle_inhouse_templates(output_path_final, templates_atom_dir, templates_atom_dir, dest_template_csv,is_tmsearch=True)
            else:
                handle_inhouse_templates(output_path_final, templates_atom_dir, templates_atom_dir, dest_template_csv)

        # return(output_path_final)

def prepare_template_csv_for_af3(input_file, output_csv, is_tmsearch=False):
    ext = input_file.split(".")[-1].lower()

    if ext == "hhr":
        parse_hhr.run_hhr_parser(hhr_file=input_file, output=output_csv)
        return output_csv

    elif ext == "csv":
        if is_tmsearch:
            parse_seq_temp_tmsearch.run_seq_temp_parser(template_file=input_file, output=output_csv)
        else:
            parse_seq_temp.run_seq_temp_parser(template_file=input_file, output=output_csv)
        return output_csv

    else:
        print("Invalid template format\nRunning Template free")
        return 0


def handle_inhouse_templates(input_csv, template_atom_path, template_cif_path, dest_template_csv,is_tmsearch=False):
    template_path = input_csv
    if not os.path.exists(template_path):
        return []
    if is_tmsearch:
        # print
        output_csv_path = prepare_template_csv_for_af3(template_path, dest_template_csv,is_tmsearch=True)
    else:
        output_csv_path = prepare_template_csv_for_af3(template_path, dest_template_csv)

    template_cif_destination_path = template_cif_path
    makedir_if_not_exists(template_cif_destination_path)
    inhouse_template_dataset_dir = template_atom_path

    if output_csv_path:
        df = pd.read_csv(output_csv_path)
        for template_name in df["template_name"]:
            atom_input = os.path.join(inhouse_template_dataset_dir, template_name + ".atom")
            atom_output = os.path.join(template_cif_destination_path, template_name + ".atom")

            cif_output = atom_output[:-4] + "cif"
            if os.path.exists(cif_output):
                continue
            if os.path.exists(atom_input):
                # print("Found:", atom_input)
                # shutil.copy(atom_input, atom_output)
                try:
                    if is_tmsearch:
                        repaired_atom = atom_output.replace(".atom",".atom_repaired")
                        repair_atom(atom_output,repaired_atom)
                        atom_output = repaired_atom

                    atom2cif(atom_output, cif_output)
                except Exception as e:
                    print(e)
                    continue

    return return_template_info(dest_template_csv, "inhouse_templates", template_cif_destination_path)

def handle_foldseek_template(predictor_directory,output_dir):
    template_cif_destination_path = os.path.join(output_dir,"templates")
    inhouse_template_dataset_dir = os.path.join(predictor_directory,"templates")
    makedir_if_not_exists(template_cif_destination_path)

    for f in os.listdir(predictor_directory):
        if f.endswith(".top50"):
            # print(f)
            destination_template_csv = os.path.join(output_dir,f.replace(".top50","_templates.csv"))
            parse_seq_temp_foldseek.run_seq_temp_parser(template_file=os.path.join(predictor_directory,f),output=destination_template_csv)
            df = pd.read_csv(destination_template_csv)
            for template_name in df["template_name"]:
                template_name = template_name.replace(".atom.gz","")
                template_name = template_name.replace(".pdb","")
                # print(template_name)
                atom_input = os.path.join(inhouse_template_dataset_dir, template_name + ".atom")
                # print(atom_input)
                atom_output = os.path.join(template_cif_destination_path, template_name + ".atom")

                cif_output = atom_output[:-4] + "cif"
                if os.path.exists(cif_output):
                    continue
                if os.path.exists(atom_input):
                    # print("yes",atom_input)
                    try:
                        shutil.copy(atom_input,atom_output)
                        atom2cif(atom_output, cif_output)
                    except Exception as e:
                        print(e)
                        continue

# handle_foldseek_template("/bmlfast/casp16_ts_qs_qa/TS_run/valid/T1235_human/N6_multimer_structure_generation/folds_iter_1","/bmlfast/bml_casp17/installation_test/MULTICOM5/test_targets/H1202_redo_template_2/T1235/input_files/templates/foldseek")


def handle_pdb_templates(N2_path, predictor, destination_path, dest_template_csv):
    # template_path = os.path.join(predictor_path, "msas", "pdb_hits.hhr")
    template_path = os.path.join(N2_path, "output.hhr")
    output_csv_path = prepare_template_csv_for_af3(template_path, dest_template_csv)

    template_cif_destination_path = os.path.join(destination_path, "templates")
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
                # print("Found:", input_cif)
                shutil.copy(input_cif, temp_cif_path)
                filter_cif_by_chain(temp_cif_path, template_name[4:], final_cif_path)

    return return_template_info(dest_template_csv, "pdb_templates", template_cif_destination_path)

def handle_pdb_seqres_templates(destination_path, dest_template_csv):
    template_cif_destination_path = os.path.join(destination_path, "templates")
    makedir_if_not_exists(template_cif_destination_path)

    if dest_template_csv:
        df = pd.read_csv(dest_template_csv)
        for template_name in df["template_name"]:
            input_cif = os.path.join(pdb_mmcif_database_dir, template_name[:4].lower() + ".cif")
            temp_cif_path = os.path.join(template_cif_destination_path, template_name[:4] + ".cif")
            final_cif_path = os.path.join(template_cif_destination_path, template_name + ".cif")

            if  os.path.exists(final_cif_path):
                continue

            if os.path.exists(input_cif):
                shutil.copy(input_cif, temp_cif_path)
                filter_cif_by_chain(temp_cif_path, template_name.split("_")[-1], final_cif_path)

    return return_template_info(dest_template_csv, "pdb_templates", template_cif_destination_path)

def handle_pdb70_templates(input_csv, predictor, destination_path, dest_template_csv):
    # template_path = os.path.join(predictor_path, "msas", "pdb_hits.hhr")
    template_path = input_csv
    output_csv_path = prepare_template_csv_for_af3(template_path, dest_template_csv)

    template_cif_destination_path = os.path.join(destination_path, "templates")
    makedir_if_not_exists(template_cif_destination_path)

    if output_csv_path:
        df = pd.read_csv(output_csv_path)
        for template_name in df["template_name"]:
            template_name = template_name.replace("_","")
            input_cif = os.path.join(pdb_mmcif_database_dir, template_name[:4].lower() + ".cif")
            temp_cif_path = os.path.join(template_cif_destination_path, template_name[:4] + ".cif")
            final_cif_path = os.path.join(template_cif_destination_path, template_name + ".cif")

            if  os.path.exists(final_cif_path):
                continue

            if os.path.exists(input_cif):
                shutil.copy(input_cif, temp_cif_path)
                filter_cif_by_chain(temp_cif_path, template_name[4:], final_cif_path)

    return return_template_info(dest_template_csv, "pdb_templates", template_cif_destination_path)

def process_target_templates(root_dir, output_root_dir):
    """
    Process monomer templates (N2) and concatenated templates (N5) for a given target.

    Args:
        root_dir (str): Root directory for the target (e.g., ".../H1245_human")
        output_root_dir (str): Output directory where all results will be saved

    Returns:
        None
    """
    # from your_module import makedir_if_not_exists, handle_pdb_templates, split_chains_csv  # replace with actual imports

    target_name = os.path.basename(root_dir)

    # -------------------------------
    # Step 1: Handle pdb_seqres templates

    n1_path = os.path.join(root_dir, "N1_monomer_alignments_generation")
    if not os.path.isdir(n1_path):
        print(f"[WARNING] N1 path not found: {n1_path}")
    else:
        for chain in sorted(os.listdir(n1_path)):
            chain_path = os.path.join(n1_path, chain)
            if not os.path.isdir(chain_path):
                continue  # skip non-directories

            destination_path = os.path.join(output_root_dir, "pdb_seqres")
            makedir_if_not_exists(destination_path)

            fasta_path = os.path.join(chain_path,f"{chain}.fasta")
            target,sequences = read_fasta(fasta_path)
            uniref_sto_path = os.path.join(chain_path,f"{chain}_uniref90.sto")

            dest_template_csv = os.path.join(destination_path, f"{chain}_templates.csv")
            # _process_single_chain(sequences[0],uniref_sto_path,destination_path,chain,dest_template_csv)
            # _process_single_chain(sequences[0],uniref_sto_path,destination_path,chain,dest_template_csv)
            # handle_pdb_seqres_templates(destination_path,dest_template_csv)
            try:
                # handle_pdb_templates(chain_path, "", destination_path, dest_template_csv)
                _process_single_chain(sequences[0],uniref_sto_path,destination_path,chain,dest_template_csv)
                handle_pdb_seqres_templates(destination_path,dest_template_csv)
                print(f"[INFO] Processed PDB templates for chain {chain}")
            except Exception as e:
                print(f"[ERROR] Failed to handle_pdb_templates for chain {chain}: {e}")


    # -------------------------------
    # Step 2: Handle foldseek templates
    # -------------------------------
    print("Handling foldseek templates")
    N6_directory = os.path.join(root_dir, "N6_multimer_structure_generation")
    for predictor in os.listdir(N6_directory):
        if predictor not in ["folds_iter_1","folds_iter_2","folds_iter_o_1","folds_iter_o_2","folds_iter_nop_1","folds_iter_nop_2"]:
            continue
        predictor_dir = os.path.join(N6_directory,predictor)
        output_dir = os.path.join(output_root_dir, predictor)

        handle_foldseek_template(predictor_dir,output_dir)


    # return
    # -------------------------------
    # Step 3: Handle concatenated (N5) templates
    # -------------------------------
    template_concat_dir = os.path.join(root_dir, "N5_monomer_templates_concatenation")
    if not os.path.isdir(template_concat_dir):
        print(f"[WARNING] N5 path not found: {template_concat_dir}")
        return

    for template_type in os.listdir(template_concat_dir):
        # if template_type != "tmsearch":
        #     continue
        template_type_dir = os.path.join(template_concat_dir, template_type)
        if not os.path.isdir(template_type_dir):
            continue

        # Find the relevant *_templates.csv
        input_csv = None
        for fl in os.listdir(template_type_dir):
            if fl.endswith("_templates.csv"):
                input_csv = os.path.join(template_type_dir, fl)
                break

        if not input_csv or not os.path.isfile(input_csv):
            print(f"[WARNING] No valid *_templates.csv found in {template_type_dir}")
            continue

        output_dir = os.path.join(output_root_dir, template_type)
        makedir_if_not_exists(output_dir)

        try:
            if template_type == "pdb70_seq":
                _ = split_chains_csv(input_csv, output_dir, template_type_dir,template_type, is_pdb70_seq=True)
            else:
                # print("############3")
                # print(input_csv, output_dir, template_type_dir,template_type)
                _ = split_chains_csv(input_csv, output_dir, template_type_dir,template_type)
            print(f"[INFO] Processed chain-split for template_type: {template_type}")
        except Exception as e:
            print(f"[ERROR] Failed to split_chains_csv for {template_type}: {e}")


# process_target_templates(
#     root_dir="/bmlfast/casp16_ts_qs_qa/TS_run/valid/H1202_human",
#     output_root_dir="/bmlfast/bml_casp17/installation_test/MULTICOM5/test_targets/H1215_redo/H1215"
# )
