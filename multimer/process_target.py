import os, csv, json
import argparse
from templates import copy_templates,prepare_template_csv_for_af3,return_template_info,handle_pdb_templates,handle_inhouse_templates
from utils.utils import makedir_if_not_exists,is_dir
import shutil
from utils.generate_jsons import generate_json,generate_default_json
import pandas as pd
from utils.convert_to_cif import atom2cif
from utils.extract_cif_chains import filter_cif_by_chain
from msa import handle_msa




def process_each_predictor(target_root_path, predictor_path, predictor, output_path_target,chain):
    monomer_msa_info,paired_msa_info = handle_msa(predictor_path, predictor, chain,output_path_target)


    if monomer_msa_info and paired_msa_info:

        return {
            "unpairedMsaPath": f"msas/{os.path.basename(monomer_msa_info)}",
            "pairedMsaPath": f"msas/{os.path.basename(paired_msa_info)}"
        }
    elif monomer_msa_info and not paired_msa_info:
        return {
            "unpairedMsaPath": f"msas/{os.path.basename(monomer_msa_info)}",
            "pairedMsaPath": ""
        }


def process_target(input_fasta, target_root_path, output_path_target, num_models,chain=None):

    af3_inputs_path = os.path.join(output_path_target, "input_files")
    
    makedir_if_not_exists(af3_inputs_path)
    generate_default_json(
            fasta_path=input_fasta,
            predictor="default_af3",
            num_models=num_models,
            output_path=af3_inputs_path,
            chain=chain
        )

    predictor_generation_dir = os.path.join(target_root_path, "N6_multimer_structure_generation")
    if not os.path.isdir(predictor_generation_dir):
        print(f"[Warning] Predictor directory missing: {predictor_generation_dir}")
        return None


    for predictor in os.listdir(predictor_generation_dir):
        
        if "afsample" in predictor or "drop" in predictor or "ptm" in predictor:
            continue
        predictor_path = os.path.join(predictor_generation_dir, predictor)
        ranking_debug_path = os.path.join(predictor_path, "ranking_debug.json")

        if not os.path.isfile(ranking_debug_path):
            continue  # Skip if ranking_debug.json doesn't exist

        result = process_each_predictor(target_root_path, predictor_path, predictor, af3_inputs_path,chain)
        if not result:
            continue  # Skip if processing fails

        generate_json(
            fasta_path=input_fasta,
            predictor=predictor,
            num_models=num_models,
            output_path=af3_inputs_path,
            unpaired_msa_path=result["unpairedMsaPath"],
            paired_msa_path=result["pairedMsaPath"],
            # templates=result["templates"],
            chain=chain
        )

    return af3_inputs_path




