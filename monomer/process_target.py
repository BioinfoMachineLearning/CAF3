import argparse
import shutil
import os, csv, json
import pandas as pd
from common.config import MONOMER_CONFIG
from utils.generate_jsons import generate_json
from utils.utils import makedir_if_not_exists,is_dir
from utils.convert_to_cif import atom2cif
from utils.extract_cif_chains import filter_cif_by_chain
from msa import handle_msa
from templates import copy_templates,prepare_template_csv_for_af3,return_template_info,handle_pdb_templates,handle_inhouse_templates

common_template_source = MONOMER_CONFIG.common_config.template_source

# pdb_mmcif_database_dir = "/bmlfast/bml_casp16/tools/alphafold_databases_multicom3/pdb_mmcif//mmcif_files/"

predictor_template_sources = {}
for predictor_name, config in MONOMER_CONFIG.predictors.items():
    template_source = config.get('template_source', common_template_source)
    predictor_template_sources[predictor_name] = template_source


def process_each_predictor(target_workdir, predictor_path, predictor, output_path_target,chain):
    msa_info = handle_msa(predictor_path, predictor, output_path_target)
    if msa_info is None:
        # print(f"No msa found for {predictor}")
        return 0

    destination_template_dir = os.path.join(output_path_target, "templates")
    makedir_if_not_exists(destination_template_dir)
    destination_template_path = os.path.join(destination_template_dir, predictor + "_templates.csv")

    predictor_template_source = predictor_template_sources.get(predictor, None)
    if predictor_template_source in ["pdb70", "pdb70_newest"]:
        templates = handle_pdb_templates(predictor_path, predictor, output_path_target, destination_template_path)
    elif predictor_template_source == "pdb_sort90":
        templates = handle_inhouse_templates(target_workdir,chain, predictor_path, predictor, output_path_target, destination_template_path)
    else:
        templates = []

    return {
        "unpairedMsaPath": f"msas/{os.path.basename(msa_info)}",
        "templates": templates
    }

# import os

def process_target(input_fasta, target_workdir, output_path_target, num_models,chain=None):
    af3_inputs_path = os.path.join(output_path_target, "input_files")
    
    makedir_if_not_exists(af3_inputs_path)

    if chain:
        predictor_generation_dir = os.path.join(target_workdir, "N3_monomer_structure_generation",chain.upper())
    else:
        predictor_generation_dir = os.path.join(target_workdir, "N3_monomer_structure_generation")
    if not os.path.isdir(predictor_generation_dir):
        print(f"[Warning] Predictor directory missing: {predictor_generation_dir}")
        return None

    for predictor in os.listdir(predictor_generation_dir):
        if "afsample" in predictor or "drop" in predictor or "ptm" in predictor or "def_esm_msa_ckpt5" in predictor or "dom" in predictor or "ori" in predictor or "pdb70" in predictor:
            continue
        predictor_path = os.path.join(predictor_generation_dir, predictor)
        ranking_debug_path = os.path.join(predictor_path, "ranking_debug.json")

        if not os.path.isfile(ranking_debug_path):
            continue  # Skip if ranking_debug.json doesn't exist

        result = process_each_predictor(target_workdir, predictor_path, predictor, af3_inputs_path,chain)
        if not result:
            continue  # Skip if processing fails

        # print(f"[Processed] {predictor} -> {result}")
        generate_json(
            fasta_path=input_fasta,
            predictor=predictor,
            num_models=num_models,
            output_path=af3_inputs_path,
            msa_path=result["unpairedMsaPath"],
            templates=result["templates"],
            chain=chain
        )

    return af3_inputs_path




