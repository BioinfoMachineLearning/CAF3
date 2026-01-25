import os
import re
import json
from itertools import product
from common.config import HETEROMULTIMER_CONFIG
import os
import argparse


common_template_source = HETEROMULTIMER_CONFIG.common_config.template_source


predictor_template_sources = {}
for predictor_name, config in HETEROMULTIMER_CONFIG.predictors.items():
    template_source = config.get('template_source', common_template_source)
    predictor_template_sources[predictor_name] = template_source

def read_fasta(file_path):
    target_name = os.path.splitext(os.path.basename(file_path))[0].lower()
    sequences = []
    with open(file_path, "r") as fasta_file:
        for line in fasta_file:
            if line.startswith(">"):
                continue
            sequence = line.strip()
            if sequence not in sequences:
                sequences.append(sequence)
            else:
                print("Contains duplicates")
    return target_name, sequences

def generate_json(fasta_path,predictor, num_models, output_path,unpaired_msa_path,paired_msa_path,chain):
    target_name,sequences = read_fasta(fasta_path)
    sequence = sequences[0]
    num_seeds = int(int(num_models)/5)
    used_alphabets = 0 
    json_skeleton = {
        "name": f"{predictor}",
        "sequences": [],
        "modelSeeds": [i for i in range(num_seeds)],
        "dialect": "alphafold3",
        "version": 1
    }

    predictor_template_source = predictor_template_sources.get(predictor, None)
    if predictor_template_source in ["notemplate"]:
        json_skeleton["sequences"].append(
        {
            "protein": {
                "id": chain,
                "sequence": sequences[0],
                "unpairedMsaPath":unpaired_msa_path,
                "pairedMsaPath":paired_msa_path,
                "templates":"notemplate"
            }
        }
    )

    else:
        json_skeleton["sequences"].append(
            {
                "protein": {
                    "id": chain,
                    "sequence": sequences[0],
                    "unpairedMsaPath":unpaired_msa_path,
                    "pairedMsaPath":paired_msa_path,
                }
            }
        )
    json_path = f'{output_path}/{predictor}.json'
    with open(json_path, 'w') as f:
        json.dump(json_skeleton, f,indent=4)

    return json_skeleton

def generate_default_json(fasta_path,predictor, num_models, output_path,chain):
    target_name,sequences = read_fasta(fasta_path)
    sequence = sequences[0]
    num_seeds = int(int(num_models)/5)
    used_alphabets = 0 
    json_skeleton = {
        "name": f"{predictor}",
        "sequences": [],
        "modelSeeds": [i for i in range(num_seeds)],
        "dialect": "alphafold3",
        "version": 1
    }


    json_skeleton["sequences"].append(
    {
        "protein": {
            "id": chain,
            "sequence": sequences[0]
        }
    })
    

    json_path = f'{output_path}/{predictor}.json'
    with open(json_path, 'w') as f:
        json.dump(json_skeleton, f,indent=4)

    return json_skeleton