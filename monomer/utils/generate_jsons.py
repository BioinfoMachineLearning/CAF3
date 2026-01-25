import os
import re
import json
from itertools import product
import random

import os
import argparse

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

def generate_json(fasta_path,predictor, num_models, output_path,msa_path,templates,chain):
    target_name,sequences = read_fasta(fasta_path)
    sequence = sequences[0]
    # num_seeds = int(int(num_models)/5)
    used_alphabets = 0 
    json_skeleton = {
        "name": f"{predictor}",
        "sequences": [],
        "modelSeeds": random.sample(range(1, 100000), int(int(num_models)/5)),
        "dialect": "alphafold3",
        "version": 1
    }
    if not chain:
        chain="A"
    json_skeleton["sequences"].append(
        {
            "protein": {
                "id": chain,
                "sequence": sequences[0],
                "unpairedMsaPath":msa_path,
                "pairedMsa":"",
                "templates":templates
            }
        }
    )
    json_path = f'{output_path}/{predictor}.json'
    with open(json_path, 'w') as f:
        json.dump(json_skeleton, f,indent=4)
    # print("Created : ",json_path)

    return json_skeleton