
import copy
import json
import subprocess
import re
import os
import time
import csv
import argparse
from utils.utils import makedir_if_not_exists,is_dir,is_file
from utils.generate_jsons import read_fasta
from process_target import process_target
from common.config import af3_program_path,af3_params_path,af3_db_path,num_models_subunit




def run_docker(each_json, input_jsons_path,target_output_path):
    original_cwd = os.getcwd()
    docker_command = [
        "docker", "run", "--rm",
        "--volume", f"{input_jsons_path}:/root/af_input",
        "--volume", f"{target_output_path}:/root/af_output",
        "--volume", f"{af3_params_path}:/root/models",
        "--volume", f"{af3_db_path}:/root/public_databases",
        "--gpus", "all",
        "alphafold3",
        "python", "run_alphafold.py",
        f"--json_path=/root/af_input/{each_json}",
        "--model_dir=/root/models",
        "--output_dir=/root/af_output"
    ]

    return(" ".join(docker_command))



    
def run_subunits(input_fasta,target_root_path,output_path,chain):

    target_name,sequences = read_fasta(input_fasta)
    output_path_target = os.path.join(output_path,chain)
    makedir_if_not_exists(output_path_target)

    output_path_predictions = os.path.join(output_path_target,"outputs")
    makedir_if_not_exists(output_path_predictions)

    output_path_inference_scripts = os.path.join(output_path_target,"inference_scripts")
    makedir_if_not_exists(output_path_inference_scripts)

    input_jsons_path = process_target(input_fasta,target_root_path,output_path_target,num_models_subunit,chain)
    all_docker_commands = []
    for each_json in os.listdir(input_jsons_path):
        json_path = os.path.join(input_jsons_path, each_json)
        if not each_json.split(".")[-1]=="json":
            continue

        docker_cmd = run_docker(each_json, input_jsons_path, output_path_predictions)
        all_docker_commands.append(docker_cmd)

        script_name = os.path.splitext(each_json)[0] + ".sh"
        script_path = os.path.join(output_path_inference_scripts, script_name)
        with open(script_path, 'w') as f:
            f.write("#!/bin/bash\n")
            f.write(docker_cmd + "\n")
        os.chmod(script_path, 0o755)  # make it executable

    # print("One can run this command directly.\nOr run the sh scripts in inference_scripts/ in output path")
    # print(" && ".join(all_docker_commands))


