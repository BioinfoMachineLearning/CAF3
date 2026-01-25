
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
from common.config import af3_program_path,af3_params_path,af3_db_path




def run_docker(each_json, input_jsons_path,target_output_path):
    original_cwd = os.getcwd()

    docker_command = [
    "docker", "run", "--rm",
    "-u", f"{os.getuid()}:{os.getgid()}",
    "--volume", f"{input_jsons_path}:/af_input",
    "--volume", f"{target_output_path}:/af_output",
    "--volume", f"{af3_params_path}:/models",
    "--volume", f"{af3_db_path}:/public_databases",
    "--gpus", "all",
    "alphafold3",
    "python", "run_alphafold.py",
    f"--json_path=/af_input/{each_json}",
    "--model_dir=/models",
    "--output_dir=/af_output"
]
    # print("Running Docker Command: ", " ".join(docker_command))

    return(" ".join(docker_command))
    # print("Running Docker Command: ", " ".join(docker_command))


        

if __name__ == '__main__':
    """

    """
    parser = argparse.ArgumentParser()
    parser.add_argument('--input_fasta',type=is_file,required=True)
    parser.add_argument('--target_workdir', required=True)
    parser.add_argument('--num_models', required=True)
    parser.add_argument('--output_path', type=is_dir, required=True)
    parser.add_argument('--chain', type=str, required=False,default=None)
    args = parser.parse_args()

    # so. input_fasta will get a single fasta. and output_path will get subunits/ dir. and target root_path will get like H1236/
    target_name,sequences = read_fasta(args.input_fasta)
    #for a subunit
    # chain = target_name.upper()
    chain = args.chain
    if chain:
        chain = chain.upper()
    output_path_target = os.path.join(args.output_path,os.path.basename((target_name)).upper())
    # output_path_target = os.path.join(args.target_workdir,"N3_monomer_structure_generation","alphafold3_workdir")
    # makedir_if_not_exists(output_path_target)

    output_path_predictions = os.path.join(output_path_target,"outputs")
    makedir_if_not_exists(output_path_predictions)

    output_path_inference_scripts = os.path.join(output_path_target,"inference_scripts")
    makedir_if_not_exists(output_path_inference_scripts)

    input_jsons_path = process_target(args.input_fasta,args.target_workdir,output_path_target,args.num_models,chain)
    all_docker_commands = []
    for each_json in os.listdir(input_jsons_path):
        json_path = os.path.join(input_jsons_path, each_json)
        if not each_json.split(".")[-1]=="json":
            continue

        docker_cmd = run_docker(each_json, input_jsons_path, output_path_predictions)
        all_docker_commands.append(docker_cmd)

        # Write individual bash script
        script_name = os.path.splitext(each_json)[0] + ".sh"
        script_path = os.path.join(output_path_inference_scripts, script_name)
        with open(script_path, 'w') as f:
            f.write("#!/bin/bash\n")
            f.write(docker_cmd + "\n")
        os.chmod(script_path, 0o755)  # make it executable

    print("One can run this command directly.\nOr run the sh scripts in inference_scripts/ in output path")
    print(" && ".join(all_docker_commands))


