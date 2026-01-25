
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
from run_subunits import run_subunits

from combine_jsons import combine_jsons
from process_templates import process_target_templates


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

    return(" ".join(docker_command))


        

if __name__ == '__main__':
    """
    --input_fasta : *.fasta file location 
    --msa_path : Comma separated stoichiometries. Example --stoichiometries A2,A3,A4
    --num_models : Number of models to generate for each stoichiometry (in the multiple of 5). Example --num_model 25  
    --output_path : Directory where output is to be saved. Example --output_path --/home/user/examples/
    """
    parser = argparse.ArgumentParser()
    parser.add_argument('--input_fasta',type=is_file,required=True)
    parser.add_argument('--target_workdir', required=True)
    parser.add_argument('--num_models', required=True)
    parser.add_argument('--output_path', type=is_dir, required=True)
    # parser.add_argument('--chain', type=str, required=False,default=None)
    args = parser.parse_args()

    target_name,sequences = read_fasta(args.input_fasta)
    is_homomultimer = False
    if len(set(sequences))==1:
        is_homomultimer = True


    output_path_target_complex = os.path.join(args.output_path,os.path.basename((target_name)).upper())
    # output_path_target_complex = os.path.join(args.target_workdir,"N6_multimer_structure_generation","alphafold3_workdir")
    makedir_if_not_exists(output_path_target_complex)

    output_path_predictions = os.path.join(output_path_target_complex,"outputs")
    makedir_if_not_exists(output_path_predictions)

    input_files_path_predictions = os.path.join(output_path_target_complex,"input_files")
    makedir_if_not_exists(input_files_path_predictions)

    complex_subunits_data_dir = os.path.join(input_files_path_predictions,"subunits")
    makedir_if_not_exists(complex_subunits_data_dir)

    complex_templates_data_dir = os.path.join(input_files_path_predictions,"templates")
    makedir_if_not_exists(complex_templates_data_dir)

    # First prepare all the templates. 

    process_target_templates(args.target_workdir,complex_templates_data_dir)

    # Then prepare msas for subunits.

    for fl in os.listdir(args.target_workdir): 
        if fl.split(".")[-1]=="fasta":
            subunit_fasta_path = os.path.join(args.target_workdir,fl)
            chain = fl.split(".")[0]
            run_subunits(subunit_fasta_path,args.target_workdir,complex_subunits_data_dir,chain)
    
    # Then combine all the information and form a json for complex.
    combine_jsons(input_files_path_predictions,args.num_models,is_homomultimer)
    output_path_inference_scripts = os.path.join(output_path_target_complex,"inference_scripts")
    makedir_if_not_exists(output_path_inference_scripts)
    input_jsons_path = input_files_path_predictions

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




