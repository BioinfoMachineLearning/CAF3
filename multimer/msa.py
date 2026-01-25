import os
import shutil
from utils.utils import makedir_if_not_exists




def handle_msa_old(predictor_path, predictor,chain, destination_path):
    monomer_msa_path = os.path.join(predictor_path, "msas",chain, "monomer_final.a3m")
    multimer_msa_path = os.path.join(predictor_path, "msas", f"{chain}.paired.a3m")
    destination_monomer_msa_path,destination_multimer_msa_path = None,None
    destination_msa_dir = os.path.join(destination_path, "msas")
    makedir_if_not_exists(destination_msa_dir)
    
    if os.path.exists(monomer_msa_path):
        destination_monomer_msa_path = os.path.join(destination_msa_dir, predictor + "_monomer_alignments.a3m")
        shutil.copy(monomer_msa_path, destination_monomer_msa_path)
    if os.path.exists(multimer_msa_path):
        destination_multimer_msa_path = os.path.join(destination_msa_dir, predictor + "_multimer_alignments.a3m")
        shutil.copy(multimer_msa_path, destination_multimer_msa_path)
    
    return destination_monomer_msa_path,destination_multimer_msa_path

def handle_msa(predictor_path, predictor,chain, destination_path):
    monomer_msa_path = os.path.join(predictor_path, "msas",chain, "monomer_final.a3m")
    multimer_msa_path = os.path.join(predictor_path, "msas",chain, "multimer_final.a3m")
    destination_monomer_msa_path,destination_multimer_msa_path = None,None
    destination_msa_dir = os.path.join(destination_path, "msas")
    makedir_if_not_exists(destination_msa_dir)
    
    if os.path.exists(monomer_msa_path):
        destination_monomer_msa_path = os.path.join(destination_msa_dir, predictor + "_monomer_alignments.a3m")
        shutil.copy(monomer_msa_path, destination_monomer_msa_path)
    if os.path.exists(multimer_msa_path):
        destination_multimer_msa_path = os.path.join(destination_msa_dir, predictor + "_multimer_alignments.a3m")
        shutil.copy(multimer_msa_path, destination_multimer_msa_path)
    
    return destination_monomer_msa_path,destination_multimer_msa_path
