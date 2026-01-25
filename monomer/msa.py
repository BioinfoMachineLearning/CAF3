import os
import shutil
from utils.utils import makedir_if_not_exists


def handle_msa(predictor_path, predictor, destination_path):
    msa_path = os.path.join(predictor_path, "msas", "monomer_final.a3m")
    if not os.path.exists(msa_path):
        return None

    destination_msa_dir = os.path.join(destination_path, "msas")
    makedir_if_not_exists(destination_msa_dir)
    destination_msa_path = os.path.join(destination_msa_dir, predictor + "_alignments.a3m")
    shutil.copy(msa_path, destination_msa_path)
    return destination_msa_path
