import json
import sys
from pathlib import Path
from collections import defaultdict
import os
from templates import return_template_info
import random
import re
from common.config import MONOMER_CONFIG,HOMOMULTIMER_CONFIG,HETEROMULTIMER_CONFIG


def get_template_source(is_homomultimer):
    predictor_template_sources = {}
    if is_homomultimer:
        common_template_source = HOMOMULTIMER_CONFIG.common_config.template_source
        
        for predictor_name, config in HOMOMULTIMER_CONFIG.predictors.items():
            template_source = config.get('template_source', common_template_source)
            predictor_template_sources[predictor_name.lower()] = template_source.lower()
    else:
        common_template_source = HETEROMULTIMER_CONFIG.common_config.template_source
        for predictor_name, config in HETEROMULTIMER_CONFIG.predictors.items():
            template_source = config.get('template_source', common_template_source)
            predictor_template_sources[predictor_name.lower()] = template_source.lower()
    return predictor_template_sources



def get_template_info_from_source(ip_dir,template_source,representative_id,predictor):
    if template_source == "notemplate":
        return []
    if template_source == "pdb_seqres":
        return return_template_info(os.path.join(ip_dir,"templates","pdb_seqres",f"{representative_id}_templates.csv"),"templates/pdb_seqres/templates",os.path.join(ip_dir,"templates","pdb_seqres","templates"))
    if template_source == "sequence_based_template_pdb70":
        return return_template_info(os.path.join(ip_dir,"templates","pdb70_seq",f"{representative_id}_templates.csv"),"templates/pdb70_seq/templates",os.path.join(ip_dir,"templates","pdb70_seq","templates"))
    if template_source == "foldseek_structure_based_template":
        return return_template_info(os.path.join(ip_dir,"templates","struct_temp",f"{representative_id}_templates.csv"),"templates/struct_temp/templates",os.path.join(ip_dir,"templates","struct_temp","templates"))
    if template_source == "tmsearch_structure_based_template":
        return return_template_info(os.path.join(ip_dir,"templates","tmsearch",f"{representative_id}_templates.csv"),"templates/tmsearch/templates",os.path.join(ip_dir,"templates","tmsearch","templates"))
    if template_source == "sequence_based_template_pdb_sort90":
        return return_template_info(os.path.join(ip_dir,"templates","pdb_seq",f"{representative_id}_templates.csv"),"templates/pdb_seq/templates",os.path.join(ip_dir,"templates","pdb_seq","templates"))
    if template_source == "sequence_based_template_pdb_complex":
        return return_template_info(os.path.join(ip_dir,"templates","complex_pdb_seq",f"{representative_id}_templates.csv"),"templates/complex_pdb_seq/templates",os.path.join(ip_dir,"templates","complex_pdb_seq","templates"))
    if template_source == "foldseek":
        return return_template_info(os.path.join(ip_dir,"templates",predictor,f"{representative_id}_templates.csv"),f"templates/{predictor}/templates",os.path.join(ip_dir,"templates",predictor,"templates"))
    else:
        return []

def extract_predictor_name(predictor: str) -> str:
    """
    Extract base predictor name by removing numeric suffix after last underscore (if present).
    
    Examples:
        "deepmsa2_1"    → "deepmsa2"
        "folds_iter_12" → "folds_iter"
        "some.name_foo" → "some.name_foo"  (unchanged)
    """
    predictor = predictor.split(".")[0]  # Remove file extension if any (e.g., .pkl, .json)

    # Check if the last part after "_" is numeric
    if "_" in predictor:
        parts = predictor.rsplit("_", 1)
        if parts[1].isdigit():
            return parts[0]  # Remove numeric suffix
    return predictor


def find_common_jsons(subunits_dir):
    subunits = sorted([p for p in Path(subunits_dir).iterdir() if p.is_dir()])

    if not subunits:
        print(f"No subunits found in {subunits_dir}")
        return []

    common_files = None

    for subunit in subunits:
        input_dir = subunit / "input_files"
        if not input_dir.exists():
            print(f"⚠️  Missing input_files in {subunit.name}")
            return []

        json_files = set(f.name for f in input_dir.glob("*.json"))

        if common_files is None:
            common_files = json_files
        else:
            common_files &= json_files  # intersect

    return sorted(common_files)
def combine(ip_dir,json_paths, output_path,predictor,num_models,is_homomultimer,is_default_af3=False):
    """
    Combine multiple JSONs by grouping on identical `sequence`.
    Also rewrites paths with representative chain id.
    """
    grouped = defaultdict(lambda: {
        "ids": [],
        "data": None  # first occurrence
    })

    for path in json_paths:
        with open(path) as f:
            j = json.load(f)

        protein = j["sequences"][0]["protein"]
        seq = protein["sequence"]
        pid = protein["id"]

        if grouped[seq]["data"] is None:
            grouped[seq]["data"] = j  # save first JSON

        grouped[seq]["ids"].append(pid)

    # Build final JSON
    output = {
        "name": list(grouped.values())[0]["data"]["name"],  # assume all have same name
        "sequences": [],
        "modelSeeds": random.sample(range(1, 100000), int(int(num_models)/5)),
        "dialect": "alphafold3",
        "version": 1
        
    }

    for seq, info in grouped.items():
        base_json = info["data"]
        protein = base_json["sequences"][0]["protein"]

        representative_id = info["ids"][0]  # e.g., "A"

        # Adjust paths
        unpaired_path = protein.get("unpairedMsaPath", "")
        if unpaired_path:
            unpaired_path = f"subunits/{representative_id}/input_files/{unpaired_path}"
        
                # Adjust paths
        paired_path = protein.get("pairedMsaPath", "")
        if paired_path:
            paired_path = f"subunits/{representative_id}/input_files/{paired_path}"

        predictor_name_for_template_source = extract_predictor_name(predictor.split(".")[0].lower())


        
        if is_default_af3:
            combined_protein = {
                "id": info["ids"],
                "sequence": seq
            }
        
        else:
            predictor_template_sources = get_template_source(is_homomultimer)
            template_source = predictor_template_sources[predictor_name_for_template_source]
            template_info = get_template_info_from_source(ip_dir,template_source,representative_id,predictor.split(".")[0])
        
            combined_protein = {
                "id": info["ids"],
                "sequence": seq,
                "unpairedMsaPath": unpaired_path,
                # "pairedMsa": "",
                "templates": template_info,
            }
            if paired_path == "":
                combined_protein["pairedMsa"] = ""
            else:
                combined_protein["pairedMsaPath"] = paired_path

        output["sequences"].append({
            "protein": combined_protein
        })

    with open(output_path, "w") as f:
        json.dump(output, f, indent=2)



def combine_jsons(ip_dir,num_models,is_homomultimer):
    subunits_dir = os.path.join(ip_dir,"subunits")
    common_jsons = find_common_jsons(subunits_dir)
    subunits = sorted(os.listdir(subunits_dir))

    for predictor in common_jsons:
        if "afsample" in predictor or "drop" in predictor or "ptm" in predictor:
            continue
        ip_jsons = [f"{subunits_dir}/{sub}/input_files/{predictor}" for sub in subunits]
        op_json = f"{ip_dir}/{predictor}"
        combine(ip_dir,ip_jsons, op_json,predictor,num_models,is_homomultimer,is_default_af3=(predictor=="default_af3.json"))


