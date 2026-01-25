# Copyright 2021 DeepMind Technologies Limited
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#      http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""Functions for building the features for the AlphaFold multimer model."""

import collections
import contextlib
import copy
import dataclasses
import json
import os
import tempfile
from typing import Mapping, MutableMapping, Sequence

from absl import logging
import numpy as np
import pickle
import pandas as pd
from template_search import hmmsearch
from template_search import templates
from template_search import parsers
import subprocess
from template_search.parse_stockholm_template import parse_sto_with_indices
from utils.utils import makedir_if_not_exists
from common.config import env_dir,database_dir,max_template_date,template_max_hits


template_searcher = hmmsearch.Hmmsearch(
    binary_path=os.path.join(env_dir, 'hmmsearch'),
    hmmbuild_binary_path=os.path.join(env_dir, 'hmmbuild'),
    database_path=os.path.join(database_dir, 'pdb_seqres/pdb_seqres.txt'))

template_featurizer = templates.HmmsearchHitFeaturizer(
    mmcif_dir=os.path.join(database_dir, 'pdb_mmcif/mmcif_files'),
    max_template_date=max_template_date,
    max_hits=template_max_hits,
    kalign_binary_path=os.path.join(env_dir, 'kalign'),
    release_dates_path=None,
    obsolete_pdbs_path=os.path.join(database_dir, 'pdb_mmcif/obsolete.dat'))

def hmm_runner(ip_sto,op_hmm,op_sto):
    # op_hmm_path = ip_sto.replace(".sto",".hmm")
    hmmbuild_cmd = [
    f"{env_dir}/hmmbuild",
    "--hand",
    "--amino",
    op_hmm,
    ip_sto,
    ]
    # print(f"Launching: {hmmbuild_cmd}")

    result = subprocess.run(
        hmmbuild_cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        check=True
    )

    # print("STDOUT:\n", result.stdout)

    hmmsearch_cmd = [
    f"{env_dir}/hmmsearch",
    "--noali",
    "--cpu", "8",
    "--F1", "0.1",
    "--F2", "0.1",
    "--F3", "0.1",
    "--incE", "100",
    "-E", "100",
    "--domE", "100",
    "--incdomE", "100",
    "-A", op_sto,
    op_hmm,
    f"{database_dir}/pdb_seqres/pdb_seqres.txt",
    ]
    result = subprocess.run(
    hmmsearch_cmd,
    stdout=subprocess.PIPE,
    stderr=subprocess.PIPE,
    text=True,
    check=True
    )
    # print(f"Launching: {hmmsearch_cmd}")
    # print("STDOUT:\n", result.stdout)

from difflib import SequenceMatcher

def normalize_seq(seq):
    """Remove gaps/dots for sequence comparison."""
    return seq.replace("-", "").replace(".", "")

def sequence_similarity(seq1, seq2):
    """Return similarity between 0 and 1 ignoring gaps/dots."""
    seq1 = normalize_seq(seq1)
    seq2 = normalize_seq(seq2)
    return SequenceMatcher(None, seq1, seq2).ratio()

def extract_top_templates(template_dict, templates_from_sto, min_similarity=0.9):
    """
    Match top templates:
      - If pdb_id is unique in template_dict, match on pdb_id
      - Else, match based on sequence similarity
    Returns a list of dicts with:
      template_name, query_indices, template_indices, template_sequence
    """
    # Check if pdb_id are unique
    pdb_ids = list(template_dict.keys())
    unique_pdb = len(pdb_ids) == len(set(pdb_ids))

    top_templates = []

    if unique_pdb:
        # Match by pdb_id only
        top_pdbs = set(template_dict.keys())
        for t in templates_from_sto:
            if t["pdb_id"] in top_pdbs:
                top_templates.append({
                    "template_name": t["pdb_id"],
                    "query_indices": t["queryIndices"],
                    "template_indices": t["templateIndices"],
                    # "template_sequence": t["template_sequence"],
                })
    else:
        # Multiple entries per pdb_id → match by sequence similarity
        for t in templates_from_sto:
            top_seq = template_dict.get(t["pdb_id"])
            if top_seq:
                sim = sequence_similarity(top_seq, t["template_sequence"])
                if sim >= min_similarity:
                    top_templates.append({
                        "template_name": t["pdb_id"],
                        "query_indices": t["queryIndices"],
                        "template_indices": t["templateIndices"],
                        # "template_sequence": t["template_sequence"],
                        # "similarity": sim,
                    })

    return top_templates


def _process_single_chain(sequence,chain_template_sto,output_path,chain,destination_template_csv):
    output_path_temp = os.path.join(output_path,"temp")
    makedir_if_not_exists(output_path_temp)
    op_sto = os.path.join(output_path_temp,f"{chain}_pdb_hits.sto")
    op_hmm = os.path.join(output_path_temp,f"{chain}_query.hmm")

    print(f"Running pdb_seqres template search for chain: {chain}")
    hmm_runner(chain_template_sto,op_hmm,op_sto)



    msa_for_templates = parsers.truncate_stockholm_msa(chain_template_sto,
                                                    max_sequences=50000)
    msa_for_templates = parsers.deduplicate_stockholm_msa(msa_for_templates)
    msa_for_templates = parsers.remove_empty_columns_from_stockholm_msa(msa_for_templates)

    pdb_templates_result = template_searcher.query(msa_for_templates)


    pdb_template_hits = template_searcher.get_template_hits(
        output_string=pdb_templates_result, input_sequence=sequence)


    chain_template_results = template_featurizer.get_templates(
        query_sequence=sequence,
        hits=pdb_template_hits)

    templates_pdb_id = [item.decode() for item in chain_template_results.features["template_domain_names"]]

    template_sequences = [item.decode() for item in chain_template_results.features["template_sequence"]]

    template_dict = {}
    for template,sequence in zip(templates_pdb_id,template_sequences):
        template_dict[template]=sequence

    templates_from_sto = parse_sto_with_indices(op_sto)

    results = extract_top_templates(template_dict,templates_from_sto)
    df = pd.DataFrame(results)
    df.to_csv(destination_template_csv, index=False)
