def parse_sto_with_indices(sto_path, one_based_query=False):
    """
    Parse a Stockholm (.sto) file and extract template alignment indices
    compatible with AlphaFold3.
    Rules:
      - '.'  : ignored completely (no advance)
      - '-'  : template gap (query advances)
      - lowercase letters : treated as template gaps (query advances)
      - uppercase letters : aligned residues (advance both)
    Returns:
      - pdb_id
      - template_range (0-based, inclusive)
      - queryIndices
      - templateIndices
      - template_sequence (ungapped, aligned residues only)
    """
    results = {}  # key = "pdb_id/start-end"
    
    with open(sto_path) as f:
        for line in f:
            line = line.rstrip()
            if not line or line.startswith("#"):
                continue
            
            parts = line.split()
            if len(parts) < 2 or "/" not in parts[0]:
                continue
            
            pdb_range = parts[0]      # e.g. 1abc_A/45-120
            aligned_seq = parts[1]
            
            pdb_id, range_str = pdb_range.split("/")
            t_start, t_end = map(int, range_str.split("-"))
            
            # convert to 0-based template numbering
            t_start -= 1
            t_end -= 1
            
            key = f"{pdb_id}/{t_start}-{t_end}"
            
            if key not in results:
                results[key] = {
                    "pdb_id": pdb_id,
                    "template_range": (t_start, t_end),
                    "queryIndices": [],
                    "templateIndices": [],
                    "template_sequence": [],
                    "template_pos": t_start,
                    "query_pos": 0,
                }
            
            entry = results[key]
            
            # Clean template: remove '.' and lowercase letters
            cleaned_template = "".join(c for c in aligned_seq if c != "." and not c.islower())
            
            # Now align: template has only '-' and uppercase letters
            template_pos = entry["template_pos"]
            query_pos = entry["query_pos"]
            
            for c in cleaned_template:
                if c == "-":
                    # Template gap - query advances
                    query_pos += 1
                else:
                    # Uppercase letter - aligned residue, both advance
                    q_idx = query_pos + (1 if one_based_query else 0)
                    entry["queryIndices"].append(q_idx)
                    entry["templateIndices"].append(template_pos)
                    entry["template_sequence"].append(c)
                    query_pos += 1
                    template_pos += 1
            
            entry["template_pos"] = template_pos
            entry["query_pos"] = query_pos
    
    # finalize output
    final_results = []
    for entry in results.values():
        assert len(entry["queryIndices"]) == len(entry["templateIndices"])
        assert len(entry["queryIndices"]) == len(entry["template_sequence"])
        
        final_results.append({
            "pdb_id": entry["pdb_id"],
            "template_range": entry["template_range"],
            "queryIndices": entry["queryIndices"],
            "templateIndices": entry["templateIndices"],
            "template_sequence": "".join(entry["template_sequence"]),
        })
    
    return final_results

