import sys

def filter_cif_by_chain(input_path, chain_id, output_path):
    atom_loop_started = False
    atom_headers = []
    atom_data_start = False

    with open(input_path, 'r') as infile, open(output_path, 'w') as outfile:
        for line in infile:
            # Skip HETATM lines
            if line.startswith("HETATM"):
                continue

            # Detect beginning of atom_site loop
            if line.strip() == "loop_":
                atom_loop_started = True
                atom_headers = []
                atom_data_start = False
                outfile.write(line)
                continue

            # If in atom loop, collect headers
            if atom_loop_started and line.strip().startswith("_atom_site."):
                atom_headers.append(line.strip())
                outfile.write(line)
                continue

            # End of header block = start of data
            if atom_loop_started and not line.strip().startswith("_") and not atom_data_start:
                atom_data_start = True

            # If in atom_site data
            if atom_data_start:
                tokens = line.strip().split()
                if len(tokens) != len(atom_headers):
                    outfile.write(line)  # malformed or non-atom line
                    continue

                # Find chain ID column dynamically
                try:
                    # chain_col = [h for h in atom_headers if "_atom_site.label_asym_id" in h][0]
                    chain_col = [h for h in atom_headers if "_atom_site.auth_asym_id" in h][0]
                    chain_idx = atom_headers.index(chain_col)
                    if tokens[chain_idx] == chain_id:
                        outfile.write(line)
                except IndexError:
                    continue  # Skip if chain index is out of range
                continue

            # Default: write all non-ATOM lines
            outfile.write(line)

