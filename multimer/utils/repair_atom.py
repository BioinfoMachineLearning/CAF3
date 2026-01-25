def repair_atom(input_atom, output_atom):
    """
    Repair malformed PDB file by adding missing chain IDs, occupancy, and B-factors
    """
    with open(input_atom, 'r') as infile, open(output_atom, 'w') as outfile:
        for line in infile:
            if line.startswith('ATOM') or line.startswith('HETATM'):
                # Extract parts
                record = line[0:6]
                serial = line[6:11]
                atom_name = line[12:16]
                alt_loc = line[16:17]
                res_name = line[17:20]
                # Chain ID at position 21 (0-indexed)
                original_chain = line[21:22].strip()
                chain = original_chain if original_chain else "A"
                
                res_seq = line[22:26]
                icode = line[26:27]
                x = line[30:38]
                y = line[38:46]
                z = line[46:54]
                
                # Default values if missing
                occupancy = line[54:60].strip() if len(line) > 60 and line[54:60].strip() else '1.00'
                bfactor = line[60:66].strip() if len(line) > 66 and line[60:66].strip() else '0.00'
                
                # Element symbol
                element = line[76:78].strip() if len(line) > 78 else atom_name.strip()[0]
                
                # Reconstruct line
                new_line = f"{record}{serial} {atom_name}{alt_loc}{res_name} {chain}{res_seq}{icode}   {x}{y}{z}{occupancy:>6}{bfactor:>6}          {element:>2}\n"
                outfile.write(new_line)
            else:
                outfile.write(line)