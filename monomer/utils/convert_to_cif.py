from Bio import SeqIO, pairwise2
from Bio.PDB.MMCIFParser import MMCIFParser
from Bio.PDB.PDBParser import PDBParser
from Bio.SeqUtils import seq1
import json
import os

from Bio.PDB import PDBParser, MMCIFIO
import os


no_temp_seq_list = list()

from Bio import PDB

restype_3to1 = {
    'ALA': 'A',
    'ARG': 'R',
    'ASN': 'N',
    'ASP': 'D',
    'CYS': 'C',
    'GLN': 'Q',
    'GLU': 'E',
    'GLY': 'G',
    'HIS': 'H',
    'ILE': 'I',
    'LEU': 'L',
    'LYS': 'K',
    'MET': 'M',
    'PHE': 'F',
    'PRO': 'P',
    'SER': 'S',
    'THR': 'T',
    'TRP': 'W',
    'TYR': 'Y',
    'VAL': 'V',
    'UNK' : 'U',
}





def atom2cif(template_pdb, template_cif):

    # Load the PDB file
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('structure', template_pdb)

    # Write the initial mmCIF file
    io = MMCIFIO()
    io.set_structure(structure)
    io.save(template_cif)

    release_date = "1999-01-01"  # change release date. this date is safe and AF3 includes it.
    
    # Need below for AF3 to find the release date. The cutoff is ~ 2021
    metadata_entries = f"""
    #
    _pdbx_database_status.recvd_initial_deposition_date {release_date}
    _struct_ref_seq.seq_release_date {release_date}
    _pdbx_audit_revision_history.revision_date {release_date}
    #
    """

    # Append the release date to the mmCIF file
    with open(template_cif, "a") as cif_file:
        cif_file.write(metadata_entries)


