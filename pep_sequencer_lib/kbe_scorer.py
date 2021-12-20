from Bio.PDB import PDBParser


kb_entropy_dict = {("ALA", "ALA"): 0.08, ("ALA", "ILE"): 0.11, ("ALA", "LEU"): 0.10, ("ALA", "VAL"): 0.10, ("ALA", "MET"): 0.10, ("ALA", "PHE"): 0.08, ("ALA", "TRP"): 0.09, ("ALA", "GLY"): 0.11, ("ALA", "PRO"): 0.13, ("ALA", "CYS"): 0.11, ("ALA", "ASN"): 0.12, ("ALA", "GLN"): 0.13, ("ALA", "SER"): 0.12, ("ALA", "THR"): 0.12, ("ALA", "TYR"): 0.09, ("ALA", "ASP"): 0.13, ("ALA", "GLU"): 0.14, ("ALA", "ARG"): 0.12, ("ALA", "HIS"): 0.09, ("ALA", "LYS"): 0.14, ("ILE", "ILE"): 0.13, ("ILE", "LEU"): 0.11, ("ILE", "VAL"): 0.09, ("ILE", "MET"): 0.11, ("ILE", "PHE"): 0.10, ("ILE", "TRP"): 0.10, ("ILE", "GLY"): 0.11, ("ILE", "PRO"): 0.12, ("ILE", "CYS"): 0.11, ("ILE", "ASN"): 0.14, ("ILE", "GLN"): 0.11, ("ILE", "SER"): 0.12, ("ILE", "THR"): 0.13, ("ILE", "TYR"): 0.10, ("ILE", "ASP"): 0.14, ("ILE", "GLU"): 0.13, ("ILE", "ARG"): 0.12, ("ILE", "HIS"): 0.13, ("ILE", "LYS"): 0.15, ("LEU", "LEU"): 0.12, ("LEU", "VAL"): 0.11, ("LEU", "MET"): 0.12, ("LEU", "PHE"): 0.10, ("LEU", "TRP"): 0.10, ("LEU", "GLY"): 0.10, ("LEU", "PRO"): 0.13, ("LEU", "CYS"): 0.13, ("LEU", "ASN"): 0.12, ("LEU", "GLN"): 0.12, ("LEU", "SER"): 0.14, ("LEU", "THR"): 0.13, ("LEU", "TYR"): 0.10, ("LEU", "ASP"): 0.14, ("LEU", "GLU"): 0.15, ("LEU", "ARG"): 0.11, ("LEU", "HIS"): 0.13, ("LEU", "LYS"): 0.15, ("VAL", "VAL"): 0.10, ("VAL", "MET"): 0.11, ("VAL", "PHE"): 0.09, ("VAL", "TRP"): 0.07, ("VAL", "GLY"): 0.10, ("VAL", "PRO"): 0.14, ("VAL", "CYS"): 0.10, ("VAL", "ASN"): 0.12, ("VAL", "GLN"): 0.12, ("VAL", "SER"): 0.12, ("VAL", "THR"): 0.11, ("VAL", "TYR"): 0.09, ("VAL", "ASP"): 0.13, ("VAL", "GLU"): 0.11, ("VAL", "ARG"): 0.11, ("VAL", "HIS"): 0.12, ("VAL", "LYS"): 0.12, ("MET", "MET"): 0.09, ("MET", "PHE"): 0.12, ("MET", "TRP"): 0.07, ("MET", "GLY"): 0.14, ("MET", "PRO"): 0.14, ("MET", "CYS"): 0.11, ("MET", "ASN"): 0.14, ("MET", "GLN"): 0.13, ("MET", "SER"): 0.12, ("MET", "THR"): 0.15, ("MET", "TYR"): 0.10, ("MET", "ASP"): 0.17, ("MET", "GLU"): 0.14, ("MET", "ARG"): 0.12, ("MET", "HIS"): 0.09, ("MET", "LYS"): 0.17, ("PHE", "PHE"): 0.09, ("PHE", "TRP"): 0.09, ("PHE", "GLY"): 0.11, ("PHE", "PRO"): 0.12, ("PHE", "CYS"): 0.11, ("PHE", "ASN"): 0.14, ("PHE", "GLN"): 0.12, ("PHE", "SER"): 0.11, ("PHE", "THR"): 0.11, ("PHE", "TYR"): 0.10, ("PHE", "ASP"): 0.14, ("PHE", "GLU"): 0.12, ("PHE", "ARG"): 0.13, ("PHE", "HIS"): 0.10, ("PHE", "LYS"): 0.15, ("TRP", "TRP"): 0.09, ("TRP", "GLY"): 0.13, ("TRP", "PRO"): 0.08, ("TRP", "CYS"): 0.11, ("TRP", "ASN"): 0.07, ("TRP", "GLN"): 0.12, ("TRP", "SER"): 0.14, ("TRP", "THR"): 0.11, ("TRP", "TYR"): 0.09, ("TRP", "ASP"): 0.14, ("TRP", "GLU"): 0.13, ("TRP", "ARG"): 0.13, ("TRP", "HIS"): 0.10, ("TRP", "LYS"): 0.12, ("GLY", "GLY"): 0.14, ("GLY", "PRO"): 0.15, ("GLY", "CYS"): 0.15, ("GLY", "ASN"): 0.15, ("GLY", "GLN"): 0.16, ("GLY", "SER"): 0.15, ("GLY", "THR"): 0.16, ("GLY", "TYR"): 0.13, ("GLY", "ASP"): 0.16, ("GLY", "GLU"): 0.16, ("GLY", "ARG"): 0.17, ("GLY", "HIS"): 0.10, ("GLY", "LYS"): 0.18, ("PRO", "PRO"): 0.15, ("PRO", "CYS"): 0.18, ("PRO", "ASN"): 0.13, ("PRO", "GLN"): 0.14, ("PRO", "SER"): 0.16, ("PRO", "THR"): 0.17, ("PRO", "TYR"): 0.12, ("PRO", "ASP"): 0.17, ("PRO", "GLU"): 0.14, ("PRO", "ARG"): 0.16, ("PRO", "HIS"): 0.15, ("PRO", "LYS"): 0.20, ("CYS", "CYS"): 0.06, ("CYS", "ASN"): 0.12, ("CYS", "GLN"): 0.15, ("CYS", "SER"): 0.16, ("CYS", "THR"): 0.16, ("CYS", "TYR"): 0.11, ("CYS", "ASP"): 0.17, ("CYS", "GLU"): 0.14, ("CYS", "ARG"): 0.17, ("CYS", "HIS"): 0.13, ("CYS", "LYS"): 0.14, ("ASN", "ASN"): 0.10, ("ASN", "GLN"): 0.13, ("ASN", "SER"): 0.14, ("ASN", "THR"): 0.17, ("ASN", "TYR"): 0.13, ("ASN", "ASP"): 0.14, ("ASN", "GLU"): 0.16, ("ASN", "ARG"): 0.17, ("ASN", "HIS"): 0.14, ("ASN", "LYS"): 0.17, ("GLN", "GLN"): 0.18, ("GLN", "SER"): 0.14, ("GLN", "THR"): 0.13, ("GLN", "TYR"): 0.12, ("GLN", "ASP"): 0.15, ("GLN", "GLU"): 0.18, ("GLN", "ARG"): 0.17, ("GLN", "HIS"): 0.15, ("GLN", "LYS"): 0.17, ("SER", "SER"): 0.15, ("SER", "THR"): 0.16, ("SER", "TYR"): 0.14, ("SER", "ASP"): 0.16, ("SER", "GLU"): 0.14, ("SER", "ARG"): 0.16, ("SER", "HIS"): 0.12, ("SER", "LYS"): 0.18, ("THR", "THR"): 0.16, ("THR", "TYR"): 0.13, ("THR", "ASP"): 0.18, ("THR", "GLU"): 0.17, ("THR", "ARG"): 0.17, ("THR", "HIS"): 0.14, ("THR", "LYS"): 0.17, ("TYR", "TYR"): 0.08, ("TYR", "ASP"): 0.13, ("TYR", "GLU"): 0.14, ("TYR", "ARG"): 0.13, ("TYR", "HIS"): 0.13, ("TYR", "LYS"): 0.15, ("ASP", "ASP"): 0.20, ("ASP", "GLU"): 0.19, ("ASP", "ARG"): 0.13, ("ASP", "HIS"): 0.12, ("ASP", "LYS"): 0.15, ("GLU", "GLU"): 0.18, ("GLU", "ARG"): 0.14, ("GLU", "HIS"): 0.13, ("GLU", "LYS"): 0.15, ("ARG", "ARG"): 0.19, ("ARG", "HIS"): 0.16, ("ARG", "LYS"): 0.20, ("HIS", "HIS"): 0.08, ("HIS", "LYS"): 0.18, ("LYS", "LYS"): 0.19}
for aa_i, aa_j in list(kb_entropy_dict.keys()):
    if aa_i != aa_j:
        kb_entropy_dict[(aa_j, aa_i)] = kb_entropy_dict[(aa_i, aa_j)]


def residues_are_in_contact_kbe(res_i, res_j, threshold=4.5):
    contact = False
    for atm_i in res_i.child_list:
        for atm_j in res_j.child_list:
            if get_atoms_distance(atm_i, atm_j) <= threshold:
                contact = True
                break
        if contact:
            break
    return contact

def get_atoms_distance(atm_i, atm_j):
    return atm_i - atm_j


def get_normalized_score(score, count):
    if count != 0:
        return score/float(count)
    else:
        return 0.0


def score_interface_kbe(pdb_filepath, receptor_chain, peptide_chain):

    structure = PDBParser(QUIET=True).get_structure("complex", pdb_filepath)
    model = structure[0]

    tot_kbe = 0
    tot_kbe_count = 0
    for pep_res in model[peptide_chain]:
        for rec_res in model[receptor_chain]:
            d_ca = pep_res["CA"] - rec_res["CA"] # Get the carbon alpha distance in order to speed up calculations (pairs too far apart will be pre-filtered).
            if d_ca <= 16.0 and residues_are_in_contact_kbe(pep_res, rec_res):
                tot_kbe += kb_entropy_dict.get((pep_res.get_resname(), rec_res.get_resname()))
                tot_kbe_count += 1
    return get_normalized_score(tot_kbe, tot_kbe_count)
