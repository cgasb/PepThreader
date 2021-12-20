from pep_sequencer_lib import three_to_one

def save_template(template_filepath, output_filepath, receptor_chain, peptide_chain, verbose=False):


    # Parses the template filepath to extract the receptor and peptide atomic coordinates.
    receptor_lines = []
    peptide_lines = []

    template_read = open(template_filepath, 'r')

    for line in template_read.readlines():

        if line.startswith(('ATOM', 'HETATM')):
            res_name = line[17:20]
            chain = line[21]
            res_symbol = three_to_one.get(res_name)
            if res_symbol:
                if chain == receptor_chain:
                    receptor_lines.append(line)
                elif chain == peptide_chain:
                    peptide_lines.append(line)
            else:
                if verbose:
                    print("- Unknown res:", res_name)

    template_read.close()

    # Checks if the input PDB file contains the right chains.
    if len(receptor_lines) == 0:
        raise KeyError("The template '%s' does not have a receptor chain named as '%s'." % (template_filepath, receptor_chain))
    if len(peptide_lines) == 0:
        raise KeyError("The template '%s' does not have a peptide chain named as '%s'." % (template_filepath, peptide_chain))


    # Writes a new template file with only the receptor and peptide atoms.
    # Also extracts the receptor and peptide sequences from their atomic coordinates.
    new_lines = []
    atom_count = 0
    template_seqs_dict = {receptor_chain: "", peptide_chain: ""}
    template_pep_seq = ""

    old_res_num = None
    for line in receptor_lines + peptide_lines:
        res_name = line[17:20]
        chain = line[21]
        res_num = line[22:26]

        new_line = line[0:6] + str(atom_count + 1).rjust(5, " ") + line[11:]
        new_lines.append(line)

        if old_res_num:
            if old_res_num != res_num:
                template_seqs_dict[chain] += three_to_one.get(res_name)
            old_res_num = res_num
        else:
            template_seqs_dict[chain] += three_to_one.get(res_name)
            old_res_num = res_num

        atom_count += 1

    if output_filepath != None:
        t_fh = open(output_filepath, "w")
        t_fh.writelines(new_lines)
        t_fh.close()

    # Returns the receptor and peptide sequences of the template.
    return template_seqs_dict[receptor_chain], template_seqs_dict[peptide_chain]


def get_seq(query_filepath):
    """
    Get a sequence from a single-entry FASTA file.
    """
    f = open(query_filepath,'r')
    lines = f.readlines()
    lines = [row.strip() for row in lines] # restituisce una copia della stringa con i caratteri iniziali e finali rimossi per le righe in fila
    lines = [row.replace('\n', '') for row in lines if not row.startswith('>')]

    return ''.join(lines) # unisci il tutto e stampa


def seq_separation(template_sequence, peptide_length):
    p = 0
    fragment_list = []
    for aa_id in template_sequence:
        fragment = (template_sequence [p:p+peptide_length])

        if len(fragment) == peptide_length:
            fragment_list.append(fragment)
        elif len(fragment) < peptide_length:
            if peptide_length-len(fragment) <= 3:
                fragment_list.append(fragment + '-'*(peptide_length-len(fragment)))
        else:
            continue
        p = p+1

    return fragment_list
