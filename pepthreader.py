import os
import sys
import argparse
import shutil
from Bio.PDB import *
import numpy as np
import pandas as pd
import time
import subprocess
import pickle
import multiprocessing as mp
import re

try:
    from tqdm.auto import tqdm
    tqdm_module = True
except:
    print("tqdm module not found")
    tqdm_module = False


from pep_sequencer_lib import three_to_one
from pep_sequencer_lib.pep_seq_functions import get_seq, save_template, seq_separation
from pep_sequencer_lib.scoring_matrices import *
from pep_sequencer_lib.model_functions import *

from modeller import environ, selection
use_soap_pep = True
from modeller import soap_pp
if use_soap_pep:
    from modeller import soap_peptide
    #check_score_peptide_file()

##############################################################################################################################
################################### PEP SEQUENCER BUILDER ###################################################################
##############################################################################################################################

########## Take input infos ###########################

parser = argparse.ArgumentParser()
parser.add_argument("-t", "--template_filepath", help="Filepath of the template PDB file.", type=str, required=True,)
parser.add_argument("-q", "--query_filepath", help="Filepath of the target FASTA file.", type=str, required=True,)
parser.add_argument("-r", "--receptor_chain", help="PDB chain of the receptor in the template.", type=str, required=True,)
parser.add_argument("-p", "--peptide_chain", help="PDB chain of the peptide in the template.", type=str, required=True,)
parser.add_argument("-o", "--out_dirpath", help="Output directory.", type=str, default=os.getcwd())
parser.add_argument('-S', '--sequence_only', help="Only shows the threaded peptide sequence aligned to the template one.", action='store_const', const=True)
parser.add_argument('-f', '--force', action='store_const', const=True)
parser.add_argument('-mx_top', '--matrix_top_value', default = 50)
parser.add_argument('-rr', '--re_ranking', action = 'store_true')
parser.add_argument('-n', "--n_decoys", help="Number of decoys for peptide-receptor complex.", type=int, default=8)
#parser.add_argument('-M', "--mode", help="TODO.", type=str, default="threading", choices=["threading", "panel"])
# parser.add_argument("-a", "--alignment_filepath", help="Filepath of a FASTA file with the template peptide aligned to the peptide residues.", type=str, default=None)
parser.add_argument("--n_jobs", help="Number of parallel jobs.", type=int, default=1)
parser.add_argument("--n_mod_jobs", help="Number of parallel jobs for building decoys with MODELLER.", type=int, default=1)
parser.add_argument("--parallel_backend", help="Parallel backend.", type=str, default="multiprocessing", choices=["multiprocessing", ])
# parser.add_argument('--use_parallel_soap_pep', help="TODO.", action='store_const', const=True)
parser.add_argument('--pt_modeling', help="Use a automodel class specific for pep_threader.", action='store_const', const=True)
parser.add_argument('-v', '--verbose', action='store_const', const=True)
parser.add_argument('-so', '--score_only', default=False)


cmd = parser.parse_args()
main_dirpath = os.getcwd()


hp_contact_forming_res_list = ['ALA', 'VAL', 'LEU', 'ISO', 'MET', 'PHE', 'TRP', 'PRO', 'TYR']
ch_contact_forming_res_list = ["GLU", "LYS", "ASP", "ARG"]
contact_forming_res_list = hp_contact_forming_res_list + ch_contact_forming_res_list


############### Initialize MODELLER score_peptide extra file ##########################


def get_path_sep():

    if sys.platform == "win32":
        path_sep = "\\"

    elif sys.platform == "linux":
        path_sep = "/"

    elif sys.platform == "darwin":
        path_sep = "/"

    return path_sep


def check_score_peptide_file():

    if sys.platform == "win32":
        lib_name = "Library"
        lib_name2 = "'Program Files'"

    elif sys.platform == "linux":
        lib_name = "lib"
        lib_name2 = "lib"

    elif sys.platform == "darwin":
        lib_name = "lib"
        lib_name2 = "lib"

    path_sep = get_path_sep()

    from pathlib import Path
    import modeller

    # Full path to modeller.py
    path_to_modeller_init = os.path.abspath(modeller.__file__)
    path_to_modeller_init_obj = Path(os.path.abspath(modeller.__file__))

    # Path length (in other terms: the number of directories)
    path_len = len(path_to_modeller_init.split(path_sep))

    # Parent directory of full path
    parent = path_to_modeller_init_obj.parent

    try:

        # Iterate over the path, as many times as the length of the path
        for dr in range(path_len):
            # Split the path.
            a = (str(parent).split(path_sep))

            if a[-1] == 'lib':
                path_to_user_lib = parent

            elif a[-1] == 'modlib':
                path_to_mod_lib = parent
                path_to_user_lib = path_to_mod_lib.parent.parent

            elif a[-1] == "Library":
                path_to_user_lib = parent

            elif a[-1] == "Program Files":
                path_to_user_lib = parent

            # else continue with the parent directory
            else:
                parent = parent.parent

        for libb in os.listdir(path_to_user_lib):
            temp = re.search("modeller", libb, re.IGNORECASE)

            if temp and os.path.isdir(os.path.join(path_to_user_lib, libb)):
                modeller_version = libb

        path_to_score_peptide_file = Path(os.path.join(path_to_user_lib, modeller_version, "modlib", "soap_peptide.hdf5"))

        if path_to_score_peptide_file.is_file():
            print("\n\n\n'soap_peptide.hdf5' file exists")
        else:
            print("\nInitializing Modeller Scorer\nCopying 'soap_peptide.hdf5' file")
            copy_soap_peptide_file(path_to_score_peptide_file)

    except Exception as e:

        print(e)
        print("\nMany attemps of Initializing the Modeller Scorer failed. \nThe 'soap_peptide.hdf5' file must be manually copied in modeller 'modlib' directory to continue with the analysis. \nThe file is located in the PEPthreader main directory. \nAfter having initialized the Modeller Scorer, it is suggested to restart the analysis with the '-so True' option")


def copy_soap_peptide_file(dst_path):

    if sys.platform == "win32":
        path_to_pepthread = os.path.join(os.path.dirname(__file__), "soap_peptide.hdf5")
        shutil.copyfile(path_to_pepthread, dst_path)

    else:
        path_to_pepthread = os.path.join(os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir)), "soap_peptide.hdf5")
        shutil.copyfile(path_to_pepthread, dst_path)



###########################################################
#### Create Matrix for Blosum Modified####################
##############################################################

def matrix_generator(pdb, peptide_chain, receptor_chain, hp_contact_threshold=5.0, ch_contact_threshold=5.0):

    path_sep = get_path_sep()

    print(str("\nGenerating Matrix for: " + "\n  Receptor: " + pdb.split(path_sep)[-1] + "\n  Peptide Chain: " + peptide_chain + "\n  Receptor Chain: " + receptor_chain))

    # Prepare the contact threshold dictionary.
    contact_threshold_dict = {}
    contact_threshold_dict.update(dict([(aa, hp_contact_threshold) for aa in hp_contact_forming_res_list]))
    contact_threshold_dict.update(dict([(aa, ch_contact_threshold) for aa in ch_contact_forming_res_list]))

    #Load PDB structure
    parser = PDBParser()
    structure = parser.get_structure('structure', pdb)
    model = structure[0]

    #Define receptor and peptide chain
    pepchain = model[peptide_chain]
    receptchain = model[receptor_chain]


    pep_interactions = []

    residue_count = 0
    #iterates on the peptide residues
    for residue in pepchain:
        res_id = residue.get_id()

        # exclude het atoms
        if res_id[0] == 'W' or res_id[0].startswith('H_'):
            continue

        residue_count += 1

        # Number of contacts of the peptide residue with the receptor chain.
        interactions_count = 0

        for chain_res in receptchain:

            if not chain_res.resname in contact_forming_res_list:
                continue

            interaction_found = False
            contact_threshold = contact_threshold_dict[chain_res.resname]
            for atom in residue:

                for chain_atom in chain_res:
                    distance = atom - chain_atom

                    if distance <= contact_threshold:
                        interactions_count += 1
                        interaction_found = True
                        break

                if interaction_found:
                    break

        pep_interactions.append(interactions_count)

    # Writes a file with the weights for each residue in the peptide.
    final_matrix = pep_interactions

    matrix_filepath = os.path.join(cmd.out_dirpath, 'matrix.txt')

    matrix_file = open(matrix_filepath, 'w')
    matrix_file.write(str(final_matrix))
    matrix_file.close()


######################################################################################
####### Build working space ###############################
######################################################################################

# Create output dir
results_dirpath = cmd.out_dirpath

bold_start = '\033[1m\033[4m'
bold_end   = '\033[0m\033[24m'

a = """
Prediction and modeling of protein-peptide interactions using a template-based approach

To carry out the analysis, the structure of a protein-peptide complex and a protein sequence are needed.

For further information visit: --github link--
"""""

print("\n", bold_start, "PepThreader", bold_end)
print(a)

print(str("\nOutput directory: " + results_dirpath))


# Remove file from dir if already exsists
if cmd.force:
    shutil.rmtree(results_dirpath)
if not os.path.exists(results_dirpath):
    os.makedirs(results_dirpath)

template_name = 'template_complex'
template_complex_filepath = os.path.join(cmd.out_dirpath, "%s.pdb" % template_name)


#Build the template complex with only receptor chain and pep chain and matrix for Blosum Modified
t1 = time.time()
template_receptor_sequence, template_peptide_sequence = save_template(cmd.template_filepath,
                                                                      template_complex_filepath,
                                                                      cmd.receptor_chain,
                                                                      cmd.peptide_chain,
                                                                      )
matrix_generator(cmd.template_filepath, cmd.peptide_chain, cmd.receptor_chain)
print("\nTo prepare input files, it took", str(round(time.time()-t1, 3)) + " s")


# Divide fasta sequence in n fragments
peptide_name = 'fragment'
peptide_length = len(template_peptide_sequence)
query_sequence = get_seq(cmd.query_filepath)
query_peptides = seq_separation(query_sequence, peptide_length)


#############################################################################################################
#################################  SCORER ###################################################################
############################################################################################################à

scores_dict = {'FRAGMENT': [],'PEPTIDE_SEQ': [], 'TEMPLATE_SEQ': []}
os.chdir(results_dirpath)

# iterates matrices from scoring_matrices.py and scores peptide
for matrix_name in matrices_names:
    scores_dict[matrix_name] = []
    scores_dict[matrix_name + '_MODIFIED'] = []
    scores_dict[matrix_name + '_GALAXY'] = []


    for pep_id, pep_seq in enumerate(query_peptides, 0):

        scores_dict[matrix_name].append(compute_alignment_score(template_peptide_sequence, pep_seq, matrix_name=matrix_name, mode = 'normal'))
        scores_dict[matrix_name + '_MODIFIED'].append(compute_alignment_score(template_peptide_sequence, pep_seq, matrix_name=matrix_name, mode = 'modified', mod_type = 'simple'))
        scores_dict[matrix_name + '_GALAXY'].append(compute_alignment_score(template_peptide_sequence, pep_seq, matrix_name=matrix_name, mode = 'modified', mod_type = 'galaxy'))

for pep_id, pep_seq in enumerate(query_peptides, 0):

    scores_dict['PEPTIDE_SEQ'].append(pep_seq)
    scores_dict['TEMPLATE_SEQ'].append(template_peptide_sequence)
    scores_dict['FRAGMENT'].append(pep_id+1)
    #scores_dict['BLOSUM_MODIFIED'].append(blosum_modified(pep_seq, template_peptide_sequence))

scores_df = pd.DataFrame(scores_dict)
scores_df.to_csv('Scoring_results.csv', index=False)


####### MATRIX DATA ANALYSIS#############

dict_scoring = {}

matrices = []

for matrix in matrices_names:
    matrices.append(matrix)

    matrices.append(matrix + '_MODIFIED')
    matrices.append(matrix + '_GALAXY')

for column in matrices:
    values = scores_df[column]
    fragments_ids = scores_df['FRAGMENT']
    peptides_seqs = scores_df['PEPTIDE_SEQ']
    dict_scoring['FRAGMENT'] = fragments_ids
    dict_scoring['PEPTIDE_SEQ'] = peptides_seqs

    rank_values = values.rank(ascending = False, method = 'first')

    if column in dict_scoring:
        dict_scoring[column].append(rank_values.tolist())
    else:
        dict_scoring[column] = rank_values.tolist()


rank_df = pd.DataFrame.from_dict(dict_scoring)
matrix_ranks = rank_df.to_csv('matrix_ranks.csv')

### Displays the results
peptide_index_dict = {}
for column in matrices:
    #print(column)
    #print(rank_df[column].sort_values(ascending = True)[:cmd.matrix_top_value].to_string())


    if column == 'BLOSUM62_MODIFIED':
        peptide_index_list = rank_df[column].sort_values(ascending = True)[:int(cmd.matrix_top_value)].index.to_list()

        #print(rank_df[column].sort_values(ascending = True)[:cmd.matrix_top_value].index +1)
        peptide_index_dict[column] = peptide_index_list
    else:
        continue

if not cmd.re_ranking:
    print("\nEND of ANALYSIS\n\nThe results can be found at " + cmd.out_dirpath + "\n\nTo continue with re-ranking, please enable -rr option")

##########     RE RANKING #####################

python_exec = sys.executable
build_models_script_filepath = os.path.join(main_dirpath, "_build_peptide.py")




def launch_external_modeling(peptide_sequence, peptide_directory, peptide_name):

    t1 = time.time()

    input_filepath = os.path.join(peptide_directory, peptide_name + ".pkl")
    i_fh = open(input_filepath, "wb")
    pickle.dump({"peptide_directory": peptide_directory,
                 "peptide_name": peptide_name,
                 "peptide_sequence": peptide_sequence,

                 "template_name": 'template_complex',
                 "supress_output": not cmd.verbose,
                 "n_decoys": cmd.n_decoys,
                 "n_mod_jobs": cmd.n_mod_jobs,
                 "pt_modeling": cmd.pt_modeling,}, i_fh)


    i_fh.close()

    try:
        subprocess.check_call([python_exec, build_models_script_filepath, input_filepath])

    except subprocess.CalledProcessError as e:
        print(e.output)



def launch_external_modeling_parallel(peptide_tuple):

    launch_external_modeling(peptide_tuple[0], peptide_tuple[1], peptide_tuple[2])



def launch_external_scoring(peptide_dirnames, peptides_count):

    t1 = time.time()

    # print("\n# Scoring %s models of %s (batch %s)." % (len(peptide_dirnames), len(dirnames_list), peptides_count))

    input_filename = "score_batch_%s.pkl" % peptides_count
    i_fh = open(input_filename, "wb")
    pickle.dump({"peptide_dirpaths": [os.path.join(matrix_dirpath, peptide_dirname) for peptide_dirname in peptide_dirnames],
                 "n_score_jobs": cmd.n_jobs,
                 "supress_output": not cmd.verbose,}, i_fh)
    i_fh.close()

    try:
        subprocess.check_call([python_exec, score_models_script_filepath, input_filename])

    except subprocess.CalledProcessError as e:
        print(e.output)

    os.remove(input_filename)
    # print("- It took %s." % (time.time()-t1))



def launch_external_scoring_parallel(peptide_tuple):

    launch_external_scoring(peptide_tuple[1], peptide_tuple[0])



def split_seq(seq, n_jobs):
    """
    Splits a 'seq' in 'n_jobs' groups of items.
    """

    newseq = []
    splitsize = int(len(seq)/n_jobs)
    for i in range(0, n_jobs):
        newseq.append(seq[i*splitsize:i*splitsize+splitsize])

    newseq_count = sum(len(s) for s in newseq)
    if newseq_count < len(seq):
        # print("Error %s %s" % (newseq_count, len(seq)))
        add_count = 0
        for element in seq[newseq_count:]:
            newseq[add_count].append(element)
            add_count += 1
            if add_count == n_jobs:
                add_count = 0

    return newseq


if cmd.re_ranking:
    #modeling = True

    if not cmd.score_only:

        if '.' in cmd.out_dirpath:
            raise Exception('full path needed for output directory')

        peptides_dict = {}
        #create dir for matrix models

        print("\nStarting Modeling of top peptides\n")

        for key in peptide_index_dict.keys():

            matrix_dir = os.path.join(cmd.out_dirpath, key + '_MODELS')

            if not os.path.isdir(matrix_dir):
                os.mkdir(matrix_dir)

            peptide_list = []
            for rank, pos in enumerate(peptide_index_dict[key]):
                peptide_seq = rank_df[rank_df['FRAGMENT'] == pos+1]['PEPTIDE_SEQ'].values[0]
                peptide_rank = rank+1
                peptide_list.append((peptide_seq, peptide_rank, pos+1))

            # Iter through top peptides.
            peptide_jobs = []
            for rank_idx in range(int(cmd.matrix_top_value)):

                rank = rank_idx + 1

                #create dir for every top value
                model_dirname = 'matrix_rank_' + str(rank)
                model_dir = os.path.join(matrix_dir, model_dirname)
                if not os.path.isdir(model_dir):
                    os.mkdir(model_dir)

                #copy the template complex structure in model dir
                source = template_complex_filepath # cmd.template_filepath
                destination = os.path.join(model_dir, "%s.pdb" % template_name)
                shutil.copy(source, destination)


                pep_seq, pep_rank, pep_seq_pos = peptide_list[rank_idx]

                generate_ali_file(model_dir, 'template_complex', peptide_name,
                                template_receptor_sequence, template_peptide_sequence, pep_seq, cmd.receptor_chain, cmd.peptide_chain)


                peptide_jobs.append([pep_seq, model_dir, peptide_name])

            # Actually runs MODELLER for top peptides.
            if cmd.n_jobs == 1:

                iteration = 1

                if tqdm_module:
                    bar = tqdm(total=len(peptide_jobs), position=1, dynamic_ncols=True, leave=True, unit='file', desc="Progress", bar_format='{l_bar}{bar}| {n_fmt}/{total_fmt} [{elapsed}]')

                for peptide_seq, model_dir, peptide_name in peptide_jobs:

                    launch_external_modeling(peptide_seq, model_dir, peptide_name)

                    if tqdm_module:
                        bar.update()

            else:

                if tqdm_module:

                    pool = mp.Pool(cmd.n_jobs)

                    results = []

                    for result in tqdm(pool.imap_unordered(launch_external_modeling_parallel, peptide_jobs), total=len(peptide_jobs), position=1, dynamic_ncols=True, leave=True, unit='file', desc="Progress", bar_format='{l_bar}{bar}| {n_fmt}/{total_fmt} [{elapsed}]'):
                        results.append(result)

                    pool.close()
                    pool.join()

                else:

                    pool = mp.Pool(cmd.n_jobs)

                    results = pool.map(launch_external_modeling_parallel, peptide_jobs)
                    pool.close()



    ###############################################################
    ################ Scoring #####################################
    ##############################################################

    python_exec = sys.executable

    # Initializes SOAP scores.
    sp_pp = soap_pp.PairScorer()
    sp_atm = soap_pp.AtomScorer()

    if use_soap_pep:
        sp_pep = soap_peptide.Scorer()

    score_models_script_filepath = os.path.join(main_dirpath, "_score_peptide.py")

    peptide_name = 'matrix_rank'

    check_score_peptide_file()


    scoring = True

    for matrix_dir in os.listdir(results_dirpath):
        frames = []

        if os.path.isdir(matrix_dir):
            matrix_dirpath = os.path.join(results_dirpath, matrix_dir)
            dirnames_list = []


            for top_value_dir in os.listdir(matrix_dirpath):
                if top_value_dir.startswith('matrix_rank'):
                    top_value_dirpath = os.path.join(matrix_dirpath, top_value_dir)
                    #dirnames_list.append(top_value_dir)
                    dirnames_list.append(top_value_dir)
                else:
                    continue


            dirnames_list = list(sorted(dirnames_list, key=lambda fn: int(fn.split("_")[2])))


            dirnames_batches = split_seq(dirnames_list, n_jobs=1)


            peptides_count = 0


            for peptide_dirnames in dirnames_batches:

                launch_external_scoring(peptide_dirnames, peptides_count)
                peptides_count += 1


            for peptide_dirname in dirnames_list:

                rank_list = []
                rank_list.extend([peptide_dirname for i in range(cmd.n_decoys)])

                #print("* Reading in %s." % peptide_dirname)

                df = pd.read_csv(os.path.join(matrix_dirpath, peptide_dirname, 'results.score.csv'))
                df['RANK'] = rank_list

                frames.append(df)

            finalscores_df = pd.concat(frames)
            finalscores_df.to_csv(os.path.join(matrix_dirpath, 'final_scores.csv'), index=False, index_label = False)

            ###########################################################################################################
            #################### MODEL RANKING #######################################################################
            ##########################################################################################################à

            mx_to_remove = [x for x in matrices_names if x != matrix_dir[:-7]]

            scoring_path = os.path.join(matrix_dirpath, 'final_scores.csv')
            scoring_data = pd.read_csv(scoring_path)
            #scoring_data.groupby('RANK')

            mean_score_path = os.path.join(matrix_dirpath, 'mean_scores.csv')

            df_list = []
            for name, group in scoring_data.groupby('RANK'):
                #for column in group.columns:
                    #if column in ['PEPTIDE_NAME','DECOY_NAME','DECOY_NUM','BUILT', 'PEPTIDE_SEQ', 'RANK']:
                        #group_dropped = group.drop(column, axis = 1)
                        #continue

                mean_df = pd.DataFrame(group.mean()).transpose()
                #print(group.columns)
                #print(group.columns)
                mean_df['RANK'] = group.iloc[0,11] #change according to the number of scoring function that you are using (TEST)
                mean_df['PEPTIDE_SEQ'] = group.iloc[0,4]

                df_list.append(mean_df)
            concat = pd.concat(df_list)
            concat = concat.drop('BUILT', axis = 1)
            concat = concat.drop('DECOY_NUM', axis = 1)
            concat = concat.drop(mx_to_remove, axis = 1)

            concat.to_csv(mean_score_path, index = False)

            for column in concat.columns:
                if column in ['RANK', 'PEPTIDE_SEQ']:
                    continue
                #if column in ['VOROMQA_TOT', 'DOPE_INT']:
                #    concat[column + '_RANK'] = concat[column].rank(method='first', ascending=False)
                else:
                    concat[column + '_RANK'] = concat[column].rank(method='first', ascending=True)

            final_rank_path = os.path.join(matrix_dirpath, 'rank_scores.csv')
            concat.to_csv(final_rank_path, index = False)

            print("\nEND of ANALYSIS\n The results can be found at:\n    " + final_rank_path)
