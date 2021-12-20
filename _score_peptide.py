import os
import sys
import pickle
import random
import subprocess
import multiprocessing

import pandas as pd

try:
    from tqdm.auto import tqdm
    tqdm_module = True
except:
    print("tqdm module not found")
    tqdm_module = False

from modeller import environ, selection
from modeller.scripts import complete_pdb

#from pep_sequencer_lib.cb_scoring_function import cb_scoring_function

import numpy as np
import ast


use_soap_pep = True
from modeller import soap_pp

if use_soap_pep:
    from os.path import exists
    from modeller import soap_peptide

from pep_sequencer_lib.scoring_matrices import compute_alignment_score, matrices_names
#from pep_sequencer_lib.kbe_scorer import score_interface_kbe

voronota_filepath = None

#-------------
# Functions. -
#-------------

def goap(pdb_filepaths, goap_filepath, return_all_scores=False, permissive=False):

    # TODO: change its behaviour so that all complexes are scored in a single batch.
    if not len(set([os.path.dirname(fp) for fp in pdb_filepaths])) == 1:
        raise Exception("Not all PDB files are in the same directory.")

    pdb_dirpath = os.path.dirname(pdb_filepaths[0])
    pdb_filenames = [os.path.basename(fp) for fp in pdb_filepaths]
    original_cwd = os.getcwd()
    if pdb_dirpath != "":
        os.chdir(pdb_dirpath)

    random_id = "".join([random.choice("1234567890qwertyuiopasdfghjklzxcvbnm_") for i in range(0, 10)])
    temp_filename = "__goap_temp_%s__.txt" % random_id
    temp_fh = open(temp_filename, "w")
    temp_fh.write("%s\n" % os.path.dirname(goap_filepath))
    temp_fh.write("\n".join(pdb_filenames).rstrip())
    temp_fh.close()

    successful = None
    try:
        # out = subprocess.check_output("%s<%s" % (goap_filepath, temp_filename), shell=True)
        temp_fh = open(temp_filename, "r")
        out = subprocess.check_output(goap_filepath, stdin=temp_fh, stderr=subprocess.STDOUT)
        successful = True
    except subprocess.CalledProcessError as e:
        successful = False

    temp_fh.close()
    os.remove(temp_filename)
    if pdb_dirpath != "":
        os.chdir(original_cwd)

    if not successful:
        if not permissive:
            raise Exception("Encountered an error when running GOAP:", e, e.returncode)

    if successful:
        goap_scores = [float(line.split()[-3]) for line in out.split("\n")[0:-1]]
        dfire_scores = [float(line.split()[-2]) for line in out.split("\n")[0:-1]]
        goap_ag_scores = [float(line.split()[-1]) for line in out.split("\n")[0:-1]]
    else:
        goap_scores = [0]*len(pdb_filepaths)
        dfire_scores = [0]*len(pdb_filepaths)
        goap_ag_scores = [0]*len(pdb_filepaths)

    if return_all_scores:
        return goap_scores, dfire_scores, goap_ag_scores
    else:
        return goap_scores


def assess_voromqa(pdb_filepaths, voronota_filepath=voronota_filepath, n_jobs=1, verbose=False):

    voronota_dirpath = os.path.dirname(voronota_filepath)
    voromqa_script_filepath = os.path.join(voronota_dirpath, "voronota-voromqa")
    scores = []

    if n_jobs == 1:
        for filepath in pdb_filepaths:
            scores.append(_launch_voromqa(voromqa_script_filepath, filepath, verbose=verbose))

    else:
        pool = multiprocessing.Pool(n_jobs)
        scores = pool.map(_launch_voromqa_par, [(voromqa_script_filepath, fp, verbose) for fp in pdb_filepaths])
        pool.close()

    return scores

def _launch_voromqa(voromqa_script_filepath, pdb_filepath, verbose=False):
    if verbose:
        print ("- Assessing with VoroMQA: %s." % pdb_filepath)
    out = subprocess.check_output([voromqa_script_filepath, "-i", pdb_filepath])

    if verbose:
        print (out.split())

    return float(out.split()[-3])

def _launch_voromqa_par(t):
    return _launch_voromqa(t[0], t[1])

#################################################################
# Blosum modified alignment
###############################################################

blosum62 = {
              ('W', 'F'): 1, ('L', 'R'): -2, ('S', 'P'): -1, ('V', 'T'): 0,
              ('Q', 'Q'): 5, ('N', 'A'): -2, ('Z', 'Y'): -2, ('W', 'R'): -3,
              ('Q', 'A'): -1, ('S', 'D'): 0, ('H', 'H'): 8, ('S', 'H'): -1,
              ('H', 'D'): -1, ('L', 'N'): -3, ('W', 'A'): -3, ('Y', 'M'): -1,
              ('G', 'R'): -2, ('Y', 'I'): -1, ('Y', 'E'): -2, ('B', 'Y'): -3,
              ('Y', 'A'): -2, ('V', 'D'): -3, ('B', 'S'): 0, ('Y', 'Y'): 7,
              ('G', 'N'): 0, ('E', 'C'): -4, ('Y', 'Q'): -1, ('Z', 'Z'): 4,
              ('V', 'A'): 0, ('C', 'C'): 9, ('M', 'R'): -1, ('V', 'E'): -2,
              ('T', 'N'): 0, ('P', 'P'): 7, ('V', 'I'): 3, ('V', 'S'): -2,
              ('Z', 'P'): -1, ('V', 'M'): 1, ('T', 'F'): -2, ('V', 'Q'): -2,
              ('K', 'K'): 5, ('P', 'D'): -1, ('I', 'H'): -3, ('I', 'D'): -3,
              ('T', 'R'): -1, ('P', 'L'): -3, ('K', 'G'): -2, ('M', 'N'): -2,
              ('P', 'H'): -2, ('F', 'Q'): -3, ('Z', 'G'): -2, ('X', 'L'): -1,
              ('T', 'M'): -1, ('Z', 'C'): -3, ('X', 'H'): -1, ('D', 'R'): -2,
              ('B', 'W'): -4, ('X', 'D'): -1, ('Z', 'K'): 1, ('F', 'A'): -2,
              ('Z', 'W'): -3, ('F', 'E'): -3, ('D', 'N'): 1, ('B', 'K'): 0,
              ('X', 'X'): -1, ('F', 'I'): 0, ('B', 'G'): -1, ('X', 'T'): 0,
              ('F', 'M'): 0, ('B', 'C'): -3, ('Z', 'I'): -3, ('Z', 'V'): -2,
              ('S', 'S'): 4, ('L', 'Q'): -2, ('W', 'E'): -3, ('Q', 'R'): 1,
              ('N', 'N'): 6, ('W', 'M'): -1, ('Q', 'C'): -3, ('W', 'I'): -3,
              ('S', 'C'): -1, ('L', 'A'): -1, ('S', 'G'): 0, ('L', 'E'): -3,
              ('W', 'Q'): -2, ('H', 'G'): -2, ('S', 'K'): 0, ('Q', 'N'): 0,
              ('N', 'R'): 0, ('H', 'C'): -3, ('Y', 'N'): -2, ('G', 'Q'): -2,
              ('Y', 'F'): 3, ('C', 'A'): 0, ('V', 'L'): 1, ('G', 'E'): -2,
              ('G', 'A'): 0, ('K', 'R'): 2, ('E', 'D'): 2, ('Y', 'R'): -2,
              ('M', 'Q'): 0, ('T', 'I'): -1, ('C', 'D'): -3, ('V', 'F'): -1,
              ('T', 'A'): 0, ('T', 'P'): -1, ('B', 'P'): -2, ('T', 'E'): -1,
              ('V', 'N'): -3, ('P', 'G'): -2, ('M', 'A'): -1, ('K', 'H'): -1,
              ('V', 'R'): -3, ('P', 'C'): -3, ('M', 'E'): -2, ('K', 'L'): -2,
              ('V', 'V'): 4, ('M', 'I'): 1, ('T', 'Q'): -1, ('I', 'G'): -4,
              ('P', 'K'): -1, ('M', 'M'): 5, ('K', 'D'): -1, ('I', 'C'): -1,
              ('Z', 'D'): 1, ('F', 'R'): -3, ('X', 'K'): -1, ('Q', 'D'): 0,
              ('X', 'G'): -1, ('Z', 'L'): -3, ('X', 'C'): -2, ('Z', 'H'): 0,
              ('B', 'L'): -4, ('B', 'H'): 0, ('F', 'F'): 6, ('X', 'W'): -2,
              ('B', 'D'): 4, ('D', 'A'): -2, ('S', 'L'): -2, ('X', 'S'): 0,
              ('F', 'N'): -3, ('S', 'R'): -1, ('W', 'D'): -4, ('V', 'Y'): -1,
              ('W', 'L'): -2, ('H', 'R'): 0, ('W', 'H'): -2, ('H', 'N'): 1,
              ('W', 'T'): -2, ('T', 'T'): 5, ('S', 'F'): -2, ('W', 'P'): -4,
              ('L', 'D'): -4, ('B', 'I'): -3, ('L', 'H'): -3, ('S', 'N'): 1,
              ('B', 'T'): -1, ('L', 'L'): 4, ('Y', 'K'): -2, ('E', 'Q'): 2,
              ('Y', 'G'): -3, ('Z', 'S'): 0, ('Y', 'C'): -2, ('G', 'D'): -1,
              ('B', 'V'): -3, ('E', 'A'): -1, ('Y', 'W'): 2, ('E', 'E'): 5,
              ('Y', 'S'): -2, ('C', 'N'): -3, ('V', 'C'): -1, ('T', 'H'): -2,
              ('P', 'R'): -2, ('V', 'G'): -3, ('T', 'L'): -1, ('V', 'K'): -2,
              ('K', 'Q'): 1, ('R', 'A'): -1, ('I', 'R'): -3, ('T', 'D'): -1,
              ('P', 'F'): -4, ('I', 'N'): -3, ('K', 'I'): -3, ('M', 'D'): -3,
              ('V', 'W'): -3, ('W', 'W'): 11, ('M', 'H'): -2, ('P', 'N'): -2,
              ('K', 'A'): -1, ('M', 'L'): 2, ('K', 'E'): 1, ('Z', 'E'): 4,
              ('X', 'N'): -1, ('Z', 'A'): -1, ('Z', 'M'): -1, ('X', 'F'): -1,
              ('K', 'C'): -3, ('B', 'Q'): 0, ('X', 'B'): -1, ('B', 'M'): -3,
              ('F', 'C'): -2, ('Z', 'Q'): 3, ('X', 'Z'): -1, ('F', 'G'): -3,
              ('B', 'E'): 1, ('X', 'V'): -1, ('F', 'K'): -3, ('B', 'A'): -2,
              ('X', 'R'): -1, ('D', 'D'): 6, ('W', 'G'): -2, ('Z', 'F'): -3,
              ('S', 'Q'): 0, ('W', 'C'): -2, ('W', 'K'): -3, ('H', 'Q'): 0,
              ('L', 'C'): -1, ('W', 'N'): -4, ('S', 'A'): 1, ('L', 'G'): -4,
              ('W', 'S'): -3, ('S', 'E'): 0, ('H', 'E'): 0, ('S', 'I'): -2,
              ('H', 'A'): -2, ('S', 'M'): -1, ('Y', 'L'): -1, ('Y', 'H'): 2,
              ('Y', 'D'): -3, ('E', 'R'): 0, ('X', 'P'): -2, ('G', 'G'): 6,
              ('G', 'C'): -3, ('E', 'N'): 0, ('Y', 'T'): -2, ('Y', 'P'): -3,
              ('T', 'K'): -1, ('A', 'A'): 4, ('P', 'Q'): -1, ('T', 'C'): -1,
              ('V', 'H'): -3, ('T', 'G'): -2, ('I', 'Q'): -3, ('Z', 'T'): -1,
              ('C', 'R'): -3, ('V', 'P'): -2, ('P', 'E'): -1, ('M', 'C'): -1,
              ('K', 'N'): 0, ('I', 'I'): 4, ('P', 'A'): -1, ('M', 'G'): -3,
              ('T', 'S'): 1, ('I', 'E'): -3, ('P', 'M'): -2, ('M', 'K'): -1,
              ('I', 'A'): -1, ('P', 'I'): -3, ('R', 'R'): 5, ('X', 'M'): -1,
              ('L', 'I'): 2, ('X', 'I'): -1, ('Z', 'B'): 1, ('X', 'E'): -1,
              ('Z', 'N'): 0, ('X', 'A'): 0, ('B', 'R'): -1, ('B', 'N'): 3,
              ('F', 'D'): -3, ('X', 'Y'): -1, ('Z', 'R'): 0, ('F', 'H'): -1,
              ('B', 'F'): -3, ('F', 'L'): 0, ('X', 'Q'): -1, ('B', 'B'): 4
              }



def blosum_modified( seq_1, seq_2, matrix_name ='blosum62'):
        #Prepares the amino acid couples. Each couple is a pair of aligned residues.
    pep_division_0 = list(seq_1)
    pep_division_1 = list(seq_2)

    list_coupled_aa = []

    for (aminoacid_0, aminoacid_1) in zip(pep_division_0, pep_division_1):
        list_coupled_aa.append((aminoacid_0, aminoacid_1))

        # Actually computes the scoring matrices values.
    list_value_aa = [] # Contains the scoring matrices values for each couple.


    selected_matrix = blosum62

    for element in list_coupled_aa:
        if element in selected_matrix:
            list_value_aa.append(int(selected_matrix[element]))
        else:
            list_value_aa.append(int(selected_matrix[(element[1], element[0])]))

    dir = os.path.dirname(os.getcwd())

    path = os.path.dirname(dir)

    matrix = open(os.path.join(path,'matrix.txt'), 'r')
    str_matrix = ast.literal_eval(matrix.read().splitlines()[0])



    final_matrix = np.array([int(i) for i in str_matrix])


    #matrix = matrix_generator('3RWG.pdb', 'C', 'A')
    return sum(list_value_aa*final_matrix)

#------------------------------
# Takes the input parameters. -
#------------------------------

input_filename = sys.argv[-1]
i_fh = open(input_filename, "rb")
input_dict = pickle.load(i_fh)
i_fh.close()

peptide_dirpaths = input_dict["peptide_dirpaths"]
n_score_jobs = input_dict["n_score_jobs"]
supress_output = input_dict["supress_output"]

# Suppresses the output.
if supress_output:
    original_stdout = sys.stdout
    sys.stdout = open(os.devnull, "w")

# Initialize MODELLER scorers.
env = environ()
env.libs.topology.read(file='$(LIB)/top_heav.lib')
env.libs.parameters.read(file='$(LIB)/par.lib')

# edat = env.edat
# edat.contact_shell = 8.0

# Initializes SOAP scores.
sp_pp = soap_pp.PairScorer()
sp_atm = soap_pp.AtomScorer()

if use_soap_pep:
    sp_pep = soap_peptide.Scorer()

if supress_output:
    sys.stdout = original_stdout


#-----------------------------------
#  Writing results on a .csv file  -
#-----------------------------------

if tqdm_module:
    bar = tqdm(total=len(peptide_dirpaths), position=1, dynamic_ncols=True, leave=True, unit='file', desc="Progress", bar_format='{l_bar}{bar}| {n_fmt}/{total_fmt} [{elapsed}]')

else:
    print("\nScoring - Please Wait")
for peptide_dirpath_i, peptide_dirpath in enumerate(peptide_dirpaths):

    #print("* Scoring in %s (%s of %s)." % (peptide_dirpath, peptide_dirpath_i+1, len(peptide_dirpaths)))
    original_dirpath = os.getcwd()

    os.chdir(peptide_dirpath)
    df = pd.read_csv('results.csv')

    scores_dict = {# "DOPE": [], # TODO: remove...
                   #"SOAP_PP": [],
                   #"SOAP_ATM": [],
                   #"DOPE_SEL": [],
                   #"DOPE_INT": [],
                   #"KBE_INT": [],
                   #"CB_SCORE_1": [],
                   #"CB_SCORE_2": [],
                   #"CB_SCORE_3": [],
                   #"BLOSUM_MODIFIED": [],
                   #"VOROMQA_TOT": []
                   }

    if use_soap_pep:
        scores_dict["SOAP_PEP"] = []

    # Scoring matrices.
    for matrix_name in matrices_names:
        scores_dict[matrix_name] = []

    # Initializes the columns with the scoring function values.
    for k in sorted(scores_dict.keys()):
        df[k] = None

    decoy_ids = df[df['BUILT'] == 1].index.values
    decoy_names = df.loc[decoy_ids]['DECOY_NAME'].values
    # print("- Decoy names:", decoy_names, decoy_ids)

    for decoy_name in decoy_names:

        #print("- Scoring %s." % decoy_name)

        # Suppresses the output.
        if supress_output:
            original_stdout = sys.stdout
            sys.stdout = open(os.devnull, "w")

        mdl = complete_pdb(env, str(decoy_name))
        atmsel = selection(mdl)
        atmsel_local = selection(mdl.chains['B']).select_sphere(9).by_residue()
        atmsel_receptor = selection(mdl.chains['A'])
        atmsel_peptide = selection(mdl.chains['B'])

        #scores_dict["DOPE_INT"].append(atmsel.assess_dope() - (atmsel_receptor.assess_dope() + atmsel_peptide.assess_dope()))

        # Score with DOPE.
        # scores_dict["DOPE"].append(atmsel.assess_dope())
        #scores_dict["DOPE_SEL"].append(atmsel_local.assess_dope())

        # Score with SOAP-PP.
        #scores_dict["SOAP_PP"].append(atmsel.assess(sp_pp))
        #scores_dict["SOAP_ATM"].append(atmsel.assess(sp_atm))

        # Score with SOAP-PEP.
        if use_soap_pep:
            scores_dict["SOAP_PEP"].append(atmsel.assess(sp_pep))

        # Restores the output.
        if supress_output:
            sys.stdout = original_stdout


        #---------------------------------------------------
        # Compute other knowledge-based scoring functions. -
        #---------------------------------------------------

        #scores_dict["KBE_INT"].append(score_interface_kbe(decoy_name, receptor_chain="A", peptide_chain="B"))

        #---------------------------------------------------
        # Compute CB scoring
        #---------------------------------------------------

        # scores_dict['CB_SCORE'].append(cb_scoring_function(os.path.join(peptide_dirpath, decoy_name), 'A', 'B'))
        #cb_score_1, cb_score_2, cb_score_3 = cb_scoring_function(decoy_name, receptor_chain='A', pep_chain='B')
        #scores_dict['CB_SCORE_1'].append(cb_score_1)
        #scores_dict['CB_SCORE_2'].append(cb_score_2)
        #scores_dict['CB_SCORE_3'].append(cb_score_3)



        #-------------------------------------
        # Computing scoring matrices values. -
        #-------------------------------------

        # Reads the alignment file.
        alignment = open('alignment.ali', 'r')
        alignment_lines = alignment.readlines()

        # Gets the template
        tem_seq = alignment_lines[2].split("/")[1].replace('*', '').rstrip()
        que_seq = alignment_lines[6].split("/")[1].replace('*', '').rstrip()

        for matrix_name in matrices_names:
            scores_dict[matrix_name] = compute_alignment_score(tem_seq, que_seq, matrix_name=matrix_name)

        # Compute BLOSUM_MODIFIED


        #scores_dict['BLOSUM_MODIFIED'].append(blosum_modified(tem_seq, que_seq, matrix_name ='blosum62'))

    #scores_dict["VOROMQA_TOT"] = assess_voromqa(decoy_names, n_jobs=n_score_jobs)
    # scores_dict["GOAP_TOT"] = goap(df['DECOY_NAME'].values)

    for k in sorted(scores_dict.keys()):
        df.at[decoy_ids, k] = scores_dict[k]


    df.to_csv('results.score.csv', index=False, index_label=False)
    #os.chdir(original_dirpath)

    if tqdm_module:
        bar.update()
