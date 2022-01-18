import os
import sys
import pickle

import pandas as pd

from modeller import *
from modeller.automodel import *
from modeller.parallel import *

from pep_sequencer_lib.pt_automodel import PT_automodel


DEVELOP = False


def launch_modeling(peptide_directory, peptide_name, peptide_sequence, template_name, n_decoys=2, n_jobs=1, pt_modeling=False, supress_output=True):
    """
    This function calls MODELLER and build a series of models for a certain peptide-receptor complex.
    """

    original_dirpath = os.getcwd()
    os.chdir(peptide_directory)

    # Suppressing the MODELLER output.
    if supress_output:
        original_stdout = sys.stdout
        sys.stdout = open("log.txt", "w") # open(os.devnull, "w")

    log.verbose()
    env = environ()

    # env.io.atom_files_directory = []
    # env.io.atom_files_directory.append('.')
    # env.io.hetatm = True

    # Use multiple CPUs in a parallel job on this machine.
    if n_jobs > 1:
        j = job()
        for i in range(0, n_jobs):
            j.append(local_slave())

    if pt_modeling:
        automodel_class = PT_automodel
    else:
        automodel_class = automodel

    a = automodel_class(env, alnfile='alignment.ali',
                        knowns=template_name, sequence=peptide_name,)
    a.starting_model= 1
    a.ending_model = n_decoys

    # Do not perform molecular dynamics to speed up (testing only).
    if DEVELOP:
        a.md_level = None

    if n_jobs > 1:
        a.use_parallel_job(j)

    a.make()

    if supress_output:
        sys.stdout = original_stdout


    # Write a results file for this receptor-peptide complex.
    model_results = []
    for single_model in a.outputs:
        model_dict = {"PEPTIDE_NAME": peptide_name,
                      "DECOY_NAME": single_model["name"],
                      "DECOY_NUM": single_model["num"],
                      "BUILT": int(single_model["failure"] == None),
                      "PEPTIDE_SEQ": peptide_sequence,}

        if single_model["failure"] == None:
            model_dict.update({# Total objective function of MODELLER.
                               "MOLPDF": single_model["molpdf"],
                               # Side chain dihedral angles terms of the objective function.
                               "SC_DH_OBJ": single_model["pdfterms"][physical.chi1_dihedral]+single_model["pdfterms"][physical.chi2_dihedral]+single_model["pdfterms"][physical.chi3_dihedral]+single_model["pdfterms"][physical.chi4_dihedral],
                               # Main chain dihedral angles terms of the objective function.
                               "MC_DH_OBJ": single_model["pdfterms"][physical.phi_psi_dihedral],
                               # Steric repulsion terms of the objective function.
                               "SOFT_SPHERE_OBJ": single_model["pdfterms"][physical.soft_sphere]
                               })

        model_results.append(model_dict)


    model_df = pd.DataFrame(model_results)
    model_df.to_csv("results.csv", index=False, index_label=False)


# Get the input parameters from a pickle file.
input_filename = sys.argv[-1]
i_fh = open(input_filename, "rb")
input_dict = pickle.load(i_fh)
i_fh.close()

peptide_directory = input_dict["peptide_directory"]
peptide_name = input_dict["peptide_name"]
peptide_sequence = input_dict["peptide_sequence"]

template_name = input_dict["template_name"]
supress_output = input_dict["supress_output"]
n_decoys = input_dict["n_decoys"]
n_mod_jobs = input_dict["n_mod_jobs"]
pt_modeling = input_dict["pt_modeling"]


launch_modeling(peptide_directory=peptide_directory,
                peptide_name=peptide_name,
                peptide_sequence=peptide_sequence,
                template_name='template_complex',
                n_decoys=n_decoys, n_jobs=n_mod_jobs,
                pt_modeling=pt_modeling,
                supress_output=supress_output)
