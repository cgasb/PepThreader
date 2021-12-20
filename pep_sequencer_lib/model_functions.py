import os
import time
import subprocess
import pickle

voronota_filepath = None

def generate_ali_file(model_directory, complex_name, model_name,  template_sequence, template_peptide_sequence, peptide_sequence, receptor_chain, peptide_chain):

        aln_dir = open(os.path.join(model_directory, "alignment.ali"), "w")
        structure_info = "structure:" + complex_name

        tem_string = ">P1;" + complex_name + "\n" + structure_info + ":FIRST:" + receptor_chain + ":END:" + peptide_chain + "::::" + "\n" + template_sequence + "/" + template_peptide_sequence + "*"
        mod_string = ">P1;" + model_name + "\n" + "sequence" + ":" + model_name + ":FIRST:.:LAST:.::::" + "\n" + template_sequence + "/" + str(peptide_sequence) + "*"
        aln_dir.write(tem_string + "\n" + "\n" + mod_string)
        aln_dir.close()

def launch_external_modeling(peptide_sequence, peptide_directory, peptide_name):

    t1 = time.time()

    input_filename = peptide_name + ".pkl"
    i_fh = open(input_filename, "wb")
    pickle.dump({"peptide_directory": peptide_directory,
                 "peptide_name": peptide_name,
                 "peptide_sequence": peptide_sequence,

                 "template_name": 'template_complex',
                 "supress_output": not cmd.verbose,
                 "n_decoys": cmd.n_decoys,
                 "n_mod_jobs": cmd.n_mod_jobs,
                 "pt_modeling": cmd.pt_modeling,}, i_fh)
    i_fh.close()
    #subprocess.check_call([python_exec, build_models_script_filepath, input_filename])
    os.remove(input_filename)
    print("- It took %s." % (time.time()-t1))

def launch_external_modeling_parallel(peptide_tuple):
    print("function of model_functions.py")
    launch_external_modeling(peptide_tuple[1][0], peptide_tuple[1][1], peptide_tuple[1][2], peptide_tuple[0])

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

def launch_external_scoring(peptide_dirnames, peptides_count):

    t1 = time.time()

    print("\n# Scoring %s models of %s (batch %s)." % (len(peptide_dirnames), len(dirnames_list), peptides_count))

    input_filename = "score_batch_%s.pkl" % peptides_count
    i_fh = open(input_filename, "wb")
    pickle.dump({"peptide_dirpaths": [os.path.join(cmd.results_dirpath, peptide_dirname) for peptide_dirname in peptide_dirnames],
                 "n_score_jobs": cmd.n_score_jobs,
                 "supress_output": not cmd.verbose,}, i_fh)
    i_fh.close()
    subprocess.check_call([python_exec, score_models_script_filepath, input_filename])
    os.remove(input_filename)
    print("- It took %s." % (time.time()-t1))

def launch_external_scoring_parallel(peptide_tuple):
    launch_external_scoring(peptide_tuple[1], peptide_tuple[0])
