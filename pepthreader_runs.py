import os
import pandas as pd
import subprocess
import shutil
import argparse
import multiprocessing

from pep_sequencer_lib.scoring_matrices import *

parser = argparse.ArgumentParser()
parser.add_argument("-s", "--scoring", action="store_const", const=True)
#parser.add_argument("-m", "--in_top_value", type=int, default=10)
parser.add_argument("--n_jobs", type=int, default=1)

cmd = parser.parse_args()

df = pd.read_csv('dataset.csv')
decoder_condition = df[df['PEP_SEQ_DECODER'] == 1]

def thread_one_query(row):

    #extract infos from dataset
    pair_id = row[1]['PAIR_ID']
    pdb_id = row[1]['PDB_ID']

    receptor_chain = row[1]['REC_CHAIN']
    pep_chain = row[1]['PEP_CHAIN']
    pep_thread_position = row[1]['THREAD_POSITION']
    benchmark = row[1]['PAIR_LETTER']

    #takes row with opposite benchmark
    opposite_row = df[(df['PAIR_ID']== pair_id) & (df['PAIR_LETTER'] != benchmark)]
    ligand_uniprot = opposite_row['PEP_UNIPROT'].iloc[0]
    opposite_benchmark = opposite_row['PAIR_LETTER'].iloc[0]

    complex_dir = os.path.join(runs_dir, str(pair_id))

    os.chdir(runs_dir)

    if not os.path.exists(complex_dir):
        os.mkdir(complex_dir)
    #create input and output dirs
    #receptor_dir = os.path.join(complex_dir, str(pair_id))
    pdb_filepath = os.path.join(pdb_data_dir, pdb_id + '_het.pdb')
    fasta_filepath = os.path.join(fasta_data_dir, ligand_uniprot + '.fasta')
    output_dir = os.path.join(complex_dir,str(pair_id) + '_' + benchmark )


    os.chdir('/mnt/linux_data/chiara/Tesi/THREADER_MODE/PEP_SEQUENCER/')
    if cmd.scoring:
        #run processes
        print('\n# Scoring ' + str(pair_id) + benchmark + ' ' + ligand_uniprot + ' ' + pdb_id)

        args = ['python', 'pep_sequencer3.py', '-t', str(pdb_filepath), '-q', str(fasta_filepath),
                '-r', str(receptor_chain), '-p', str(pep_chain), '-o' , str(output_dir),'--n_jobs', str(8), '-rr' ]
        if os.path.exists(output_dir):
            print('- forced process')
            args.append("-f")

        subprocess.call(args)


if cmd.n_jobs == 1:
    for row in decoder_condition.iterrows():
        thread_one_query(row)
else:
    pool = multiprocessing.Pool(cmd.n_jobs)
    pool.map(thread_one_query, [row for row in decoder_condition.iterrows()])
    pool.close()
