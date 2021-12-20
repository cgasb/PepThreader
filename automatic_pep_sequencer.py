import os
import pandas as pd
import subprocess
import shutil
import argparse
import multiprocessing

from pep_sequencer_lib.scoring_matrices import *


parser = argparse.ArgumentParser()
parser.add_argument("-s", "--scoring", action="store_const", const=True)
parser.add_argument("-m", "--in_top_value", type=int, default=10)
parser.add_argument("--n_jobs", type=int, default=1)

cmd = parser.parse_args()


#----------------------------
# 1/2 Calculate the scores. -
#----------------------------

df = pd.read_csv('dataset.csv')
decoder_condition = df[df['PEP_SEQ_DECODER'] == 1]

'''

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


    os.chdir('/home/chiara/Scrivania/Tesi//PEP_THREADER_ULTIMATE/PEP_SEQUENCER')
    if cmd.scoring:
        #run processes
        print('\n# Scoring ' + str(pair_id) + benchmark + ' ' + ligand_uniprot + ' ' + pdb_id)

        args = ['python', 'pep_sequencer2.py', '-t', str(pdb_filepath), '-q', str(fasta_filepath),
                '-r', str(receptor_chain), '-p', str(pep_chain), '-o' , str(output_dir)]
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
'''

#------------------------
# 2/2 Analyse the data. -
#-----------------------

#######################################################################
# Extract Results files
#######################################################################

matrices = []

for matrix in matrices_names:
    matrices.append(matrix)

    matrices.append(matrix + '_MODIFIED')
    matrices.append(matrix + '_GALAXY')

analysis = True
scores_dir = '/home/chiara/Scrivania/Tesi/PEP_THREADER_ULTIMATE/PEP_SEQUENCER/scores_results'

if not os.path.exists(scores_dir):
    os.mkdir(scores_dir)


#extract scores
if analysis:
    for dir in os.listdir(runs_dir):
        complex_dir = os.path.join(runs_dir, dir)
        for benchmark_dir in os.listdir(complex_dir):

        #if dir == 'scores_results':
            #continue
            out_dirpath = os.path.join(complex_dir, benchmark_dir)
            scores_filepath = os.path.join(out_dirpath, 'Scoring_results.csv')
            renamed_scores_filepath = os.path.join(out_dirpath, benchmark_dir + '.csv')


            if os.path.exists(scores_filepath):

                os.rename(scores_filepath, renamed_scores_filepath)

            scores_finalpath = os.path.join(scores_dir, benchmark_dir + '.csv')

            shutil.copy(renamed_scores_filepath, scores_dir)


    ###########################################################################
    # Analyse Data
    ###########################################################################
    os.chdir(scores_dir)
    dict_scoring = {"COMPLEX_NAME": []} #"CUSTOM": []
    custom_scoring_count = 0
    file_counter = 0


    for file in os.listdir(scores_dir):



        #not counting documents not necessary to the analysis
        if not ("_p" in file or "_q" in file) or "lock" in file:
            continue
        dict_scoring["COMPLEX_NAME"].append(file[:-4])
        file_counter += 1

        #read the results file
        protein_df = pd.read_csv(file)
        #extract useful infos for analysis
        pair_id = file.split('_')[0]


        benchmark = file[-5]


        protein_data = df[(df['PAIR_ID'] == int(pair_id)) & (df['PAIR_LETTER'] == benchmark)]
        position = np.array(protein_data['THREAD_POSITION']).astype(int)

        #display results name and right position
        #rint('File name is:', os.path.basename(file))
        #print('True position is =', int(position))

        #add to dict scoring rank value of true position

        for column in matrices:

            values = protein_df[column]


            rank_value = values.rank(ascending = False)[position-1]



            if column in dict_scoring:
                dict_scoring[column].append(int(rank_value))
            else:
                dict_scoring[column] = [int(rank_value)]



    #add dict to DF and save csv file
    sdf = pd.DataFrame(dict_scoring)
    sdf.to_csv(os.path.join(scores_dir,'data.csv'))


    #create dict with n of true position with rank under m
    dict_val = {}

    for column in sdf.columns:
        if column == "COMPLEX_NAME":
            continue

        score = sum(sdf[column].values <= cmd.in_top_value)
        dict_val[column] = score

    print('total cases are = ', file_counter)
    #print ordered matrices with as key n of correct prediction + probability over all calculations
    for i in sorted(dict_val.items(), key=lambda x: x[1]):
        print(i, i[1]/file_counter)
