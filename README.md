# PepThreader
Protein-peptide interactions (PPIs) are a subset of the overall protein-protein interaction network in the living cell and are pivotal for the majority of cell processes and functions. High-throughput methods to detect such interactions usually require time and costs that are not always affordable. Therefore, reliable in silico predictions represent a valid and effective alternative. 
Here a new algorithm is described, implemented in a freely available tool, i.e., ‘**PepThreader**’, to carry out PPIs predictions and analysis. PepThreader threads multiple fragments derived from a protein sequence (or from a sequence library) onto a second template peptide, in complex with a protein target. The hypothetical binding peptides are then identified and ranked according to a sequence-based and structure-based threading score. The threading algorithm first makes use of a scoring function that is based on peptides sequence similarity. Then, a re-rank of the initial hits is performed, according to structure-based scoring functions. PepThreader has been benchmarked on a dataset of 293 protein-peptide complexes that were collected from existing databases of protein-peptide interactions. An accuracy of 80%, when considering the top predicted 25 hits, was achieved, which performs in a comparable way with the other state-of-art tools in PPI modeling. Therefore, PepThreader adds to the already available tools supporting the experimental PPI identification and characterization. 
 



                REQUIREMENTS


PEPthread has been tested and runs on Windows, macOS and Linux.

The following packages are required: 

 - Biopython 
 - pandas 
 - tqdm

 - MODELLER (*)

(*) Before starting modeling, the presence of 'soap_peptide.hdf5' file (MODELLER SCORER) is checked. 
If the file is not present, it is automatically copied in the proper location to be used by MODELLER.
Because of non-standard installations or denial of writing permission, MODELLER SCORER initialization may fail. In these cases the 'soap_peptide.hdf5' file, which is located in PepThreader main directory must be manually copied in modeller 'modlib' directory to continue with the analysis. 




    DOWNLOAD 
                
                
PepThreader is freely-available at:

https://github.com/cgasb/PepThreader.git



    HOW TO USE PepThreader


PepThreader has a command-line usage. 
Commands and options have to be typed on the command-line as reported next. 
The "pepthreader.py" script must be in the same directory from which it is run. 
The output directory must already exists


COMMAND EXAMPLE:
python pepthreader.py -t ~/THREADER_MODE/PEP_THREADER/DATA/PDBHET_DIR/4QRQ_het.pdb -q ~/THREADER_MODE/PEP_THREADER/DATA/PEP_UNIPROT/P03169.fasta -r A -p C -o ~/THREADER_MODE/PEP_THREADER/test/ --n_jobs 4 -rr 



    OPTIONS:
 
 
-t = path to receptor protein 

-q = path to sequence in FASTA format

-r = to specify the receptor's chain (e.g. A)

-p = to specify the peptide's chain (e.g. C)

-o = output directory (it must already exists) 

-f = force the creation of new directories when generating models

-mx_top = threshold of top-peptides with which carry out the re-rank

-n = Number of decoys for peptide-receptor complex (default = 8)

--n_jobs = number of parallel jobs (default = 1)

-rr = enable re-ranking

-v = verbose

-so = score only. If True, skips modeling phase. In this case models must have been already generated throught PepThreader.  



 CONTACTS
 e-mail: serena.rosignoli@uniroma1.it; alessandro.paiardini@uniroma1.it



 LICENCE
