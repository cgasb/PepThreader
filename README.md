# PepThreader
A tool for predicting protein-peptide interactions with a template-based approach. 



                REQUIREMENTS

PEPthread has been tested and runs on Windows, macOS and Linux.

The following packages are required: 

 - Biopython 
 - pandas 
 - tqdm

 - MODELLER
   Before starting modelization, the presence of 'soap_peptide.hdf5' file (MODELLER SCORER) is checked. 
   If the file is not present, it is automatically copied in the proper location to be used by MODELLER.
   
   Because of non-standard installations or denial of writing permission, MODELLER SCORER initialization may fail. 
   In these cases the 'soap_peptide.hdf5' file, which is located in PepThreader main directory must be manually copied in modeller 'modlib' directory to continue with the analysis. 



                DOWNLOAD 
                  
PepThreader is freely-available at:

https://github.com/cgasb/PepThreader.git



              HOW TO USE PEPthreader

PEPthread has a command-line usage. Commands and options have to be typed on the command-line as reported next. The "pepthreader.py" script must be in the same directory from which it is run. 


COMMAND EXAMPLE:
python pepthreader.py -t ~/THREADER_MODE/PEP_THREADER/DATA/PDBHET_DIR/4QRQ_het.pdb -q ~/THREADER_MODE/PEP_THREADER/DATA/PEP_UNIPROT/P03169.fasta -r A -p C -o ~/THREADER_MODE/PEP_THREADER/test/ --n_jobs 4 -rr 


OPTIONS:
 
-t = path to receptor protein 
-q = path to sequence in FASTA format
-r = to specify the receptor's chain 
-p = to specify the peptide's chain 
-o = output directory (it must already exists) 
-f = force the creation of new directories when generating models
-mx_top = threshold of top-peptides with which carry out the re-rank
-n = Number of decoys for peptide-receptor complex (default = 8)
--n_jobs = number of parallel jobs (default = 1)
-rr = enable re-ranking
-v = verbose
-so = score only. If True, skips modeling phase. In this case models must have been already generated throught PEPthreader.  


 CONTACTS
 
 e-mail: 


 LICENCE
