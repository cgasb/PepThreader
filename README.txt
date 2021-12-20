                 PEPthread

--------introduction ----------------

                REQUIREMENTS

PEPthread has been tested and runs on macOS and Linux.

The following packages are required: 

 - Biopython 
 - pandas 
 - tqdm

 - MODELLER
   Before starting modelization, the presence of 'soap_peptide.hdf5' file is checked. 
   If the file is not present, it is copied in the proper directory to be used by MODELLER.

                  DOWNLOAD 

--------- link per download ----

                 INSTALLATION



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
  

############################################################################################4QRQ_het.pdb
