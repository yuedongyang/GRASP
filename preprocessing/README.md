# TH-RSS: Data Preprocessing.

Yaobin Ke, Jiahua Rao 2019-01-15

## Introduction
Labeling and Encoding of RNA sequences.

## Input data  
Three datasets are selected for training and testing:  

### Raw data:
* [PARS-Yeast (Kertesz, et al., 2010)](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE22393)

* [PARS-Human (Rouskin, et al., 2014)](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE50676)  

* [SS_PDB (Andronescu, et al., 2008)]()
   
 
*Two of datasets are from PARS method, including PARS-Yeast and PARS-Human. Another dataset, SS_PDB, containing structures experimentally validated with NMR and X-ray, was downloaded from RNAstrand, a wide-used database.*
  
### File Formats:  
* fasta: containing gene_name, RNA sequences.  
    * example of fasta file:

    gene_name|sequences
    ---|---
    YKL152C|ATATTACAATAATGCC......GTAGCTATATATAGTCAA
* tab: containing gene_name, experimental scores for each nucleotide on RNA sequences.
    * example of tab file:  

    gene_name|sequence_length|experimental scores
    ---|---|---
    YKL152C|872|2.18; 0.03; -0.94; ...; 0; 0;
* ss: containing gene_name, sequences, position of corresponding paired nucleotides
    * example of ss file:  

    gene_name|sequences|position
    ---|---|---
    PDB_00001|CGCGAAUUAGCGCGCGAAUUAGCG|24 23 22 21 20 19 18 17 16 15 14 13 12 11 10 9 8 7 6 5 4 3 2 1

*Two of raw datasets have been loaded, including PARS-Yeast and SS_PDB. However, the file of PARS-Human is too large to upload, and you can find it in [here](https://pan.baidu.com/s/1CsqmX-zz0wMuHz7B-cGrhA "PARS-Human").*


## Output  
labeling and encoding RNA sequences, so that it can be input to Machine Learning Algorithms.

### output file:
* nocode.csv: contain gene_name, sequences, position, label
    * example of nocode.csv: such as py_nocode_37.csv  (*37 is window_size*)

    gene_name|f1|...|f37|label|position
    ---|---|---|---|---|---
    YKL152C|T| ... |C|1|99


* encode.csv: containing gene_name, sequences, encoded features, label
    * example of encode.csv: such as py_encode_37.csv

    gene_name|sequences|f1|...|f37|label
    ---|---|---|---|---|---
    YKL152C|TTGATGTTAAATTGTCTGCCAAGGGTCAACAAGAAGC|0 1 0 0| ... |0 0 1 0|1


## Usage
* **PARS-Yeast**:  
    Command line: 
    ```
    usage: get_input_py.py [-h] [--input_fasta INPUT_FASTA]
                        [--input_score INPUT_SCORE] [--window_size WINDOW_SIZE]
                        [--dataset DATASET]

    optional arguments:
    -h, --help            show this help message and exit
    --input_fasta INPUT_FASTA
                            FASTA format file, containing RNA sequences.
    --input_score INPUT_SCORE
                            The experimental scores of RNA structural profile.
    --window_size WINDOW_SIZE
                            The window size when truncating RNA sequences.
    --dataset DATASET     
                            The type of dataset: PARS_human:ph, PARS_yeast:py or NMR_X-Ray:pdb
    ```
    To preprocess the data, you can run:  
    > python3 get_input_py.py    

    Also, you can run:  
    > python3 get_input_py.py --input_fasta='./raw_data/PARS_yeast/sce_genes.fasta' --input_score='./raw_data/PARS_yeast/sce_Score.tab' --window_size=37 --dataset='py'

* **PARS-Human**:  
    Command line: 
    ```
    usage: get_input_ph.py [-h] [--input_fasta INPUT_FASTA]
                       [--input_score INPUT_SCORE] [--window_size WINDOW_SIZE]
                       [--dataset DATASET]

    optional arguments:
    -h, --help            show this help message and exit
    --input_fasta INPUT_FASTA
                        FASTA format file, containing RNA sequences.
    --input_score INPUT_SCORE
                        The experimental scores of RNA structural profile.
    --window_size WINDOW_SIZE
                        The window size when truncating RNA sequences.
    --dataset DATASET     The type of dataset: PARS_human:ph, PARS_yeast:py or NMR_X-Ray:pdb
    ```
    To preprocess the data, you can run:  
    > python3 get_input_ph .py    

    Also, you can run:  
    > python3 get_input_py.py --input_fasta='./raw_data/PARS_human/pars_h_gene_haveseq.fasta' --input_score='./raw_data/PARS_human/GM12878_score_no+5score_haveSeq.tab' --window_size=37 --dataset='ph'

* **SS_PDB**:  
    Command line: 
    ```
    usage: get_input_pdb.py [-h] [--input_fasta INPUT_FASTA] [--input_ss INPUT_SS]
                            [--window_size WINDOW_SIZE] [--dataset DATASET]

    optional arguments:
    -h, --help            show this help message and exit
    --input_fasta INPUT_FASTA
                            FASTA format file, containing RNA sequences.
    --input_ss INPUT_SS   The RNA sequences structural profile.
    --window_size WINDOW_SIZE
                            The window size when truncating RNA sequences.
    --dataset DATASET     The type of dataset: PARS_human:ph, PARS_yeast:py or
                            NMR_X-Ray:pdb
    ```
    To preprocess the data, you can run:  
    > python3 get_input_pdb.py    

    Also, you can run:  
    > python3 get_input_pdb.py --input_fasta='./raw_data/NMR_X-Ray/new75.fa' --input_ss='./raw_data/NMR_X-Ray/pdball.SS' --window_size=37 --dataset='pdb'
