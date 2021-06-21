# master_thesis
This repository contains all the scripts and functions used in Blas Chaves-Urbano master's thesis made in the Macintyre Lab at CNIO

This repository is divided in different folder which contains the scripts and files neccesary to reproduce undertaken in this master thesis:

- get_results: This folder contains the code for getting the copy-number signatures results and figures (get_results.R). This script is compiled from BRITOC CNsignatures which contains all the code necessary to reproduce the analyses for the accompaning manuscript "Copy-number signatures and mutational processes in ovarian carcinoma" (https://www.biorxiv.org/content/early/2017/09/04/174201). Details on how to compile the document can be found in the repository README: https://bitbucket.org/britroc/cnsignatures. Much of the code for the signature analysis can be found in the main_functions.R and helper_functions.R files in the base directory of the previous repository. This folder also contains the file cdks_30_kb_14_30kb_ds_absCopyNumber.rds. This file has the absolute copy number profiles of the different experimental cell lines generated using a modified implementation of the QDNAseq package in R designed by Dr. Geoff Macintyre (unplished).
- 
