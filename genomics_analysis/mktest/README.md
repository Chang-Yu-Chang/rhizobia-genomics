Follow the steps below to do the MK test

Requirements
- Multifasta alignment for each gene. This is the core gene alignment from Panaroo
- Outgroup genomes. MK test requires a output for each test. Use S. meliloti EM1021 and S. medicae WSM419
- The outgroup genomes have to be annotated (using Prokka) in order to get the gene sequences

1. `prepare_msa.sh` prepares the MSA file with a outgroup sequence
2. `count_nt.sh` uses the python script `sfs_from_fasta_2.py` to count the contingency table for each gene
3. `mktest.R` takes the div and daf tables to perform MK test. Output an aggregated MK test result csv `mktests.csv`
4. `plot_mktest.R` plots the results
