# De novo assembly

1. `assess_reads.sh` performs quality control on the raw nanopore long reads
2. `denovo_assembly.sh` assembles the genomes from raw reads following the plasmidsaurus pipeline
3. `assess_assemblies.sh` performs quality control on the raw nanopore long reads
4. `consolidate_genomes.sh` consolidates genome fasta files into one folder
5. `extract_contigs.sh` splits each genome into contigs. This includes 32 genomes of our isolates and 5 NCBI reference genomes 
