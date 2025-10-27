# De novo assembly

1. `assess_reads.sh` performs quality control on the raw nanopore long reads. Also filter out bad quality reads
2. `denovo_assembly.sh` assembles the genomes from raw reads following the plasmidsaurus pipeline
3. `assess_assemblies.sh` performs quality control on the assemblies
4. `manual_concat.sh` concatenates the contigs based on the blast result
5. `consolidate_qcs.R` aggregates the qc results of assemblies (flye, quast, checkm, busco)
