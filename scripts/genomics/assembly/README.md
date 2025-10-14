# De novo assembly

1. `assess_reads.sh` performs quality control on the raw nanopore long reads. Also filter out bad quality reads
2. `denovo_assembly.sh` assembles the genomes from raw reads following the plasmidsaurus pipeline
3. `assess_assemblies.sh` performs quality control on the raw nanopore long reads
4. `manual_concat.sh` concatenates the contigs based on the blast result

To reproduce it, do the following

```
zsh assess_reads.sh
zsh denovo_assembly.sh
zsh assess_assemblies.sh
zsh manual_concat.sh
```

# deprecated
6. `extract_contigs.sh` splits each genome into contigs. This includes 32 genomes of our isolates and 5 NCBI reference genomes 
