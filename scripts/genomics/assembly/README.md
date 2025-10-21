# De novo assembly

1. `assess_reads.sh` performs quality control on the raw nanopore long reads. Also filter out bad quality reads
2. `denovo_assembly.sh` assembles the genomes from raw reads following the plasmidsaurus pipeline
3. `assess_assemblies.sh` performs quality control on the assemblies
4. `manual_concat.sh` concatenates the contigs based on the blast result
5. `aggregate_qc.R` aggregates the qc results of assemblies


### DEPRECATED
1. `raw_reads`: examine raw reads
4. `taxonomy`: aggregates the blast results
2. `genomes` and `contigs`: examine assembled genome size and aggregate contig data
3. `gene_content`: examine gene content on panaroo outputs
5. `fst`, `gcv_fst`: compute fst for SNPs and GCV (gene content variation)
6. `dxy`, `gcv_dxy`: compute dxy for SNPs and GCV
7. `go`, `gcv_go`: perform GO enrichment analysis for SNPs and GCV
