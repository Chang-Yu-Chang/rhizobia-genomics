# Pangenome analysis and tree inference

- `pangenome.sh` performs the pangenome using panaroo
- `clean_gpa.R` cleans and consolidates the panaroo output tables
- `concatenate_alignment.sh` concatenates the single-copy core genes
- `compute_trees1.sh` uses iqtree to compute single-copy core-gene trees
- `compute_trees2.R` computes trees for gene content variation
- `curate_trees.R` combines the tree objects into a Rdata file `trees.rdata`
- `rf_tree.R` computes the RF distance between trees
- `recA_dnaE2_analysis.sh`
- `recA_dnaE2_analysis.R`
