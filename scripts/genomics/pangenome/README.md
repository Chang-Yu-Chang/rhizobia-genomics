# Pangenome analysis and tree inference

- `pangenome.sh` performs the pangenome using panaroo
- `clean_gpa.R` clean and consolidate the panaroo output tables
- `concatenate_alignment.sh` concatenates the single-copy core genes
- `compute_trees1.sh` uses iqtree to compute genome-level single copy core gene trees
- `compute_trees2.R` computes trees for gene content variation
- `curate_trees.R` combines the tree objects into a Rdata file `trees.rdata`

To reproduce the intermediate files without invoking the exploratory plots

```
zsh compute_trees1.sh
cd ../../
Rscript -e "renv::activate('.'); source('phylogenomics_analysis/trees/compute_trees2.R'); source('phylogenomics_analysis/trees/curate_trees.R')"
```
