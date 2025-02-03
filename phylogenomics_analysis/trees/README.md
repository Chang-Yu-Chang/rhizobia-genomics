Compute genome wide tree

1. `concatenate_alignment.sh` concatenates the single-copy core-gene alignment into one alignment file. It also copies the aligment file to the folder tree/
    1a. `concatenate_alignment.py` is the command line tool for 1.
2. `compute_trees1.sh` uses iqtree to compute genome-level single copy core gene trees
3. `compute_trees2.R` computes trees for gpa, ani, and kmers
4. `curate_trees.R` combines the tree objects into a Rdata file `trees.rdata`
5. `plot.trees.R` plots trees

To reproduce the intermediate files without invoking the exploratory plots

```
zsh concatenate_alignment.sh
zsh compute_trees1.sh
cd ../../
Rscript -e "renv::activate('.'); source('phylogenomics_analysis/trees/compute_trees2.R'); source('phylogenomics_analysis/trees/curate_trees.R')"
```
