Compute replicon-level tree

1. `concatenate_alignment.sh` concatenates the single-copy core-gene alignment into one alignment file
    1a. `concatenate_alignment.py` is a command line tool that concatenates alignments
    1b. `curate_replicon_sccg.R` is a command line tool that curates the list of single copy core gene for each replicon
2. `compute_trees1.sh` uses iqtree to compute genome-level single copy core gene trees
3. `compute_trees2.R` computes trees for gpa, ani, and kmers
4. `curate_trees.R` combines the tree objects into a Rdata file `trees.rdata`

To reproduce the intermediate files without invoking the exploratory plots

```
zsh concatenate_alignment.sh
zsh compute_trees1.sh
cd ../../
Rscript -e "renv::activate('.'); source('phylogenomics_analysis/replicon_trees/compute_trees2.R'); source('phylogenomics_analysis/replicon_trees/curate_trees.R')"
```
