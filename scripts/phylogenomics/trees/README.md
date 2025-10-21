Compute genome wide tree

2. `compute_trees1.sh` uses iqtree to compute genome-level single copy core gene trees
3. `compute_trees2.R` computes trees for gpa, ani, and kmers
4. `curate_trees.R` combines the tree objects into a Rdata file `trees.rdata`

To reproduce the intermediate files without invoking the exploratory plots

```
zsh compute_trees1.sh
cd ../../
Rscript -e "renv::activate('.'); source('phylogenomics_analysis/trees/compute_trees2.R'); source('phylogenomics_analysis/trees/curate_trees.R')"
```
