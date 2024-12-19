Compute genome wide tree

1. `concatenate_alignment.sh` concatenates the single-copy core-gene alignment into one alignment file
2. `concatenate_alignment.py` is the command line tool for 1.
3. `compute_trees1.sh` uses iqtree to compute genome-level single copy core gene trees
4. `compute_trees2.R` computes trees for gpa
5. `curate_trees.R` combines the tree objects into a Rdata file `trees.rdata`
