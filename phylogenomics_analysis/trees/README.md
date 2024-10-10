# Maximum likelihood tree

0. `concatenate_alignment.sh` concatenate the single copy core gene alignment into one alignment file
1. `iqtree_core.sh` builds ML tree based on multiple sequence alignment of core genes
2. `binarize_gpa.R` converts the gene presence absence table to a binary file
3. `iqtree_gpa.sh` builds ML tree based on gene presence absence pattern (binary matrix)
4. `iqtree_syn.sh` builds ML tree based on synteny blocks presecen absence (binary matrix)
