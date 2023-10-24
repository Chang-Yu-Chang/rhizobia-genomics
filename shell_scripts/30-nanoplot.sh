cd
source ~/.zshrc

conda activate
mamba activate nanoplot
mamba env list

raw_reads=$1
nanoplot_folder=$2

NanoPlot -t 2 --verbose  --outdir $nanoplot_folder --fastq $raw_reads --maxlength 100000 --plots kde -c crimson  -cm RdBu --N50
# `-t` thread
# `--verbose` print  Write log messages also to terminal
# `--outdir` output directory
# `--fastq` data is in one or more fastq files
# `--maxlength`  Hide reads longer than length specified
# `--plots [{kde,hex,dot} ...]` Specify which bivariate plots have to be made.
# `-c`  Specify a valid matplotlib color for the plots
# `-cm`  Specify a valid matplotlib colormap for the heatmap
# `--N50` Show the N50 mark in the read length histogram
