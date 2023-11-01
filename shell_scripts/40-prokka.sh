cd
source ~/.zshrc


# Check
mamba activate prokka
mamba env list

medaka_consensus=$1
prokka_folder=$2

prokka --force --outdir $prokka_folder --kingdom Bacteria --locustag $medaka_consensus --prefix annotated --gcode 11 $medaka_consensus
# `--force` force overwriting existing output folder
# `--outdir` o utput folder
# `--kingdom`
# `--locustag` locus tax prefix
# `--prefix` filname output prefix
# `--gcode` genetic code / translation table (set if --kingdom is set)
