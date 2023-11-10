# THis classifies the contigs into plasmids

mamba create -y -n plasforest python=3.8 # scikit-learn==0.22.2 requires python <= 3.8
mamba activate plasforest


folder_data="/Users/cychang/Dropbox/lab/local-adaptation/data"
plasforest="/Users/cychang/bioinformatics/plasforest/PlasForest.py"

python3 $plasforest -b -f \
    -i "$folder_data.fasta" \
    -o "$folder_data.outputfile.csv"


