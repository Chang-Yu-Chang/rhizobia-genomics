#!/usr/bin/env zsh

folder_data="/Users/cychang/Dropbox/lab/local-adaptation/data"
folder_anvio="$folder_data/temp/anvio"

# 6. Index the Reference Genome
# 6.1 For long reads, tools like Minimap2 or NGMLR are often preferred
mamba activate minimap2
mkdir -p "$folder_anvio/06-alignment/sam"
ref_genome="$folder_anvio/06-alignment/reference_genome.mmi"
minimap2 -d $ref_genome "$folder_anvio/01-fasta/em1021.fasta"

# 6.2 Align Each Assembled Genome to the Reference
for i in Chang_Q5C_{1..19} em1021 em1022 wsm419
do
    minimap2 -ax map-ont $ref_genome "$folder_anvio/01-fasta/$i.fasta" > "$folder_anvio/06-alignment/sam/$i.sam"
done

# 6.3 Cover sam to bam
mamba activate samtools
mkdir -p "$folder_anvio/06-alignment/bam"
for i in Chang_Q5C_{1..19} em1021 em1022 wsm419
do
    # Convert SAM to BAM: For each alignment use samtools to convert the SAM file to BAM and sort it.
    samtools view -bS "$folder_anvio/06-alignment/sam/$i.sam" | samtools sort -o "$folder_anvio/06-alignment/bam/$i.bam"
    # Index the BAM file: Index the resulting BAM file
    samtools index "$folder_anvio/06-alignment/bam/$i.bam"
done

# i=Chang_Q5C_1
# samtools mpileup "$folder_anvio/bam/$i.bam"

# 6.4. move hdf file from medaka
mkdir -p "$folder_anvio/06-alignment/hdf"
for i in Chang_Q5C_{1..19}
do
    cp "$folder_data/temp/plasmidsaurus/$i/04-medaka/consensus_probs.hdf" "$folder_anvio/06-alignment/hdf/$i.hdf"
done


# 7. call SNPs variation
folder_data="/Users/cychang/Dropbox/lab/local-adaptation/data"

# 7.1 use snippy for each genome
# For snippy, I dont need to index my fasta. I can use the fasta file directly
mamba activate snippy
mkdir -p "$folder_anvio/06-alignment/snippy"
i=Chang_Q5C_1
for i in Chang_Q5C_{1..19} em1021 em1022 wsm419
do
    snippy \
        --ref "$folder_anvio/01-fasta/em1021.fasta" \
        --ctgs "$folder_anvio/01-fasta/$i.fasta" \
        --outdir "$folder_anvio/06-alignment/snippy/$i"
done

# 7.2 merge thee vcf
mkdir -p "$folder_anvio/06-alignment/snippy/core"
cd "$folder_anvio/06-alignment/snippy"
snippy-core --ref Chang_Q5C_1/ref.fa --prefix core/core  Chang_Q5C_1 Chang_Q5C_10 Chang_Q5C_11 Chang_Q5C_12 Chang_Q5C_13 Chang_Q5C_14 Chang_Q5C_15
#snippy-core --ref Chang_Q5C_1/ref.fa Chang_Q5C_1 Chang_Q5C_10 Chang_Q5C_11 Chang_Q5C_12 Chang_Q5C_13 Chang_Q5C_14 Chang_Q5C_15 Chang_Q5C_16 Chang_Q5C_17 Chang_Q5C_18 Chang_Q5C_19 Chang_Q5C_2 Chang_Q5C_3 Chang_Q5C_4 Chang_Q5C_5 Chang_Q5C_6 Chang_Q5C_7 Chang_Q5C_8 Chang_Q5C_9 em1021 em1022 wsm419

# 7.3 quick plot
cd "$folder_anvio/06-alignment/snippy/core"
snippy-clean_full_aln core.full.aln > clean.full.aln
run_gubbins.py -p gubbins clean.full.aln
snp-sites -c gubbins.filtered_polymorphic_sites.fasta > clean.core.aln
#FastTree -gtr -nt clean.core.aln > clean.core.tree

# 8. Handle vcf files
mkdir -p "$folder_anvio/08-vcf"
cp "$folder_anvio/06-alignment/snippy/core/core.vcf" "$folder_anvio/08-vcf/core.vcf"

# 8.1 Filter variants
mamba activate bcftools
cd "$folder_anvio/08-vcf"
# Apply filters based on variant quality (QUAL, DP, QD, etc.)
bcftools filter -i 'QUAL > 30 && DP > 10 && QD > 2' core.vcf -o filtered.vcf



# 7.2 use snippy for all genomes at once -> the shell script is the same as the previous section
# for i in Chang_Q5C_{1..19} em1021 em1022 wsm419; do; echo "$i\t$folder_anvio/01-fasta/$i.fasta" >> "$folder_anvio/06-alignment/snippy/input.txt"; done
# snippy-multi "$folder_anvio/06-alignment/snippy/input.txt" --ref "$folder_anvio/01-fasta/em1021.fasta" --cpus 16 > runme.sh

#mamba activate snippy
# i=Chang_Q5C_1
# cd "$folder_anvio/06-alignment/snippy/$i"
# snippy-vcf_report --cpus 8 --auto > snps.report.txt



# 7.1 use medaka to call SNPs
#mamba activate medaka
#cd "$folder_anvio/06-alignment/hdf"
#medaka variant "$folder_anvio/01-fasta/em1021.fasta" $folder_anvio/06-alignment/hdf/Chang_Q5C_2.hdf  "$folder_anvio/06-alignment/medaka.vcf"
#Chang_Q5C_2,Chang_Q5C_3,Chang_Q5C_4,Chang_Q5C_5,Chang_Q5C_6,Chang_Q5C_8,Chang_Q5C_9,Chang_Q5C_10,Chang_Q5C_11,Chang_Q5C_13,Chang_Q5C_15,Chang_Q5C_16,Chang_Q5C_17,Chang_Q5C_19 \

# 7.1 Profile a BAM file. For EACH SAMPLE
# mamba activate anvio-8
# mkdir -p "$folder_anvio/07-single_profiles"
# mkdir -p "$folder_anvio/07-single_profiles/profile_db"
# i=Chang_Q5C_1
# for i in Chang_Q5C_{1..19} em1021 em1022 wsm419
# do
#     anvi-profile \
#         -i "$folder_anvio/06-alignment/bam/$i.bam" \
#         -c "$folder_anvio/03-contigs_db/$i.db" \
#         -o "$folder_anvio/07-single_profiles/$i" \
#         --sample-name $i \
#         --overwrite-output-destinations
#     cp "$folder_anvio/07-single_profiles/$i/PROFILE.db" "$folder_anvio/07-single_profiles/profile_db/$i.db"
# done
#
# # 7.2 Merge all anvi’o profiles using the program anvi-merge
# cd "$folder_anvio/07-single_profiles/profile_db"
# anvi-merge \
#     -c "$folder_anvio/03-contigs_db/Chang_Q5C_1.db" \
#     Chang_Q5C_1.db Chang_Q5C_2.db \
#     -o "$folder_anvio/07-single_profiles/merged"
#
# anvi-merge "$folder_anvio/07-single_profiles/*/PROFILE.db" -o SAMPLES-MERGED -c contigs.db
#
# anvi-gen-variability-profile \
#     -p profile-db \
#     -c contigs-db \
#     -C DEFAULT \
#     -b EVERYTHING






