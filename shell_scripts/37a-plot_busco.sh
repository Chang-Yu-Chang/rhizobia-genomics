
cd
source ~/.zshrc

# Check genome quality
mamba activate busco
mamba env list

folder_temp_result="/Users/cychang/Dropbox/lab/local-adaptation/data/temp/plasmidsaurus"
mkdir -p "$folder_temp_result/summary/07-busco"

# Copy the generic busco
for i in {1..19}
do
    cp "$folder_temp_result/Chang_Q5C_$i/07-busco/results/short_summary.generic.bacteria_odb10.results.txt" \
    "$folder_temp_result/summary/busco/short_summary.generic.bacteria_odb10.g$i.txt"
done


# Copy the specigic busco
for i in {1..19}
do
    cd "$folder_temp_result"/Chang_Q5C_$i/07-busco/results/
    summ_sp=$(ls | grep "specific" | grep 'txt')
    output_txt=$folder_temp_result/summary/busco/$summ_sp

    cp $summ_sp ${output_txt/.results.txt/.g$i.txt}
done

# plot
generate_plot.py -wd "$folder_temp_result/summary/07-busco"
