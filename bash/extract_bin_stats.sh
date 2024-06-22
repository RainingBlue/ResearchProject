# #!/bin/bash

filename=$(basename "$1")
filename_proc="${filename%.*}_proc.bam"

echo "processing ${filename}"
# samtools view $1 chr{1..22} -F 260 -h -b > out/$filename
# echo "Extracted chromosome 1-22, filtered on mapq 30"
# samtools index out/$filename
bedtools intersect -abam $1 -b dukeExcludeRegions.bed -v > "$2/$filename_proc"
echo "Removed blacklisted regions"
samtools index "$2/${filename_proc}"
echo "Created temp index"
python collect_bin_stats.py "$2/${filename_proc}" $1 100Kb_sorted.bed
echo "Done collecting bin counts"
python LoessRegression.py "$2/${filename_proc}"
echo "Done correcting GC"
python scale_counts.py $1 5Mb_sorted.bed "$2/${filename_proc}_counts_gc_corrected.npy"
echo "Done scaling counts based on GC"
python ratio.py "$2/${filename_proc}" "$2/${filename_proc}_stats.npy"
python analyze_ratios.py median_profile.npy "$2/${filename_proc}_ratio.npy" "$2/${filename}_deeptools.bam_stats.npy"
# echo "Done"
# rm "$2/$filename_proc"
# rm "$2/${filename_proc}.bai"
# rm out/$filename
# rm "out/${filename}.bai"