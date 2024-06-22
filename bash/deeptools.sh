#!/bin/bash

#SBATCH --qos=short                 # Request Quality of Service. Default is 'short' (maximum run time: 4 hours)
#SBATCH --time=1:30:00              # Request run time (wall-clock). Default is 1 minute
#SBATCH --job-name=deeptools       
#SBATCH --ntasks=1                   # 1 task per job in array
#SBATCH --cpus-per-task=4            # Specify, otherwise will run into issues
#SBATCH --mem=1200MB                  # Memory needed per task
#SBATCH --array=1-75              # Array with 10 tasks
#SBATCH --output=out/slurm-%A_%a.out # Set name of output log. %A is SLURM_ARRAY_JOB_ID and %a is SLURM_ARRAY_TASK_ID
#SBATCH --error=out/slurm-%A_%a.err  # Set name of error log. %A is SLURM_ARRAY_JOB_ID and %a is SLURM_ARRAY_TASK_ID

config=luad_dt.config          # Path to config file

path=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $2}' $config)
filename=$(echo $path | sed 's:.*/::')
echo $path
echo $filename
eval "$($HOME/miniforge3/bin/conda shell.bash hook)"
conda activate env
# python deepTools/deeptools/computeGCBias.py -b $path --effectiveGenomeSize 2913022398 -g hg38.2bit -o out/${filename}_freq.txt -bl dukeExcludeRegions.bed -p 4
python deepTools/deeptools/correctGCBias.py -b $path --effectiveGenomeSize 2913022398 -g hg38.2bit -p 4 -freq out/${filename}_freq.txt -o out/${filename}_deeptools.bam