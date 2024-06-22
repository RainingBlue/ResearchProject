#!/bin/bash

#SBATCH --qos=short                 # Request Quality of Service. Default is 'short' (maximum run time: 4 hours)
#SBATCH --time=0:15:00             # Request run time (wall-clock). Default is 1 minute
#SBATCH --job-name=gc_deep       
#SBATCH --ntasks=1                   # 1 task per job in array
#SBATCH --cpus-per-task=1            # Specify, otherwise will run into issues
#SBATCH --mem=100MB                  # Memory needed per task
#SBATCH --array=1-68                # Array with 10 tasks
#SBATCH --output=out/slurm-%A_%a.out # Set name of output log. %A is SLURM_ARRAY_JOB_ID and %a is SLURM_ARRAY_TASK_ID
#SBATCH --error=out/slurm-%A_%a.err  # Set name of error log. %A is SLURM_ARRAY_JOB_ID and %a is SLURM_ARRAY_TASK_ID

config=control_dt.config          # Path to config file

path=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $2}' $config)
filename=$(echo $path | sed 's:.*/::')
echo $path
echo $filename
eval "$($HOME/miniforge3/bin/conda shell.bash hook)"
conda activate env
python get_gc_bias.py out/${filename}_deeptools.bam 100Kb_sorted.bed
