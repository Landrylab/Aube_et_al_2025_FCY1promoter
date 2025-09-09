#!/bin/bash
#SBATCH -D /project/chlandry/users/siaub8
#SBATCH -J FCY1_150_CDS_array
#SBATCH -o Logs_FCY1_25-04-07/FCY1_CDS_150_PD-%j.log
#SBATCH -c 1
#SBATCH --mail-type=ALL
#SBATCH --mail-user=simon.aube.2@ulaval.ca
#SBATCH --array=6,7,8,9,20,21,22,23,24,25,26,27
#SBATCH --time=2:00:00
#SBATCH --mem=4096
#SBATCH --partition=small

module load python/3.11
module load vsearch/2.28.1
module load pandaseq
module load emboss/6.5.7

source /project/chlandry/users/siaub8/venv_FCY1_300/bin/activate

/project/chlandry/users/siaub8/Reads_analysis.py --id $SLURM_ARRAY_TASK_ID --in /project/chlandry/users/siaub8/LANC016/Raw_modified/CDS_libs --ref ref_CDS_150_trimmed.fa --pre CDS --info Libs_info_150_final.csv --g_open 50 --g_extend 0.5 --min_len 135 --max_len 170 --out FCY1_150_CDS_PD_25-04-07 --seq_type CDS
