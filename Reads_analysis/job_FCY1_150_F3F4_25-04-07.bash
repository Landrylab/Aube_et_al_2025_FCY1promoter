#!/bin/bash
#SBATCH -D /project/chlandry/users/siaub8
#SBATCH -J FCY1_150_F3F4_array
#SBATCH -o Logs_FCY1_25-04-07/FCY1_F3F4_150-%j.log
#SBATCH -c 1
#SBATCH --mail-type=ALL
#SBATCH --mail-user=simon.aube.2@ulaval.ca
#SBATCH --array=1,12,13,20,21,3,4,5,6,7,8,9,14,15,16,17,18,19,22,23,24,25,26,27
#SBATCH --time=2:00:00
#SBATCH --mem=4096
#SBATCH --partition=small

module load python/3.11
module load vsearch/2.28.1
module load pandaseq
module load emboss/6.5.7

source /project/chlandry/users/siaub8/venv_FCY1_300/bin/activate

/project/chlandry/users/siaub8/Reads_analysis.py --id $SLURM_ARRAY_TASK_ID --in /project/chlandry/users/siaub8/LANC016/Raw_modified/F3F4_libs --ref ref_F3F4_150_trimmed.fa --pre F3F4 --info Libs_info_150_final.csv --min_len 160 --max_len 190 --out FCY1_150_F3F4_25-04-07 --min_reads 125
