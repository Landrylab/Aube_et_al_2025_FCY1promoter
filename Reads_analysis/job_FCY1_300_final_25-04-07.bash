#!/bin/bash
#SBATCH -D /project/chlandry/users/siaub8
#SBATCH -J FCY1_300_final_array
#SBATCH -o Logs_FCY1_25-04-07/FCY1_300_final-%j.log
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

/project/chlandry/users/siaub8/Reads_analysis.py --id $SLURM_ARRAY_TASK_ID --in /project/chlandry/users/siaub8/FC_300_SA --ref ref_F3F4_300_trimmed.fa --pre FCY1 --nested True --info Libs_info_300.csv --min_len 310 --max_len 345 --out FCY1_300_final_25-04-07
