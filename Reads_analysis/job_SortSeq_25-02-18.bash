#!/bin/bash
#SBATCH -D /project/chlandry/users/siaub8
#SBATCH -J SortSeq_FCY1
#SBATCH -o SortSeq_24-02-18_logs/SortSeq%j.log
#SBATCH -c 1
#SBATCH --mail-type=ALL
#SBATCH --mail-user=simon.aube.2@ulaval.ca
#SBATCH --array=1-33
#SBATCH --time=10:00
#SBATCH --mem=2048
#SBATCH --partition=small

module load python/3.11
module load vsearch/2.28.1
module load pandaseq
module load emboss/6.5.7

source /project/chlandry/users/siaub8/venv_FCY1_300/bin/activate

/project/chlandry/users/siaub8/Mut_calling.py --id $SLURM_ARRAY_TASK_ID --in /project/chlandry/users/siaub8/AVITI_sortseq_2025-02-18 --ref ref_SortSeq_25-02-18.fa --info Libs_info_SortSeq_25-02-18_SA.csv --out SortSeq_25-02-18
