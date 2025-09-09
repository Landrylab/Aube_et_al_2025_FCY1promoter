#!/bin/bash
#SBATCH -D /project/chlandry/users/siaub8
#SBATCH -J sort_seq_compile
#SBATCH -o sort_25-02-18__compile.log
#SBATCH -c 2
#SBATCH --mail-type=ALL
#SBATCH --mail-user=simon.aube.2@ulaval.ca
#SBATCH --time=1:00
#SBATCH --mem=2048
#SBATCH --partition=small

module load python/3.11

source /project/chlandry/users/siaub8/venv_FCY1_300/bin/activate

/project/chlandry/users/siaub8/compile_reads.py --in /project/chlandry/users/siaub8/SortSeq_25-02-18 --out_suffix SortSeq_25-02 --file_pre Muts_1n
