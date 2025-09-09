#!/bin/bash
#SBATCH -D /project/chlandry/users/siaub8
#SBATCH -J FCY1_150_F3F4_freqs
#SBATCH -o Compile_logs_25-04-08/FCY1_150_F3F4_freqs.log
#SBATCH -c 2
#SBATCH --mail-type=ALL
#SBATCH --mail-user=simon.aube.2@ulaval.ca
#SBATCH --time=1:00
#SBATCH --mem=2048
#SBATCH --partition=small

module load python/3.11

source /project/chlandry/users/siaub8/venv_FCY1_300/bin/activate

/project/chlandry/users/siaub8/compile_mut_freqs.py --in /project/chlandry/users/siaub8/FCY1_150_F3F4_25-04-07 --out_suffix FCY1_150_F3F4 --out Data_compiled_25-04-08 --cultures Libs_info_150_final.csv --lib_type F3F4 --n_muts 1
