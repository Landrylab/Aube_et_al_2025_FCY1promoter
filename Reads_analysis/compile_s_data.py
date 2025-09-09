#!/usr/bin/env python

# Short script to compile multiple csv files of selection coefficients (one per culture in a DMS screen)
# produced by the Reads_analysis.py pipeline

# Import statements
import argparse
import os
import pandas as pd
import datetime
import Mut_call_func

# Setting up the argparser
PARSER = argparse.ArgumentParser(description="Short script to compile selection coefficients obtained for multiple "
                                             "replicate cultures in the same screen.")
PARSER.add_argument("--in",
                    metavar="STRING", dest="input_path", type=str,
                    help="Path to the root directory of the output on the relevant Reads_analysis.py run.")

PARSER.add_argument("--out_suffix",
                    metavar="STRING", default=f'{str(datetime.date.today())}', dest="name_suffix", type=str,
                    help="Suffix to be added to the name of final file (compiled selection coefficients). "
                         "Defaults to the current date, in YYYY-MM-DD format.")

PARSER.add_argument("--out",
                    metavar="STRING", dest="out_path", type=str,
                    help="Path to the directory where output files should be saved")

PARSER.add_argument("--cultures",
                    metavar="STRING", dest="cultures_path", type=str,
                    help="Name of the csv file containing information on the different cultures included in a "
                         "given screen (biological replicates, with or without 5-FC, with or without spike-ins). "
                         "This should be given for each alias included in the info csv given to Reads_analysis.py "
                         "when computing the selection coefficients. This file should be in the starting working "
                         "directory where the script is called.")

PARSER.add_argument("--lib_type",
                    metavar="STRING", dest="lib_type", type=str,
                    help="Type of samples for which the data will be compiled. This identifier needs to exist"
                         " in the 'Lib_type' column of the cultures file.")

PARSER.add_argument("--n_muts",
                    metavar="INT", dest="mut_num", type=int,
                    help="Type of mutants for which the selection coefficients should be fetched and compiled."
                         "Set to 1 for single mutants, to 2 for double mutant, etc (if the data exists).")

PARSER.add_argument("--seq_type",
                    metavar="STRING", dest="seq_type", type=str, choices=['prom', 'CDS'], default='prom',
                    help="Type of sequencing library analyzed. Should be 'prom' for promoter and 'CDS' for "
                         "coding sequence. Defaults to 'prom'.")

PARSER.add_argument("--start_num",
                    metavar="INT", dest="start_num", type=int, default=1,
                    help="Only used with 'CDS' type libraries. Specifies the codon (numbered according to the "
                         "reference CDS) within which the first base of the sequenced fragment is found.")

PARSER.add_argument("--codon_off",
                    metavar="INTEGER", dest="codon_off", type=int, default=0,
                    help="Only used with 'CDS' type libraries. Specifies the position within the first codon at which"
                         " the reference sequence begins. If 0, the reference starts on the first position of a codon, "
                         "while 1 and 2 respectively refer to the second and third positions of the first codon. "
                         "Defaults to 0.")

ARGS = PARSER.parse_args()

in_path = ARGS.input_path
name_suffix = ARGS.name_suffix
out_dir = ARGS.out_path
cultures_info = ARGS.cultures_path
lib_type = ARGS.lib_type
mut_num = ARGS.mut_num
seq_type = ARGS.seq_type
start_num = ARGS.start_num
codon_off = ARGS.codon_off

# Keeping the starting working directory as path and making sure to get full path to the csv file
start_path = os.getcwd()
cultures_full = os.path.join(start_path, cultures_info)

# Creating the out directory if it does not already exist
if not os.path.exists(out_dir):
    try:
        os.makedirs(out_dir)
    except FileExistsError:
        pass

# A) Preparing a dataframe of culture infos without duplicates (one row per alias)

# Importing the datafame
alias_info_df = pd.read_csv(cultures_full)
alias_info_df = alias_info_df[alias_info_df['Lib_type'] == lib_type].copy().reset_index(drop=True)

# Keeping only the relevant columns
alias_info_df = alias_info_df[['Alias', 'Bio_rep', 'Tech_rep', 'With_CDS_spikeins', 'Condition']].copy()

# Removing duplicate rows
alias_info_ready = alias_info_df.drop_duplicates(ignore_index=True)

# B) Reading all the data (selection coefficients)

# The first step is to get the list of subdirectories (one per culture) and obtain df structure from one of the files
# The name of each subdirectory should be the unique name used throughout Reads_analysis.py for the corresponding
# culture
os.chdir(in_path)
in_list = os.listdir()

dirs_list = []
for entry in in_list:
    if os.path.isdir(entry):
        dirs_list += [entry]

df_first = pd.read_csv(f'{dirs_list[0]}/s_coefficients/s_all_{dirs_list[0]}_{mut_num}n.csv', nrows=5)
model_cols = list(df_first.columns) + ['Alias']
df_model = pd.DataFrame(columns=model_cols)
df_all_s = df_model.copy()

# Then all csv files of selection coefficients are iteratively read and concatenated to df_all_s
for lib_dir in dirs_list:
    lib_alias = lib_dir.split('_L')[1]
    try:
        s_df = pd.read_csv(f'{lib_dir}/s_coefficients/s_all_{lib_dir}_{mut_num}n.csv')
        s_df['Alias'] = int(lib_alias)
        df_all_s = pd.concat([df_all_s, s_df]).reset_index(drop=True)

    except FileNotFoundError:
        print(f'No selection coefficients are available for library {lib_alias}.')

# The full dataframe is then merged with the alias infos
df_all_s = pd.merge(df_all_s, alias_info_ready, on='Alias', how='left')

# Separating the WT and mutant data
df_WT = df_all_s[df_all_s['Genotype'] == 'WT'].copy().reset_index(drop=True)
df_s_mutants = df_all_s[df_all_s['Genotype'] != 'WT'].copy().reset_index(drop=True)

# For the mutants, columns with mutation position and type are added
type_dict = {'S': 'Substitution', 'I': 'Insertion', 'D': 'Deletion'}

if mut_num == 1:

    if seq_type == 'prom':
        df_s_mutants['Mutation type'] = df_s_mutants['Genotype'].apply(lambda x: type_dict[x[0]])
        df_s_mutants['Position'] = df_s_mutants['Genotype'].apply(lambda x: int(x.split('-')[1].split('_')[0]))

    elif seq_type == 'CDS':
        df_s_mutants['Mutation types'] = df_s_mutants['Genotype'].apply(lambda x: list(set([type_dict[code[0]] for code in x.split(', ')])))
        df_s_mutants['Position (codon)'] = df_s_mutants['Genotype'].apply(lambda x: int(Mut_call_func.codes_to_codons(x.split(', '), start_num, codon_off)[0]))
        # The information on WT and mutant codons will be added locally during the post-processing of s coefficients

elif mut_num > 1:

    if seq_type == 'prom':
        df_s_mutants['Mutation types'] = df_s_mutants['Genotype'].apply(lambda x: list(set([type_dict[code[0]] for code in x.split(', ')])))
        df_s_mutants['Positions'] = df_s_mutants['Genotype'].apply(lambda x: [int(code.split('-')[1].split('_')[0]) for code in x.split(', ')])

    elif seq_type == 'CDS':
        df_s_mutants['Mutation types'] = df_s_mutants['Genotype'].apply(lambda x: list(set([type_dict[code[0]] for code in x.split(', ')])))
        df_s_mutants['Positions (codons)'] = df_s_mutants['Genotype'].apply(lambda x: Mut_call_func.codes_to_codons(x.split(', '), start_num, codon_off))

    # As for the single mutant case, further processing into amino acid changes will be done locally

# Saving both dataframes
os.chdir(start_path)

df_WT.to_csv(f'{out_dir}/s_data_WT_{mut_num}n_{name_suffix}.csv', index=False)
df_s_mutants.to_csv(f'{out_dir}/s_all_{mut_num}n_{name_suffix}.csv', index=False)
