#!/usr/bin/env python

# Short script to compile multiple csv files of relative frequencies (one per culture in a DMS screen)
# produced by the Reads_analysis.py pipeline

# Import statements
import argparse
import os
import pandas as pd
import datetime
import Mut_call_func

# Setting up the argparser
PARSER = argparse.ArgumentParser(description="Short script to compile relative frequencies obtained for multiple "
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

PARSER.add_argument("--time",
                    metavar="INT", dest="n_times", type=int, default=4,
                    help="Number of sequencing timepoints (including T0) in the experiment. A default value of 4 will "
                         "be interpreted as files having suffixes 'T0', 'T1', 'T2' and 'T3'")

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
n_timepoints = ARGS.n_times
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

# B) Reading all the data (relative frequencies)

# The first step is to get the list of subdirectories (one per culture) and obtain df structure from one of the files
# The name of each subdirectory should be the unique name used throughout Reads_analysis.py for the corresponding
# culture
os.chdir(in_path)
in_list = os.listdir()

dirs_list = []
for entry in in_list:
    if os.path.isdir(entry):
        dirs_list += [entry]

# Instead of taking the columns from the dfs, the desired ones are specified directly
model_cols = ['Library', 'Mutations', 'Alias', 'Rel_wt', 'Rel_tot']
df_model = pd.DataFrame(columns=model_cols)
df_all_freqs = df_model.copy()

# Then, the files for each timepoint within each directory are iteratively read and merged, before being concatenated
# to the full df
for lib_dir in dirs_list:
    lib_alias = lib_dir.split('_L')[1]

    # The T0 file is read, which will allow to iteratively merge the other ones
    freqs_df = pd.read_csv(f'{lib_dir}/Frequencies/rel_freqs_{lib_dir}_T0_T0_{mut_num}n.csv')
    freqs_df = freqs_df[['Library', 'Mutations', 'Rel_wt', 'Rel_tot']].copy()
    freqs_df['Alias'] = int(lib_alias)
    freqs_df['Library'] = lib_dir  # Otherwise, the timepoint will also be included, which will prevent merging
    freqs_df['Timepoint'] = 'T0'

    for time_num in range(1, n_timepoints):
        try:
            freqs_time = pd.read_csv(f'{lib_dir}/Frequencies/rel_freqs_{lib_dir}_T{time_num}_T{time_num}_{mut_num}n.csv')

        except FileNotFoundError:
            print(f'No frequencies are available for T{time_num} of {lib_dir}')
            continue

        freqs_time = freqs_time[['Library', 'Mutations', 'Rel_wt', 'Rel_tot']].copy()
        freqs_time['Alias'] = int(lib_alias)
        freqs_time['Library'] = lib_dir
        freqs_time['Timepoint'] = f'T{time_num}'

        freqs_df = pd.concat([freqs_df, freqs_time]).reset_index(drop=True)

    # Once all timepoints of the current library have been assembled into a df, it is concatenated to df_all_freqs
    df_all_freqs = pd.concat([df_all_freqs, freqs_df]).reset_index(drop=True)

# The full dataframe is then merged with the alias infos
df_all_freqs = pd.merge(df_all_freqs, alias_info_ready, on='Alias', how='left')

# Separating the WT and mutant data
df_WT = df_all_freqs[df_all_freqs['Mutations'] == 'WT'].copy().reset_index(drop=True)
df_mutants = df_all_freqs[df_all_freqs['Mutations'] != 'WT'].copy().reset_index(drop=True)

# For the mutants, columns with mutation position and type are added
type_dict = {'S': 'Substitution', 'I': 'Insertion', 'D': 'Deletion'}

if mut_num == 1:

    if seq_type == 'prom':
        df_mutants['Mutation type'] = df_mutants['Mutations'].apply(lambda x: type_dict[x[0]])
        df_mutants['Position'] = df_mutants['Mutations'].apply(lambda x: int(x.split('-')[1].split('_')[0]))

    elif seq_type == 'CDS':
        df_mutants['Mutation types'] = df_mutants['Mutations'].apply(lambda x: list(set([type_dict[code[0]] for code in x.split(', ')])))
        df_mutants['Position (codon)'] = df_mutants['Mutations'].apply(lambda x: int(Mut_call_func.codes_to_codons(x.split(', '), start_num, codon_off)[0]))

        # The information on WT and mutant codons will be added locally during the post-processing of s coefficients

elif mut_num > 1:

    if seq_type == 'prom':
        df_mutants['Mutation types'] = df_mutants['Mutations'].apply(lambda x: list(set([type_dict[code[0]] for code in x.split(', ')])))
        df_mutants['Positions'] = df_mutants['Mutations'].apply(lambda x: [int(code.split('-')[1].split('_')[0]) for code in x.split(', ')])

    elif seq_type == 'CDS':
        df_mutants['Mutation types'] = df_mutants['Mutations'].apply(lambda x: list(set([type_dict[code[0]] for code in x.split(', ')])))
        df_mutants['Positions (codons)'] = df_mutants['Mutations'].apply(lambda x: Mut_call_func.codes_to_codons(x.split(', '), start_num, codon_off))

    # As for the single mutant case, further processing into amino acid changes will be done locally

# Saving both dataframes
os.chdir(start_path)

df_WT.to_csv(f'{out_dir}/freqs_WT_{mut_num}n_{name_suffix}.csv', index=False)
df_mutants.to_csv(f'{out_dir}/freqs_all_{mut_num}n_{name_suffix}.csv', index=False)
