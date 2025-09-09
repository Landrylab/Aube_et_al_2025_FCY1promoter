#!/usr/bin/env python

# Short script to compile all the read counts from multiple samples which are part of the same sort-seq experiment

# Import statements
import argparse
import os
import pandas as pd
import warnings
import datetime

PARSER = argparse.ArgumentParser(description="Short script to compile read counts obtained for each of the different "
                                             "samples in a Sort-seq experiment.")
PARSER.add_argument("--in",
                    metavar="STRING", dest="input_path", type=str,
                    help="Path to the root directory of the output on the relevant Mut_calling.py run.")

PARSER.add_argument("--out_suffix",
                    metavar="STRING", default=f'{str(datetime.date.today())}', dest="name_suffix", type=str,
                    help="Suffix to be added to the name of final file (compiled read counts). "
                         "Defaults to the current date, in YYYY-MM-DD format.")

PARSER.add_argument("--file_pre",
                    metavar="STRING", dest="name_pre", type=str,
                    help="Filename prefix common to all csv files which should be combined. This should be "
                         "followed by {sample_name}.csv, where sample_name is also the name of the sub-directory "
                         "contaning the corresponding file.")

ARGS = PARSER.parse_args()

in_path = ARGS.input_path
name_suffix = ARGS.name_suffix
files_prefix = ARGS.name_pre

# Keeping the starting working directory as path
start_path = os.getcwd()

# The first step is to get the list of subdirectories (one per culture) and obtain df structure from one of the files
# The name of each subdirectory should be the unique name used throughout Reads_analysis.py for the corresponding
# culture
os.chdir(in_path)
in_list = os.listdir()

dirs_list = []
for entry in in_list:
    if os.path.isdir(entry):
        dirs_list += [entry]

df_first = pd.read_csv(f'{dirs_list[0]}/Mutations/{files_prefix}_{dirs_list[0]}.csv', nrows=5)
df_first = df_first.drop(columns=['Mut_categories'])
model_cols = list(df_first.columns)
df_model = pd.DataFrame(columns=model_cols)
df_all_reads = df_model.copy()

# Then all csv files of read counts are iteratively read and concatenated to df_all_reads
for lib_dir in dirs_list:
    try:
        reads_df = pd.read_csv(f'{lib_dir}/Mutations/{files_prefix}_{lib_dir}.csv')
        df_all_reads = pd.concat([df_all_reads, reads_df]).reset_index(drop=True)

    except FileNotFoundError:
        print(f'No read counts are available for sample {lib_dir}.')


# Saving the dataframe
os.chdir(start_path)

df_all_reads.to_csv(f'counts_all_{name_suffix}.csv', index=False)
