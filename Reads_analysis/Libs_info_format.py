#!/usr/bin/env python

# Short script to reformat the library names of the 150 paired-end libraries,
# to make them consistent with those used in the subsequent 300 paired-end run.

# This is done both within the libs_info csv file and in the directory containing the R1 and R2 files.
# The directory containing the renamed files is then uploaded directly to cluster, where the Reads_analysis.py
# pipeline is run.

import pandas as pd
import numpy as np
import os
import glob

# 1) Reformatting the names in the libs_info csv file

libs_info_150 = pd.read_csv("C:/Users/tiger/OneDrive - Université Laval/PhD/Chapitre 1/DMS_planning"
                            "/paired_end_300_sequencing/Libs_info_150_final.csv")


# Defining a function to reformat the sample names
def rename_samp(samp_name):
    lib_alias = samp_name.split('-')[0]
    timepoint = samp_name.split('-')[1].split('_')[0]
    lib_type = samp_name.split('-')[1].split('_')[1]
    new_name = f'{lib_type}-L{lib_alias}-{timepoint}'

    return new_name


# Reformatting the names
libs_info_150['Sample_name'] = libs_info_150['Sample_name'].apply(rename_samp)

# Saving the new version of the csv
libs_info_150.to_csv("C:/Users/tiger/OneDrive - Université Laval/PhD/Chapitre 1/DMS_planning"
                     "/paired_end_300_sequencing/Libs_info_150_final.csv", index=False)

# 2) Reformatting the file names

# Moving to the directory containing the "raw modified" files (in which the headers have been reformatted to prevent
# crashes when running pandaseq) and producing a list of filenames
os.chdir("F:/Sequencing_LANC016/LANC016/Raw_modified")
files_list = glob.glob('*.fastq.gz')

# Then, renaming each file
for gz_name in files_list:
    lib_name = gz_name.split('.fastq.gz')[0]
    lib_id = lib_name.split('-')[0]
    timepoint = lib_name.split('-')[1].split('_')[0]
    lib_type = lib_name.split('-')[1].split('_')[1]
    read_type = lib_name.split('-')[1].split('_')[2]

    new_filename = f'{lib_type}_L{lib_id}_{timepoint}_{read_type}.fastq.gz'

    os.rename(gz_name, new_filename)
