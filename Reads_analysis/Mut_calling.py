#!/usr/bin/env python

# Import statements
import os
import subprocess
import pandas as pd
import numpy as np
import glob
import scipy.stats as stats
import statsmodels.api as sm
import argparse
import Mut_call_func
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import upsetplot
import seaborn as sns
from pandas.api.types import CategoricalDtype

# Setting up the argparser
PARSER = argparse.ArgumentParser(description="Script to analyze raw data from a sequencing run, call mutations and "
                                             "return counts for all genotypes within each provided sample.")

PARSER.add_argument("-n", "--id",
                    metavar="STRING", dest="samp_id", type=str,
                    help="Unique (numerical) identifier for the current sample. Should refer to the content of a "
                         "'Unique_id' column in the provided info file. This allows to map the indexes of a Slurm "
                         "job array to the different samples.")

PARSER.add_argument("--in",
                    metavar="STRING", dest="input_path", type=str,
                    help="Path to the directory containing all the input fastq data (should contain one subdirectory "
                         "for each sample). This script assumes one subdirectory for each pair of R1 and R2 files, "
                         "named according to the sample name associated with the 'Unique_id' in the info file.")

PARSER.add_argument("--ref",
                    metavar="STRING", dest="ref_path", type=str,
                    help="Path to a Fasta file containing the reference sequence (WT of the sequenced region). "
                         "Since the primers are removed prior to the alignments, this reference must EXCLUDE the "
                         "primer regions at both ends.")

PARSER.add_argument("--info",
                    metavar="STRING", dest="info_path", type=str,
                    help="Path to the csv file containing the necessary information on all sequencing libraries "
                         "(number of N, etc)")

PARSER.add_argument("--prime_l",
                    metavar="INTEGER", dest="prime_l", type=int, default=20,
                    help="Length of the forward and reverse primers, in bases. Assumed to be the same for both, so "
                         "only one value needs to be supplied. Defaults to 20 bp.")

PARSER.add_argument("--g_open",
                    metavar="INTEGER", dest="g_open", type=int, default=4,
                    help="Gap open penalty, to be supplied to Needle when aligning the reads to the reference. "
                         "Defaults to 4.")

PARSER.add_argument("--g_extend",
                    metavar="INTEGER", dest="g_extend", type=int, default=4,
                    help="Gap extend penalty, to be supplied to Needle when aligning the reads to the reference. "
                         "Defaults to 4.")

# This script will output everything generate files containing EVERY call, as well as more useful files contaning only
# the WT and all single mutants

PARSER.add_argument("--out",
                    metavar="STRING", dest="out_path", type=str,
                    help="Path to the directory where all subdirectories containing output files should be created "
                         "(if they do not already exist).")

ARGS = PARSER.parse_args()
sample_id = int(ARGS.samp_id)
in_dir = ARGS.input_path
ref_path = ARGS.ref_path
info_path = ARGS.info_path
primers_length = ARGS.prime_l
gap_open = ARGS.g_open  # Penalties of (4,4) work well for the F3F4 libraries,
gap_extend = ARGS.g_extend  # which contain both substitutions and indels
out_dir = ARGS.out_path

# Keeping the starting working directory as path
start_path = os.getcwd()
info_full = os.path.join(start_path, info_path)


# Creating the out directory, if it does not exist
if not os.path.exists(out_dir):
    try:
        os.makedirs(out_dir)
    except FileExistsError:
        pass

# Importing the Info file and obtaining the sample name corresponding to the unique_id
libs_info = pd.read_csv(info_full)
sample_name = libs_info[libs_info['Unique_id'] == sample_id]['Sample_name'].values[0]

# Creating the subdirectory for the current sample, unless it exists
path_lib = os.path.join(out_dir, f'{sample_name}')

if not os.path.exists(path_lib):
    try:
        os.makedirs(path_lib)
    except FileExistsError:
        pass

# 1) Merging the R1 and R2 reads using pandaseq
# This time, there is no need to reformat the headers before running pandaseq.
# This is done iteratively for each of the timepoints

os.chdir(path_lib)  # Moving to the output subdirectory for the current pool

# Creating the subdirectory where the merged reads will be saved
if not os.path.exists('Merged_pandaseq'):
    try:
        os.makedirs('Merged_pandaseq')
    except FileExistsError:
        pass

# Preparing the corresponding paths
names_list = [f'{sample_name}_R1.fastq.gz', f'{sample_name}_R2.fastq.gz']

R1_path = os.path.join(in_dir, sample_name, names_list[0])
R2_path = os.path.join(in_dir, sample_name, names_list[1])

merged_path = os.path.join(start_path, path_lib, 'Merged_pandaseq')
merged_out = open(f'{merged_path}/{sample_name}_merged.fasta', 'w')

# Calling pandaseq to perform the merging
panda_list = ['pandaseq', '-f', f'{R1_path}', '-r', f'{R2_path}', '-k', '4', '-B', '-N', 'T', '1']
subprocess.call(panda_list, stdout=merged_out, stderr=subprocess.DEVNULL, cwd=merged_path)

# 2a) Trimming the added Ns using vsearch
# Create the subdirectory where the trimmed reads will be stored
if not os.path.exists('Trimmed_vsearch'):
    try:
        os.makedirs('Trimmed_vsearch')
    except FileExistsError:
        pass

# Extracting the numbers of N for the current sample
N_left = libs_info[libs_info['Unique_id'] == sample_id]['N_forward'].values[0]
N_right = libs_info[libs_info['Unique_id'] == sample_id]['N_reverse'].values[0]

# Preparing the corresponding paths
merged_full = os.path.join(start_path, path_lib, 'Merged_pandaseq', f'{sample_name}_merged.fasta')
n_trimmed_out = os.path.join(start_path, path_lib, 'Trimmed_vsearch')

vsearch_list = ['vsearch', '--fastx_filter', f'{merged_full}', '--fastq_stripleft', f'{N_left}',
                '--fastq_stripright', f'{N_right}', '--fastaout', f'{n_trimmed_out}/{sample_name}_noNs.fasta']

subprocess.call(vsearch_list, stderr=subprocess.DEVNULL, cwd=n_trimmed_out)

# 2b) Trimming the primer sequences (first and last 20 bp of each read)
# It would be better to keep these sequences somewhere, but it does not seem possible to redirect them to a file
noNs_path = os.path.join(start_path, path_lib, 'Trimmed_vsearch', f'{sample_name}_noNs.fasta')

vsearch_primers = ['vsearch', '--fastx_filter', f'{noNs_path}', '--fastq_stripleft', f'{primers_length}',
                   '--fastq_stripright', f'{primers_length}', '--fastaout',
                   f'{n_trimmed_out}/{sample_name}_trimmed.fasta']

subprocess.call(vsearch_primers, stderr=subprocess.DEVNULL, cwd=n_trimmed_out)

# 3) Aggregating all identical trimmed reads
# To loop through files in the newly created 'Trimmed_vsearch' directory

# Create the subdirectory where the aggregated reads will be stored
if not os.path.exists('Aggregated_vsearch'):
    try:
        os.makedirs('Aggregated_vsearch')
    except FileExistsError:
        pass

# Preparing the corresponding paths
trimmed_path = os.path.join(start_path, path_lib, 'Trimmed_vsearch', f'{sample_name}_trimmed.fasta')
out_agg = os.path.join(start_path, path_lib, 'Aggregated_vsearch')

# Calling vsearch to perform the aggregation
agg_list = ['vsearch', '--derep_fulllength', f'{trimmed_path}', '--relabel', 'seq', '--lengthout', '--sizeout',
            '--output', f'{out_agg}/{sample_name}_agg.fasta']

subprocess.call(agg_list, stderr=subprocess.DEVNULL, cwd=out_agg)

# 4) Aligning all aggregated reads to the reference using Needle
# Create the subdirectory where the alignments will be stored
if not os.path.exists('Aligned_Needle'):
    try:
        os.makedirs('Aligned_Needle')
    except FileExistsError:
        pass

# Preparing the paths
agg_path = os.path.join(start_path, path_lib, 'Aggregated_vsearch', f'{sample_name}_agg.fasta')
ref_updated = os.path.join(start_path, ref_path)
out_needle = os.path.join(start_path, path_lib, 'Aligned_Needle')

# Calling Needle to perform the alignment
needle_list = ['needle', '-auto', '-gapopen', f'{gap_open}', '-gapextend', f'{gap_extend}', '-bsequence',
               f'{agg_path}', '-asequence', f'{ref_updated}', '-aformat3', 'markx10', '-outfile',
               f'{out_needle}/{sample_name}_aligned.needle']

subprocess.call(needle_list, stderr=subprocess.DEVNULL, cwd=out_needle)

# 5) Calling mutations in all alignments
# Using functions defined in Mut_call_func.py

# Create the subdirectory where the csv of called mutations will be stored
if not os.path.exists('Mutations'):
    try:
        os.makedirs('Mutations')
    except FileExistsError:
        pass

Mut_call_func.find_mutations('Aligned_Needle', 'Mutations')

# Generating an Upset plot of mutation categories
os.chdir(os.path.join(start_path, out_dir, f'{sample_name}', 'Mutations'))
muts_list = glob.glob('*.csv')
muts_list.sort()

os.chdir('..')

# Initializing the pdf which will contain the figures\n",
mut_cats_fig = PdfPages(f"Mutations/{sample_name}_Muts_Upset.pdf")

for muts_file in muts_list:
    muts_all = pd.read_csv(f'Mutations/{muts_file}')
    muts_cat = muts_all[['Mut_categories', 'N_reads']].copy()
    muts_cat['Mut_categories'] = muts_cat['Mut_categories'].fillna('Wild-type')
    muts_cat = muts_cat.groupby(by='Mut_categories', as_index=False).sum()

    # Formatting the data for the generation of the Upset plot
    data_upset = upsetplot.from_memberships(muts_cat['Mut_categories'].str.strip(', ').str.split(', '), data=muts_cat)

    # Constructing the Upset plot
    upset_fig = plt.figure()
    upsetplot.UpSet(data_upset, sum_over='N_reads', sort_by='cardinality').plot(fig=upset_fig)
    upset_fig.suptitle(f"{sample_name}", fontsize=24)

    upset_fig.savefig(mut_cats_fig, format='pdf')
    plt.close()

mut_cats_fig.close()

# 6) Generating a "simplified" file with only single mutants
# The muts_all file previously imported is used. This assumes that there is only one csv file in the directory
if len(muts_list) != 1:
    raise Exception('There should be exactly one csv of mutation calls in the current folder!')

# Adding 'WT' annotations for genotypes which do not contain any mutation
muts_df = muts_all[['Library', 'Ref_seq', 'Mutations', 'Mut_categories', 'N_reads']].copy()
muts_df[['Mutations', 'Mut_categories']] = muts_df[['Mutations', 'Mut_categories']].fillna('WT')
n_total = muts_df['N_reads'].sum()
muts_agg = muts_df.groupby(by=['Library', 'Ref_seq', 'Mutations', 'Mut_categories'], as_index=False).sum()
muts_agg = muts_agg.sort_values(by='N_reads', ascending=False).reset_index(drop=True)

# Separating the WT from the mutants
if muts_agg[muts_agg['Mutations'] == 'WT'].shape[0] != muts_agg[muts_agg['Mut_categories'] == 'WT'].shape[0]:
    raise Exception(f"There are unexpected 'WT' labels in the mutation data!")

wt_subset = muts_agg[muts_agg['Mutations'] == 'WT'].copy().reset_index(drop=True)
n_wt = wt_subset['N_reads'].sum()
mutants_subset = muts_agg[muts_agg['Mutations'] != 'WT'].copy().reset_index(drop=True)

# Converting the lists of mutations and mutation categories into frozen sets
mutants_subset['Mutations'] = mutants_subset['Mutations'].apply(lambda x: frozenset(x.split(', ')))
mutants_subset['Mut_categories'] = mutants_subset['Mut_categories'].apply(lambda x: frozenset(x.strip(', ').split(', ')))

# Selecting variants where the total number of mutated nucleotides equals 1
mutants_subset['n_mutated'] = mutants_subset['Mutations'].apply(lambda x: np.sum(np.array([float(mut.split('_')[1].split('nt')[0]) for mut in list(x)])))
mutants_subset['n_mutated'] = mutants_subset['n_mutated'].apply(lambda x: x if x == 1 else np.NaN)
mutants_subset = mutants_subset.dropna(subset=['n_mutated']).reset_index(drop=True)

# Any identical mutations are aggregated (should not be necessary, but just in case)
mut_count_agg = mutants_subset.copy().drop(['Mut_categories'], axis='columns')
mut_count_agg = mut_count_agg.groupby(by=['Library', 'Ref_seq', 'Mutations', 'n_mutated'],
                                      as_index=False).sum().reset_index(drop=True)

# Adding columns for mut position and mut_type
type_dict = {'S': 'Substitution', 'I': 'Insertion', 'D': 'Deletion'}
mut_count_agg = mut_count_agg.rename(columns={'Mutations': 'Genotype'})
mut_count_agg['Genotype'] = mut_count_agg['Genotype'].apply(lambda x: ', '.join(list(x)))  # To avoid keeping a frozenset
mut_count_agg['Mutation type'] = mut_count_agg['Genotype'].apply(lambda x: type_dict[x[0]])
mut_count_agg['Position'] = mut_count_agg['Genotype'].apply(lambda x: int(x.split('-')[1].split('_')[0]))

# Adding back the WT
wt_subset = wt_subset.rename(columns={'Mutations': 'Genotype'})
wt_subset['Mutation type'] = np.nan
wt_subset['Position'] = np.nan
mut_count_agg = pd.concat([mut_count_agg, wt_subset]).reset_index(drop=True)

# Adding the total number of reads and number of WT reads as columns
mut_count_agg['WT_reads'] = n_wt
mut_count_agg['Total_reads'] = n_total

# Saving the file with single mutants
mut_count_agg.to_csv(f"Mutations/Muts_1n_{sample_name}.csv", index=False)
