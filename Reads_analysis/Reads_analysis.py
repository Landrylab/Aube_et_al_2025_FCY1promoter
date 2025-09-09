#!/usr/bin/env python

# Import statements
import os
import subprocess
import pandas as pd
import numpy as np
import glob
import math
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
PARSER = argparse.ArgumentParser(description="Script to analyse raw data from a DMS screen with multiple timepoints, "
                                             "one pool at a time.")

PARSER.add_argument("-n", "--id",
                    metavar="STRING", dest="pool_id", type=str,
                    help="Unique identifier for the current pool (used for the data of all the timepoints"
                         " for the same culture)")

PARSER.add_argument("--in",
                    metavar="STRING", dest="input_path", type=str,
                    help="Path to the directory containing all the input fastq data (should contain one subdirectory "
                         "for each sample)")

PARSER.add_argument("--ref",
                    metavar="STRING", dest="ref_path", type=str,
                    help="Path to a Fasta file containing the reference sequence (WT of the sequenced region). "
                         "Since the primers are removed prior to the alignments, this reference must EXCLUDE the "
                         "primer regions at both ends.")

PARSER.add_argument("--pre",
                    metavar="STRING", dest="lib_pre", type=str,
                    help="Prefix for the sample names. Should be followed by '_L{n}_T{timepoint}'")

PARSER.add_argument("--nested",
                    metavar="BOOL", dest="nested", type=bool, default=False,
                    help="Boolean specifying whether the R1 and R2 files are in library-specific subdirectories "
                         "(True) or instead all in the same 'in' directory (False). Defaults to False.")

PARSER.add_argument("--info",
                    metavar="STRING", dest="info_path", type=str,
                    help="Path to the csv file containing the necessary information on all sequencing libraries "
                         "(number of N, etc)")

PARSER.add_argument("--prime_l",
                    metavar="INTEGER", dest="prime_l", type=int, default=20,
                    help="Length of the forward and reverse primers, in bases. Assumed to be the same for both, so "
                         "only one value needs to be supplied. Defaults to 20 bp.")

PARSER.add_argument("--g_open",
                    metavar="FLOAT", dest="g_open", type=float, default=4,
                    help="Gap open penalty, to be supplied to Needle when aligning the reads to the reference. "
                         "Defaults to 4.")

PARSER.add_argument("--g_extend",
                    metavar="FLOAT", dest="g_extend", type=float, default=4,
                    help="Gap extend penalty, to be supplied to Needle when aligning the reads to the reference. "
                         "Defaults to 4.")

PARSER.add_argument("--min_len",
                    metavar="INTEGER", dest="min_len", type=int,
                    help="Minimum sequence length threshold used to remove aberrant reads when compiling the "
                         "relative frequencies of all mutants.")

PARSER.add_argument("--max_len",
                    metavar="INTEGER", dest="max_len", type=int,
                    help="Maximum sequence length threshold used to remove aberrant reads when compiling the "
                         "relative frequencies of all mutants.")

PARSER.add_argument("--time",
                    metavar="INTEGER", dest="n_times", type=int, default=4,
                    help="Number of sequencing timepoints (including T0) in the experiment. A default value of 4 will "
                         "be interpreted as files having suffixes 'T0', 'T1', 'T2' and 'T3'")

PARSER.add_argument("--out",
                    metavar="STRING", dest="out_path", type=str,
                    help="Path to the directory where all subdirectories containing output files should be created "
                         "(if they do not already exist)")

PARSER.add_argument("--min_reads",
                    metavar="INTEGER", dest="reads_arg", type=int, default=85,
                    help="Minimum abundance (in reads) at T0 OR Tfinal necessary for selection coefficient "
                         "calculations to be performed. Only used for double mutants.")

PARSER.add_argument("--seq_type",
                    metavar="STRING", dest="seq_type", type=str, choices=['prom', 'CDS'], default='prom',
                    help="Type of sequencing library analyzed. Should be 'prom' for promoter and 'CDS' for "
                         "coding sequence. Defaults to 'prom'.")

PARSER.add_argument("--codon_off",
                    metavar="INTEGER", dest="codon_off", type=int, default=0,
                    help="Only used with 'CDS' type libraries. Specifies the position within the first codon at which"
                         " the reference sequence begins. If 0, the reference starts on the first position of a codon, "
                         "while 1 and 2 respectively refer to the second and third positions of the first codon. "
                         "Defaults to 0.")

ARGS = PARSER.parse_args()
pool_n = ARGS.pool_id
in_dir = ARGS.input_path
ref_path = ARGS.ref_path
name_pre = ARGS.lib_pre
nested_samples = ARGS.nested
info_path = ARGS.info_path
primers_length = ARGS.prime_l
gap_open = ARGS.g_open  # Penalties of (4,4) work well for the F3F4 libraries,
gap_extend = ARGS.g_extend  # which contain both substitutions and indels
min_length = ARGS.min_len  # Originally, thresholds of 350 bp and 385 bp were used for the 300 paired-end libraries
max_length = ARGS.max_len  # (when the primers were NOT removed prior to the aggregation of sequences)
n_timepoints = ARGS.n_times
out_dir = ARGS.out_path
min_reads = ARGS.reads_arg
seq_type = ARGS.seq_type
codon_off = ARGS.codon_off

# Keeping the starting working directory as path
start_path = os.getcwd()
info_full = os.path.join(start_path, info_path)

# Generating the list of timepoints
times_list = [f'T{n_time}' for n_time in range(n_timepoints)]

# Creating the out directory, if it does not exist
if not os.path.exists(out_dir):
    try:
        os.makedirs(out_dir)
    except FileExistsError:
        pass

# Creating the subdirectory for the current pool, unless it exists
path_lib = os.path.join(out_dir, f'{name_pre}_L{pool_n}')

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

# Looping through the (annoying) subfolders for each timepoint
for timepoint in times_list:
    sample_name = f'{name_pre}_L{pool_n}_{timepoint}'
    names_list = [f'{sample_name}_R1.fastq.gz', f'{sample_name}_R2.fastq.gz']

    if nested_samples:
        R1_path = os.path.join(in_dir, sample_name, names_list[0])
        R2_path = os.path.join(in_dir, sample_name, names_list[1])

    elif not nested_samples:
        R1_path = os.path.join(in_dir, names_list[0])
        R2_path = os.path.join(in_dir, names_list[1])

    out_path = os.path.join(start_path, path_lib, 'Merged_pandaseq')
    out_file = open(f'{out_path}/{sample_name}_merged.fasta', 'w')

    panda_list_v2 = ['pandaseq', '-f', f'{R1_path}', '-r', f'{R2_path}', '-k', '4', '-B', '-N', 'T', '1']

    subprocess.call(panda_list_v2, stdout=out_file, stderr=subprocess.DEVNULL, cwd=out_path)

# It might be worthwhile to trim the last 40 or so base pairs before merging the reads

# 2a) Trimming the added Ns using vsearch
libs_info = pd.read_csv(info_full)  # Importing the csv file containing the numbers of Ns added in each library

# Looping through the files in the newly created 'Merged_pandaseq' subfolder
os.chdir('Merged_pandaseq')
merged_list = glob.glob('*_merged.fasta')

os.chdir('..')  # Go back one level

# Create the subdirectory where the trimmed reads will be stored
if not os.path.exists('Trimmed_vsearch'):
    try:
        os.makedirs('Trimmed_vsearch')
    except FileExistsError:
        pass

for merged_file in merged_list:
    lib_filename = merged_file.split('_merged.fasta')[0]
    lib_code = lib_filename.split('_')
    lib_name = '-'.join(lib_code)

    # Extracting the numbers of N
    info_lib = libs_info[libs_info['Sample_name'] == lib_name].copy()
    N_left = info_lib['N_forward'].values[0]
    N_right = info_lib['N_reverse'].values[0]

    # Preparing the paths
    merged_path = os.path.join(start_path, path_lib, 'Merged_pandaseq', f'{lib_filename}_merged.fasta')
    out_fasta = os.path.join(start_path, path_lib, 'Trimmed_vsearch')

    vsearch_list = ['vsearch', '--fastx_filter', f'{merged_path}', '--fastq_stripleft', f'{N_left}',
                    '--fastq_stripright', f'{N_right}', '--fastaout', f'{out_fasta}/{lib_filename}_noNs.fasta']

    subprocess.call(vsearch_list, stderr=subprocess.DEVNULL, cwd=out_fasta)

# 2b) Trimming the primer sequences (first and last 20 bp of each read)
# It would be better to keep these sequences somewhere, but it does not seem possible to redirect them to a file
os.chdir('Trimmed_vsearch')
noNs_list = glob.glob('*_noNs.fasta')

os.chdir('..')

for noNs_file in noNs_list:
    lib_filename = noNs_file.split('_noNs.fasta')[0]
    noNs_path = os.path.join(start_path, path_lib, 'Trimmed_vsearch', f'{lib_filename}_noNs.fasta')
    out_trimmed = os.path.join(start_path, path_lib, 'Trimmed_vsearch')

    vsearch_primers = ['vsearch', '--fastx_filter', f'{noNs_path}', '--fastq_stripleft', f'{primers_length}',
                       '--fastq_stripright', f'{primers_length}', '--fastaout',
                       f'{out_trimmed}/{lib_filename}_trimmed.fasta']

    subprocess.call(vsearch_primers, stderr=subprocess.DEVNULL, cwd=out_trimmed)

# 3) Aggregating all identical trimmed reads
# To loop through files in the newly created 'Trimmed_vsearch' directory
os.chdir('Trimmed_vsearch')
trimmed_list = glob.glob('*_trimmed.fasta')

os.chdir('..')  # Go back one level

# Create the subdirectory where the aggregated reads will be stored
if not os.path.exists('Aggregated_vsearch'):
    try:
        os.makedirs('Aggregated_vsearch')
    except FileExistsError:
        pass

for trimmed_file in trimmed_list:
    lib_filename = trimmed_file.split('_trimmed.fasta')[0]

    # Preparing the paths
    trimmed_path = os.path.join(start_path, path_lib, 'Trimmed_vsearch', f'{lib_filename}_trimmed.fasta')
    out_agg = os.path.join(start_path, path_lib, 'Aggregated_vsearch')

    agg_list = ['vsearch', '--derep_fulllength', f'{trimmed_path}', '--relabel', 'seq', '--lengthout', '--sizeout',
                '--output', f'{out_agg}/{lib_filename}_agg.fasta']

    subprocess.call(agg_list, stderr=subprocess.DEVNULL, cwd=out_agg)

# 4) Aligning all aggregated reads to the reference using Needle

# To loop through files in the newly created 'Aggregated_vsearch' directory
os.chdir('Aggregated_vsearch')
agg_list = glob.glob('*_agg.fasta')

os.chdir('..')  # Go back one level

# Create the subdirectory where the alignments will be stored
if not os.path.exists('Aligned_Needle'):
    try:
        os.makedirs('Aligned_Needle')
    except FileExistsError:
        pass

for agg_file in agg_list:
    lib_filename = agg_file.split('_agg.fasta')[0]

    # Preparing the paths
    agg_path = os.path.join(start_path, path_lib, 'Aggregated_vsearch', f'{lib_filename}_agg.fasta')
    ref_updated = os.path.join(start_path, ref_path)
    out_needle = os.path.join(start_path, path_lib, 'Aligned_Needle')

    needle_list = ['needle', '-auto', '-gapopen', f'{gap_open}', '-gapextend', f'{gap_extend}', '-bsequence',
                   f'{agg_path}', '-asequence', f'{ref_updated}', '-aformat3', 'markx10', '-outfile',
                   f'{out_needle}/{lib_filename}_aligned.needle']

    subprocess.call(needle_list, stderr=subprocess.DEVNULL, cwd=out_needle)

# 5) Calling mutations in all alignments
# Using functions defined in Mut_call_func.py

# Because of the way that the find_mutations() function has been defined, we only need to run it on the subdirectory
# containing the Needle alignment files

# Create the subdirectory where the csv of called mutations will be stored
if not os.path.exists('Mutations'):
    try:
        os.makedirs('Mutations')
    except FileExistsError:
        pass

Mut_call_func.find_mutations('Aligned_Needle', 'Mutations')

# Generating Upset plots of mutation categories (one pdf with pages for T0, T1, T2 and T3)
os.chdir(os.path.join(start_path, out_dir, f'{name_pre}_L{pool_n}', 'Mutations'))
muts_list = glob.glob('*.csv')
muts_list.sort()

os.chdir('..')

# Initializing the pdf which will contain the figures\n",
mut_cats_fig = PdfPages(f"Mutations/{name_pre}_L{pool_n}_Muts_Upset.pdf")

for muts_file in muts_list:
    lib_filename = muts_file.split('_')[1]
    lib_timepoint = muts_file.split('_')[3]

    muts_all = pd.read_csv(f'Mutations/{muts_file}')
    muts_cat = muts_all[['Mut_categories', 'N_reads']].copy()
    muts_cat['Mut_categories'] = muts_cat['Mut_categories'].fillna('Wild-type')
    muts_cat = muts_cat.groupby(by='Mut_categories', as_index=False).sum()

    # Formatting the data for the generation of the Upset plot
    data_upset = upsetplot.from_memberships(muts_cat['Mut_categories'].str.strip(', ').str.split(', '), data=muts_cat)

    # Constructing the Upset plot
    upset_fig = plt.figure()
    upsetplot.UpSet(data_upset, sum_over='N_reads', sort_by='cardinality').plot(fig=upset_fig)
    upset_fig.suptitle(f"{muts_file.split('_')[2]}, {lib_timepoint}", fontsize=24)

    upset_fig.savefig(mut_cats_fig, format='pdf')
    plt.close()

mut_cats_fig.close()

# 6) Compute relative frequencies for all mutations

# Create the subdirectory where the csv of relative frequencies will be stored
if not os.path.exists('Frequencies'):
    try:
        os.makedirs('Frequencies')
    except FileExistsError:
        pass

Mut_call_func.mutants_rel('Mutations', '*.csv', 'Frequencies', 1, seq_type=seq_type,
                          offset=codon_off, min_length=min_length, max_length=max_length)  # For single mutants

Mut_call_func.mutants_rel('Mutations', '*.csv', 'Frequencies', 2, seq_type=seq_type,
                          offset=codon_off, min_length=min_length, max_length=max_length)  # For double mutants
# This step includes a length-based filtering of reads

# 7) Compute selection coefficients for all genotypes (one file each for single and double)

# Creating the subdirectory for the output selection coefficients
if not os.path.exists('s_coefficients'):
    try:
        os.makedirs('s_coefficients')
    except FileExistsError:
        pass

# Defining the 'full' dictionaries of time (in hours) and approximate numbers of generations
gen_df_path = os.path.join(start_path, 'Gen_time_df.csv')
gen_df = pd.read_csv(gen_df_path)
hours_dict = Mut_call_func.get_time_dict(gen_df, 'Library', 'time', times_list[1:])
total_dict = Mut_call_func.get_time_dict(gen_df, 'Library', 'gen_tot', times_list[1:])
ref_dict = Mut_call_func.get_time_dict(gen_df, 'Library', 'gen_ref', times_list[1:])

# Performing the selection coefficient calculations, for single and double mutants
Mut_call_func.compile_s('Frequencies', '1n.csv', 's_coefficients', times_list,
                        total_dict, ref_dict, hours_dict, dicts_type='full')

# For double mutants, a min_reads threshold is used. It is applied at T0 and T3 (needs to be reached in at least one of
# the two timepoints). This threshold has been determined empirically from the 1n selection coefficients.
Mut_call_func.compile_s('Frequencies', '2n.csv', 's_coefficients', times_list,
                        total_dict, ref_dict, hours_dict, dicts_type='full', min_reads=min_reads)

# Summary figures of the selection coefficients

# First figure: Distribution of the selection coefficient through time
os.chdir('s_coefficients')
s_list = glob.glob('*.csv')
s_list.sort()

os.chdir('..')

s_dists = PdfPages(f"s_coefficients/s_dists_{name_pre}_L{pool_n}.pdf")

for s_file in s_list:

    s_data = pd.read_csv(f's_coefficients/{s_file}')

    # Removing all "true singletons" (variants with no more than 1 read at both T0 and T3)
    s_data = s_data[(s_data['N_T0'] > 1) | (s_data['N_T3'] > 1)].copy().reset_index(drop=True)

    cols_dict = {'Time': ['s_T0toT1_Time', 's_T0toT2_Time', 's_T0toT3_Time'],
                 'Total': ['s_T0toT1_Total', 's_T0toT2_Total', 's_T0toT3_Total'],
                 'Reference': ['s_T0toT1_Ref', 's_T0toT2_Ref', 's_T0toT3_Ref']}

    fig, axs = plt.subplots(1, 3, figsize=(18, 6), sharex=True, sharey=True)
    fig.suptitle(f'For library {name_pre}_L{pool_n}')

    # Plotting the wt data
    time_data = s_data[['Library', 'Genotype'] + cols_dict['Time'] + ['N_T0', 'freq_T0']].copy()
    time_tolong = time_data[['Genotype'] + cols_dict['Time']].copy()
    time_long = pd.melt(time_tolong, id_vars=['Genotype'], value_vars=cols_dict['Time'], value_name='s_coeff',
                        var_name='Time interval')

    sns.kdeplot(data=time_long, x='s_coeff', hue='Time interval', cut=0, ax=axs[0])
    axs[0].set_title('Time-based approach', size='large')

    # Plotting the Total data
    total_data = s_data[['Library', 'Genotype'] + cols_dict['Total'] + ['N_T0', 'freq_T0']].copy()

    total_tolong = total_data[['Genotype'] + cols_dict['Total']].copy()
    total_long = pd.melt(total_tolong, id_vars=['Genotype'], value_vars=cols_dict['Total'], value_name='s_coeff',
                         var_name='Time interval')

    sns.kdeplot(data=total_long, x='s_coeff', hue='Time interval', cut=0, ax=axs[1])
    axs[1].set_title('OD-based approach', size='large')

    # Plotting the ref-based data
    ref_data = s_data[['Library', 'Genotype'] + cols_dict['Reference'] + ['N_T0', 'freq_T0']].copy()

    ref_tolong = ref_data[['Genotype'] + cols_dict['Reference']].copy()
    ref_long = pd.melt(ref_tolong, id_vars=['Genotype'], value_vars=cols_dict['Reference'], value_name='s_coeff',
                       var_name='Time interval')

    sns.kdeplot(data=ref_long, x='s_coeff', hue='Time interval', cut=0, ax=axs[2])
    axs[2].set_title('Reference-based approach', size='large')

    s_fig = plt.gcf()
    s_fig.savefig(s_dists, format='pdf', bbox_inches='tight')
    plt.close()

s_dists.close()

# Second figure: s of mutants along the sequence
# Using the T0toT2 selection coefficients
os.chdir('s_coefficients')
s_unique = glob.glob('*1n.csv')
os.chdir('..')

s_single = pd.read_csv(f's_coefficients/{s_unique[0]}')

# Filtering the data to remove the WT as well as all variants which had only one read or less at both T0 and T3
s_single = s_single[s_single['Genotype'] != 'WT'].copy().reset_index(drop=True)
s_single = s_single[(s_single['N_T0'] > 1) | (s_single['N_T3'] > 1)].copy().reset_index(drop=True)

# Extract position information for each variant
s_single['Position'] = s_single['Genotype'].apply(lambda x: int(x.split('-')[1].split('_')[0]))
s_single['Type'] = s_single['Genotype'].apply(lambda x: str(x.split('-')[0]))
s_single['Change'] = s_single['Genotype'].apply(lambda x: str(x.split('_')[2]))

loc_ready = s_single[['Library', 's_T0toT2_Time', 's_T0toT2_Total', 's_T0toT2_Ref', 'Position', 'Type',
                      'Change']].copy().infer_objects()

fig, axs = plt.subplots(3, 1, figsize=(42, 36))

sns.stripplot(data=loc_ready, x='Position', y='s_T0toT2_Time', hue='Type', ax=axs[0])
# Removing 'native_scale=True', because it is not supported by the version available on the server
sns.stripplot(data=loc_ready, x='Position', y='s_T0toT2_Total', hue='Type', ax=axs[1])
# Removing 'native_scale=True', because it is not supported by the version available on the server
sns.stripplot(data=loc_ready, x='Position', y='s_T0toT2_Ref', hue='Type', ax=axs[2])
# Removing 'native_scale=True', because it is not supported by the version available on the server

for ax in [axs[0], axs[1], axs[2]]:
    ax.axvline(x=71.5, c='red', linestyle='--')
    ax.axvline(x=225.5, c='red', linestyle='--')
    ax.axvline(x=148.5, c='grey', linewidth=3)
    ax.axhline(y=-0.05, c='black', linestyle='--')
    ax.axhline(y=0.05, c='black', linestyle='--')

    ax.set_xlim(left=1, right=373)
    ax.set_ylabel('Selection coefficient', fontsize=20)
    ax.set_xlabel('Position within the sequenced fragment (red dashed lines)', fontsize=20)

# Adding annotation for mutants with s>0.1
s_sub_time = loc_ready[abs(loc_ready['s_T0toT2_Time']) >= 0.05].copy().reset_index(drop=True)
for row in list(s_sub_time.index):
    mut_nat = s_sub_time.at[row, 'Change']
    point_xy = (s_sub_time.at[row, 'Position'], s_sub_time.at[row, 's_T0toT2_Time'])
    axs[0].annotate(f'{mut_nat}', point_xy, (2, 2), textcoords='offset points', fontsize=12)

s_sub_tot = loc_ready[abs(loc_ready['s_T0toT2_Total']) >= 0.1].copy().reset_index(drop=True)
for row in list(s_sub_tot.index):
    mut_nat = s_sub_tot.at[row, 'Change']
    point_xy = (s_sub_tot.at[row, 'Position'], s_sub_tot.at[row, 's_T0toT2_Total'])
    axs[1].annotate(f'{mut_nat}', point_xy, (2, 2), textcoords='offset points', fontsize=12)

s_sub_ref = loc_ready[abs(loc_ready['s_T0toT2_Ref']) >= 0.1].copy().reset_index(drop=True)
for row in list(s_sub_ref.index):
    mut_nat = s_sub_ref.at[row, 'Change']
    point_xy = (s_sub_ref.at[row, 'Position'], s_sub_ref.at[row, 's_T0toT2_Ref'])
    axs[2].annotate(f'{mut_nat}', point_xy, (2, 2), textcoords='offset points', fontsize=12)

# Adding titles for each plot
axs[0].set_title("Using the time-based approach", fontsize=24)
axs[1].set_title("Using the OD-based approach", fontsize=24)
axs[2].set_title("Using the reference-based approach", fontsize=24)

fig_1nt_map = plt.gcf()
fig_1nt_map.savefig(f"s_coefficients/Seq_Mapping_1n_{name_pre}_L{pool_n}.pdf", bbox_inches='tight')
