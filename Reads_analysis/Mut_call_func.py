#!/usr/bin/env python

import os
import pandas as pd
import numpy as np
import glob
import math
import scipy.stats as stats
import statsmodels.api as sm
from pandas.api.types import CategoricalDtype


def parse_needle_output(path, file_name, lib_name):
    # This function parses the needle alignments and extract the aligned reference sequence and query sequences.
    # It takes as input the path to the needle outputs, as well as the name of the current FASTA file and the unique
    # name of the corresponding library.
    # THIS FUNCTION ASSUMES THAT THE REFERENCE SEQUENCE IS WRITTEN BEFORE THE QUERY IN THE NEEDLE OUTPUT

    n_aligns = 0
    # counter for the number of alignments

    align_seqs_dict = {}
    # empty container that will hold the aligned sequences

    needle_file_path = path + '/' + file_name
    # path to the alignments

    with open(needle_file_path, 'r') as source:
        # open the alignment

        current_align = ''
        current_qseq = ''
        current_rseq = ''
        # empty container objects for data processing

        rseq_done = 0
        # counter for the number of alignments processed

        for line in source:
            # loop through the file

            if line.startswith('>>>'):
                # detect headers

                n_aligns += 1
                # increment alignment counter by one

                align_info = line.strip('>>>')
                # get alignment name, looks like this:  CDS_ref, 194 nt vs seq1;size=59000;length=214, 214 nt
                ref_id = align_info.split(',')[0]
                seq_name = align_info.split(' vs ')[1].split(';')[0]
                n_reads = align_info.split('size=')[1].split(';')[0]
                query_length = align_info.split('length=')[1].split(',')[0]

                align_name = seq_name + '_' + lib_name

                if n_aligns != 1:
                    # if this is not the first alignment

                    align_seqs_dict[current_align] = {
                        'Infos': {'N_reads': current_params[1], 'Length_query': current_params[2],
                                  'Ref_name': current_params[0], 'Align_score': current_score},
                        'Seq_query': current_qseq, 'Seq_ref': current_rseq}
                    # add the information on the previous alignment to the dict

                    current_align = align_name
                    current_params = [ref_id, n_reads, query_length]
                    # update the name for the new entry

                    current_qseq = ''
                    current_rseq = ''
                    current_score = np.NaN
                    # reset temporary variables

                    rseq_done = 0
                    # reset indicator for reference sequence extraction

                else:
                    current_align = align_name
                    current_params = [ref_id, n_reads, query_length]
                    # for the first sequence, just need to store the align name and features

            elif line.startswith('; sw_score:'):
                # extract the score of the alignment
                current_score = line.split(': ')[1]


            elif (line.startswith(';') == False) and (line.startswith('>') == False) and (
                    line.startswith('\n') == False) and (line.startswith('#') == False):
                # skip all the useless lines to process only the aligned sequences

                if rseq_done == 1:
                    current_qseq += line.strip('\n')
                    # if the reference seq is done (rseq = 1), add sequence to the query

                elif rseq_done == 0:
                    current_rseq += line.strip('\n')
                    # if the ref seq is not done, continue to update it

            elif line.startswith('#--') == True:
                align_seqs_dict[current_align] = {
                    'Infos': {'N_reads': current_params[1], 'Length_query': current_params[2],
                              'Ref_name': current_params[0], 'Align_score': current_score},
                    'Seq_query': current_qseq, 'Seq_ref': current_rseq}
                # update dict with info from the very last entry in the alignment sequence

            else:
                if rseq_done == 0 and current_rseq != '':
                    rseq_done = 1
                    # if the reference seq has been recorded, update value

    return align_seqs_dict, n_aligns


def mut_cat(muts_list):
    # Function to infer the mutation category codes according to a list of normalized mutation strings.
    # These strings are as follows: "Type-Start_Lnt_WTtoNew", eg "S-138_1nt_TtoC".

    cats_list = ""  # Initializing container for mutation categories
    number_string = ""  # And empty stings to construct it
    contig_string = ""

    totals = {'S': 0, 'D': 0, 'I': 0}  # Counters for total mutations and contiguous mutated sequences
    max_contigs = {'S': 0, 'D': 0, 'I': 0}

    for mut in muts_list:
        type_mut = mut[0]
        length_mut = int(mut.split('_')[1].split('nt')[0])

        totals[type_mut] += length_mut

        if length_mut > 1:
            max_contigs[type_mut] = max(length_mut, max_contigs[type_mut])

    # Once all mutations in the list have been processed, the list of mutation types can be produced
    names_dict = {'S': 'substitution', 'D': 'deletion', 'I': 'insertion'}

    for entry in totals.keys():
        mut_val = totals[entry]

        # First for number of mutated positions
        if mut_val == 0:
            continue

        elif mut_val > 0:

            if mut_val == 1:
                type_word = names_dict[entry]

            elif mut_val > 1:
                type_word = names_dict[entry] + 's'

            if mut_val < 5:
                n_mut = str(mut_val)

            elif mut_val >= 5:
                n_mut = "5+"

            # Assembling the string
            number_string = n_mut + ' ' + type_word

            # Next, for max length of contiguous mutated sequence
        contig_max = max_contigs[entry]

        if contig_max == 0:
            n_contig = ''
            word_contig = ''

        elif contig_max > 0:

            if contig_max >= 10:
                n_contig = "Large "
                word_contig = names_dict[entry]

            elif contig_max > 1:
                word_contig = ' contig. ' + names_dict[entry] + 's'

                if contig_max >= 4:
                    n_contig = "4+"

                elif contig_max < 4:
                    n_contig = str(contig_max)

            elif contig_max == 1:
                n_contig = "1"
                word_contig = ' contig. ' + names_dict[entry]

            contig_string = n_contig + word_contig

        # Adding the two strings to the enumeration (if needed)
        if len(number_string) > 0:
            cats_list += f'{number_string}, '

        if len(contig_string) > 0:
            cats_list += f'{contig_string}, '

        # Setting each string back to its original empty state
        number_string = ''
        contig_string = ''

    return cats_list


def find_mutations(in_path, out_path):
    # For a given path, runs the Needle alignment parser on each file.
    # The function then goes through the parsed alignments, file-by-file, to output mutations
    # information into a newly generated csv file.

    # Obtaining the list of files
    os.chdir(in_path)
    files_list = glob.glob('*.needle')

    # Going back to the initial working directory
    os.chdir('..')

    # Parsing files one at a time and outputting the mutation data
    for needle_file in files_list:
        dict_of_dfs = {}  # Empty dictionary that will contain all dfs for the current file until they are concatenated into one

        lib_name = needle_file.split('_aligned')[0]

        df_rows = []  # List which will contain a dictionary per entry in the Needle file, from which a df will be generated afterwards

        # Parsing the alignment file
        align_dict, align_count = parse_needle_output(in_path, needle_file, lib_name)

        for entry in list(align_dict.keys()):
            # loop through the parsed alignments
            n_reads = align_dict[entry]['Infos']['N_reads']
            length_q = align_dict[entry]['Infos']['Length_query']
            ref_name = align_dict[entry]['Infos']['Ref_name']
            align_score = align_dict[entry]['Infos']['Align_score'].strip('\n')

            # The alignment entry is then analysed for mutations

            # Initalization of temporary variables
            mut_type = ''
            mut_start = np.NaN
            mut_length = np.NaN
            pos_wt = ''
            pos_mutant = ''

            read_muts_list = []  # temporary holders for mutations found in the current sequence

            query_seq = align_dict[entry]['Seq_query'].upper()
            # aligned sequence of the strain

            align_ref = align_dict[entry]['Seq_ref'].upper()
            # aligned sequence of the reference

            if len(query_seq) != len(align_ref):
                raise Exception("The aligned query and reference sequences should have the same length!")

            # The position of all mismatches between the two sequences are fist identified
            query_array = np.array(list(query_seq))
            ref_array = np.array(list(align_ref))

            miss_indices = list(np.nonzero(ref_array != query_array)[0])

            if len(miss_indices) == 0:
                read_mut_list = []
                reads_muts_cats = ""

            elif len(miss_indices) > 0:

                miss_diffs = np.diff(
                    miss_indices)  # Compute differences between successive indices so that we can identify contiguous mutations

                types_list = []  # List of mismatch types

                to_skip = []  # List of alignment positions to skip

                # Mutation types are first compiled
                for position in miss_indices:

                    query_pos = query_seq[position]
                    ref_pos = align_ref[position]

                    # Testing first possibility: The muismatch is a substitution
                    if (query_pos != '-') and (ref_pos != '-'):
                        mut_type = 'S'

                    # Testing second possibility: The mismatch is a deletion
                    elif (query_pos == '-') and (ref_pos != '-'):
                        mut_type = 'D'

                    # Testing third possibility: The mismatch is an insertion
                    elif (query_pos != '-') and (ref_pos == '-'):
                        mut_type = 'I'

                    types_list += mut_type  # Adding the mutation type to the list

                # Then, mutations are enumerated
                ref_indices = miss_indices.copy()  # List of reference-based index (ajusted for insertions in the query)
                ref_adjust = 0  # Temporary variable for adjustment of reference coordinates in case of insertions

                for mut_num in range(len(miss_indices)):

                    if mut_num in to_skip:
                        continue

                    elif mut_num not in to_skip:
                        mut_type = types_list[mut_num]
                        mut_start = miss_indices[mut_num]
                        mut_length = 1
                        pos_WT = align_ref[miss_indices[mut_num]]
                        pos_mutant = query_seq[miss_indices[mut_num]]

                        # Prevent an error in case this is the last mutated position
                        if mut_num == (len(miss_indices) - 1):
                            if mut_type == 'I':
                                start_ref = ref_indices[mut_num]  # No +1 because of the adjustment for the insertion

                            elif mut_type != 'I':
                                start_ref = ref_indices[mut_num] + 1  # The +1 is added to make the coordinates 1-based

                            read_muts_list += [f'{mut_type}-{start_ref}_{mut_length}nt_{pos_WT}to{pos_mutant}']

                        elif mut_num != (len(miss_indices) - 1):
                            type_next = types_list[mut_num + 1]
                            dist_next = miss_diffs[mut_num]

                            init_while = mut_num  # Mutation number at the beginning of the while

                            # To be able to ajust reference coordinates in case of insertion
                            if mut_type == 'I':
                                ins_end = init_while
                                ins_length = 1

                            while (mut_type == type_next) and (dist_next == 1):
                                init_while += 1

                                if mut_type == 'I':
                                    ins_end = init_while
                                    ins_length += 1

                                mut_length += 1
                                to_skip += [init_while]
                                pos_WT += align_ref[miss_indices[init_while]]
                                pos_mutant += query_seq[miss_indices[init_while]]

                                # Prevent an error in case the while has already reached the final mutation
                                if init_while == (len(miss_indices) - 1):
                                    dist_next = 200

                                elif init_while < (len(miss_indices) - 1):
                                    type_next = types_list[init_while + 1]
                                    dist_next = miss_diffs[init_while]

                            # Make adjustments in case this is the last base of an insertion
                            if (mut_type == 'I') and ((type_next != 'I') or (dist_next > 1)):
                                ref_adjust += ins_length

                                ref_indices = ref_indices[: ins_end + 1] + [pos - ins_length for pos in
                                                                            ref_indices[ins_end + 1:]]

                            # Adding the mutation to the list
                            if (mut_type == 'I') and (ref_indices[mut_num] == 0):
                                start_ref = 0
                                # Special case of an insertion at the very beginning

                            elif (mut_type == 'I'):
                                start_ref = ref_indices[mut_num]
                                # Insertions are defined according to the base to the right of which they occur, thus -1 to the
                                # coordinate. Because +1 is also done to convert into 1-based, the original value is kept.

                            else:
                                start_ref = ref_indices[mut_num] + 1

                            read_muts_list += [f'{mut_type}-{start_ref}_{mut_length}nt_{pos_WT}to{pos_mutant}']

            # Generating the corresponding mutation category codes
            reads_muts_cats = mut_cat(read_muts_list)

            # Filling the list of dictionaries with data from the current alignment
            reads_mut_str = ', '.join(read_muts_list)

            current_dict = {'Name': entry, 'Library': lib_name, 'Ref_seq': ref_name, 'Variant_seq': query_seq,
                            'Mutations': reads_mut_str,
                            'Mut_categories': reads_muts_cats, 'N_reads': n_reads, 'Query_length': length_q,
                            'Align_score': align_score}
            df_rows.append(current_dict)

        # Once the current Needle file has been fully parsed, the dataframe is saved
        df_mutations = pd.DataFrame(df_rows)
        df_mutations.to_csv(f'{out_path}/mutations_{lib_name}_{ref_name}.csv', index=False)

    # No return statement


def codes_to_codons(mut_codes, seq_start, offset):
    """Small function to convert a list of mutation codes into a list of mutated codons.
    mut_codes = list of mutation codes to be converted
    seq_start = codon of the first position in the reference sequence (one-based)
    offset =  position within the first codon at which the reference sequence begins
    if offset = 0: first position of the codon
    if offset = 1: second position of the codon
    if offset = 2: third position of the codon"""
    codons_list = []

    for mut_code in mut_codes:
        begin_mut = int(mut_code.split('_')[0].split('-')[1])
        end_mut = begin_mut + (int(mut_code.split('_')[1].split('nt')[0]) -1)

        start_codon = math.ceil((begin_mut + offset)/3) + (seq_start - 1)
        end_codon = math.ceil((end_mut + offset)/3) + (seq_start - 1)

        if start_codon != end_codon:

            if end_codon - start_codon > 1:
                codons_list += list(range(start_codon, end_codon + 1))

            elif end_codon - start_codon == 1:
                codons_list += [start_codon, end_codon]

        elif start_codon == end_codon:
            codons_list += [start_codon]

    return codons_list


def mutants_rel(data_dir, file_end, out_dir, mut_n, seq_type, offset=0, min_length=350, max_length=385):
    """**Final version which compiles denominators for frequencies only once all unwanted reads have been removed**
    **New version allowing to also use the "Total" approach for the downstream calculation of selection coefficients.**
    Does not perform any abundance-based filtering of the variants. Only reads shorter than the corresponding "length"
    threshold are discarded. Function to return a dataframe compiling frequencies (relative to wt) of all n mutants
    (single, double, etc.) for all DMS libraries present in the directory. The file_end argument is used to specify
    a string specific to the end of the relevant file names.
    The mut_n parameters allows to specify which kind of variants should be compiled.
    Single mutants (at the nucleotide level are obtained with mut_n=1,
    while mut_n=2 specifies double mutants.

    The last two arguments before min and max sequence lengths allow to filter mutants differently according to
    'seq_type' (either 'prom' for promoter fragments or 'CDS' for coding sequence fragments). Argument 'offset'
    is only used in the 'CDS' case. It specifies the position within the first codon at which the reference sequence
    begins. If 0, the reference starts on the first position of a codon, while 1 and 2 respectively refer to the second
    and third positions of the first codon."""

    # Generating a model df and an empty df to receive all the data
    df_model = pd.DataFrame(columns=['Library', 'Timepoint', 'Mutations',
                                     'N_reads', 'Total_reads', 'Rel_tot', 'Rel_wt', 'WT_reads'])

    os.chdir(data_dir)
    files_list = glob.glob(f"*{file_end}")

    # Going back to the initial working directory
    os.chdir('..')

    for mut_file in files_list:
        full_df = pd.read_csv(f'{data_dir}/{mut_file}')

        # Prior to the filtering of mutations, reads are filtered according to length
        full_df = full_df[full_df['Query_length'] >= min_length].reset_index(drop=True)
        full_df = full_df[full_df['Query_length'] <= max_length].reset_index(drop=True)

        full_df = full_df[['Library', 'Ref_seq', 'Mutations', 'Mut_categories', 'N_reads']].copy()
        full_df[['Mutations', 'Mut_categories']] = full_df[['Mutations', 'Mut_categories']].fillna('WT')
        df_agg = full_df.groupby(by=['Library', 'Ref_seq', 'Mutations', 'Mut_categories'], as_index=False).sum()
        df_agg = df_agg.sort_values(by='N_reads', ascending=False).reset_index(drop=True)

        # Separating the wild-type from the mutants and obtaining wt-relative frequencies for all variants
        if df_agg[df_agg['Mutations'] == 'WT'].shape[0] != df_agg[df_agg['Mut_categories'] == 'WT'].shape[0]:
            raise Exception(f"There are unexpected 'WT' labels in the mutation data!")

        wt_subset = df_agg[df_agg['Mutations'] == 'WT'].copy().reset_index(drop=True)
        wt_total = np.sum(wt_subset['N_reads'])
        mutants_subset = df_agg[df_agg['Mutations'] != 'WT'].copy().reset_index(drop=True)

        mutants_subset['Rel_wt'] = mutants_subset['N_reads'] / wt_total
        wt_subset['Rel_wt'] = wt_subset['N_reads'] / wt_total
        wt_subset.at[0, 'Mutations'] = frozenset({'WT'})
        wt_subset = wt_subset.copy().drop(['Mut_categories'], axis='columns')

        # Reformatting the 'Mutations' and 'Mut_categories' columns into frozensets
        mutants_subset['Mutations'] = mutants_subset['Mutations'].apply(lambda x: frozenset(x.split(', ')))
        mutants_subset['Mut_categories'] = mutants_subset['Mut_categories'].apply(lambda x: frozenset(x.strip(', ').split(', ')))

        # Then, variants are filtered. This is done differently for promoter and CDS mutations
        if seq_type == 'prom':
            # Selecting variants where the total number of mutated nucleotides equals mut_n
            mutants_subset['n_mutated'] = mutants_subset['Mutations'].apply(lambda x: np.sum(np.array([float(mut.split('_')[1].split('nt')[0]) for mut in list(x)])))
            mutants_subset['n_mutated'] = mutants_subset['n_mutated'].apply(lambda x: x if x == mut_n else np.NaN)
            mutants_subset = mutants_subset.dropna(subset=['n_mutated']).reset_index(drop=True)

        elif seq_type == 'CDS':
            # Only keeping variants when the mutated positions are in mut_n distinct codons
            mutants_subset['codons_set'] = mutants_subset['Mutations'].apply(codes_to_codons, args=(1, offset))
            # Since the codon number won't be kept, an arbitrary value of 1 is passed for the seq_start argument
            mutants_subset['n_mutated'] = mutants_subset['codons_set'].apply(lambda x: len(set(x)))
            mutants_subset['n_mutated'] = mutants_subset['n_mutated'].apply(lambda x: x if x == mut_n else np.NaN)
            mutants_subset = mutants_subset.dropna(subset=['n_mutated']).reset_index(drop=True)
            mutants_subset = mutants_subset.drop(['codons_set'], axis='columns')

        else:
            raise Exception(f"The only two valid values from seq_type are 'prom' and 'CDS'.")

        # Any identical mutations are aggregated (should not be necessary, but just in case)
        mut_freqs_agg = mutants_subset.copy().drop(['Mut_categories'], axis='columns')
        mut_freqs_agg['Mutations'] = mut_freqs_agg['Mutations'].apply(lambda x: frozenset(x))
        mut_freqs_agg = mut_freqs_agg.groupby(by=['Library', 'Ref_seq', 'Mutations', 'n_mutated'], as_index=False).sum().reset_index(drop=True)

        # The WT data is added back
        mut_freqs_agg = pd.concat([mut_freqs_agg, wt_subset]).reset_index(drop=True)

        # Reformatting the library and reference annotations
        mut_freqs_agg['Timepoint'] = 'None'
        mut_freqs_agg['Timepoint'] = mut_freqs_agg['Library'].apply(lambda x: x.split('_')[2])

        # Keeping the library and the timepoint as variables
        lib_num = mut_freqs_agg.at[0, 'Library']
        file_timepoint = mut_freqs_agg.at[0, 'Timepoint']

        # Adding total number of reads, total frequency relative to total and wild-type abundance (in reads)
        mut_freqs_agg['Total_reads'] = np.sum(mut_freqs_agg['N_reads'])
        mut_freqs_agg['Rel_tot'] = mut_freqs_agg['N_reads'] / np.sum(mut_freqs_agg['N_reads'])
        mut_freqs_agg['WT_reads'] = wt_total

        # Adding the data to the df of frequencies (according to the existing columns in each case
        data_ready = mut_freqs_agg[['Library', 'Timepoint', 'Mutations', 'N_reads', 'Total_reads',
                                    'Rel_tot', 'Rel_wt', 'WT_reads']].copy()

        # All columns which are frozensets are converted to strings.
        # At this point, each row contains one of the selected mutant genotypes
        data_ready['Mutations'] = data_ready['Mutations'].apply(lambda x: ', '.join(list(x)))

        data_ready.to_csv(f"{out_dir}/rel_freqs_{lib_num}_{file_timepoint}_{mut_n}n.csv", index=False)

    # No return statement


# Functions to compute and compile selection coefficients from the frequencies
def compile_freqs(freqs_template, files_list, lib_id):
    """Function to compile the frequencies of genotypes through all timepoints for one library.
       It should receive one of the source genotype frequencies files as template (freqs_template),
       which will be used to generate a new df with the same columns. The other parameters are:
       files_list = List of files, each containing genotype frequencies at one timepoint for the same library.
                    This file list is given as a list of paths produced by glob.glob()
       lib_id = Unique identifiers for the library/pool for which the timepoints are being processed"""

    # Preparing the dataframe from the template
    freqs_all = pd.DataFrame(columns=freqs_template.columns)

    # Dictionaries (to be filled further) containing WT and Total read counts by timepoint are defined
    wt_reads = {}
    total_reads = {}

    for freqs_df in files_list:
        time_freqs = pd.read_csv(f'{freqs_df}')

        if len(list(time_freqs['Timepoint'].unique())) != 1:
            raise Exception(f'There should be only one timepoint in {freqs_df}')

        current_timepoint = time_freqs.at[0, 'Timepoint']  # Timepoint of the current file
                                                           # (which should contain data for ONLY ONE timepoint)

        # The 'Library' column of the current file is reformatted, so that the timepoint annotation is absent
        time_freqs['Library'] = lib_id

        # While each file is imported, dictionaries of WT and Total frequencies through time are built
        wt_abun = time_freqs['WT_reads'].values[0]
        wt_reads[current_timepoint] = wt_abun

        total_abun = np.sum(time_freqs['N_reads'])
        total_reads[current_timepoint] = total_abun

        # Combining the frequency data for each timepoint into one df
        freqs_all = pd.concat([freqs_all, time_freqs]).reset_index(drop=True)

    freqs_all = freqs_all.rename(columns={'Mutations': 'Genotype'})

    # Both the dataframe and the two dictionaries are returned
    return freqs_all, wt_reads, total_reads


def get_time_dict(gen_df, lib_col, time_suffix, timepoints):
    """Function to generate a dictionary of time (in hours) or generations for each timepoint of each library.
    This dict will have a library -> Timepoint -> Time/Generations structure.
    gen_df = Dataframe containing the data to be parsed
    lib_col = Name of the column containing the library identifiers
    time_suffix = Suffix of the column names containing the relevant time/generations data. All data-containing
    columns should be named as Tn_{time_suffix}. For example: T1_gen_ref
    timepoints = List of timepoints (strings, eg T1, T2...), excluding the T0"""

    time_dict = {}

    for lib_id in gen_df[lib_col].unique():
        time_dict[lib_id] = {}
        df_subset = gen_df[gen_df[lib_col] == lib_id].copy()

        time_dict[lib_id]['T0'] = 0

        for time_id in timepoints:
            time_dict[lib_id][time_id] = df_subset[f'{time_id}_{time_suffix}'].values[0]

    return time_dict


def compute_s(freqs_df, time_dict, dict_type='full', final_time_id='T3', with_n_reads=False):
    """General function for the calculation of a selection coefficient from genotype frequencies through time.
    The calculation is performed relative to the WT genotype.
    The WT-relative frequencies over the sampling points should be given as a dataframe with the following columns:
    ['Timepoint', 'Library', 'Genotype', 'N_reads', 'Rel_wt']. This dataframe should only contain the data for the
    timepoints to use in the current selection coefficient calculation.
    The time_dict argument should be a dictionary specifying the (cumulative) number of generations (or time in hours)
    for each of the successive sampling timepoints. Whether generations or actual time units are used does not
    impact the calculation - only the units of the returned selection coefficient. The dictionary can either be
    'simplified' and have only one level of timepoint-to-time mappings, or be 'full' and define all the mappings
    individually for each library (such that the values specific to the relevant library are used in the calculation.
    The 'with_n_reads' option allows to specify whether the calculation should also return the T0 and final reads
    and wt-relative frequencies."""

    # Making sure that the timepoints are ordered
    n_range = int(final_time_id.split('T')[1]) + 1
    times_list = [f'T{num_id}' for num_id in range(n_range)]
    timepoints_type = CategoricalDtype(categories=times_list, ordered=True)

    freqs_df = freqs_df.astype({'Timepoint': timepoints_type})
    freqs_df = freqs_df.sort_values(by='Timepoint', ascending=True).reset_index(drop=True)

    # Separating the T0 from the rest of the data
    data_T0 = freqs_df[freqs_df['Timepoint'] == 'T0'].copy().reset_index(drop=True)
    data_T0 = data_T0[['Library', 'Genotype', 'N_reads', 'Rel_wt']].copy()
    data_T0 = data_T0.rename({'N_reads': 'Reads_T0', 'Rel_wt': 'Rel_T0'}, axis=1)

    # Merging the data
    freqs_wt = freqs_df.copy()
    freqs_wt = freqs_wt.rename({'Rel_wt': 'Rel_timepoint'}, axis=1)
    rel_merged = pd.merge(freqs_wt, data_T0, how='inner', on=['Genotype', 'Library']).reset_index(drop=True)

    # Keeping the genotype for later
    all_genotypes = rel_merged['Genotype'].unique()
    lib_id = rel_merged.at[0, 'Library']
    if all_genotypes.shape[0] != 1:
        raise Exception(f'Less or more than one genotype detected in a selection coefficient calculation ({lib_id})')
    else:
        current_genotype = all_genotypes[0]

    # Adding the corresponding time variable (generations or hours) for each frequency measurement
    rel_merged['Sample_time'] = np.NaN

    for row in list(rel_merged.index):
        timepoint = rel_merged.at[row, 'Timepoint']
        lib_name = rel_merged.at[row, 'Library'].split('_')[1]

        if dict_type == 'simplified':
            rel_merged.at[row, 'Sample_time'] = time_dict[timepoint]

        elif dict_type == 'full':
            rel_merged.at[row, 'Sample_time'] = time_dict[lib_name][timepoint]

    # Conversion of relative frequencies to log ratios
    rel_merged = rel_merged.infer_objects()
    rel_merged['Rel_T0'] = np.log(rel_merged['Rel_T0'])
    rel_merged['Rel_timepoint'] = np.log(rel_merged['Rel_timepoint'])

    # Extracting the vectors for the OLS linear regression
    Y_current = rel_merged['Rel_timepoint'].copy()
    X_current = rel_merged['Sample_time'].copy()
    X_current = sm.add_constant(X_current)

    # Performing the regression
    try:
        model_current = sm.OLS(Y_current, X_current)
        fit_current = model_current.fit()
        s_wt = fit_current.params['Sample_time']
        intercept_fit = fit_current.params['const']

    except ValueError:
        # In case a regression is attempted from a dataframe which does not contain any observations
        print(f'ValueError raised when defining the OLS model for genotype {current_genotype} in {lib_id}')
        s_wt = np.NaN
        intercept_fit = np.NaN

    # Computing the other version of the selection coefficient, from single-timepoint log frequencies
    if final_time_id in rel_merged['Timepoint'].unique():
        final_index = rel_merged[rel_merged['Timepoint'] == final_time_id].index.values[0]
        time_diff = rel_merged.at[final_index, 'Rel_timepoint'] - rel_merged.at[final_index, 'Rel_T0']
        ratio_s = time_diff / rel_merged.at[final_index, 'Sample_time']

    elif final_time_id not in rel_merged['Timepoint'].unique():
        ratio_s = np.NaN

    # Preparing the list to be returned
    if not with_n_reads:
        return_list = [s_wt, intercept_fit, ratio_s]

    elif with_n_reads:
        T0_reads = data_T0.at[0, 'Reads_T0']
        T0_freq = data_T0.at[0, 'Rel_T0']

        # Test if the final timepoint exists (to avoid breaking the pipeline when a sequencing library is missing)
        # This happens for library 24, for which the T3 PCR di not work and was not included in the final run
        try:
            Tend_reads = freqs_df[freqs_df['Timepoint'] == final_time_id].copy().reset_index(drop=True).at[0, 'N_reads']
            Tend_freq = freqs_df[freqs_df['Timepoint'] == final_time_id].copy().reset_index(drop=True).at[0, 'Rel_wt']

        except KeyError:
            print(f'{final_time_id} data does not exist for {lib_name}')
            Tend_reads = np.NaN
            Tend_freq = np.NaN

        return_list = [s_wt, intercept_fit, T0_reads, T0_freq, Tend_reads, Tend_freq, ratio_s]

    return return_list


def loop_s_calc(freqs_all, wt_dict, lib_id, timepoints, time_dict, dict_type='full', with_n_reads=False,
                s_suffix='Total'):
    """Function to loop the calculation of selection coefficients (compute_s() function) over variants. Also
    handles to computation over any possible T0-to-Tn time intervals. This function needs a dataframe of genotype
    frequencies (produced by the compile_freqs() function), all the arguments for the compute_s() function (except
    final_time_id, which is inferred within this function) and an additional 's_suffix' specifying which identifier
    to add to the newly produced column of selection coefficients (among those used in compile_s when generating the
    template dataframe of selection coefficients). A list of sampled timepoints (given as 'timepoints') is also
    needed. In addition, the 'wt_dict' argument is a dictionary of WT read counts (with timepoints as keys) generated
    by the compile_freqs() function"""

    # Generating the dictionary of time intervals from the list of timepoints, as within the compile_s() function
    timepoints_dict = {}
    time_ids_list = []  # list of T0toTn id, in order
    for t_final in range(1, len(timepoints)):
        time_string = f'{timepoints[0]}to{timepoints[t_final]}'
        time_ids_list += [time_string]
        times_list = timepoints[0:t_final+1]
        timepoints_dict[time_string] = times_list

    final_time = timepoints[t_final]  # Keeping the last t_final, for later

    # From the freqs_all dataframe combining the frequencies through all timepoints for all variants in the current
    # library, all variants detected at T0 are first extracted
    var_T0 = list(freqs_all[freqs_all['Timepoint'] == 'T0']['Genotype'].unique())

    # Constructing a dataframe which will be filled with selection coefficients
    s_variables = [f's_{time_id}_{s_suffix}' for time_id in time_ids_list]
    int_variables = [f'intercept_{time_id}_{s_suffix}' for time_id in time_ids_list]
    ratio_variables = [f'ratio_{time_id}_{s_suffix}' for time_id in time_ids_list]

    if with_n_reads:
        df_columns = (['Library', 'Genotype'] + s_variables + int_variables + ratio_variables +
                      ['N_T0', 'freq_T0', f'N_{final_time}', f'freq_{final_time}'])
        loop_df = pd.DataFrame(columns=df_columns)

    elif not with_n_reads:
        df_columns = ['Library', 'Genotype'] + s_variables + int_variables + ratio_variables
        loop_df = pd.DataFrame(columns=df_columns)

    loop_df['Genotype'] = var_T0
    loop_df.index = var_T0
    loop_df['Library'] = lib_id

    # Then, all the selection coefficient calculations (one per time_id) are performed for each variant
    for var_code in var_T0:
        var_subset = freqs_all[freqs_all['Genotype'] == var_code].copy().reset_index(drop=True)

        # For each of the T0toTn time intervals considered, the data for the current variant is extracted - adding
        # a one when the variant is not detected -, and the selection coefficient calculation is performed
        for time_id in time_ids_list:
            to_keep = timepoints_dict[time_id]  # Time points to keep

            # Selecting only the necessary data
            freqs_current = var_subset[var_subset['Timepoint'].isin(to_keep)].copy().reset_index(drop=True)

            # For any variant which does not exist at a later timepoint than T0, a read count of 1 is added
            count_1_added = False  # So that 1 is only added for the first missing observation

            for timepoint in to_keep:
                points_n = freqs_current[freqs_current['Timepoint'] == timepoint]['N_reads'].values.shape[0]
                # Zero if the variant has not been detected at the current timepoint
                if points_n == 0 and not count_1_added:
                    df_toadd = pd.DataFrame(columns=freqs_current.columns)

                    # This assumes a dataframe with the following columns:
                    # ['Library', 'Timepoint' 'Genotype', 'N_reads', 'Rel_wt']
                    # The WT may not exist for the current timepoint. If it is the case, an error message is printed
                    # instead. This should only happen in rare instances where a timepoint is entirely missing from the
                    # sequencing run.
                    try:
                        df_toadd.loc[0, :] = [lib_id, timepoint, var_code, 1, 1/wt_dict[timepoint]]

                    except KeyError:
                        print(f'WT not found at {timepoint} in library {lib_id}')

                    freqs_current = pd.concat([freqs_current, df_toadd]).reset_index(drop=True)

            # If the function has been called with with_n_reads=True AND the current T0toTn interval is the last one,
            # the selection coefficient calculation is performed with with_n_reads=True
            if with_n_reads and time_id.split('to')[1] == final_time:
                s_current, int_current, T0_reads, T0_freq, Tend_reads, Tend_freq, ratio_s = compute_s(freqs_current,
                                                                                                      time_dict,
                                                                                                      dict_type,
                                                                                                      final_time_id=final_time,
                                                                                                      with_n_reads=True)
                loop_df.at[var_code, f's_{time_id}_{s_suffix}'] = s_current
                loop_df.at[var_code, f'intercept_{time_id}_{s_suffix}'] = int_current
                loop_df.at[var_code, f'ratio_{time_id}_{s_suffix}'] = ratio_s
                loop_df.at[var_code, 'N_T0'] = T0_reads
                loop_df.at[var_code, 'freq_T0'] = T0_freq
                loop_df.at[var_code, f'N_{final_time}'] = Tend_reads
                loop_df.at[var_code, f'freq_{final_time}'] = Tend_freq

            else:
                s_current, int_current, ratio_s = compute_s(freqs_current, time_dict, dict_type,
                                                            final_time_id=time_id.split('to')[1], with_n_reads=False)

                loop_df.at[var_code, f's_{time_id}_{s_suffix}'] = s_current
                loop_df.at[var_code, f'intercept_{time_id}_{s_suffix}'] = int_current
                loop_df.at[var_code, f'ratio_{time_id}_{s_suffix}'] = ratio_s

    return loop_df


def compile_s(data_dir, file_end, out_dir, timepoints, total_dict, ref_dict, time_dict, dicts_type, min_reads=None):
    """Function to handle the computation/compilation of selection coefficients over all possible time intervals
       for one library sampled at multiple timepoints.
       All selection coefficients are computed in three ways: 1) Using the total generations (from OD),
       2) Using the duration (in hours) of each passage and 3) Using reference-based generations
       (computed from the duration of the passages)."""

    # Generating the files list by looking into the data_dir
    files_list = glob.glob(f"{data_dir}/*{file_end}")

    # Reading the beginning of the first file, only to generate a model df with the same columns
    freqs_template = pd.read_csv(f'{files_list[0]}', nrows=5)

    # Keeping lib_id and mut_n info
    lib_id = freqs_template.at[0, 'Library'][:-3]  # The timepoint annotation is removed ('_T0')
    mut_n = file_end.split('.csv')[0]

    # Combining the frequencies for the timepoints, while also compiling dictionaries of
    # WT and Summed read counts through the sampling points
    freqs_all, wt_reads, summed_reads = compile_freqs(freqs_template, files_list, lib_id)
    freqs_all = freqs_all.drop(['Total_reads', 'Rel_tot', 'WT_reads'], axis='columns')

    # If a min_reads threshold has been specified, the dataframe is filtered accordingly before the loop_s_calc calls
    # are made. Because frozensets were used to aggregate within the mutants_rel() function, there should be only one
    # possible "Genotype" string for each genotype. We can thus use strings here, without converting back to frozensets
    if min_reads:
        df_T0 = freqs_all[freqs_all['Timepoint'] == timepoints[0]].copy().reset_index(drop=True)
        T0_tokeep = df_T0[df_T0['N_reads'] >= min_reads]['Genotype'].to_list()
        df_Tfinal = freqs_all[freqs_all['Timepoint'] == timepoints[-1]].copy().reset_index(drop=True)
        Tfinal_tokeep = df_Tfinal[df_Tfinal['N_reads'] >= min_reads]['Genotype'].to_list()
        muts_to_calc = set(T0_tokeep) | set(Tfinal_tokeep)
        freqs_all = freqs_all[freqs_all['Genotype'].isin(muts_to_calc)].copy().reset_index(drop=True)

    # Calling loop_s_calc to perform the three types of selection coefficient calculations on all variants in the
    # current library. The with_n_reads=True option only needs to be used once, as the added data will be the same.
    s_total_df = loop_s_calc(freqs_all, wt_reads, lib_id, timepoints, total_dict, dict_type=dicts_type,
                             with_n_reads=True, s_suffix='Total')
    s_ref_df = loop_s_calc(freqs_all, wt_reads, lib_id, timepoints, ref_dict, dict_type=dicts_type,
                           with_n_reads=False, s_suffix='Ref')
    s_time_df = loop_s_calc(freqs_all, wt_reads, lib_id, timepoints, time_dict, dict_type=dicts_type,
                            with_n_reads=False, s_suffix='Time')

    # Merging the three dataframes into one
    merge_1 = pd.merge(s_total_df, s_ref_df, on=['Library', 'Genotype'], how='outer')
    merge_final = pd.merge(merge_1, s_time_df, on=['Library', 'Genotype'], how='outer')

    # Saving the data
    merge_final.to_csv(f"{out_dir}/s_all_{lib_id}_{mut_n}.csv", index=False)

    # No return statement
