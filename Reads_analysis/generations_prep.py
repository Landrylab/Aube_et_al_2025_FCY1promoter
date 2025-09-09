#!/usr/bin/env python

# Import statements
import pandas as pd
import numpy as np
import scipy.stats as stats

# Short script to prepare the generations data (both total-based and reference-based)
gen_df = pd.read_csv("C:/Users/tiger/OneDrive - Université Laval/PhD/Chapitre 1/DMS_planning/paired_end_300_sequencing/OD_time_df.csv")

# Computing total (OD-based) generations
gen_df['T1_gen_tot'] = np.log2(gen_df['T1_OD'] / 0.01)
gen_df['T2_to_tot'] = np.log2(gen_df['T2_OD'] / 0.01)
gen_df['T2_gen_tot'] = gen_df['T2_to_tot'] + gen_df['T1_gen_tot']
gen_df['T3_to_tot'] = np.log2(gen_df['T3_OD'] / 0.01)
gen_df['T3_gen_tot'] = gen_df['T3_to_tot'] + gen_df['T2_gen_tot']

# Computing reference-based generations
# For 12 uM 5-FC, a growth rate of 0.178904 h^-1 is used (median measurements for AKD1123 in these conditions (Fig5)
# For the controls without 5-FC, we can consider a WT growth rate of 0,36 h^-1 (Does it even matter?)
# From this, the WT generation times are computed
gen_time_5FC = np.log(2) / 0.178904
gen_time_No_5FC = np.log(2) / 0.36

# Splitting the dataframe
no_5FC_subset = gen_df[gen_df['Library'] <= 11].copy().reset_index(drop=True)
with_5FC_subset = gen_df[gen_df['Library'] > 11].copy().reset_index(drop=True)

# Computing reference-based number of generations within each df
no_5FC_subset['T1_gen_ref'] = no_5FC_subset['T1_time'] / gen_time_No_5FC
no_5FC_subset['T2_gen_ref'] = no_5FC_subset['T2_time'] / gen_time_No_5FC
no_5FC_subset['T3_gen_ref'] = no_5FC_subset['T3_time'] / gen_time_No_5FC

with_5FC_subset['T1_gen_ref'] = with_5FC_subset['T1_time'] / gen_time_5FC
with_5FC_subset['T2_gen_ref'] = with_5FC_subset['T2_time'] / gen_time_5FC
with_5FC_subset['T3_gen_ref'] = with_5FC_subset['T3_time'] / gen_time_5FC

# Concatenating the two subset dfs
gen_all_df = pd.concat([no_5FC_subset, with_5FC_subset]).reset_index(drop=True)

# Preparing and exporting the final df
gen_df_ready = gen_all_df[['Library', 'T1_time', 'T2_time', 'T3_time', 'T1_gen_tot',
                           'T2_gen_tot', 'T3_gen_tot', 'T1_gen_ref', 'T2_gen_ref',
                           'T3_gen_ref']].copy()

gen_df_ready['Library'] = gen_df_ready['Library'].apply(lambda x: 'L' + str(x))

gen_df_ready.to_csv("C:/Users/tiger/OneDrive - Université Laval/PhD/Chapitre 1/DMS_planning/paired_end_300_sequencing/Gen_time_df.csv",
                    index=False)

# How to generate the dictionaries to be used during the calculation of selection coefficient?
# (needs to be copied in the Reads_analysis.py script

# We want a dict with structure Library -> Timepoint -> Generations or time
# This quick loop is defined as a function in Mut_call_func.py, and then used in Reads_analysis.py
time_dict = {}
for lib_id in gen_df_ready['Library'].unique():
    time_dict[lib_id] = {}
    df_subset = gen_df_ready[gen_df_ready['Library'] == lib_id].copy()

    time_dict[lib_id]['T0'] = 0
    time_dict[lib_id]['T1'] = df_subset['T1_time'].values[0]
    time_dict[lib_id]['T2'] = df_subset['T2_time'].values[0]
    time_dict[lib_id]['T3'] = df_subset['T3_time'].values[0]
