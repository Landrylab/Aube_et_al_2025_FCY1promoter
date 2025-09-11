# Aube_et_al_2025_FCY1promoter
All code used to analyze bulk competition and sort-seq data in Aubé et al. (2025).

The repository is split into two levels: 1) Reads_analysis and 2)Data_analysis.
The first one contains the scripts which were used to obtain read counts and/or selection coefficients from raw sequencing data. This was performed on a server.
The second one contains all jupyter notebooks which were run locally to perform additional data analysis and generate figures used in the manuscript.

**All data files necessary to run the scripts are provided in the corresponding Zenodo repository (https://doi.org/10.5281/zenodo.17088129), also split in "Reads_analysis" and "Data_analysis".**

## A) Reads_analysis

### Common dependencies
* create_venv_FCY1_300.txt -- Command lines which were used to create the virtual environment within which all Python script from this section were run
* requirements_FCY1_300.txt -- Requirement file listing the version of Python modules loaded within this virtual environment
* Mut_call_func.py -- Custom Python functions which are called within the different analysis scripts

### A1) Bulk competition

#### For preparation (performed locally before setting up the analysis pipeline on the server)
* Libs_info_format.py -- Short script used to reformat the names of 150 paired-end sequencing libraries
* generations_prep.py -- Short script to compute total and (OD-based) reference-based generations for each sampling timepoint of each replicate culture

#### For analysis of sequencing data
* Reads_analysis.py -- Main Python script to compute selection coefficients for all single-nucleotide variants from sequencing data. This script is designed to be run on one series of sequencing libraries (all timepoints T0 to T3 of the same culture).
* compile_s_data.py -- Short script to compile selection coefficients obtained from many cultures into a single file
* compile_mut_freqs.py -- Sort script to compile WT-relative frequencies (for each variant) obtained from many cultures into a single file

#### To run the above analysis across all samples (exact command lines used)
* job_FCY1_300_final_25-04-07.bash -- To run Reads_analysis.py on all samples of the 300 paired-end sequencing run
* job_FCY1_150_F3F4_25-04-07.bash -- To run Reads_analysis.py on all F3F4 (promoter) sequencing libraries from the 150 paired-end sequencing run.
* job_FCY1_150_CDS_PD_25-04-07.bash -- To run Reads_analysis.py on all CDS (FCY1) sequencing libraries from the 150 paired-end sequencing run. This only includes a subset of cultures, in which a set of previously characterized FCY1 mutations had been added as controls.
* compile_FCY1_300_25-04-08.bash -- To compile selection coefficients (compile_s_data.py) from the 300 paired-end run in one single file
* compile_FCY1_150_F3F4_25-04-08.bash -- To compile selection coefficients (compile_s_data.py) from the 150 paired-end F3F4 run in one single file 
* compile_FCY1_150_CDS_PD_25-04-08.bash -- To compile selection coefficients (compile_s_data.py) from the 150 paired-end CDS run in one single file
* freqs_FCY1_300_25-04-08.bash -- To compile WT-relative mutant frequencies (compile_mut_freqs.py) from the 300 paired-end run in one single file
* freqs_FCY1_150_F3F4_25-04-08.bash -- To compile WT-relative mutant frequencies (compile_mut_freqs.py) from the 150 paired-end F3F4 run in one single file
* freqs_FCY1_150_CDS_PD_25-04-08.bash -- To compile WT-relative mutant frequencies (compile_mut_freqs.py) from the 150 paired-end CDS run in one single file

### A2) Sort-seq

#### For analysis of sequencing data
* Mut_calling.py -- Main Python scripts, which computed read counts for all genotypes within each sequencing sample (bins from FACS)
* compile_reads.py -- Short Python script used to combined all counts into a single file

#### To run the above analysis across all samples (exact command lines used)
* job_SortSeq_25-02-18.bash -- To run Mut_calling.py on all sequencing libraries from the sort-seq experiment (300 paired-end)
* compile_sort_25-02-18.bash -- To compile per-bin read counts of genotypes (compile_reads.py) across all sequencing libraries

## B) Data_analysis
The different notebooks are listed in the order that they should be executed. Main, Extended Data and Supplementary Figures from the manuscript are listed in parenthesis next to the notebook where they are generated.

All notebooks were executed within Jupyter Lab, using the same "FCY1prom" kernel. File Creating_FCY1prom_venv.txt describe how to define this kernel. **The corresponding versions of all Python modules are listed in requirements_FCY1prom.txt.**

All paths within these analysis notebooks are relative, such that the content of the corresponding "Data_anlysis" folder of the Zenodo repository only needs to be copied within the same directory (where the notebooks are saved and Jupyter Lab is executed).

1. **Final_s_coefficients.ipynb** -- Abundance-based filtering of genotypes and identifications of samples to keep for fitness estimates shown in the paper.
2. **Aggregate_s_coeffs.ipynb** (Extended Data Fig 4; Supplementary Figs 9-10) -- Removal of outliers in selection coefficients and statistical analyses to identify significant fitness differences.
3. **Despres_2022_comparisons.ipynb** -- Comparisons of FCY1 spike-ins with the full FCY1 DMS from Després et al. (2022). Also compilation of variant effects from all single-nucleotide substitutions in the FCY1 gene.
4. **Fig4.ipynb** (Fig 3; Extended Data Fig 3) -- Generation of Fig. 3 from the manuscript.
5. **SortSeq_Analysis.ipynb** -- Abundance-based filtering of sort-seq data and calculation of corresponding Delta_Ref scores. Also included the selection of validation mutants and the generation of oligonucleotides for fusion PCR.
6. **Sort_calibration.ipynb** (Extended Data Figs 6,7,8) -- Validation and calibration of sort-seq using cytometry measurements for individually reconstructed validation mutants.
7. **Fig5.ipynb** (Fig 4) -- Generation of Fig. 4 from the manuscript.
8. **FCY1_fit_function_with_valid.ipynb** (Extended Data Fig 9) -- Fitting of expression-fitness curves and selection coefficients calculations for reconstructed validation mutants.
9. **Fig6.ipynb** (Fig 5) -- Generation of Fig 5 from the manuscript.
10. **Sim_muts_resist.ipynb** (Supplementary Figs 5,8) -- Simulations of competition experiments with "invisible" resistant mutants.
11. **Fig3.ipynb** (Supplementary Fig 1) -- Generation of Supplementary Fig 1 from the manuscript.
12. **S1Fig.ipynb** (Extended Data Fig 2) -- Generation of Extended Data Fig 1 from the manuscript.
13. **S2Fig.ipynb** (Supplementary Fig 2) -- Generation of Supplementary Fig 2 from the manuscript.
14. **S14Fig.ipynb** (Extended Data Fig 5) -- Generation of Extended Data Fig 5 from the manuscript.
15. **S3Fig.ipynb** (Supplementary Fig 3) -- Generation of Supplementary Fig 3 from the manuscript.
16. **S5Fig.ipynb** (Supplementary Fig 4) -- Generation of Supplementary Fig 4 from the manuscript.
17. **S7S8Figs.ipynb** (Supplementary Figs 6,7) -- Generation of Supplementary Figs 6 and 7 from the manuscript.
18. **FCY1_CSM_selection.ipynb** (Extended Data Fig 2) -- Script which was used to choose codon positions from the FCY1 DMS (Després et al., 2022) to include as controls in the bulk competition experiment.
