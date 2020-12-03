# HamiltonRuleMicrobiome_gitRepos
 
This repository holds data, output and analysis scripts for "Kin selection explains the evolution of cooperation in the gut microbiota" (Simonet & McNally, 2020).


## Data

**./metagenomes**:
- "manifest files": from HMP portal HMP portal, accessed April 2020, under Project > HMP, Body Site > feces, Studies > WGS-PP1, File Type > WGS raw sequences set, File format > FASTQ.
- MOCAT.cfg: for reads processing with MOCAT, using defaults settings (used by script 1.MIDAS_per_sample.sh)
- ``.metagenomes/metagenomes`` directory receives metagenomes and processed output when running samples processing script

**./patric**:
- receives amino-acids fasta files and features files from PATRIC database when running secretome and GO analyses (scripts 5_secretome.sh, 6_identify_social_GO.sh)

**./species_info_files**:
- genomes_info.txt, species_info.txt, species_tree.newick are directly taken from MIDAS db. Used by scripts during analysis.
- midas_tree_renamed.newick is species_tree.newick with species names made complete
- gram_profiles_db.txt is manually assembled table of species and gram profiles. Gram profiles determined from Bergey's Manual of Systematic Bacteriology or original publication describing the species
- node_genus.txt and species_plot_names.txt are used for figures production

**./sporulation_genes**:
- material taken from Brown et al (2016), used to compute sporulation scores in script 8_sporulation_scores.sh
- sporulation_genes_aa_fasta: aa fasta sequences retreived from NCBI for the 66 genes listed by Browne et al (2016).

## Output

Outputs from the scripts as detailed below. Outputs too large to be kept on github are not tracked. Repository structure can be used to re-run analysis in cloned repository anywhere.
- **analyses** contains all models output from R
- **tables** contains all tables output required to re-run analysis. "ANALYSIS_DATA_ASSEMBLED.txt" assembles all the other tables and is the table used in the statistical analysis scripts.
- **figures** are all figures parts used for main MS figures. Directly output from R scripts detailed below
- Other directories: receive output from early steps of the pipeline (samples processing, MIDAS's output, raw AA fasta and features files, blast output for computing sporulation scores etc.)

## Scripts

The Rmd scripts extracts results from the RData output to produce a pdf with all key results.

The scripts to re-run the analysis are in ./scripts/analysis

- Scripts 1-3 require working installation of MOCAT and MIDAS.
- Scripts 0-10 are for re-running analysis from scratch. They output data tablesin ./output/tables
- Script 11 is the statistical analysis in R.
	- Works with  dataframe ./output/tables/ANALYSIS_DATA_ASSEMBLED.txt
	- Model outputs saved in ./output/analyses (.RData files)
- Scripts 12-13 produce figures and tables as presented in the manuscript


Scripts work with user defined environmental variables:
- $programs_install_dir: pointing to those installations directory
- $user_dir and $local_project_dir: path where repository HamiltonRuleMicrobiome_gitRepos is cloned


### ./analysis

**0_samples_data.R**:
- parses the HMP portal manifest file to extarct unique sample/subject metagenome access link. This gives 251 samples. Out of these, 12 were filtered at different stage of the analysis (to see which ones and why, see sheet 1 of the SI dataset in ./output/tables, excel file)
- Output: list of samples and access links: ``./output/tables/HMP_first_visit_smaples.txt``

**1_midas_per_sample.sh**:
- downloada HMP metagenomes, processes with MOCAT, runs MIDAS. Task array jobs, using ``$SGE_TASK_ID`` env variable to process the ``$SGE_TASK_ID`` th sample. Run as many task as there are samples listed in ``./data/HMP_samples_first_visit_available`` (251). Requires user defined environemental variables ``$programs_install_dir`` and ``$user_dir``.
- run command on SGE will look like: ``qsub -v user_dir='path/where/gitrepo/is/cloned',programs_install_dir='/path/to/install/dir/of/MOCAT_and_MIDAS' -t 1-251 ./scripts/analysis/1_midas_per_sample.sh``
- Ouput:
	- in ``./data/metagenomes/metagenomes``, MOCAT quality statistics for each sample
	- in ``./output/midas/per-sample``, MIDAS's relative abundance tables for each sample


**2_midas_merge.sh**:
- merges species abundance files across samples (``merge_midas.py species``) to output relative abundance matrix and table of species prevalence. Then perform ``pooled-sample`` core-genome SNP calling across all core genomic sites for each species (``merge_midas.py snps``). Both commands take directory ``./output/midas/per_sample/`` as input. In SNP calling, ``--sample_depth (5.0)`` and ``--site_ratio (2.0)`` left to default. But set: ``--site_depth 5`` for better accuracy in allele frequency estimates, ``--snp_type any`` to take into account all variants even if observed at low frequency (very low frequency variants will simply have minimal impact on relatedness calculation), ``--site_prev 0.90`` for more relaxed core sites identification setting (95% prevalence is too stringent for species detected in a small number of hosts).
- run command on SGE will look like: ``qsub -v user_dir='path/where/gitrepo/is/cloned',programs_install_dir='/path/to/install/dir/of/MOCAT_and_MIDAS' ./scripts/analysis/2_midas_merge.sh``
- Output:
	- in ``./output/midas/merge/species_merged``, merged read counts, coverage, relative abundance and species prevalence tables
	- in ``./output/midas/merge/snps_merged``, merged allele frequencies tables


**3_midas_diversity.sh**:
- run command on SGE will look like:``qsub -v user_dir='path/where/gitrepo/is/cloned',programs_install_dir='/path/to/install/dir/of/MOCAT_and_MIDAS' -t1-XXX ./scripts/analysis/3_midas_diversity.sh``
- Output: in ``./output/midas/mdiversity``, within and across diversity estimates tables for each species



**4_relatedness.R**:
- computes relatedness from the diversity.within and diversity.between output of MIDAS. Also gather various diversity patterns statistics to compile an additional table fo SI dataset (SI table 7).
- Output:
	- relatedness table: ``/output/tables/relatedness.txt``
	- diversity table (for SI dataset table 7): ``/output/tables/diversity_patterns.txt``



**5_secretome.sh**:
- for the 101 species included in the analysis, downloads the representative genome (as defined in MIDASdb) AA fasta file fromt PATRIC, run PSORTb on them, compute secretome size. Code to install psortb (based on docker image) included. Part of the script runs in R.
- Output:
	- fasta and feature files in ``./data/patric``
	- table of secretome sizes in ``./output/tables/secretome.txt``



**6_identify_social_GO.sh**:
- creates a bacterial GO slim by annotated all representative genomes of MIDASdb with GO using Pannzer. Then use a list of keywords of bacterial social behaviour to browse this bacteria GO slim and retreive GO terms capturing cooperative behaviour categorised in 5 types of cooperation. The last part browsing the assembled GO  slim runs in R. Requires a working installation of Pannzer, specify path to installation directory with environemental variable programs_install_dir. The Pannzer step code is written to run iteratively through a loop working through all genomes one by one because large ammount of parrallels works launched on Pannzer at once (e.g. with task array job) seem to dramatically slow pannzer server (priority issue?) so the jobs never finish. When a single task is launched to Pannzer, it takes about 20min to process a genome. GO slim is simply the assembly of all unique GO terms identified across all >5900 genomes (taking all ARGOT ranks).
NOTE: the GO retreived when browsing the GO slim may differ as the GO.db gets updated. This may give slightly different social_go_list_wide.txt output. Current table (used for analysis) ran August 2019.
- Output:
	- All MIDAS species fasta files saved in ``./data/patric/fasta_all_midas`` (not git tracked)
	- All pannzer output saved in ``./output/pannzer`` (not git tracked)
	- Assembled GO slim: ``./output/tables/bacteria_go_slim.txt``
	- list of identified bacteria sociality keywords, identified across 10 top reviews in web of science search. Search string for Web of Science search  included in script. List of keywors in: ``./output/tables/bacteria_social_keywords.txt``
	- Table of all potential bacterial sociality GOs identified by raw keyword matching  + inclusion of direct parent and all child terms: ``./output/tables/social_go_list_wide.txt``
	- Manual curation of that table, with comments of decisions: ``./output/tables/social_go_list_curation.xls``
	- Final retained set of bacteria social GO terms: ``./output/tables/social_go_list_final.xls``


**7_go_categories.R**:
- computes the measure of cooperation based on GO annotation for each species, with the different cooperation split between the 5 cooperation categories. Runs in R. Uses the social_go_list_final.xls and the GO annotations (pannzer output) as input.
Script update: following review, script now also output table with per-gene break down details of GO and secretome annotations.
- Output:
	- table of measure of cooperation based on GO categories ``/output/tables/go_cooperation_categories.txt`` 
	- per-gene break down of GO categories + secretome: ``/output/tables/per_gene_annotation.txt``


**8_sporulation_scores.sh**:
- computes sporulation scores for the 101 species following method from Browne et al. We gathered (from NCBI) the AA fasta sequences of the 66 sporulation genes listed by Browne et al. Used this as target for a blastp with each species genes AA fasta files as query. Then compute sporulation score following Browne et al. Script include some exploration and visual output to check sensitivity to evalue. 
- Output:
	- blastp output: ``./output/sporulation_genes_blast/`` (not git tracked)
	- Some additional figures (not included in MS): ``./output/figures/additional_figs``
	- Sporulation scores table: ``./output/tables/sporulation_scores.txt``


**9_relative_abundance.R**:
- extracts relative abundance per specie per host from the MIDAS per-sample species outputs.
- Output: table of relative abundance: ``./output/tables/relative_abundance.txt``

10_assemble_data.R: assemble all the tables in ./output/table to create a dataframe used in statistical analyses.
Output: ./output/tables/ANALYSIS_DATA_ASSEMBLED.txt


**11_statistical_analyses.R**:
- all codes for statistical analyses in R.
- Output: models RData output in ``./output/analyses/``

**12_figures_main.R** and **13_figures_SI.R**
- codes for main manuscript figures and supplementary datasets figures
- Output: all figures in ``./output/figures``

