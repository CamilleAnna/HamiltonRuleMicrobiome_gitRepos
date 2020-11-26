# HamiltonRuleMicrobiome_gitRepos
 
This repository holds data, output and analysis scripts for "Kin selection explains the evolution of cooperation in the gut microbiota" (Simonet & McNally, 2020).


## data

### metagenomes

- "manifest files": from HMP portal HMP portal, accessed April 2020, under Project > HMP, Body Site > feces, Studies > WGS-PP1, File Type > WGS raw sequences set, File format > FASTQ.
- HMP_samples_first_visit_available: list from manifest file filtered to keep earliest time point per host.
- MOCAT.cfg: for reads processing with MOCAT, using defaults settings (used by script 1.MIDAS_per_sample.sh)



## scripts

Script 1 requiress working installation of MOCAT and MIDAS. Scripts 2,3 require a working installation of MOCAT and MIDAS.
Scripts 1-3 are for re-running analysis from scratch, starting from metagenomic samples
Scripts 4-X are for re-running analysis from MIDAS output.

Before executing scripts, define an environmental variables
$programs_install_dir: pointing to those installations directory
$user_dir: path where repository HamiltonRuleMicrobiome_gitRepos is cloned

Scripts 1-3 requires environmental variables:

[...] user defined paths variables --> also local_project_dir

### analysis


0_samples_data.R: parse the HMP portal manifest file to extarct unique sample/subject metagenome access link. This gives 251 samples. Out of these, 12 were filtered at different stage of the analysis. Hence, in the final analysis, 239 samples are included. Specifically:
HMP_SRS013215_454	Download failed. Not a WGS file
HMP_SRS015794_454	Download failed. Not a WGS file
HMP_SRS016989_454	Download failed. Not a WGS file
HMP_SRS057049	Quality filtering (MOCAT) left no read
HMP_SRS011084	Large file cause pipeline crash
HMP_SRS014235	Large file cause pipeline crash
HMP_SRS015890	Large file cause pipeline crash
HMP_SRS017433	Large file cause pipeline crash
HMP_SRS018313	Large file cause pipeline crash
HMP_SRS021219	Large file cause pipeline crash
HMP_SRS019286	No species met minimum read depths for SNP calling (MIDAS)
HMP_SRS045324	Contains only species unique to that host, hence filtered at relatedness computation step. Given the species found in that host, likely there were a mistake in the HMP metadata, those were all species commonly found in oral cavity. This must be an oral sample, not feces

Output:
- List of samples and access links: ./output/tables/HMP_first_visit_smaples.txt

1_midas_per_sample.sh: download HMP metagenomes, process with MOCAT, run MIDAS. Task array jobs, using $SGE_TASK_ID env variable to process the $SGE_TASK_ID th sample. Run as many task as there are samples listed in ./data/HMP_samples_first_visit_available (251). Requires user defined environemental variables $programs_install_dir and $user_dir.

run as: 
qsub -v user_dir='path/where/gitrepo/is/cloned',programs_install_dir='/path/to/install/dir/of/MOCAT_and_MIDAS' -t 1-251 ./scripts/analysis/1_midas_per_sample.sh


2_midas_merge.sh: merge species abundance files across samples (merge_midas.py species) to output relative abundance matrix and table of species prevalence. Then perform pooled-sample core-genome SNP calling across all core genomic sites for each species (merge_midas.py snps). Both commands take directory ./output/midas/per_sample/ as input. In SNP calling, --sample_depth (5.0) and --site_ratio (2.0) left to default. But set: --site_depth 5 for better accuracy in allele frequency estimates, --snp_type any to take into account all variants even if observed at low frequency (very low frequency variants will simply have minimal impact on relatedness calculation), --site_prev 0.90 for more relaxed core sites identification setting (95% prevalence is too stringent for species detected in a small number of hosts).

run as:
qsub -v user_dir='path/where/gitrepo/is/cloned',programs_install_dir='/path/to/install/dir/of/MOCAT_and_MIDAS' ./scripts/analysis/2_midas_merge.sh


run as:
qsub -v user_dir='path/where/gitrepo/is/cloned',programs_install_dir='/path/to/install/dir/of/MOCAT_and_MIDAS' -t1-XXX ./scripts/analysis/3_midas_diversity.sh

3_midas_diversity.sh



4_relatedness.R: compute relatedness from the diversity.within and diversity.between output of MIDAS.
Output: - relatedness table: /output/tables/relatedness.txt



5_secretome.sh: for the 101 species included in the analysis, download the representative genome (as defined in MIDASdb) AA fasta file fromt PATRIC, run PSORTb on them, compute secretome size. Code to install psortb (based on docker image) included. Part of the script runs in R. Saves fasta and feature files in ./data/patric. Saves table of secretome sizes in ./output/tables/secretome.txt



6_identify_social_GO.sh: creates a bacterial GO slim by annotated all representative genomes of MIDASdb with GO using Pannzer. Then use a list of keywords of bacterial social behaviour to browse this bacteria GO slim and retreive GO terms capturing cooperative behaviour categorised in 5 types of cooperation. The last part browsing the assembled GO  slim runs in R. Requires a working installation of Pannzer, specify path to installation directory with environemental variable programs_install_dir. The Pannzer step code is written to run iteratively through a loop working through all genomes one by one because large ammount of parrallels works launched on Pannzer at once (e.g. with task array job) seem to dramatically slow pannzer server (priority issue?) so the jobs never finish. When a single task is launched to Pannzer, it takes about 20min to process a genome. GO slim is simply the assembly of all unique GO terms identified across all >5900 genomes (taking all ARGOT ranks).
NOTE: the GO retreived when browsing the GO slim may differ as the GO.db gets updated. This may give slightly different social_go_list_wide.txt output. Current table (used for analysis) ran August 2019.
Output:
- All MIDAS species fasta files saved in ./data/patric/fasta_all_midas (not git tracked)
- All pannzer output saved in ./output/pannzer (not git tracked)
- Assembled GO slim: ./output/tables/bacteria_go_slim.txt
- list of identified bacteria sociality keywords, identified across 10 top reviews in web of science search. Search string for Web of Science search  included in script. List of keywors in: ./output/tables/bacteria_social_keywords.txt
- Table of all potential bacterial sociality GOs identified by raw keyword matching  + inclusion of direct parent and all child terms: ./output/tables/social_go_list_wide.txt
- Manual curation of that table, with comments of decisions: ./output/tables/social_go_list_curation.xls
- Final retained set of bacteria social GO terms: ./output/tables/social_go_list_final.xls


7_go_categories.R: computes the measure of cooperation based on GO annotation for each species, with the different cooperation split between the 5 cooperation categories. Runs in R. Uses the social_go_list_final.xls and the GO annotations (pannzer output) as input.
output:
- table of measure of cooperation based on GO categories /output/tables/go_cooperation_categories.txt 


8_sporulation_scores.sh: computes sporulation scores for the 101 species following method from Browne et al. We gathered (from NCBI) the AA fasta sequences of the 66 sporulation genes listed by Browne et al. Used this as target for a blastp with each species genes AA fasta files as query. Then compute sporulation score following Browne et al. Script include some exploration and visual output to check sensitivity to evalue. 
Output:
- blastp output: ./output/sporulation_genes_blast/ (not git tracked)
- Some additional figures (not included in MS): ./output/figures/additional_figs 
- Sporulation scores table: ./output/tables/sporulation_scores.txt


9_relative_abundance.R: extract relative abundance per specie per host from the MIDAS per-sample species outputs.
Output:
- Table of relative abundance: ./output/tables/relative_abundance.txt

10_assemble_data.R: assemble all the tables in ./output/table to create a dataframe used in statistical analyses.
Output: ./output/tables/ANALYSIS_DATA_ASSEMBLED.txt


