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


### analysis

1_midas_per_sample.sh: download HMP metagenomes, process with MOCAT, run MIDAS. Task array jobs, using $SGE_TASK_ID env variable to process the $SGE_TASK_ID th sample. Run as many task as there are samples listed in ./data/HMP_samples_first_visit_available (251). Requires user defined environemental variables $programs_install_dir and $user_dir.

run as: 
qsub -v user_dir='path/where/gitrepo/is/cloned',programs_install_dir='/path/to/install/dir/of/MOCAT_and_MIDAS' -t 1-251 ./scripts/analysis/1_midas_per_sample.sh


2_midas_merge.sh: merge species abundance files across samples (merge_midas.py species) to output relative abundance matrix and table of species prevalence. Then perform pooled-sample core-genome SNP calling across all core genomic sites for each species (merge_midas.py snps). Both commands take directory ./output/midas/per_sample/ as input. In SNP calling, --sample_depth (5.0) and --site_ratio (2.0) left to default. But set: --site_depth 5 for better accuracy in allele frequency estimates, --snp_type any to take into account all variants even if observed at low frequency (very low frequency variants will simply have minimal impact on relatedness calculation), --site_prev 0.90 for more relaxed core sites identification setting (95% prevalence is too stringent for species detected in a small number of hosts).

run as:
qsub -v user_dir='path/where/gitrepo/is/cloned',programs_install_dir='/path/to/install/dir/of/MOCAT_and_MIDAS' ./scripts/analysis/2_midas_merge.sh



run as:
qsub -v user_dir='path/where/gitrepo/is/cloned',programs_install_dir='/path/to/install/dir/of/MOCAT_and_MIDAS' -t1-XXX ./scripts/analysis/3_midas_diversity.sh








