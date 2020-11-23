# HamiltonRuleMicrobiome_gitRepos
 
This repository holds data, output and analysis scripts for "Kin selection explains the evolution of cooperation in the gut microbiota" (Simonet & McNally, 2020).

## data

### metagenomes

- "manifest files": from HMP portal HMP portal, accessed April 2020, under Project > HMP, Body Site > feces, Studies > WGS-PP1, File Type > WGS raw sequences set, File format > FASTQ.
- HMP_samples_first_visit_available: list from manifest file filtered to keep earliest time point per host.
- MOCAT.cfg: for reads processing with MOCAT, using defaults settings (used by script 1.MIDAS_per_sample.sh)



## scripts

Require a working installation of MOCAT, and MIDAS
Before executing scripts, define an environmental variables
$programs_install_dir: pointing to those installations directory
$user_dir: path where repository HamiltonRuleMicrobiome_gitRepos is cloned


### analysis

1.MIDAS_per_sample.sh: download HMP metagenomes, process with MOCAT, run MIDAS. Task array jobs, using $SGE_TASK_ID env variable to process the $SGE_TASK_ID th sample. Run as many task as there are samples listed in ./data/HMP_samples_first_visit_available (251)
