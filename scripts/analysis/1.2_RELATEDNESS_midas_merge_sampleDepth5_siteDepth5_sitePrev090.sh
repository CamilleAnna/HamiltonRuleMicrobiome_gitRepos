#!/bin/sh 
#$ -o /exports/eddie/scratch/s1687811/HamiltonRuleMicrobiome/logs
#$ -e /exports/eddie/scratch/s1687811/HamiltonRuleMicrobiome/logs
#$ -N midas_merge_sampleDepth5_siteDepth5_sitePrev090
#$ -l h_rt=12:00:00
#$ -l h_vmem=15G
#$ -pe sharedmem 8
#$ -P bio_MPGS_csimonet


# Initialise the environment modules
. /etc/profile.d/modules.sh

module load anaconda
source activate mypythonMIDAS


# Update environment
export PYTHONPATH=$PYTHONPATH:/exports/csce/eddie/biology/groups/mcnally/camille/programs/MIDAS
export PATH=$PATH:/exports/csce/eddie/biology/groups/mcnally/camille/programs/MIDAS/scripts
export MIDAS_DB=/exports/csce/eddie/biology/groups/mcnally/camille/programs/MIDAS/midas_db_v1.2


cd /exports/eddie/scratch/s1687811/HamiltonRuleMicrobiome/data

mkdir -p /exports/eddie/scratch/s1687811/HamiltonRuleMicrobiome/data/midas_merged/species_merged_sampleDepth5_siteDepth5_sitePrev090
merge_midas.py species ./midas_merged/species_merged_sampleDepth5_siteDepth5_sitePrev090 -i ./midas_output -t dir


mkdir -p /exports/eddie/scratch/s1687811/HamiltonRuleMicrobiome/data/midas_merged/snps_merged_sampleDepth5_siteDepth5_sitePrev090
merge_midas.py snps ./midas_merged/snps_merged_sampleDepth5_siteDepth5_sitePrev090 -i ./midas_output -t dir --sample_depth 5 --snp_type any --site_depth 5 --site_ratio 2.0 --site_prev 0.90 --threads 8 #--max_sites 1000


mkdir -p /exports/eddie/scratch/s1687811/HamiltonRuleMicrobiome/data/midas_merged/diversity_sampleDepth5_siteDepth5_sitePrev090
