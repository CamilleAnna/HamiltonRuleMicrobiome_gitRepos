#!/bin/sh 
#$ -o /exports/eddie/scratch/s1687811/HamiltonRuleMicrobiome/logs
#$ -e /exports/eddie/scratch/s1687811/HamiltonRuleMicrobiome/logs
#$ -N midas_diversity
#$ -l h_rt=24:00:00
#$ -l h_vmem=15G
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



SPECIES=$(ls ./midas_merged/snps_merged_sampleDepth5_siteDepth5_sitePrev090/ | awk "NR==$SGE_TASK_ID")

# quantify within sample heterogeneity
# input path is path to output of `merge_midas.py snps` for one species
snp_diversity.py /exports/eddie/scratch/s1687811/HamiltonRuleMicrobiome/data/midas_merged/snps_merged_sampleDepth5_siteDepth5_sitePrev090/$SPECIES --genomic_type genome-wide --sample_type per-sample --out /exports/eddie/scratch/s1687811/HamiltonRuleMicrobiome/data/midas_merged/diversity_sampleDepth5_siteDepth5_sitePrev090/$SPECIES.diversity.within
# quantify pooled samples heterogeneity
snp_diversity.py /exports/eddie/scratch/s1687811/HamiltonRuleMicrobiome/data/midas_merged/snps_merged_sampleDepth5_siteDepth5_sitePrev090/$SPECIES --genomic_type genome-wide --sample_type pooled-samples --out /exports/eddie/scratch/s1687811/HamiltonRuleMicrobiome/data/midas_merged/diversity_sampleDepth5_siteDepth5_sitePrev090/$SPECIES.diversity.between

