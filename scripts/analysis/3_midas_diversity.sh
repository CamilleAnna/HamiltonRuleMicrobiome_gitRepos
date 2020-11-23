#!/bin/sh 
#$ -o ./logs
#$ -e ./logs
#$ -N midas_diversity
#$ -cwd
#$ -l h_rt=24:00:00
#$ -l h_vmem=15G

# Array job, each task processes one species
# run as many tasks as there are species for which snps were computed
# species listed in snps_merged directory (output of merge_midas.py snps command in previous step script '2_midas_merge.sh')

# Initialise the environment modules
. /etc/profile.d/modules.sh

module load anaconda
source activate mypythonMIDAS


# Update environment
export PYTHONPATH=$PYTHONPATH:$programs_install_dir/MIDAS
export PATH=$PATH:$programs_install_dir/MIDAS/scripts
export MIDAS_DB=$programs_install_dir/MIDAS/midas_db_v1.2


cd $user_dir/HamiltonRuleMicrobiome_gitRepos


SPECIES=$(ls ./output/midas/merge/snps_merged/ | awk "NR==$SGE_TASK_ID")

# quantify within sample heterogeneity
# input path is path to output of `merge_midas.py snps` for the focal $SPECIES
snp_diversity.py ./output/midas/merge/snps_merged/$SPECIES --genomic_type genome-wide --sample_type per-sample --out ./output/midas/diversity/$SPECIES.diversity.within
# quantify pooled samples heterogeneity
snp_diversity.py ./output/midas/merge/snps_merged/$SPECIES --genomic_type genome-wide --sample_type pooled-samples --out ./output/midas/diversity/$SPECIES.diversity.between

