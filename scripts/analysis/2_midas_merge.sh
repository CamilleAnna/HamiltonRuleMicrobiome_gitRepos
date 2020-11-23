#!/bin/sh 
#$ -o ./logs
#$ -e ./logs
#$ -N midas_merge
#$ -cwd
#$ -l h_rt=12:00:00
#$ -l h_vmem=15G
#$ -pe sharedmem 8


# Initialise the environment modules
. /etc/profile.d/modules.sh

module load anaconda
source activate mypythonMIDAS


# Update environment
export PYTHONPATH=$PYTHONPATH:$programs_install_dir/MIDAS
export PATH=$PATH:$programs_install_dir/MIDAS/scripts
export MIDAS_DB=$programs_install_dir/MIDAS/midas_db_v1.2


#cd /exports/eddie/scratch/s1687811/HamiltonRuleMicrobiome/data
cd $user_dir/HamiltonRuleMicrobiome_gitRepos

merge_midas.py species ./output/midas/merge/species_merged/ -i ./output/midas/per_sample/ -t dir
merge_midas.py snps ./output/midas/merge/snps_merged/ -i ./output/midas/per_sample/ -t dir --snp_type any --site_depth 5 --site_prev 0.90 --threads 8












