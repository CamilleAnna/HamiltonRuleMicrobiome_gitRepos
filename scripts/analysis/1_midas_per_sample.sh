#!/bin/sh 
#$ -o ./logs
#$ -e ./logs
#$ -N download_process_runMIDAS
#$ -cwd
#$ -l h_rt=24:00:00
#$ -l h_vmem=15G



# IMPORT & MOVE IN REPOSITORY BEFORE REUNNING SCRIPT, LAUNCH SCRIPT FROM ROOT OF REPOSITORY
# git clone https://github.com/CamilleAnna/HamiltonRuleMicrobiome_gitRepos.git
# cd ./HamiltonRuleMicrobiome_gitRepos/
# Have a working installation of MOCAT, and MIDAS
# Before executing scripts, define an environmental variables
# $programs_install_dir: pointing to those installations directory
# $user_dir: path where repository HamiltonRuleMicrobiome_gitRepos is cloned


############################################# PART 1: DOWNLOAD METAGENOMES #############################################
echo 'PART 1: DOWNLOADING FILES:'

# Available HMP stool metagenomes at HMP portal, accessed April 2020, under:
# Project > HMP, Body Site > feces, Studies > WGS-PP1, File Type > WGS raw sequences set, File format > FASTQ
# Get "manifest file" and corresponding metadata
# filtered using below command to keep earliest time point collected per hosts
# list of final samples: HMP_samples_first_visit_available
# launch this script as a task array job to process each $SGE_TASK_ID th sample

# cat ./data/metagenomes/HMP_manifest_metadata_1955206147.tsv | sed  1d | sort -k 2,5 |  awk '!seen[$2]++' > ./data/metagenomes/HMP_samples_first_visit_available


# Download metagenome under a host-specific directory per host
SAMPLE=$(cut -f 1 $user_dir/HamiltonRuleMicrobiome_gitRepos/data/metagenomes/HMP_samples_first_visit_available | awk "NR==$SGE_TASK_ID")
LINK=$(grep $SAMPLE $user_dir/HamiltonRuleMicrobiome_gitRepos/data/metagenomes/HMP_manifest_68fb5dcb39.tsv | cut -f 4 | cut -f 1 -d',')
HOST=$(echo HMP_$(echo $LINK | rev | cut -d '/' -f 1 | rev | cut -d '.' -f 1))
mkdir -p $user_dir/HamiltonRuleMicrobiome_gitRepos/data/metagenomes/metagenomes/$HOST
cd $user_dir/HamiltonRuleMicrobiome_gitRepos/data/metagenomes/metagenomes/$HOST
echo Downloading: $LINK
wget -q $LINK


# Housekeeping for MOCAT processing.
# MOCAT requires reads1 and reads2 in separate files, with names ending by .1.fq and .2.fq
tar xjf *.tar.bz2
rm *.tar.bz2
dir1=$(ls)
find ./$dir1 -mindepth 1 -type f -exec mv -t . -n '{}' +
rm -r $dir1	
rm *singleton*
mv *.1.f* $HOST.1.fq
mv *.2.f* $HOST.2.fq
gzip *

# checkpoint for log file
if [ -f $HOST.1.fq.gz ] && [ -f $HOST.2.fq.gz ]; then echo 'ALL DONE WITH DOWNLOAD AND FILE PREP'; else echo 'At least one file missing after download & file prep step, a problem occured. Exiting'; exit; fi


############################################# PART 2: PROCESS #############################################
echo 'PART 2: READS PROCESSING'
# MOCAT launches a sh script.
# Hence can't run basic MOCAT wrapper command within this script, because SGE does not let launching an sh script from within a script
# So run MOCAT "manually", directly executing here the command that is generated by the MOCAT wrapper.
# human reads and adapters already removed from available HMP files --> only run a readtrimfilter to ensure uniform quality filtering

cp $user_dir/HamiltonRuleMicrobiome_gitRepos/data/metagenomes/MOCAT.cfg .
mkdir -p logs/readtrimfilter/samples
exec 2>> ./logs/readtrimfilter/samples/MOCATJob_readtrimfilter.$HOST.log && mkdir -p temp && $programs_install_dir/MOCAT//src/MOCATReadTrimFilter_aux.pl -sample $HOST -trim_5prime_end yes -src_dir $programs_install_dir/MOCAT//src -paired_end_data yes -file_formats_array ss -length_cutoff 45 -qual_cutoff 20 -solexaqa_or_fastx solexaqa -bin_dir $programs_install_dir/MOCAT/bin -cwd ../ -use3files 0 -use_5prime_file no -temp_dir ../ -zcat "gunzip -c" 2>> ./logs/readtrimfilter/samples/MOCATJob_readtrimfilter.$HOST.log >> ./logs/readtrimfilter/samples/MOCATJob_readtrimfilter.$HOST.log
size_rtf=$(du -sh ./reads.processed.solexaqa/*.pair.1.fq.gz | cut -f 1)
if [ $size_rtf != 0 ]; then echo 'read trim filter OK'; else echo 'read trim filter left no read, a problem occured. Exiting'; exit; fi


# Housekeeping
rm MOCAT.cfg
rm -r logs
rm -r temp

#mkdir raw.reads
#mv *.fq.gz ./raw.reads
rm *.fq.gz # don't need to keep raw reads

############################################# PART 3: RUN MIDAS #############################################
echo "PART 3: RUN MIDAS ON $HOST"

# Set environment modules
. /etc/profile.d/modules.sh

module load anaconda
source activate mypythonMIDAS # anaconda environment with all pre-requirements for MIDAS
export PYTHONPATH=$PYTHONPATH:$programs_install_dir/MIDAS
export PATH=$PATH:$programs_install_dir/MIDAS/scripts
export MIDAS_DB=$programs_install_dir/MIDAS/midas_db_v1.2


cd $user_dir/HamiltonRuleMicrobiome_gitRepos/output/midas/per_sample
mkdir $HOST

# run midas on the processed files
echo "run MIDAS species $HOST" # for species abundance data + prevalence info
run_midas.py species ./$HOST -1 $user_dir/HamiltonRuleMicrobiome_gitRepos/data/metagenomes/metagenomes/$HOST/reads.processed.solexaqa/*pair.1.fq.gz -2 $user_dir/HamiltonRuleMicrobiome_gitRepos/data/metagenomes/metagenomes/$HOST/reads.processed.solexaqa/*pair.2.fq.gz #-n 1000000


echo "run MIDAS snps $HOST" # for diversity analysis
run_midas.py snps ./$HOST -1 $user_dir/HamiltonRuleMicrobiome_gitRepos/data/metagenomes/metagenomes/$HOST/reads.processed.solexaqa/*pair.1.fq.gz -2 $user_dir/HamiltonRuleMicrobiome_gitRepos/data/metagenomes/metagenomes/$HOST/reads.processed.solexaqa/*pair.2.fq.gz --remove_temp #-n 50000


rm $user_dir/HamiltonRuleMicrobiome_gitRepos/data/metagenomes/metagenomes/$HOST/reads.processed.solexaqa/*.fq.gz
# keep qual_stats and stats directory for records, delete actual reads for space issues

if [ -f ./$HOST/snps/summary.txt ]; then echo 'MIDAS ran fine until the end'; else echo 'There might have been an issue with midas or no species sastisfied filtering criteria'; fi
echo 'ALL DONE! :) '







