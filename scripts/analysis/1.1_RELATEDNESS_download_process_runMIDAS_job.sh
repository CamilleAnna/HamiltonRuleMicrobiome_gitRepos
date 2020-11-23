#!/bin/sh 
#$ -o ./logs
#$ -e ./logs
#$ -N download_process_runMIDAS
#$ -cwd
#$ -l h_rt=24:00:00
#$ -l h_vmem=15G
#$ -P bio_MPGS_csimonet



# IMPORT & MOVE IN REPOSITORY BEFORE REUNNING SCRIPT, LAUNCH SCRIPT FROM ROOT OF REPOSITORY
# git clone https://github.com/CamilleAnna/HamiltonRuleMicrobiome_gitRepos.git
# cd ./HamiltonRuleMicrobiome_gitRepos/



############################################# PART 1: DOWNLOAD METAGENOMES #############################################
echo 'PART 1: DOWNLOADING FILES:'

# Extract samples corresponding to first (available) visit of each host
cat ./data/metagenomes/HMP_manifest_metadata_1955206147.tsv | sed  1d | sort -k 2,5 |  awk '!seen[$2]++' > ./data/metagenomes/HMP_samples_first_visit_available


# Download metagenome, create one directory per host for next step (MOCAT) to run in parrals task array job
SAMPLE=$(cut -f 1 ./data/metagenomes/HMP_samples_first_visit_available | awk "NR==$SGE_TASK_ID")
LINK=$(grep $SAMPLE ./data/metagenomes/HMP_manifest_68fb5dcb39.tsv | cut -f 4 | cut -f 1 -d',')
HOST=$(echo HMP_$(echo $LINK | rev | cut -d '/' -f 1 | rev | cut -d '.' -f 1))
mkdir -p ./data/metagenomes/metagenomes/$HOST
cd ./data/metagenomes/metagenomes/$HOST
echo Downloading: $LINK
wget -q $LINK


# Set directory structure and files names the right way for MOCAT processing
tar xjf *.tar.bz2
rm *.tar.bz2
dir1=$(ls)
find ./$dir1 -mindepth 1 -type f -exec mv -t . -n '{}' +
rm -r $dir1	
rm *singleton*
mv *.1.f* $HOST.1.fq
mv *.2.f* $HOST.2.fq

gzip *


if [ -f $HOST.1.fq.gz ] && [ -f $HOST.2.fq.gz ]; then echo 'ALL DONE WITH DOWNLOAD AND FILE PREP'; else echo 'At least one file missing after download & file prep step, a problem occured. Exiting'; exit; fi
#echo "at $(date) , $(du -sh /exports/eddie/scratch/s1687811)"


############################################# PART 2: PROCESS #############################################
echo 'PART 2: READS PROCESSING'
cp /exports/eddie/scratch/s1687811/HamiltonRuleMicrobiome/data/metagenomes/MOCAT.cfg .

# RUN MOCAT MANUALLY
mkdir -p logs/other
mkdir -p logs/status
mkdir -p logs/readtrimfilter/commands
mkdir -p logs/readtrimfilter/jobs
mkdir -p logs/readtrimfilter/samples
mkdir -p logs/readtrimfilter/startstop
mkdir -p logs/screen/commands
mkdir -p logs/screen/jobs
mkdir -p logs/screen/samples
mkdir -p logs/screen/startstop
mkdir -p logs/screen_fasta_file/commands
mkdir -p logs/screen_fasta_file/jobs
mkdir -p logs/screen_fasta_file/samples
mkdir -p logs/screen_fasta_file/startstop

# MOCAT COMMAND (based on generated command script from MOCAT wrapper)
exec 2>> ./logs/readtrimfilter/samples/MOCATJob_readtrimfilter.$HOST.log && mkdir -p temp && /exports/csce/eddie/biology/groups/mcnally/camille/programs/MOCAT//src/MOCATReadTrimFilter_aux.pl -sample $HOST -trim_5prime_end yes -src_dir /exports/csce/eddie/biology/groups/mcnally/camille/programs/MOCAT//src -paired_end_data yes -file_formats_array ss -length_cutoff 45 -qual_cutoff 20 -solexaqa_or_fastx solexaqa -bin_dir /exports/csce/eddie/biology/groups/mcnally/camille/programs/MOCAT//bin -cwd /exports/eddie/scratch/s1687811/HamiltonRuleMicrobiome/data/metagenomes/metagenomes/ -use3files 0 -use_5prime_file no -temp_dir /exports/eddie/scratch/s1687811/HamiltonRuleMicrobiome/data/metagenomes/metagenomes/ -zcat "gunzip -c" 2>> /exports/eddie/scratch/s1687811/HamiltonRuleMicrobiome/data/metagenomes/metagenomes/$HOST/logs/readtrimfilter/samples/MOCATJob_readtrimfilter.$HOST.log >> /exports/eddie/scratch/s1687811/HamiltonRuleMicrobiome/data/metagenomes/metagenomes/$HOST/logs/readtrimfilter/samples/MOCATJob_readtrimfilter.$HOST.log

size_rtf=$(du -sh ./reads.processed.solexaqa/*.pair.1.fq.gz | cut -f 1)
if [ $size_rtf != 0 ]; then echo 'read trim filter OK'; else echo 'read trim filter left no read, a problem occured. Exiting'; exit; fi

# Housekeeping
rm MOCAT.cfg
rm -r logs
rm -r temp

mkdir raw.reads
mv *.fq.gz ./raw.reads


############################################# PART 3: RUN MIDAS #############################################
echo "PART 3: RUN MIDAS ON $HOST"

# Set environment modules
. /etc/profile.d/modules.sh

module load anaconda
source activate mypythonMIDAS
export PYTHONPATH=$PYTHONPATH:/exports/csce/eddie/biology/groups/mcnally/camille/programs/MIDAS
export PATH=$PATH:/exports/csce/eddie/biology/groups/mcnally/camille/programs/MIDAS/scripts
export MIDAS_DB=/exports/csce/eddie/biology/groups/mcnally/camille/programs/MIDAS/midas_db_v1.2


cd /exports/eddie/scratch/s1687811/HamiltonRuleMicrobiome/data/midas_output
mkdir $HOST

# run midas on the processed files
echo "run MIDAS species $HOST"
run_midas.py species ./$HOST -1 /exports/eddie/scratch/s1687811/HamiltonRuleMicrobiome/data/metagenomes/metagenomes/$HOST/reads.processed.solexaqa/*pair.1.fq.gz -2 /exports/eddie/scratch/s1687811/HamiltonRuleMicrobiome/data/metagenomes/metagenomes/$HOST/reads.processed.solexaqa/*pair.2.fq.gz #-n 1000000


echo "run MIDAS snps $HOST"
run_midas.py snps ./$HOST -1 /exports/eddie/scratch/s1687811/HamiltonRuleMicrobiome/data/metagenomes/metagenomes/$HOST/reads.processed.solexaqa/*pair.1.fq.gz -2 /exports/eddie/scratch/s1687811/HamiltonRuleMicrobiome/data/metagenomes/metagenomes/$HOST/reads.processed.solexaqa/*pair.2.fq.gz --remove_temp #-n 50000


if [ -f ./$HOST/snps/summary.txt ]; then echo 'MIDAS ran fine until the end'; else echo 'There might have been an issue with midas or no sample was suitable for snps search'; fi
echo 'ALL DONE! :) '







