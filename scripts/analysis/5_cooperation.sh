# other data assembly

# We processed the 251 HMP stool metagenomic samples
# We computed measures of relatedness of 101 species
# Now, measure cooperativity for those 101 species, listed in ./data/species_info_files/species_list.txt
# Work with representative genome of that species, as listed in MIDASdb

# Procedure is:
# 1) Download genomes and features from PATRICdb
# 2) Install psortb (easier on local machine, unless you've got admin rights)


### 1) Download genomes and features from PATRICdb

cd $local_project_dir/HamiltonRuleMicrobiome_gitRepos/data/patric/

cat $local_project_dir/HamiltonRuleMicrobiome_gitRepos/output/species_list.txt | sed '1d' | while read line
do
SP_MIDAS=$line
SP_PATRIC=$(grep $line $local_project_dir/HamiltonRuleMicrobiome_gitRepos/data/species_info_files/species_info.txt | cut -f 2)
wget -N "ftp://ftp.patricbrc.org/genomes/$SP_PATRIC/$SP_PATRIC.PATRIC.faa";
mv $SP_PATRIC.PATRIC.faa $SP_MIDAS.fasta
wget -N "ftp://ftp.patricbrc.org/genomes/$SP_PATRIC/$SP_PATRIC.PATRIC.features.tab";
mv $SP_PATRIC.PATRIC.features.tab $SP_MIDAS.features
done

# housekeeping
mv *.fasta ./fasta
mv *.features ./features
cd ./fasta/ # shorten sequences names to avoid clash with psortb
ls | while read line
do
sed 's/   .*//' $line > temp.txt && mv temp.txt $line.edited && rm $line && mv $line.edited $line
done


### 2) Install PSORTb on local machine

/usr/bin/ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)"
brew cask install docker
# Have docker open and then run:
sudo docker pull brinkmanlab/psortb_commandline:1.0.0
brew install wget
wget https://raw.githubusercontent.com/brinkmanlab/psortb_commandline_docker/master/psortb
chmod +x psortb

### 3) Generate psortb commands













