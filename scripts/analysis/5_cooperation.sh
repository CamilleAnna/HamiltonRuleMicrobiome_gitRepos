# other data assembly

# We processed the 251 HMP stool metagenomic samples
# We computed measures of relatedness of 101 species
# Now, measure cooperativity for those 101 species, listed in ./data/species_info_files/species_list.txt
# Work with representative genome of that species, as listed in MIDASdb

# Procedure is:
# 1) Browse PATRIC (download genomes AA fasta and feature files)
# 2) Run PSORTb (install, generate commands, run)



# ~~~~~~~~~~~~~~~~~~~ #
#  1) Browse PATRIC   #
# ~~~~~~~~~~~~~~~~~~~ #


cd $local_project_dir/HamiltonRuleMicrobiome_gitRepos/data/patric/

cat $local_project_dir/HamiltonRuleMicrobiome_gitRepos/output/species_list.txt | sed '1d' | while read line
do
SP_MIDAS=$line
SP_PATRIC=$(grep $line $local_project_dir/HamiltonRuleMicrobiome_gitRepos/data/species_info_files/species_info.txt | cut -f 2)
wget -N "ftp://ftp.patricbrc.org/genomes/$SP_PATRIC/$SP_PATRIC.PATRIC.faa";
wget -N "ftp://ftp.patricbrc.org/genomes/$SP_PATRIC/$SP_PATRIC.PATRIC.features.tab";
mv $SP_PATRIC.PATRIC.faa $SP_MIDAS.fasta
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


# ~~~~~~~~~~~~~~~~~~~ #
#    2) Run PSORTb    #
# ~~~~~~~~~~~~~~~~~~~ #

# Install (easier on local machine, unless you've got admin rights)
/usr/bin/ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)"
brew cask install docker
# Have docker open and then run:
sudo docker pull brinkmanlab/psortb_commandline:1.0.0
brew install wget
wget https://raw.githubusercontent.com/brinkmanlab/psortb_commandline_docker/master/psortb
chmod +x psortb


# Generate psortb commands
psortb_install_dir='/Applications/psortDB/psortb'
cd $local_project_dir/HamiltonRuleMicrobiome_gitRepos/output/psortb
cat $local_project_dir/HamiltonRuleMicrobiome_gitRepos/output/species_list.txt | sed '1d' | while read line
do
FILE=$line
GRAM=$(cat $local_project_dir/HamiltonRuleMicrobiome_gitRepos/data/species_info_files/gram_profiles_db.txt | grep $line | cut -f 2)
if [[ $GRAM == *n* ]];   then echo "$psortb_install_dir -i $local_project_dir/HamiltonRuleMicrobiome_gitRepos/data/patric/fasta/$FILE.fasta -r $local_project_dir/HamiltonRuleMicrobiome_gitRepos/output/psortb/$FILE -n -o long && rm $local_project_dir/HamiltonRuleMicrobiome_gitRepos/output/psortb/$FILE/*\.fasta && mv $local_project_dir/HamiltonRuleMicrobiome_gitRepos/output/psortb/$FILE/*.txt $local_project_dir/HamiltonRuleMicrobiome_gitRepos/output/psortb/$FILE/$FILE\.psortb.out && echo 'done with $FILE' >> $local_project_dir/HamiltonRuleMicrobiome_gitRepos/logs/psortb.log" >> psortb_commands.sh;
elif [[ $GRAM == *p* ]]; then echo "$psortb_install_dir -i $local_project_dir/HamiltonRuleMicrobiome_gitRepos/data/patric/fasta/$FILE.fasta -r $local_project_dir/HamiltonRuleMicrobiome_gitRepos/output/psortb/$FILE -p -o long && rm $local_project_dir/HamiltonRuleMicrobiome_gitRepos/output/psortb/$FILE/*\.fasta && mv $local_project_dir/HamiltonRuleMicrobiome_gitRepos/output/psortb/$FILE/*.txt $local_project_dir/HamiltonRuleMicrobiome_gitRepos/output/psortb/$FILE/$FILE\.psortb.out && echo 'done with $FILE' >> $local_project_dir/HamiltonRuleMicrobiome_gitRepos/logs/psortb.log" >> psortb_commands.sh;
else echo "For $FILE gram not determined, can't run psortb" >> noGram.txt;
fi
done


# Run
cd $local_project_dir/HamiltonRuleMicrobiome_gitRepos/output/psortb
cat $local_project_dir/HamiltonRuleMicrobiome_gitRepos/output/species_list.txt | sed '1d' | while read line
do
FILE=$(echo $line | cut -f 2 -d' ')
mkdir $FILE
done

sudo su
sh psortb_commands.sh
exit



# Generate psortb commands
psortb_install_dir='/Applications/psortDB/psortb'
cd $local_project_dir/HamiltonRuleMicrobiome_gitRepos/output/psortb
cat $local_project_dir/HamiltonRuleMicrobiome_gitRepos/output/species_list.txt | sed '1d' | while read line
do
FILE=$line
GRAM=$(cat $local_project_dir/HamiltonRuleMicrobiome_gitRepos/data/species_info_files/gram_profiles_db.txt | grep $line | cut -f 2)
if [[ $GRAM == *n* ]];   then echo "mkdir $FILE" >> psortb_commands.sh && echo "$psortb_install_dir -i $local_project_dir/HamiltonRuleMicrobiome_gitRepos/data/patric/fasta/$FILE.fasta -r $local_project_dir/HamiltonRuleMicrobiome_gitRepos/output/psortb/$FILE -n -o long && rm $local_project_dir/HamiltonRuleMicrobiome_gitRepos/output/psortb/$FILE/*\.fasta && mv $local_project_dir/HamiltonRuleMicrobiome_gitRepos/output/psortb/$FILE/*.txt $local_project_dir/HamiltonRuleMicrobiome_gitRepos/output/psortb/$FILE/$FILE\.psortb.out && echo 'done with $FILE' >> $local_project_dir/HamiltonRuleMicrobiome_gitRepos/logs/psortb.log" >> psortb_commands.sh;
elif [[ $GRAM == *p* ]]; then echo "mkdir $FILE" >> psortb_commands.sh && echo "$psortb_install_dir -i $local_project_dir/HamiltonRuleMicrobiome_gitRepos/data/patric/fasta/$FILE.fasta -r $local_project_dir/HamiltonRuleMicrobiome_gitRepos/output/psortb/$FILE -p -o long && rm $local_project_dir/HamiltonRuleMicrobiome_gitRepos/output/psortb/$FILE/*\.fasta && mv $local_project_dir/HamiltonRuleMicrobiome_gitRepos/output/psortb/$FILE/*.txt $local_project_dir/HamiltonRuleMicrobiome_gitRepos/output/psortb/$FILE/$FILE\.psortb.out && echo 'done with $FILE' >> $local_project_dir/HamiltonRuleMicrobiome_gitRepos/logs/psortb.log" >> psortb_commands.sh;
else echo "For $FILE gram not determined, can't run psortb" >> noGram.txt;
fi
done
sudo su # required to run docker image of  psortb if you don't have administrative rights
sh psortb_commands.sh
exit







