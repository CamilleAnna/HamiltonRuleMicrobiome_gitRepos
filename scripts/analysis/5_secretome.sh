# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
#                    Simonet & McNally 2020                     #
#                   Computing secretome size                    #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


# For the species included in the analysis (for which we computed relatedness in previous step)
# We measure secretome size
# That's 101 species, listed in ./data/species_info_files/species_list.txt
# We work with the representative genome of those species, as listed in MIDASdb

# STEPS ARE:
# 1) Browse PATRIC (download genomes AA fasta and feature files)
# 2) Run PSORTb (install, generate commands, run)
# 3) compute secretome size (code to execute in R)


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
if [[ -f psortb_commands.sh ]]; then echo 'deleting existing psortb_commands.sh file' && rm psortb_commands.sh; fi

cat $local_project_dir/HamiltonRuleMicrobiome_gitRepos/output/species_list.txt | sed '1d' | while read line
do
FILE=$line
GRAM=$(cat $local_project_dir/HamiltonRuleMicrobiome_gitRepos/data/species_info_files/gram_profiles_db.txt | grep $line | cut -f 2)
if [[ $GRAM == *n* ]];   then echo "mkdir $FILE" >> psortb_commands.sh && echo "$psortb_install_dir -i $local_project_dir/HamiltonRuleMicrobiome_gitRepos/data/patric/fasta/$FILE.fasta -r $local_project_dir/HamiltonRuleMicrobiome_gitRepos/output/psortb/$FILE -n -o long && rm $local_project_dir/HamiltonRuleMicrobiome_gitRepos/output/psortb/$FILE/*\.fasta && mv $local_project_dir/HamiltonRuleMicrobiome_gitRepos/output/psortb/$FILE/*.txt $local_project_dir/HamiltonRuleMicrobiome_gitRepos/output/psortb/$FILE/$FILE\.psortb.out && echo 'done with $FILE' >> $local_project_dir/HamiltonRuleMicrobiome_gitRepos/logs/psortb.log" >> psortb_commands.sh;
elif [[ $GRAM == *p* ]]; then echo "mkdir $FILE" >> psortb_commands.sh && echo "$psortb_install_dir -i $local_project_dir/HamiltonRuleMicrobiome_gitRepos/data/patric/fasta/$FILE.fasta -r $local_project_dir/HamiltonRuleMicrobiome_gitRepos/output/psortb/$FILE -p -o long && rm $local_project_dir/HamiltonRuleMicrobiome_gitRepos/output/psortb/$FILE/*\.fasta && mv $local_project_dir/HamiltonRuleMicrobiome_gitRepos/output/psortb/$FILE/*.txt $local_project_dir/HamiltonRuleMicrobiome_gitRepos/output/psortb/$FILE/$FILE\.psortb.out && echo 'done with $FILE' >> $local_project_dir/HamiltonRuleMicrobiome_gitRepos/logs/psortb.log" >> psortb_commands.sh;
else echo "For $FILE gram not determined, can't run psortb" >> noGram.txt;
fi
done

# run commands
# you may run into a 'Permission denied' error when running the following, related to github tracking
# either run in a separate local directory, or unlink directory from github by running
# rm -rf .git*
sudo su # required to run docker image of  psortb if you don't have administrative rights
sh psortb_commands.sh
exit

# housekeeping
mkdir psortb_output
mv ./*/*.psortb.out ./psortb_output/
rmdir * # removes NOW EMPTY directories



# ~~~~~~~~~~~~~~~~~~~~~~ #
#  3) Compute Secretome  #
# ~~~~~~~~~~~~~~~~~~~~~~ #
# ** code to execute in R **


# local_project_dir='/path/to/cloned/repository/'
setwd(local_project_dir)
library(dplyr)
library(tidyr)

# Gathering PSORTb output
path.to.psotb.output<- paste0(local_project_dir, '/HamiltonRuleMicrobiome_gitRepos/output/psortb/psortb_output/')
files<- paste0(path.to.psotb.output, list.files(path.to.psotb.output))
sps<- list.files(path.to.psotb.output)
loc<- list()
compartments<- c("Cytoplasmic_Score", "CytoplasmicMembrane_Score", "Periplasmic_Score", "OuterMembrane_Score", "Extracellular_Score", "Cellwall_Score", "Final_Localization")
for(i in 1:length(files)){
  
  
  coop<- read.csv(files[i], sep='\t', row.names = NULL)
  
  if(nrow(coop) == 0){
    loc[[i]]<- data.frame(species = sub(".psortb.out", "", sps[i]),
                          location = c('Cytoplasmic', 'CytoplasmicMembrane', 'Extracellular', 'OuterMembrane', 'Periplasmic', 'Unknown', 'Cellwall'),
                          freq = NA,
                          nb_sequences = NA,
                          nb_sequences_known = NA)
  }else{
    colnames(coop)<- c(colnames(coop)[-1], 'foo')
    coop2<- coop[,which(colnames(coop) %in% compartments)] 
    
    loc[[i]]<- data.frame(table(coop2$Final_Localization)) %>%
      rename(location = Var1, freq = Freq) %>%
      mutate(species = sub(".psortb.out", "", sps[i]),
             nb_sequences = nrow(coop),
             nb_sequences_known = nrow(coop[which(coop$Final_Localization != 'Unknown'),])) %>% 
      select(species, location, freq, nb_sequences, nb_sequences_known)
  }
  
  print(i)
}
locfl<- do.call('rbind', loc)

# Compute total number/localisation, keep Extracellular location only,
locfl<- locfl %>%
  mutate(other_all = locfl$nb_sequences - locfl$freq,
         other_known = locfl$nb_sequences_known - locfl$freq)

locfl<- locfl[which(locfl$location == 'Extracellular'),] %>%
  select(species, freq, nb_sequences, nb_sequences_known) %>%
  rename(nb_extracellular = freq)
colnames(locfl)[1]<- 'species_id'

# Write secretome_size table
write.table(locfl, paste0(local_project_dir, '/HamiltonRuleMicrobiome_gitRepos/output/secretome.txt'), col.names = TRUE, row.names = FALSE, quote = FALSE, sep = '\t')




# ~~~~~~~~~~~~~~~~~ END OF SCRIPT 5_secretome ~~~~~~~~~~~~~~~~~~~~ #







