# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
#                    Simonet & McNally 2020                     #
#                      Identify social GOs                      #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

# For the species included in the analysis (for which we computed relatedness in previous step)
# We now measure cooperativity split in 5 GO categories
# To start, we must identify GO terms tha capture bacterial social behaviour
# to do so, we:
#  - annotate all representative genomes of MIDASdb species with GO using Pannzer
#  - use a web of science search to identify bacterial social behaviour keywords
#  - look through all GO terms identified in bacteria for those matching keywords


# STEPS ARE:
# 1) Browse PATRIC (get all genomes AA fasta of MIDAS db representative genomes, n = 5952)
# 2) Run Pannzer on all
# 3) Assemble GO slim
# 4) Web of science search of social keywords
# 5) Identify social GOs


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
#  1) Browse PATRIC for all MIDAS representative genomes    #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

cd $local_project_dir/HamiltonRuleMicrobiome_gitRepos/data/patric/

cat $local_project_dir/HamiltonRuleMicrobiome_gitRepos/data/species_info_files/species_info.txt | sed '1d' | while read line
do
SP_MIDAS=$(echo $line | cut -f 1 -d' ')
SP_PATRIC=$(echo $line | cut -f 2 -d' ')
wget -N "ftp://ftp.patricbrc.org/genomes/$SP_PATRIC/$SP_PATRIC.PATRIC.faa";
mv $SP_PATRIC.PATRIC.faa $SP_MIDAS.fasta
done
# housekeeping
mv *.fasta ./fasta_all_midas
cd ./fasta_all_midas # shorten sequences names to avoid clash with psortb
ls | while read line
do
sed 's/   .*//' $line > temp.txt && mv temp.txt $line.edited && rm $line && mv $line.edited $line
done


# ~~~~~~~~~~~~~~~~~~~~~~~~ #
#  2) Run Pannzer on all   #
# ~~~~~~~~~~~~~~~~~~~~~~~~ #

# Have a working installation of Pannzer in $programs_install_dir
# see installation instructions at: http://ekhidna2.biocenter.helsinki.fi/sanspanz/


# Initialise the environment modules
. /etc/profile.d/modules.sh

module load anaconda
source activate mypythonMIDAS #anaconda environement with requirements for Pannzer

cd $local_project_dir/HamiltonRuleMicrobiome_gitRepos/output/pannzer

# run with NR >= 1 && NR <= 5952 or smaller range for quick test
cut -f 1 $local_project_dir/HamiltonRuleMicrobiome_gitRepos/data/species_info_files/species_info.txt | sed '1d'| awk 'NR >= 4 && NR <= 5' | while read line ; do 
    cd $local_project_dir/HamiltonRuleMicrobiome_gitRepos/output/pannzer # at each round of loop, get back to overall pannzer output directory
    SP=$line
    mkdir $SP
    cd ./$SP
  
    if ! python $programs_install_dir/SANSPANZ.3/runsanspanz.py -R -m Pannzer -i $local_project_dir/HamiltonRuleMicrobiome_gitRepos/data/patric/fasta_all_midas/$SP\.fasta -o ',,GO.out,'
    then
        echo "Failed to run panzzer on $line"
        continue
    fi
    echo "done running $line"
    mv GO.out $SP.GO.out
done



# ~~~~~~~~~~~~~~~~~~~~~~~~ #
#  3) Assemble GO slim     #
# ~~~~~~~~~~~~~~~~~~~~~~~~ #

cd $local_project_dir/HamiltonRuleMicrobiome_gitRepos/output/pannzer
cat ./*/* | grep -v 'ontology' | cut -f 2-4 | sort | uniq > $local_project_dir/HamiltonRuleMicrobiome_gitRepos/output/tables/bacteria_go_slim.txt




# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
#  4) Web of science search of social keywords    #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

TI=((microb* OR bacter* OR microorganis* OR micro-organis*) AND (coop* OR social*)

# selecting English
# Reviews
# All field
# filter out reviews that were not included but were not specifically on microbial cooperation. 

# --> Write a table in /HamiltonRuleMicrobiome_gitRepos/output/tables/bacteria_social_keywords.txt


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
#  5) Identify bacterial social GOs   #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

# **  code to execute in R **
# *** Note exact output may depend on version of GO.db used ***

# local_project_dir='/path/to/cloned/repository/'
setwd(local_project_dir)

if (!requireNamespace("BiocManager", quietly = TRUE))
   install.packages("BiocManager")
BiocManager::install("GO.db")

library(GO.db)
library(data.table)


# Some functions to retreive GOs
GO_child <- function(node, ontology) {
  if (ontology == "BP") res <- GOBPOFFSPRING[[node]]
  if (ontology == "CC") res <- GOCCOFFSPRING[[node]]
  if (ontology == "MF") res <- GOMFOFFSPRING[[node]]
  return(res)
}
GO_parent <- function(node, ontology) {
  if (ontology == "BP") res <- GOBPPARENTS[[node]]
  if (ontology == "CC") res <- GOCCPARENTS[[node]]
  if (ontology == "MF") res <- GOMFPARENTS[[node]]
  return(res)
}
GO_ancestors <- function(node, ontology) {
  if (ontology == "BP") res <- GOBPANCESTOR[[node]]
  if (ontology == "CC") res <- GOCCANCESTOR[[node]]
  if (ontology == "MF") res <- GOMFANCESTOR[[node]]
  return(res)
}

# Function to search GOs matching a keyword and retreive all child terms and direct parent term
getGO<- function(keyword, allgos){
  
  # Get all matches of a given keyword
  keyword_retreive<- allgos[allgos$description %like% keyword,] # retreive all GOs that contain the keyword
  
  if(nrow(keyword_retreive) > 0){ # For each keyword match, do:
  keyword_retreive$why = 'keyword match'
  keyword_retreive$original_keyword = keyword
  
  # Get all the children below those matches
  look_for_children<- vector('list', length = nrow(keyword_retreive))
  for(i in 1:length(look_for_children)){
    
    look_for_children[[i]]<- allgos[allgos$GO_id %in% GO_child(keyword_retreive[i,1], ontology = keyword_retreive[i,3]), ]
    if(nrow(look_for_children[[i]]) > 0){
      look_for_children[[i]]$why = paste0('child of ', keyword_retreive[i,1])
      look_for_children[[i]]$original_keyword = keyword
    }
  }
  children<- do.call('rbind', look_for_children)
  
  # Add the direct parent of each keyword match term
  look_for_directparent<- vector('list', length = nrow(keyword_retreive))
  for(i in 1:length(look_for_directparent)){
    
    look_for_directparent[[i]]<- allgos[allgos$GO_id %in% GO_parent(keyword_retreive[i,1], ontology = keyword_retreive[i,3]), ]
    if(nrow(look_for_directparent[[i]]) > 0){
      look_for_directparent[[i]]$why = paste0('parent of ', keyword_retreive[i,1])
      look_for_directparent[[i]]$original_keyword = keyword
    }
  }
  directparent<- do.call('rbind', look_for_directparent)
  
  # Assemble the direct matches, the children and the parent
  # match to same GO_id can occur between direct match children terms or between different terms childs, we keep first occurence of it, which keeps the keyword match instead of the child term hit
  d<- rbind(keyword_retreive, children)
  d2<- d[order(d$GO_id),]
  d3<- d2[as.factor(duplicated(d2$GO_id))==FALSE,]
  
  # Add the ancestors for context
  d3<- rbind(d3, directparent)
  }else{
    d3<- data.frame(GO_id = NA, description = NA, ontology = NA, why = 'no match with keyword', original_keyword = keyword)
  }
  return(d3)
}


# Get the bacteria GO slim
bacteria_go_slim_new<- read.csv(paste0(local_project_dir, '/HamiltonRuleMicrobiome_gitRepos/output/tables/bacteria_go_slim.txt'), header=FALSE, sep = '\t', colClasses = rep('character', 3))
colnames(bacteria_go_slim_new)<- c('ontology', 'GO_id', 'description')
bacteria_go_slim_new$GO_id<- paste0('GO:', bacteria_go_slim_new$GO_id)
bacteria_go_slim_new<- bacteria_go_slim_new[,c(2,3,1)]


# Get the list of social keywords identified with Web of science search
keywords<- read.table(paste0(local_project_dir, '/HamiltonRuleMicrobiome_gitRepos/output/tables/bacteria_social_keywords.txt'), header=TRUE, sep = '\t')


# Fetch the GOs
golist2_UPDATED<- vector('list', length = nrow(keywords))
for(j in 1:nrow(keywords)){
  print(paste0(j, ': ', keywords[j,'keyword']))
  golist2_UPDATED[[j]]<- getGO(keywords[j,'keyword'], bacteria_go_slim_new)
  golist2_UPDATED[[j]]$behaviour<-keywords[j, 'behaviour']
}
golistfinal2_UPDATED<- do.call('rbind', golist2_UPDATED)
golistfinal2_UPDATED<- golistfinal2_UPDATED[order(golistfinal2_UPDATED$GO_id),]
golistfinal_trim_goslim_UPDATED<- golistfinal2_UPDATED[which(duplicated(paste0(golistfinal2_UPDATED$GO_id, golistfinal2_UPDATED$behaviour)) == FALSE),]


write.table(x = golistfinal_trim_goslim_UPDATED, paste0(local_project_dir, '/HamiltonRuleMicrobiome_gitRepos/output/tables/social_go_list_wide.txt'), col.names = TRUE, row.names = FALSE, sep = '\t')


# This table of 'potential' social GO is then manually curated
# details of curation decisions are in 'social_go_list_curation.xls'
# final list is of social GO used in analysis is saved as in 'social_go_list_final.xls'









