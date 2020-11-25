# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
#                    Simonet & McNally 2020                     #
#                     Relative abundance                        #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

# The species relative abundance is included in the phylogenetic mixed model
# We extract these directly from MIDAS species_profiles output
#local_project_dir='/path/to/where/repo/is/cloned'

setwd(paste0(local_project_dir, '/HamiltonRuleMicrobiome_gitRepos/output/midas/per_sample/'))

hosts<- list.files()
tabs<- vector('list', length = length(hosts))
for(i in 1:length(hosts)){
  
  tabs[[i]]<- read.table(paste0(hosts[i], '/species/species_profile.txt'), header=TRUE)
  tabs[[i]]$host<- hosts[i]
}
d2<- do.call('rbind', tabs)
d2<- d2[,c(1,5,4)]
colnames(d2)<- c('species', 'host', 'relative_abundance')

write.table(d2, file = paste0(local_project_dir, '/HamiltonRuleMicrobiome_gitRepos/output/relative_abundance.txt'), sep = '\t', col.names = TRUE, row.names = FALSE)

