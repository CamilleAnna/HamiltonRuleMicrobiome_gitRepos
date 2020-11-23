# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# COMPUTE RELATEDNESS FROM DIVERSITY ESTIMATES #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

# Set path to directory containing output of script 1.3_RELATEDNESS_midas_snp_diversity
# We use the within host and between host diversity measures to compute 1 - pi which is within and between hosts genomic similarity
# And from this relatedness = (sim_within - sim_between)/(1 - sim_between)

# Load within and between diversity tables
setwd("~/Documents/PhD/Research/HamiltonRuleMicrobiome/HamiltonRuleMicrobiome_gitRepos/output/diversity_sampleDepth5_siteDepth5_sitePrev090//")
source("~/Documents/PhD/Research/background_scripts/basic_packages.R")

pi.within.files<- list.files(pattern = 'within')
pi.between.files<- list.files(pattern = 'between')
pi.within.ls<- vector('list', length = length(pi.within.files))
pi.between.ls<- vector('list', length = length(pi.between.files))

for(i in 1:length(pi.within.files)){
  pi.within.ls[[i]]<- read.table(pi.within.files[i], header=TRUE, stringsAsFactors = FALSE) %>%
    mutate(species = strsplit(pi.within.files[i], '.', fixed = TRUE)[[1]][1])
  
  pi.between.ls[[i]]<- read.table(pi.between.files[i], header=TRUE, stringsAsFactors = FALSE) %>%
    mutate(species = strsplit(pi.between.files[i], '.', fixed = TRUE)[[1]][1])
}

pi.within<- do.call('rbind', pi.within.ls)
pi.between<- do.call('rbind', pi.between.ls)




# Compute within and between genomic similariti
sim.within.MM<- pi.within %>%
  mutate(sim_within = 1-pi_bp) %>%
  select(species, sample_id, sim_within, sites) %>%
  rename(host = sample_id,
         nb_site_within = sites)

sim.between.MM<- pi.between %>%
  mutate(sim_between = 1-pi_bp) %>%
  select(species, sim_between, sites, samples) %>%
  rename(nb_host = samples,
         nb_site_between = sites)

# Compute relatedness
sims.HMPnew.MM<- left_join(sim.within.MM, sim.between.MM, by = 'species') %>%
  mutate(within_host_relatedness = (sim_within - sim_between)/(1 - sim_between)) %>% # point estimate of relatedness (i.e. within sample)
  group_by(species) %>%
  mutate(mean_relatedness = mean(within_host_relatedness)) %>% # mean relatedness
  select(species, host, sim_within, nb_site_within, sim_between, nb_site_between, within_host_relatedness, mean_relatedness, nb_host) %>%
  as.data.frame() %>%
  rename(species_id = species) %>%
  filter(
    nb_site_within > 1000), # species for which no core-genome site could be identified
    nb_host > 1 # need at least 2 hosts to estimate a between host diversity!)
   ) 


length(unique(sims.HMPnew.MM$species_id)) # 141 species
length(unique(sims.HMPnew.MM$species_id)) # 140 after filtering for number of core genomic sites
length(unique(sims.HMPnew.MM$species_id)) # 101 after further filtering for at least two hosts


length(unique(sims.HMPnew.MM$species_id))


write.table(sims.HMPnew.MM, '~/Documents/GitHub/HamiltonRuleMicrobiome/output/RELATEDNESS.txt', col.names = TRUE, row.names = FALSE)

# + write a "species_list.txt" file with list of species

d<- read.table('~/Documents/PhD/Research/HamiltonRuleMicrobiome/HamiltonRuleMicrobiome_gitRepos/output/RELATEDNESS.txt', header=TRUE, sep = ' ')
d<- data.frame(unique(d$species_id))
colnames(d)<- 'species_id'
write.table(d, '~/Documents/PhD/Research/HamiltonRuleMicrobiome/HamiltonRuleMicrobiome_gitRepos/output/species_list.txt', col.names = TRUE, row.names = FALSE, quote = FALSE, sep = '\t')


