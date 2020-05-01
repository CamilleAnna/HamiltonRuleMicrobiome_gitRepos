# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# GET RELATIVE ABUNDANCE FROM MIDASE MERGE OUTPUT #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

setwd("~/Documents/GitHub/HamiltonRuleMicrobiome/output/species_merged/")

# RELATEDNESS
relatedness<- read.table('~/Documents/GitHub/HamiltonRuleMicrobiome/output/RELATEDNESS.txt', header=TRUE, stringsAsFactors = FALSE) %>% mutate(species.host = paste0(species_id, '.', host))

mysp<- unique(relatedness$species_id)

# RELATIVE ABUNDANCE
abundance<- read.table('~/Documents/GitHub/HamiltonRuleMicrobiome/output/species_merged/relative_abundance.txt', header=TRUE, stringsAsFactors = FALSE) %>%
  filter(species_id %in% mysp) %>%
  gather('host', 'within_host_relative_abundance', 2:241)  %>%
  mutate(species.host = paste0(species_id, '.', host)) %>%
  filter(species.host %in% relatedness$species.host) %>%
  select(species.host, within_host_relative_abundance)


# SPORULATION SCORES ----
spo.scores<- read.table('~/Documents/GitHub/HamiltonRuleMicrobiome/output/SPORULATION_SCORES_all_midas_species.txt', header=TRUE, sep = '\t', colClasses = c('character', 'numeric')) %>%
  rename(species_id = species) %>%
  filter(species_id %in% mysp)

# GO TERMS TRAIT ----
gos<- read.table('~/Documents/GitHub/HamiltonRuleMicrobiome/output/GO_traits_quantification.txt', header=TRUE, sep = '\t', colClasses = c('character', rep('numeric', 7))) %>%
  rename(species_id = species) %>%
  filter(species_id %in% mysp)


# SECRETOME SIZE ----
ss<- read.table('~/Documents/GitHub/HamiltonRuleMicrobiome/output/SECRETOME_quantification.txt', header=TRUE, sep = '\t', colClasses = c('character', 'numeric')) %>%
  rename(species_id = species) %>%
  filter(species_id %in% mysp) # Those not in my sp are because they don't have gram_profile available

# GRAM PROFILES ----
grams<- read.table('~/Documents/GitHub/HamiltonRuleMicrobiome/data/species_info_files/gram_profiles_db.txt', colClasses = 'character', stringsAsFactors = FALSE, col.names = c('species_id', 'gram_profile')) %>%
  filter(species_id %in% mysp)


# ASSEMBLE ALL ----

# Assemble data for focus on new measures of relatedness
sp.focus<- data.frame(species_id = mysp)

dat<- sp.focus %>%
  left_join(gos, by = 'species_id') %>%
  left_join(grams, by = 'species_id') %>%
  left_join(ss, by = 'species_id') %>%
  left_join(spo.scores, by = 'species_id') %>%
  left_join(relatedness, by = 'species_id') %>%
  left_join(abundance, by = 'species.host')


write.table(dat, '/Users/s1687811/Documents/GitHub/HamiltonRuleMicrobiome/output/ANALYSIS_DATA_ASSEMBLED.txt', col.names = TRUE, row.names = FALSE, quote = FALSE, sep = '\t')


