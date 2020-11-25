# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
#                    Simonet & McNally 2020                     #
#                      Sporulation scores                       #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

# Computing the sporulation scores for the 101 species of the analysis

# STEPS ARE:
# 1) Run blastp on all the AA fasta files
# 2) compute sporulation scores following method describedd by Browne et al


# ~~~~~~~~~~~~~~~~ #
#  1) run blastp   #
# ~~~~~~~~~~~~~~~~ #

# Have a working installation of blastp located at blastp_install_dir
# uses file ./data/sporulation_genes/sporulation_genes_aa.fasta, containing sequences of the sporulation genes listed by Browne et al. Sequences retreived from NCBI

cd $local_project_dir/HamiltonRuleMicrobiome_gitRepos/output/sporulation_genes_blast

ls $local_project_dir/HamiltonRuleMicrobiome_gitRepos/data/patric/fasta | while read line
do
# sporulation genes are target, aa fasta sequences are query
$blastp_install_dir/blastp -query $local_project_dir/HamiltonRuleMicrobiome_gitRepos/data/patric/fasta/$line -out blastp.output.$(echo $line | cut -f 1 -d'.') -evalue 10 -subject $local_project_dir/HamiltonRuleMicrobiome_gitRepos/data/sporulation_genes/sporulation_genes_aa.fasta -outfmt 6
echo $line
done


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
#  2) Compute sporulation scores  #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

# **  code to execute in R **


#local_project_dir='/path/to/where/repo/is/cloned'

setwd(paste0(local_project_dir, '/HamiltonRuleMicrobiome_gitRepos/output/sporulation_genes_blast'))
source(paste0(local_project_dir, '/HamiltonRuleMicrobiome_gitRepos/scripts/analysis/0_sourced_packages_and_themes.R'))


# DEFINING A FUNCTION TO COMPUTE SPORULATION SCORE FOR A GIVEN SPECIES
# Unclear from H. Browne what did they do. According to their paper they weight the matches with this weighing they established with a machine learning algorithm. But according to his last email, they just count the number of sporulation genes present (ignoring duplicate genes) 

# The function fetches the blasp output for a given species, as well as the table with the list of sporulation genes. Trims the blast table to keep only hits below the evalue cutoff chosen. If there are several hits on a single sporulation genes, keep only the best hit. The score is then simple the number of sporulation genes that had a hit. This number is divided by 66, the maximum possible score. All scores therefore scale between 0 and 1.

compute_score_new<- function(sp, eval_cutoff){
  dat<- read.csv(paste0('blastp.output.',sp), header=FALSE, sep = '\t')
  colnames(dat)<- c('query_peg','target_spo_gene','perc_id','length_al','nb_mistmatch','gaps','start_query','end_query','start_subject','end_subject','eval','bitScore')
  dat$query_peg<- as.character(dat$query_peg)
  dat$target_spo_gene<- as.character(dat$target_spo_gene)
  dat$species<- do.call('rbind', strsplit(dat$query_peg, ';'))[,1]
  
  dat$peg<- dat$query_peg
  #dat$peg<- do.call('rbind', strsplit(dat$query_peg, ';'))[,2]
  #dat$peg<- paste0(do.call('rbind', strsplit(dat$peg, '[|]'))[,1],'|',do.call('rbind', strsplit(dat$peg, '[|]'))[,2])
  dat2<- dat[,c(13,14,2,3,4,11,12)]
  
  # Add on table what the spo genes is
  spoinfo<- read.csv(paste0(local_project_dir, '/HamiltonRuleMicrobiome_gitRepos/data/sporulation_genes/sporulation_genes_info.txt'), header=TRUE, sep = '\t')
  spoinfo<- spoinfo[,c(3,2,4,1)]
  colnames(spoinfo)<- c('target_spo_gene', 'gene_length', 'weighting', 'gene_symbol')
  
  
  spoinfo$target_spo_gene<- as.character(spoinfo$target_spo_gene)
  dat3<- left_join(dat2, spoinfo, 'target_spo_gene')
  
  # Apply cutoff and select best hit
  dat3_trim<- dat3[dat3$eval < eval_cutoff,]

  # If the same peg hits to several genes, keep only the best hit (smallest evalue)
  dat3_trim<- as.data.frame(dat3_trim %>% group_by(peg) %>% filter(eval == min(eval)))
  
  
  #### COMPUTE SCORE 
  # Simple count of uniaue sporulation genes match ('raw')
  rawnb_scaled_score<- length(unique(dat3_trim$target_spo_gene))/66 # if a gene present twice, count as one

  return(rawnb_scaled_score)
  
}


# Computing scores for various evalue cutoff
scores<- data.frame(species = gsub('blastp.output.', '', list.files()),
                    score_30_raw = NA,
                    score_25_raw = NA,
                    score_20_raw = NA,
                    score_10_raw = NA)

for(i in 1:nrow(scores)){
scores$score_30_raw[i] = compute_score_new(scores$species[i], 1e-30)
scores$score_25_raw[i] = compute_score_new(scores$species[i], 1e-25)
scores$score_20_raw[i] = compute_score_new(scores$species[i], 1e-20)
scores$score_10_raw[i] = compute_score_new(scores$species[i], 1e-10)
print(i)
}



# Check sensitivity to evalue cutoff

# We cannot expect all genes to be a sporulation gene obvisouly, or all those species to all have the 66 sporulation genes.
# After some exploration of the data, see that when evalue cutoff set to 1e-20 or 1e-30, get reasonnable number of 10-30 of the 66 sporulation genes being present in a species.

pc2<- ggplot(scores, aes(x = score_10_raw, y = score_30_raw))+geom_point()+ ggtitle('evalue = 1e-10') +mytheme
pc3<- ggplot(scores, aes(x = score_20_raw, y = score_30_raw))+geom_point()+ ggtitle('evalue = 1e-20') +mytheme
pc4<- ggplot(scores, aes(x = score_25_raw, y = score_30_raw))+geom_point()+ ggtitle('evalue = 1e-25') +mytheme


# Comparision with scores from nature paper

nat<- read.csv(paste0(local_project_dir, '/HamiltonRuleMicrobiome_gitRepos/data/sporulation_genes/browne_scores_sp_in_common.csv'), header=TRUE, colClasses = c('character', 'numeric', 'character', 'numeric'))

scores2<- scores %>% right_join(nat, 'species')
head(scores2)


p1<- ggplot(scores2, aes(x = scores2$score_30_raw))+
  geom_histogram()+ggtitle('My scores\n (scores in common only)')+
  mytheme

p2<- ggplot(scores2, aes(x = browne_score))+
  geom_histogram()+ggtitle('Browne et al. scores\n (scores in common only)')+
  geom_vline(xintercept =  0.5, linetype="dotted", col = 'dodgerblue')+
  mytheme

p3<- ggplot(scores2, aes(x = score_30_raw, y = browne_score))+
  geom_point()+xlab('\n My scores')+ylab('\n Brown et al. scores')+
  scale_colour_manual(values = c('gold', 'firebrick', 'grey'))+
  mytheme_leg

p4<- ggplot(scores, aes(x = score_30_raw))+
  geom_histogram()+ggtitle('(scores of all 101 species)')+
  mytheme

pdf(file = paste0(local_project_dir, '/HamiltonRuleMicrobiome_gitRepos/output/figures/additional_figs/Sporulation_scores_comparison.pdf'), width = 10, height = 10)
grid.arrange(p1, p2, p3, p4, ncol = 2)
dev.off()

pdf(file = paste0(local_project_dir, '/HamiltonRuleMicrobiome_gitRepos/output/figures/additional_figs/Sporulation_scores_evalue_sensitivity.pdf'), width = 10, height = 3)
grid.arrange(pc2, pc3, pc4, ncol = 3)
dev.off()


# WRITE TABLE
scores_final<- scores %>%
  select(species, score_30_raw) %>%
  rename(sporulation_score = score_30_raw)

write.table(x = scores, file = paste0(local_project_dir, '/HamiltonRuleMicrobiome_gitRepos/output/sporulation_scores.txt'), col.names = TRUE, row.names = FALSE, quote = FALSE, sep = '\t')




