
Description of output files and file formats from 'merge_midas.py species'

Output files
############
count_reads.txt  
  number of reads mapped to 15 marker genes per species
coverage.txt  
  average read-depth of 15 marker genes per species (total bp of mapped reads/total bp of 15 marker-genes)
relative_abundance.txt  
  values from coverage.txt scaled to sum to 1.0 across species per sample
species_prevalence.txt
  summary stats across species

Output formats
############
count_reads.txt, coverage.txt, relative_abundance.txt
  tab-delimited matrix files
  field names are sample ids
  row names are species ids
species_prevalence.txt
  species_id: species identifier
  mean_coverage: average read-depth of marker-genes for species across samples
  median_coverage: median read-depth of marker-genes for species across samples
  mean_abundance: average relative abundance of marker-genes for species across samples
  median_abundance: median relative abundance of marker-genes for species across samples
  prevalence: proportion of samples where species occured with at least 1.0 read-depth

Additional information for each species can be found in the reference database:
 /exports/csce/eddie/biology/groups/mcnally/camille/programs/MIDAS/midas_db_v1.2
