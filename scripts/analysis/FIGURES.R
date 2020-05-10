# SETUP ----
# installing/loading the latest installr package:
#install.packages("installr"); library(installr)
#updateR() # updating R.
# install.packages('dplyr')
# install.packages('ggplot2')
# install.packages('bindrcpp')
# install.packages('units')
# install.packages('ggraph')
# install.packages('igraph')
# install.packages('tidyverse')
# install.packages('RColorBrewer')
# install.packages('phytools')
# install.packages("devtools")
# install.packages('ape')
# install.packages('scales')
# install.packages("ggsci")
# install.packages("metafor")
# install.packages('DescTools')
# install.packages('jpeg')

# install.packages("BiocManager")
# BiocManager::install("ggtree", version = "3.8")
#devtools::install_github("eliocamp/ggnewscale@v0.1.0")

#library(ggnewscale)
library(ggsci)
library(dplyr)
library(ggraph)
library(igraph)
library(tidyverse)
library(RColorBrewer)
library(ape) 
library(ggtree)
library(phytools)
library(ggplot2)
library(DescTools)
library(ggthemes)
library(scales)
library(gridExtra)
library(metafor)
library(jpeg)
library(readxl)
library("ape")
library('MCMCglmm')
library('modeest')
library('dplyr')
library('tidyr')
library('readxl')

# LOADING REQUIRED DATA ----
source('/Users/s1687811/Documents/PhD/Research/Hamilton_project/MPGS_December18/scripts_clean/0_useful_functions.R')

setwd('/Users/s1687811/Documents/PhD/Research/Hamilton_project/MPGS_December18/output/Figures/Figs_parts')


setwd('/Users/s1687811/Documents/PhD/Research/Hamilton_project/MPGS_December18/')

# d: full dataset, to use for models using within host values, i.e. (2) uncertainty model and (3) drivers of relatedness
# d_109: same as d but to use when having secretome size as response. We remove the 5 species for which gram profile could not be determined and thus psortB could not be run.
# d_mean: relatedness and abundance averaged by species, for model (1) and (4)
# d_mean_109: same as d_mean but to use for the model where secretome size is the response



toremove<- c('Guyana_massiliensis_60772', 'Lachnospiraceae_bacterium_51870', 'Lachnospiraceae_bacterium_56833', 'Clostridiales_bacterium_56470', 'Clostridiales_bacterium_61057')

d<- as.data.frame(read.csv('data/PROCESSED_DATA_TABLES/data_for_analysis.csv') %>%
                    mutate(species = as.character(species)))

d_109<- as.data.frame(d %>%
                        filter(!is.element(species, toremove)))

d_mean<- as.data.frame(d %>%
                         group_by(species) %>% 
                         mutate(mean_relative_abundance = mean(relative_abundance_within_host)) %>% 
                         select(-host, -within_host_relatedness, -relative_abundance_within_host) %>% 
                         subset(first == TRUE))

d_mean_109<- as.data.frame(d_mean %>%
                             filter(!is.element(species, toremove)))





# Phylogeny
mytree<- read.tree('/Users/s1687811/Documents/PhD/Research/Hamilton_project/MPGS_December18/data/PHYLOGENY/my_tree_MPGS_new')

types<- c('character', 'character', 'numeric', 'numeric', 'character', 'character', 'character', 'character', 'character', 'numeric', 'numeric', 'numeric', 'numeric', 'numeric', 'numeric', 'numeric', 'numeric', 'numeric', 'numeric', 'logical') 

dat<- read.csv('/Users/s1687811/Documents/PhD/Research/Hamilton_project/MPGS_December18/data/PROCESSED_DATA_TABLES/data_for_analysis.csv', colClasses = types)

wideall<- dat[dat$first == TRUE,]


t_col <- function(color, percent = 50, name = NULL) {
  #	  color = color name
  #	percent = % transparency
  #	   name = an optional name for the color
  ## Get RGB values for named color
  rgb.val <- col2rgb(color)
  ## Make new color using input color as base and alpha set by transparency
  t.col <- rgb(rgb.val[1], rgb.val[2], rgb.val[3],
               max = 255,
               alpha = (100-percent)*255/100,
               names = name)
  ## Save the color
  invisible(t.col)
  
}


# FIGURE 1 ----


setwd('/Users/s1687811/Documents/PhD/Research/Hamilton_project/MPGS_December18/output/Figures/Figs_parts/')

# Fig 1A: Manual cartoon ----
# Fig 1B: Manual cartoon ----
# Fig 1C: relatedness values ----


# To add number of host on the plot
dat$n<- paste0('n = ', as.character(dat$nb_host))
dat$n[duplicated(paste0(dat$species,dat$n))]<- ''

# resort dat table to have them sorted by mean relatedness on the plot (plotting code loops through rows of data frame)
means_r<- wideall[,c('species', 'mean_relatedness')]
means_r<- means_r[order(means_r$mean_relatedness),]

dat_fig1<- data.frame(c(species = character(), host = character(), within_host_relatedness = numeric(), relative_abundance_within_host = numeric(), patric_genome_id = character(), patric_genome_name = character(), species_id_plot = character(), species_id_plot_short = character(), gram_profile = character(), sporulation_score = numeric(), total_cds = numeric(), mean_relatedness = numeric(), biofilm = numeric(), ab_degradation = numeric(), quorum_sensing = numeric(), siderophores = numeric(), secretion_system_no4 = numeric(), nb_extracellular = numeric(), nb_host = numeric(), first = character()))

for(s in 1:nrow(means_r)){
  dat_fig1<- rbind(dat_fig1, dat[dat$species == means_r$species[s],])
}


pdf('Fig1C.pdf', width = 4 , height = 9)
par(mar = c(1.5,0,0,0))

axes_lwd = 0.8
joint_lwd = 0.8
values_lwd = 0.5
point_cex = 0.3


ticks_text_cex = 0.5
axis_label_cex = 0.5
nb_host_text_cex = 0.2
sp_names_cex = 0.2


axis_label_font = 2
axis_label_line = -0.6
sp_names_adjust = -8
joint_lines_adjust = -12
joint_lines_type = 3


ylim = c(0, 115)
xlim = c(-2, 2)

plot(0, xaxt = 'n', yaxt = 'n', bty = 'n', pch = '', ylab = '', xlab = '', xlim = xlim, ylim = ylim)


# Axes
range_x_axis = c(-1.1, 1)
range_y_axis = c(ylim[1], ylim[2])


for(y in 1:114){
  lines(y = rep(y, 2), x = c(range_x_axis[1], -(-range_x_axis[1] + (1/115)*(range_x_axis[2] - range_x_axis[1]))))
  recol<- ifelse(y %% 2 == 0, 'gray98', 'gray93')
  rect(xleft = range_x_axis[1], ybottom = y-0.5, xright = range_x_axis[2], ytop = y+0.5, border = NA, col = recol)
}

lines(x = range_x_axis, y = rep(range_y_axis[1], 2))
lines(x = rep(range_x_axis[1], 2), y = range_y_axis)

# Ticks
#x1 = -1.5
x2 = -1
x3 = -0.5
x4 = 0
x5 = 0.5
x6 = 1

#lines(x = rep(x1, 2), y = c(0, -1))
lines(x = rep(x2, 2), y = c(0, -0.5))
lines(x = rep(x3, 2), y = c(0, -0.5))
lines(x = rep(x4, 2), y = c(0, -0.5))
lines(x = rep(x5, 2), y = c(0, -0.5))
lines(x = rep(x6, 2), y = c(0, -0.5))


# Ticks labels
text(y = -3, x = c(x2, x3, x4, x5, x6), labels = c(x2, x3, x4, x5, x6), cex = 0.5)

for(y in 1:114){
  
  spfocus = unique(dat_fig1$species)[y]
  sp = unique(dat_fig1[which(dat_fig1$species == spfocus), 'species_id_plot_short'])
  text(x = -(-range_x_axis[1] + (-3.5/115)*(range_x_axis[2] - range_x_axis[1])) , font = 3, y = y, labels = sp, cex = 0.4, pos = 2)
  
  focustab<- dat_fig1[which(dat_fig1$species == spfocus), ]
  
  q1<- quantile(focustab$within_host_relatedness, 0.25)
  q3<- quantile(focustab$within_host_relatedness, 0.75)
  
  focustab1<- focustab[which(focustab$within_host_relatedness >= q1 & focustab$within_host_relatedness <= q3),]
  focustab2<- focustab[which(focustab$within_host_relatedness < q1 | focustab$within_host_relatedness > q3),]
  
  if(nrow(focustab1) == 0 | nrow(focustab2) == 0){
    for(i in 1:nrow(focustab)){
      lines(x = rep(focustab$within_host_relatedness[i],2) , y = c(y+0.3, y-0.3), col = 'dodgerblue', lwd = values_lwd)
      
    }
  }else{
    # vertical ticks at each within host measure
    for(i in 1:nrow(focustab1)){
      lines(x = rep(focustab1$within_host_relatedness[i],2) , y = c(y+0.3, y-0.3), col = 'dodgerblue', lwd = values_lwd)
      
    }
    
    for(i in 1:nrow(focustab2)){
      lines(x = rep(focustab2$within_host_relatedness[i],2) , y = c(y+0.3, y-0.3), col = 'gray66', lwd = values_lwd)
      
    }
  }
  
  # Point at mean of the distribution
  points(x = mean(focustab$within_host_relatedness), y = y, cex = 0.5, pch = 16)
  
  # Indicate number of hosts
  text(x = 0.93 , font = 1, y = y, labels = focustab[1,'n'], cex = 0.4, pos = 4)
  
}


# Axes labels
mtext('Relatedness', side = 1, cex = 0.5, line = -0.3)

dev.off()



# FIGURE 2 ----

#pdf(file = 'Fig2.pdf', width = 7.2, height = 7.2)

plot_character_2<- function(s, d, l, adj, col, nbtips, trans){
  a<- ((2*pi)/nbtips)
  alpha = s*a
  eps = (d/(tan(asin(d))))-1
  
  if(i <= nbtips/4){
    xt = cos(alpha)
    yt = sin(alpha)
    
    x2 = cos(asin(d) + alpha) * (1 + eps)
    y2 = sin(asin(d) + alpha) * (1 + eps)
    
    x3 = cos(alpha) * (1 + l)
    y3 = sin(alpha) * (1 + l)
    
    x4 = -xt + x2 + x3
    y4 = y2 - yt + y3
    
    x5 = xt + abs(xt - x2)
    y5 = yt - abs(y2 - yt)
    
    x6 = x3 + abs(x5 - xt)
    y6 = y3 - abs(y5 - yt)
    
    
    tx = abs(abs(x4) - abs(x2))
    ty = abs(abs(y4) - abs(y2))
    
    x2 = x2+adj
    x4 = x4+adj
    x6 = x6+adj
    x5 = x5+adj
    
    y2 = y2+adj
    y4 = y4+adj
    y6 = y6+adj
    y5 = y5+adj
    
    
    polygon(x = c(x2 + trans*tx,
                  x4 + trans*tx,
                  x6 + trans*tx,
                  x5 + trans*tx),
            
            y = c( y2 + trans*ty,
                   y4 + trans*ty,
                   y6 + trans*ty,
                   y5 + trans*ty),
            
            col = col, border = NA)
    
  }else if (i >= nbtips/4 && i < nbtips/2){
    xt = cos(alpha)
    yt = sin(alpha)
    
    x2 = cos(asin(d) + alpha) * (1 + eps)
    y2 = sin(asin(d) + alpha) * (1 + eps)
    
    x3 = cos(alpha) * (1 + l)
    y3 = sin(alpha) * (1 + l)
    
    x4 =  x3-abs(abs(xt) - abs(x2))
    y4 =  y3-abs(abs(yt) - abs(y2))
    
    x5 = xt + abs(abs(xt) - abs(x2))
    y5 = yt + abs(abs(y2) - abs(yt))
    
    x6 = x3 + abs(abs(x5) - abs(xt))
    y6 = y3 + abs(abs(y5) - abs(yt))
    
    
    tx = abs(abs(x6) - abs(x5))
    ty = abs(abs(y6) - abs(y5))
    
    x2 = x2-adj
    x4 = x4-adj
    x6 = x6-adj
    x5 = x5-adj
    
    y2 = y2+adj
    y4 = y4+adj
    y6 = y6+adj
    y5 = y5+adj
    
    
    polygon(x = c(x2 - trans*tx,
                  x4 - trans*tx,
                  x6 - trans*tx,
                  x5 - trans*tx),
            
            y = c( y2 + trans*ty,
                   y4 + trans*ty,
                   y6 + trans*ty,
                   y5 + trans*ty),
            
            col = col, border = NA)
    
    
  }else if (i >= nbtips/2 && i < nbtips*0.75){
    xt = cos(alpha)
    yt = sin(alpha)
    
    x2 = cos(asin(d) + alpha) * (1 + eps)
    y2 = sin(asin(d) + alpha) * (1 + eps)
    
    x3 = cos(alpha) * (1 + l)
    y3 = sin(alpha) * (1 + l)
    
    x4 =  x3+abs(abs(xt) - abs(x2))
    y4 =  y3-abs(abs(yt) - abs(y2))
    
    x5 = xt - abs(abs(xt) - abs(x2))
    y5 = yt + abs(abs(y2) - abs(yt))
    
    x6 = x3 - abs(abs(x5) - abs(xt))
    y6 = y3 + abs(abs(y5) - abs(yt))
    
    tx = abs(abs(x6) - abs(x5))
    ty = abs(abs(y6) - abs(y5))
    
    x2 = x2-adj
    x4 = x4-adj
    x6 = x6-adj
    x5 = x5-adj
    
    y2 = y2-adj
    y4 = y4-adj
    y6 = y6-adj
    y5 = y5-adj
    
    
    polygon(x = c(x2 - trans*tx,
                  x4 - trans*tx,
                  x6 - trans*tx,
                  x5 - trans*tx),
            
            y = c( y2 - trans*ty,
                   y4 - trans*ty,
                   y6 - trans*ty,
                   y5 - trans*ty),
            
            col = col, border = NA)
    
    
  }else if (i >= nbtips*0.75 && i <= nbtips){
    xt = cos(alpha)
    yt = sin(alpha)
    
    x2 = cos(asin(d) + alpha) * (1 + eps)
    y2 = sin(asin(d) + alpha) * (1 + eps)
    
    x3 = cos(alpha) * (1 + l)
    y3 = sin(alpha) * (1 + l)
    
    x4 =  x3+abs(abs(xt) - abs(x2))
    y4 =  y3+abs(abs(yt) - abs(y2))
    
    x5 = xt - abs(abs(xt) - abs(x2))
    y5 = yt - abs(abs(y2) - abs(yt))
    
    x6 = x3 - abs(abs(x5) - abs(xt))
    y6 = y3 - abs(abs(y5) - abs(yt))
    
    tx = abs(abs(x4) - abs(x2))
    ty = abs(abs(y4) - abs(y2))
    
    x2 = x2+adj
    x4 = x4+adj
    x6 = x6+adj
    x5 = x5+adj
    
    y2 = y2-adj
    y4 = y4-adj
    y6 = y6-adj
    y5 = y5-adj
    
    
    polygon(x = c(x2 + trans*tx,
                  x4 + trans*tx,
                  x6 + trans*tx,
                  x5 + trans*tx),
            
            y = c( y2 - trans*ty,
                   y4 - trans*ty,
                   y6 - trans*ty,
                   y5 - trans*ty),
            
            col = col, border = NA)
    
  }else{
    print('It seems you are attempting to plot more species than there is in your tree')
  }
  
  
  
  
}


# get the order of species in the same as tree from ggtree
g<- ggtree(mytree)
spl<- as.data.frame(g$data)
spl<- spl[spl$isTip==TRUE,]


# build a color map from it
names(spl)<- c(names(spl)[1:3], 'species', names(spl)[5:9])
spl2<- left_join(spl, wideall, by = 'species')


rvals<- spl2[,c('species', 'mean_relatedness')]
rvals<-setNames(rvals[,2],rvals[,1])
obj<-contMap(mytree,rvals,plot=FALSE)
mycols<- rep('black', 1001)
names(mycols)<- 0:1000
obj$cols[]<-mycols
obj<-setMap(obj,invert=TRUE)



colmin = 'snow'
colmax = 'navy'
bias = 0.01
trait = 'secretion_system_no4'

map2color2<- function(colmin, colmax, bias, trait){
  colortest3<- colorRampPalette(c(colmin, colmax), bias = bias)(length(unique(spl2[,which(colnames(spl2)==trait)][!is.na(spl2[,which(colnames(spl2)==trait)])])))
  test2<- spl2[,c(which(colnames(spl2)== 'species'), which(colnames(spl2)== trait))]
  
  mapping<- data.frame(values = sort(unique(test2[,2])), colors = colortest3)
  colnames(mapping)<- c(trait, 'colors')
  
  
  yo<- data.frame(value = NA, colors = 'black')
  colnames(yo)<- c(trait, 'colors')
  mapping<- rbind(mapping, yo)
  
  
  test2<- left_join(test2, mapping, trait)
  test2$colors<- as.character(test2$colors)
  
  return(list(test2$colors, mapping))
}

my.arc.cladelabel<- function(plotted_obj, phylogeny, node, arc.offset, lab.offset, text_color, arc.lwd, marknode, clade_label, cex, text.offset){
  
  g<- ggtree(phylogeny)
  gdata<- as.data.frame(g$data)
  tree = plotted_obj$tree
  
  objtree <- get("last_plot.phylo", envir = .PlotPhyloEnv)
  
  
  # Add the arc
  arc.cladelabels(tree=obj$tree,' ', node, ln.offset = arc.offset, lwd = 2, mark.node=marknode)
  
  # Add name of clade
  
  h <- max(sqrt(objtree$xx^2 + objtree$yy^2))
  d <- getDescendants(tree, node)
  d <- sort(d[d <= Ntip(tree)])
  deg <- atan(objtree$yy[d]/objtree$xx[d]) * 180/pi
  ii <- intersect(which(objtree$yy[d] >= 0), which(objtree$xx[d] < 0))
  deg[ii] <- 180 + deg[ii]
  ii <- intersect(which(objtree$yy[d] < 0), which(objtree$xx[d] < 0))
  deg[ii] <- 180 + deg[ii]
  ii <- intersect(which(objtree$yy[d] < 0), which(objtree$xx[d] >=0))
  deg[ii] <- 360 + deg[ii]
  
  x0 <- lab.offset * cos(median(deg) * pi/180) * h
  y0 <- lab.offset * sin(median(deg) * pi/180) * h
  
  #text.offset=1+(nchar(clade_label)*0.015)
  text(x = text.offset*x0, y = text.offset*y0, clade_label, col = text_color, srt = atan(y0/x0)*(180/pi), cex = cex, font = 3)
  
}

legendbar<- function(colormap, x, y, ytop, digits, label, cex_ticks, cex_label, font_label, vert_bar_adj, tick_length, label_adjust, ticks_labels_adjust){
  
  # Continuous legend colorbar for relatedness
  colmap<- colormap[[2]][-(nrow(colormap[[2]])),]
  colmap$colors<- as.character(colmap$colors)
  # clbr2<-matrix(ncol=nrow(colmap),nrow=2)
  # clbr2[1,]<-seq(0,1,length.out=nrow(colmap))
  # clbr2[2,]<-seq(0,1,length.out=nrow(colmap))
  
  xb2 = x
  yb2 = y
  yb2top = ytop
  
  
  # image(z=clbr2,y=seq(yb2, yb2+yb2top,length.out=nrow(colmap)),x=c(xb2,xb2+0.05),col=colmap$colors,zlim=c(0,1), yaxt="n",xlab="",ylab="",useRaster=T,cex.lab=1.5,add=T)
  
  
  polygon(c(xb2+vert_bar_adj,xb2+vert_bar_adj), c(yb2, yb2+yb2top), lwd = 1.2)
  polygon(c(xb2+vert_bar_adj,xb2+vert_bar_adj+tick_length), c(yb2,yb2), lwd = 1.2)
  polygon(c(xb2+vert_bar_adj,xb2+vert_bar_adj+tick_length), c(yb2+0.25*yb2top,yb2+0.25*yb2top), lwd = 1.2)
  polygon(c(xb2+vert_bar_adj,xb2+vert_bar_adj+tick_length), c(yb2+0.5*yb2top,yb2+0.5*yb2top), lwd = 1.2)
  polygon(c(xb2+vert_bar_adj,xb2+vert_bar_adj+tick_length), c(yb2+0.75*yb2top,yb2+0.75*yb2top), lwd = 1.2)
  polygon(c(xb2+vert_bar_adj,xb2+vert_bar_adj+tick_length), c(yb2+yb2top,yb2+yb2top), lwd = 1.2)
  
  
  labs<- c(round(min(colmap[,1]),digits),
           round(quantile(colmap[,1], 0.25),digits),
           round(quantile(colmap[,1], 0.5),digits),
           round(quantile(colmap[,1], 0.75),digits),
           round(max(colmap[,1]),digits))
  
  text(y = c(yb2, yb2+0.25*yb2top, yb2+0.5*yb2top, yb2+0.75*yb2top, yb2+yb2top), x = xb2+ticks_labels_adjust, pos = 4, labels = labs, font = 1, cex = cex_ticks)
  
  
  text(y = yb2+0.5*yb2top, x = xb2-label_adjust, labels = label, font = font_label, cex = cex_label, srt = 90)
  
  
}

coltoplot1<- map2color2('snow', 'darkred', 0.01, 'mean_relatedness')
coltoplot2<- map2color2('snow', 'darkorchid4', 0.01, 'nb_extracellular')
coltoplot3<- map2color2('snow', 'navy', 0.01, 'secretion_system_no4')
coltoplot4<- map2color2('snow', 'darkorange', 0.01, 'siderophores')
coltoplot5<- map2color2('snow', 'yellow', 0.01, 'quorum_sensing')
coltoplot6<- map2color2('snow', 'forestgreen', 0.01, 'biofilm')
coltoplot7<- map2color2('snow', 'magenta', 0.01, 'ab_degradation')


spl2_mock<- spl2[c(1,2),]
spl2_mock[c(1,2),]<- NA
spl2_mock[,'nb_extracellular']<- c(0,1)
spl2_mock[,'mean_relatedness']<- c(0,1)

spl2<- rbind(spl2, spl2_mock)

coltoplot1<- map2color2('snow', 'darkred', 0.01, 'mean_relatedness')
coltoplot2<- map2color2('snow', 'darkorchid4', 0.01, 'nb_extracellular')

# Fig 2: tree and clades ----
pdf(file = 'Fig2_tree_and_clades.pdf', width = 10, height = 10)
plot(obj,type="fan",ftype="off",lwd=c(2,6),outline=FALSE,
     xlim=c(-1,1),
     ylim = c(-2.5,1.8), legend=FALSE)


for(i in 1:114){
  plot_character_2(i, 0.02, 0.03, 0.01, coltoplot1[[1]][i] , 114, 0)
  plot_character_2(i, 0.02, 0.03, 0.01, coltoplot2[[1]][i] , 114, 1.2)
  plot_character_2(i, 0.02, 0.03, 0.01, coltoplot3[[1]][i] , 114, 2.4)
  plot_character_2(i, 0.02, 0.03, 0.01, coltoplot4[[1]][i] , 114, 3.6)
  plot_character_2(i, 0.02, 0.03, 0.01, coltoplot5[[1]][i] , 114, 4.8)
  plot_character_2(i, 0.02, 0.03, 0.01, coltoplot6[[1]][i] , 114, 6.0)
  plot_character_2(i, 0.02, 0.03, 0.01, coltoplot7[[1]][i] , 114, 7.2)
  
}
  
  plotted_obj = obj
  phylogeny = mytree
  arc.offset = 1.3
  lab.offset = 1.31
  text_color = 'black'
  arc.lwd = 2
  marknode = FALSE
  cex = 0.8
  
  
  nodes_and_genus<- read.csv('/Users/s1687811/Documents/PhD/Research/Hamilton_project/MPGS_December18/data/PROCESSED_DATA_TABLES/Nodes_genus.csv', colClasses = c('numeric', 'character', 'numeric'))
  nodes_and_genus$text.offset<- 1+(nchar(nodes_and_genus$genus)*0.015)
  
  
  for(i in 1:nrow(nodes_and_genus)){
    my.arc.cladelabel(plotted_obj, phylogeny, node = nodes_and_genus$node[i], arc.offset, lab.offset, 'black', arc.lwd, marknode, clade_label = nodes_and_genus$genus[i], cex, text.offset = nodes_and_genus$text.offset[i])
  }

dev.off()


# Fig 2: legends ----
colors_legend<- data.frame(colors_1 = colorRampPalette(c('snow', 'darkred'), bias = 0.01)(114),
                           colors_2 = colorRampPalette(c('snow', 'darkorchid4'), bias = 0.01)(114),
                           colors_3 = colorRampPalette(c('snow', 'navy'), bias = 0.01)(114),
                           colors_4 = colorRampPalette(c('snow', 'darkorange'), bias = 0.01)(114),
                           colors_5 = colorRampPalette(c('snow', 'yellow'), bias = 0.01)(114),
                           colors_6 = colorRampPalette(c('snow', 'forestgreen'), bias = 0.01)(114),
                           colors_7 = colorRampPalette(c('snow', 'magenta'), bias = 0.01)(114)
                           )


y = -2.5
ytop = 0.65


colors_legend$colors_1<- as.character(colors_legend$colors_1)
colors_legend$colors_2<- as.character(colors_legend$colors_2)
colors_legend$colors_3<- as.character(colors_legend$colors_3)
colors_legend$colors_4<- as.character(colors_legend$colors_4)
colors_legend$colors_5<- as.character(colors_legend$colors_5)
colors_legend$colors_6<- as.character(colors_legend$colors_6)
colors_legend$colors_7<- as.character(colors_legend$colors_7)

clbr2<-matrix(ncol=nrow(colors_legend),nrow=2)
clbr2[1,]<-seq(0,1,length.out=nrow(colors_legend))
clbr2[2,]<-seq(0,1,length.out=nrow(colors_legend))



x = 0.5
y = 0
ytop = 1
digits = 1
cex_ticks = 1.8
cex_label = 1.8
vert_bar_adj = 0.15
tick_length = 0.02
label_adjust = 0.43
ticks_labels_adjust = 0.13
font_label = 1

pdf('Fig2_bar_1.pdf', width = 1.5, height = 4)
par(mar = c(0,0,0,0))
plot(0, xaxt = 'n', yaxt = 'n', bty = 'n', pch = '', ylab = '', xlab = '', xlim = c(0,1), ylim = c(0, 1))
image(z=clbr2,y=seq(0, 1,length.out=nrow(colors_legend)),x=c(0.3,0.5),col=colors_legend$colors_1, zlim=c(0,1), yaxt="n",xlab="",ylab="",useRaster=T,cex.lab=1.5,add=T)
legendbar(coltoplot1, x, y, ytop, digits, label = 'Relatedness', cex_ticks, cex_label, font_label, vert_bar_adj, tick_length, label_adjust, ticks_labels_adjust)
dev.off()


digits = 0

pdf('Fig2_bar_2.pdf', width = 1.5, height = 4)
par(mar = c(0,0,0,0))
plot(0, xaxt = 'n', yaxt = 'n', bty = 'n', pch = '', ylab = '', xlab = '', xlim = c(0,1), ylim = c(0, 1))
image(z=clbr2,y=seq(0, 1,length.out=nrow(colors_legend)),x=c(0.3,0.5),col=colors_legend$colors_2, zlim=c(0,1), yaxt="n",xlab="",ylab="",useRaster=T,cex.lab=1.5,add=T)
legendbar(coltoplot2, x, y, ytop, digits, label = 'Secretome size', cex_ticks, cex_label, font_label, vert_bar_adj, tick_length, label_adjust, ticks_labels_adjust)
dev.off()


pdf('Fig2_bar_3.pdf', width = 1.5, height = 4)
par(mar = c(0,0,0,0))
plot(0, xaxt = 'n', yaxt = 'n', bty = 'n', pch = '', ylab = '', xlab = '', xlim = c(0,1), ylim = c(0, 1))
image(z=clbr2,y=seq(0, 1,length.out=nrow(colors_legend)),x=c(0.3,0.5),col=colors_legend$colors_3, zlim=c(0,1), yaxt="n",xlab="",ylab="",useRaster=T,cex.lab=1.5,add=T)
legendbar(coltoplot3, x, y, ytop, digits, label = 'Secretion systems', cex_ticks, cex_label, font_label, vert_bar_adj, tick_length, label_adjust, ticks_labels_adjust)
dev.off()


pdf('Fig2_bar_4.pdf', width = 1.5, height = 4)
par(mar = c(0,0,0,0))
plot(0, xaxt = 'n', yaxt = 'n', bty = 'n', pch = '', ylab = '', xlab = '', xlim = c(0,1), ylim = c(0, 1))
image(z=clbr2,y=seq(0, 1,length.out=nrow(colors_legend)),x=c(0.3,0.5),col=colors_legend$colors_4, zlim=c(0,1), yaxt="n",xlab="",ylab="",useRaster=T,cex.lab=1.5,add=T)
legendbar(coltoplot4, x, y, ytop, digits, label = 'Siderophores', cex_ticks, cex_label, font_label, vert_bar_adj, tick_length, label_adjust, ticks_labels_adjust)
dev.off()


pdf('Fig2_bar_5.pdf', width = 1.5, height = 4)
par(mar = c(0,0,0,0))
plot(0, xaxt = 'n', yaxt = 'n', bty = 'n', pch = '', ylab = '', xlab = '', xlim = c(0,1), ylim = c(0, 1))
image(z=clbr2,y=seq(0, 1,length.out=nrow(colors_legend)),x=c(0.3,0.5),col=colors_legend$colors_5, zlim=c(0,1), yaxt="n",xlab="",ylab="",useRaster=T,cex.lab=1.5,add=T)
legendbar(coltoplot5, x, y, ytop, digits, label = 'Quorum sensing', cex_ticks, cex_label, font_label, vert_bar_adj, tick_length, label_adjust, ticks_labels_adjust)
dev.off()



pdf('Fig2_bar_6.pdf', width = 1.5, height = 4)
par(mar = c(0,0,0,0))
plot(0, xaxt = 'n', yaxt = 'n', bty = 'n', pch = '', ylab = '', xlab = '', xlim = c(0,1), ylim = c(0, 1))
image(z=clbr2,y=seq(0, 1,length.out=nrow(colors_legend)),x=c(0.3,0.5),col=colors_legend$colors_6, zlim=c(0,1), yaxt="n",xlab="",ylab="",useRaster=T,cex.lab=1.5,add=T)
legendbar(coltoplot6, x, y, ytop, digits, label = 'Biofilm', cex_ticks, cex_label, font_label, vert_bar_adj, tick_length, label_adjust, ticks_labels_adjust)
dev.off()

pdf('Fig2_bar_7.pdf', width = 1.5, height = 4)
par(mar = c(0,0,0,0))
plot(0, xaxt = 'n', yaxt = 'n', bty = 'n', pch = '', ylab = '', xlab = '', xlim = c(0,1), ylim = c(0, 1))
image(z=clbr2,y=seq(0, 1,length.out=nrow(colors_legend)),x=c(0.3,0.5),col=colors_legend$colors_7, zlim=c(0,1), yaxt="n",xlab="",ylab="",useRaster=T,cex.lab=1.5,add=T)
legendbar(coltoplot7, x, y, ytop, digits, label = 'Antibiotic degradation', cex_ticks, cex_label, font_label, vert_bar_adj, tick_length, label_adjust, ticks_labels_adjust)
dev.off()


# FIGURE 3 ----

types = c('character', 'character', 'numeric', 'numeric', 'numeric', 'character', 'character',
          'character', 'character', 'numeric', 'numeric', 'numeric', 'numeric', 
          'numeric', 'numeric', 'numeric', 'numeric', 'numeric', 'numeric', 
          'logical')

dat<- read.csv('/Users/s1687811/Documents/PhD/Research/Hamilton_project/MPGS_December18/data/PROCESSED_DATA_TABLES/data_for_analysis.csv')#, colClasses = types)

dat2<- read.csv('/Users/s1687811/Documents/PhD/Research/Hamilton_project/MPGS_December18/data/PROCESSED_DATA_TABLES/dat_summary_models_new.csv')
dat2<- dat2[dat2$test == 'relatedness only',]




cex.coop = 0.8
cex.main.axes = 0.5
cex.y.subs = 0.4
cex.x.A = 0.8
cex.points.subs = 0.7
cex.points.A = 1


pdf('/Users/s1687811/Documents/PhD/Research/Hamilton_project/MPGS_December18/output/Figures/Figs_parts/Fig3_raw.pdf', width=7.2,height=3.4)


nf<-layout(mat=matrix(c(1,1,1,1,2,5,3,6,4,7),nrow=2,ncol=5))
par(mar=c(2,3,3,2),oma=c(4,4,1,1))
plot(NA,NA,bty="l",xlab="",ylab="",cex=1.5,xlim=c(-1.5,3.8),ylim=c(0.5,7.5), yaxt = 'n', xaxt = 'n')

axis(side=1,at=c(-1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5), cex.axis = cex.x.A)

axis(side=2,at=seq(0,8,by=1),labels=rev(c("",'Secretome', 'Biofilm', 'Quorum \n sensing', 'Secretion \n systems', 'Siderophores', "Antibiotic \n degradation",'Overall effect',"")),las=1,cex.axis=cex.coop,srt=45)


cols = rev(c('darkorchid4', 'forestgreen', 'yellow', 'navy', 'darkorange', 'magenta', "black"))


mtext("a",side=3,adj=-0.2/2,line=1,font=2,cex=1)
mtext("Estimated regression coefficient",side=1,padj=2.8,cex=cex.x.A)
abline(v=0,lty=2)
for (i in 1:7){
  #x<-rev(seq(1,7))[i]
  
  x<-seq(1,7)[i]
  
  lines(c(dat2$hpd_lower[i], dat2$hpd_higher[i]),c(x,x),lwd=1,col=c(cols,"black")[i])
 
   lines(c(dat2$hpd_lower[i], dat2$hpd_lower[i]),c(x-0.1,x+0.1),lwd=1,col=c(cols,"black")[i])
  
  lines(c(dat2$hpd_higher[i], dat2$hpd_higher[i]),c(x-0.1,x+0.1),lwd=1,col=c(cols,"black")[i])
}

points(dat2$effect[1:7], seq(1,7),col=cols,cex=cex.points.A,pch=16)


dat<- read.csv('/Users/s1687811/Documents/PhD/Research/Hamilton_project/MPGS_December18/data/PROCESSED_DATA_TABLES/data_for_analysis.csv')#, colClasses = types)


cols = rev(c('darkorchid4', 'forestgreen', 'yellow', 'navy', 'darkorange', 'magenta'))
cols = rev(cols)


D<-dat[dat$first==TRUE,]

mod<-lm(I(nb_extracellular/total_cds)~ mean_relatedness+gram_profile,data=D)

plot(I(nb_extracellular[gram_profile=="n"]/total_cds[gram_profile=="n"])~ mean_relatedness[gram_profile=="n"],data=D,col=t_col(cols[1]),ylim=c(min(I(nb_extracellular/total_cds),na.rm=T),max(I(nb_extracellular/total_cds),na.rm=T)),bty="l",xlab="",ylab="", cex = cex.points.subs, las = 1)



points(I(nb_extracellular[gram_profile=="p"]/total_cds[gram_profile=="p"])~ mean_relatedness[gram_profile=="p"],data=D,col=t_col(cols[1]),pch=16, cex=cex.points.subs)


abline(coef(mod)[1],coef(mod)[2],lwd=1,lty=2,col=rgb(r=0,g=0,b=0))
abline(coef(mod)[1]+coef(mod)[3],coef(mod)[2],lwd=1,col=rgb(r=0,g=0,b=0))

mtext("Proportion cooperative genes",side=2,padj=-2.5,adj=1.5,cex=cex.x.A,las=3)
mtext("b",side=3,adj=-0.2,line=1,font=2,cex=1)



mod<-lm(I(biofilm/total_cds)~ mean_relatedness,data=D)
plot(I(biofilm/total_cds)~ mean_relatedness,data=D,col=t_col(cols[2]),pch=16,bty="l",xlab="",ylab="",cex=cex.points.subs, las = 1)

abline(coef(mod)[1],coef(mod)[2],lwd=1)
mtext("c",side=3,adj=-0.2,line=1,font=2,cex=1)


mod<-lm(I(quorum_sensing/total_cds)~ mean_relatedness,data=D)
plot(I(quorum_sensing/total_cds)~ mean_relatedness,data=D,col=t_col(cols[3]),pch=16,bty="l",xlab="",ylab="",cex=cex.points.subs, las = 1)
abline(coef(mod)[1],coef(mod)[2],lwd=1)
mtext("d",side=3,adj=-0.2,line=1,font=2,cex=1)

mod<-lm(I(secretion_system_no4/total_cds)~ mean_relatedness,data=D)
plot(I(secretion_system_no4/total_cds)~ mean_relatedness,data=D,col=t_col(cols[4]),pch=16,bty="l",xlab="",ylab="",cex=cex.points.subs, las = 1)
abline(coef(mod)[1],coef(mod)[2],lwd=1)
mtext("e",side=3,adj=-0.2,line=1,font=2,cex=1)

mod<-lm(I(siderophores/total_cds)~ mean_relatedness,data=D)
plot(I(siderophores/total_cds)~ mean_relatedness,data=D,col=t_col(cols[5]),pch=16,bty="l",xlab="",ylab="",cex=cex.points.subs, las = 1)
abline(coef(mod)[1],coef(mod)[2],lwd=1)
mtext("Relatedness",side=1,padj=2.8,cex=1)
mtext("f",side=3,adj=-0.2,line=1,font=2,cex=1)

mod<-lm(I(ab_degradation/total_cds)~ mean_relatedness,data=D)
plot(I(ab_degradation/total_cds)~ mean_relatedness,data=D,col=t_col(cols[6]),pch=16,bty="l",xlab="",ylab="",cex=cex.points.subs, las = 1)
abline(coef(mod)[1],coef(mod)[2],lwd=1)
mtext("g",side=3,adj=-0.2,line=1,font=2,cex=1)


dev.off()




# FIGURE 4 ----


load('/Users/s1687811/Documents/PhD/Research/Hamilton_project/MPGS_December18/output/analyses_output/MODELS_OUTPUT_MEAN_RELATEDNESS_DRIVERS_2502.RData')

#types<- c('character', 'character', 'numeric', 'numeric', 'character', 'character', 'character', 'character', 'character', 'numeric', 'numeric', 'numeric', 'numeric', 'numeric', 'numeric', 'numeric', 'numeric', 'numeric', 'numeric', 'logical') 

#dat<- read.csv('/Users/s1687811/Documents/PhD/Research/Hamilton_project/MPGS_December18/data/PROCESSED_DATA_TABLES/data_for_analysis.csv', colClasses = types)


intercept = mean(m3.mean$Sol[,1])
b_sporulation = mean(m3.mean$Sol[,'sporulation_score'])
b_abundance = mean(m3.mean$Sol[,'mean_relative_abundance'])



pdf('/Users/s1687811/Documents/PhD/Research/Hamilton_project/MPGS_December18/output/Figures/Figs_parts/Fig4.a.pdf', width = 5.5, height = 5.5)
par(mfrow=c(2,2),mar=c(4.5,5,2,2))
m<-seq(0,0.05,length.out=101)
n<-50
r = 1/(n-(n-1)*(1-m)^2)
plot(NA,NA,xlim=c(min(m),max(m)),ylim=c(0,1),xlab="Migration rate",ylab="Relatedness",bty="l",cex.lab=1.2)
lines(r~m,lwd=2)
m<-0.05
n<-seq(1,50,length.out=101)
r = 1/(n-(n-1)*(1-m)^2)
mtext("a",side=3,adj=-0.4,line=0,font=2,cex=1.3)
plot(NA,NA,xlim=c(min(n),max(n)),ylim=c(0,1),xlab="Group size",ylab="Relatedness",bty="l",cex.lab= 1.2)
lines(r~n,lwd=2)
mtext("b",side=3,adj=-0.4,line=0,font=2,cex=1.3)


plot(mean_relatedness ~ sporulation_score,pch=NA,data=dm1,bty="l",xlab="Sporulation score",ylab="Relatedness",cex.lab= 1.2)

points(mean_relatedness ~ sporulation_score,col=t_col("forestgreen",percent=60),data=dm1,pch=1)

abline(intercept+mean(dm1$mean_relative_abundance)* b_abundance,b_sporulation,lwd=2)
mtext("c",side=3,adj=-0.4,line=0,font=2,cex=1.3)

plot(mean_relatedness ~ mean_relative_abundance,pch=NA,data=dm1,bty="l",xlab="Mean relative abundance",ylab="Relatedness",cex.lab= 1.2)
points(mean_relatedness ~ mean_relative_abundance,col=t_col("darkorange",percent=60),data=dm1,pch=1)
#abline(intercept+mean(dat$sporulation_score)*b_sporulation,b_abundance,lwd=2) # not adding line because not significant
mtext("d",side=3,adj=-0.4,line=0,font=2,cex=1.3)

dev.off()




# SUPPLEMENTARY FIGURES AND TABLES ----

# GO CONTRIBUTIONS ----
# See output/various_plots/GO_TERMS_CONTRIBUTIONS.pdf, output of Script 11.3_quantifying_go_traits.R


# ESTIMATES FULL MODEL ----

dat2<- read.csv('/Users/s1687811/Documents/PhD/Research/Hamilton_project/MPGS_December18/data/PROCESSED_DATA_TABLES/dat_summary_models_new.csv')

pdf('/Users/s1687811/Documents/PhD/Research/Hamilton_project/MPGS_December18/output/Figures/Figs_parts/Fig4.b.pdf', width = 5.5, height = 5.5)

d<- dat2
plot(0, pch = '', ylim = c(0.9, 7.3), xlim = c(min(d$hpd_lower), max(d$hpd_higher)), yaxt = 'n', ylab = '', xlab = '', main = '', xaxt = 'n')
rect(min(d$hpd_lower)-0.55, 1-0.3, max(d$hpd_higher)+0.55, 1+0.5, border = NA, col = 'grey90')
rect(min(d$hpd_lower)-0.55, 2-0.5, max(d$hpd_higher)+0.55, 2+0.5, border = NA, col = 'grey100')
rect(min(d$hpd_lower)-0.55, 3-0.5, max(d$hpd_higher)+0.55, 3+0.5, border = NA, col = 'grey90')
rect(min(d$hpd_lower)-0.55, 4-0.5, max(d$hpd_higher)+0.55, 4+0.5, border = NA, col = 'grey100')
rect(min(d$hpd_lower)-0.55, 5-0.5, max(d$hpd_higher)+0.55, 5+0.5, border = NA, col = 'grey90')
rect(min(d$hpd_lower)-0.55, 6-0.5, max(d$hpd_higher)+0.55, 6+0.5, border = NA, col = 'grey100')
rect(min(d$hpd_lower)-0.55, 7-0.5, max(d$hpd_higher)+0.55, 7+0.5, border = NA, col = 'grey90')
#rect(min(d$hpd_lower)-0.55, 8-0.5, max(d$hpd_higher)+0.55, 8+0.2, border = NA, col = 'grey90')


addline<- function(dat, y, col, lwd, off, trait, test, x, cex){
  xs<- c(dat[dat$trait == trait & dat$x == x & dat$test == test, 'hpd_lower'], dat[dat$trait == trait & dat$x == x & dat$test == test, 'hpd_higher'])
  
  ys<- c(y,y)+off
  lines(x = xs, y = ys, col = col, lwd = lwd)
  lines(x = rep(xs[1], 2), y = c(y+off+0.05, y+off-0.05), col = col, lwd = lwd)
  lines(x = rep(xs[2], 2), y = c(y+off+0.05, y+off-0.05), col = col, lwd = lwd)
  points(x = dat[dat$trait == trait & dat$x == x & dat$test == test, 'effect'], y = y+off, col = col, pch = 16, cex = cex)
}


# dat = dfoc
# y = 1
# i = 1
# col = 'black'
# lwd = 2
# off = 0.2
# trait = dfoc$trait[i]
# test = 'relatedness only'
# x = 'mean_relatedness'


dat2$trait<- as.character(dat2$trait)


for(i in 1:7){
  test.1 = 'relatedness only'
  x.1 = 'mean_relatedness'
  focaldat.1<- dat2[dat2$x == x.1 & dat2$test == test.1, ]
  
  test.2 = 'full model'
  x.2 = 'mean_relatedness'
  focaldat.2<- dat2[dat2$x == x.2 & dat2$test == test.2, ]
  
  test.3 = 'full model'
  x.3 = 'mean_relative_abundance'
  focaldat.3<- dat2[dat2$x == x.3 & dat2$test == test.3, ]
  
  test.4 = 'full model'
  x.4 = 'sporulation_score'
  focaldat.4<- dat2[dat2$x == x.4 & dat2$test == test.4, ]
  
  addline(focaldat.1, i, 'dodgerblue', 1.5, 0.4, trait = focaldat.1$trait[i], test = test.1, x = x.1, cex = 0.8)
  addline(focaldat.2, i, 'deepskyblue', 1.5, 0.2, trait = focaldat.2$trait[i], test = test.2, x = x.2, cex = 0.8)
  addline(focaldat.3, i, 'darkorange', 1.5, 0, trait = focaldat.3$trait[i], test = test.3, x = x.3, cex = 0.8)
  addline(focaldat.4, i, 'forestgreen', 1.5, -0.2, trait = focaldat.4$trait[i], test = test.4, x = x.4, cex = 0.8)
  
}

axis(side=2,at=seq(1,7,by=1),labels=rev(c('Secretome', 'Biofilm', 'Quorum \n sensing', 'Secretion \n systems','Siderophores', "Antibiotic \n degradation", 'Overall effect')),las=1,cex.axis=0.6,srt=45)
abline(v=0,lty=2)

axis(side = 1, at = c(-5, 0, 5, 10), cex.axis = 0.8, padj = -1)

mtext("Estimated regression coefficient",side=1,padj=2.8,cex=0.8)


dev.off()


# ED TAB 1 ----


load('./output/analyses_output/model_output_cooperation_simple_and_full.RData')

mods.R<- list(ab_degradation = pmm.new.ab_degradation.R, 
              siderophores = pmm.new.siderophores.R,
              secretion_system_no4 = pmm.new.secretion_system_no4.R,
              quorum_sensing = pmm.new.quorum_sensing.R,
              biofilm = pmm.new.biofilm.R,
              secretome = pmm.new.secretome.R)



model = pmm.new.ab_degradation.R

# Function extracts the summary table for fixed effects and variances.
# for variances replaces the posterior mean by the posterior mode (better to report this for the variance estimates)
# Asks if the model include the gram profile as predictor because if TRUE, there is an additional row in the table. Gives row names accordingly.
# Title will appear on first row of first column to record which model the table is output from.

# Second function does the same but for the full model so has two additional rows for the two additional fixed effects (sporulation score and relative abundance)

groomit.v2<- function(model, title, has.gram){
  
  fix<- as.data.frame(summary(model)$solutions) # fixed effects
  fix$pMCMC<- round(fix$pMCMC, 4)
  
  
  rand<- as.data.frame(summary(model)$Gcovariances) %>% mutate(pMCMC = c(' '), post.mean = posterior.mode(model$VCV[,'species'])) # species phylo, post.mean replaced by mode
  
  resid<- as.data.frame(summary(model)$Rcovariances) %>% mutate(pMCMC = ' ', post.mean = posterior.mode(model$VCV[,'units'])) # residuals (= non phylogenetic species variance in this case)
  
  randresid<- as.data.frame(rbind(fix, rand, resid))
  row.names(randresid)<- NULL
  randresid[,c(1:3)]<- round(randresid[,c(1:3)], 5)
  
  
  
  
  if(has.gram == TRUE){
    str1<- c('Fixed effects', ' ', ' ', ' ', 'Variance structure', ' ')
    str2<- c('Intercept', 'Mean relatedness', 'Log(genome size)', 'Gram positive', 'Phylogenetic', 'Residual (non-phylogenetic)')
  }else{
  str1<- c('Fixed effects', ' ', ' ', 'Variance structure', ' ')
  str2<- c('Intercept', 'Mean relatedness', 'Log(genome size)', 'Phylogenetic', 'Residual (non-phylogenetic)')
  }
  
  tab4<- data.frame(str1, str2, randresid)
  tab4$eff.samp<- round(tab4$eff.samp, 0)
  colnames(tab4)<- c(title, ' ', 'Posterior mean/mode', '95% CI lower', '95% CI higher', 'Effective \n sampling', 'pMCMC')
  
  return(tab4)
}
groomit.v2.fullModel<- function(model, title, has.gram){
  
  fix<- as.data.frame(summary(model)$solutions) # fixed effects
  fix$pMCMC<- round(fix$pMCMC, 4)
  
  
  rand<- as.data.frame(summary(model)$Gcovariances) %>% mutate(pMCMC = c(' '), post.mean = posterior.mode(model$VCV[,'species'])) # species phylo, post.mean replaced by mode
  
  resid<- as.data.frame(summary(model)$Rcovariances) %>% mutate(pMCMC = ' ', post.mean = posterior.mode(model$VCV[,'units'])) # residulas (= non phylogenetic species variance in this case)
  
  randresid<- as.data.frame(rbind(fix, rand, resid))
  row.names(randresid)<- NULL
  randresid[,c(1:3)]<- round(randresid[,c(1:3)], 5)
  
  
  
  
  if(has.gram == TRUE){
    str1<- c('Fixed effects', ' ', ' ', ' ', ' ', ' ', 'Variance structure', ' ')
    str2<- c('Intercept', 'Mean relatedness', 'Mean relative abundance', 'Sporulation score', 'Log(genome size)', 'Gram positive', 'Phylogenetic', 'Residual (non-phylogenetic)')
  }else{
    str1<- c('Fixed effects', ' ', ' ', ' ', ' ', 'Variance structure', ' ')
    str2<- c('Intercept', 'Mean relatedness', 'Mean relative abundance', 'Sporulation score', 'Log(genome size)', 'Phylogenetic', 'Residual (non-phylogenetic)')
  }
  
  tab4<- data.frame(str1, str2, randresid)
  tab4$eff.samp<- round(tab4$eff.samp, 0)
  colnames(tab4)<- c(title, ' ', 'Posterior mean/mode', '95% CI lower', '95% CI higher', 'Effective \n sampling', 'pMCMC')
  
  return(tab4)
}


# Table for model including only relatedness (ED table 1)
write.csv(groomit.v2(mods.R$secretome, 'Cooperative trait: secretome', has.gram = TRUE), 'output/analyses_output/ED_table_1.1.secretome.csv')
write.csv(groomit.v2(mods.R$biofilm, 'Cooperative trait: biofilm', has.gram = FALSE), 'output/analyses_output/ED_table_1.2.biofilm.csv')
write.csv(groomit.v2(mods.R$quorum_sensing, 'Cooperative trait: quorum sensing', has.gram = FALSE), 'output/analyses_output/ED_table_1.3.qs.csv')
write.csv(groomit.v2(mods.R$secretion_system_no4, 'Cooperative trait: secretion systems', has.gram = FALSE), 'output/analyses_output/ED_table_1.4.secretionSystems.csv')
write.csv(groomit.v2(mods.R$siderophores, 'Cooperative trait: siderophores', has.gram = FALSE), 'output/analyses_output/ED_table_1.5.siderophores.csv')
write.csv(groomit.v2(mods.R$ab_degradation, 'Cooperative trait: antibiotic degradation', has.gram = FALSE), 'output/analyses_output/ED_table_1.6.abDegradation.csv')


# ED TAB 2 ----

library(MCMCglmm)
pmcmc_vcv<- function(model_output){
  p2<- model_output$VCV[,2] # phy  cov
  p5<- model_output$VCV[,6] # spe  cov
  
  pmcmc<- function(post){
    p = 2*(sum(post<0)/length(post))
    return(p)
  }
  
  return(c('',pmcmc(p2),'','',pmcmc(p5),'',''))
  
}



# SANITY CHECK
# model_output = pmm_ab_degradation
# 
# # FIXED EFFCTS
# fixed = round(as.vector(summary(model_output)$solutions[,1]), 4) # MEAN
# 
# # PHYLOGENETIC EFFECTS
# p.mean.G = round(as.vector(summary(model_output)$Gcovariances[,1]), 4) # POSTERIOR MEAN
# round(as.vector(apply(model_output$VCV[,c(1:4)], 2, mean)), 4) # Should be the same
# p.mode.G = round(as.vector(posterior.mode(model_output$VCV[,c(1:4)])), 4) # POSTERIOR MODE
# 
# # RESIDUAL = NON-PHYLOGENETIC
# p.mean.R = round(as.vector(summary(model_output)$Rcovariances[,1]), 4)  # POSTERIOR MEAN
# round(as.vector(apply(model_output$VCV[,c(5:10)], 2, mean)), 4) # Should be the same
# p.mode.R = round(as.vector(posterior.mode(model_output$VCV[,c(5:10)])), 4) # POSTERIOR MODE
# 
# c(fixed, p.mode.G[1], p.mean.G[2], p.mode.G[4], p.mode.R[1], p.mean.R[2], p.mode.R[4], p.mode.R[5])
# c(fixed, p.mean.G[1], p.mean.G[2], p.mean.G[4], p.mean.R[1], p.mean.R[2], p.mean.R[4], p.mean.R[5])
# c(fixed, p.mode.G[1], p.mode.G[2], p.mode.G[4], p.mode.R[1], p.mode.R[2], p.mode.R[4], p.mode.R[5])
# summary(model_output)




tabs<- function(model_output, name){
  
  
  fix<- as.data.frame(summary(model_output)$solutions)
  row.names(fix)<- NULL
  
  phy<- as.data.frame(summary(model_output)$Gcovariances[-3,]) 
  phy[1,1]<- posterior.mode(model_output$VCV[,c(1)]) # replace posterior mean by posterior mode for the variance terms
  phy[3,1]<- posterior.mode(model_output$VCV[,c(4)])
  
  
  spe<- as.data.frame(summary(model_output)$Rcovariances[-c(3,6),])
  spe[1,1]<- posterior.mode(model_output$VCV[,c(5)]) # replace posterior mean by posterior mode for the variance terms
  spe[3,1]<- posterior.mode(model_output$VCV[,c(8)])
  spe[4,1]<- posterior.mode(model_output$VCV[,c(9)])
  
  
  physpe<- rbind(phy, spe)
  physpe$pMCMC<- pmcmc_vcv(model_output)
  
  outtab<- rbind(fix, physpe)
  return(outtab)
}


load('/Users/s1687811/Documents/PhD/Research/Hamilton_project/MPGS_December18/output/analyses_output/model_output_cooperation.RData')


library(tibble)

t.secretome<- tabs(pmm_secretome, 'secretome') %>% rownames_to_column(var = ' ')
t.secretome[,1]<- c('intercept relatedness', 'intercept trait', 'gram profile', 'Log(genome size)', 'relatedness', 'relatedness, trait', 'trait', 'relatedness', 'relatedness, trait', 'trait', 'relatedness')
t.secretome$structure<- c('Main effects', ' ', ' ', ' ', 'Phylogenetic (co)-variances', ' ', ' ', 'Species (co)-variances', ' ', ' ', 'Residual')
t.secretome<- t.secretome[,c(7, 1:6)]
colnames(t.secretome)<- c('Cooperative trait: secretome', ' ', 'Posterior mean/mode', '95% CI lower', '95% CI higher', 'Effective sampling', 'pMCMC')


t.secretome[,c(3:5)]<- round(t.secretome[,c(3:5)], 4)
t.secretome[,c(6)]<- round(t.secretome[,c(6)], 0)


write.csv(t.secretome, 'output/analyses_output/ED_table_2.1.secretome.csv')





# All other tabs are processed in the same way (only secretome size was different because has gram negative as well in main factors), use a function
groomit<- function(table, title){
  
  table<- table %>% rownames_to_column(var = ' ')
  
  table[,1]<- c('intercept relatedness', 'intercept trait', 'Log(genome size)', 'relatedness', 'relatedness, trait', 'trait', 'relatedness', 'relatedness, trait', 'trait', 'relatedness')
  
  table$structure<- c('Main effects', ' ', ' ', 'Phylogenetic (co)-variances', ' ', ' ', 'Species (co)-variances', ' ', ' ', 'Residual')
  table<- table[,c(7, 1:6)]
  
  colnames(table)<- c(title, ' ', 'Posterior mean/mode', '95% CI lower', '95% CI higher', 'Effective sampling', 'pMCMC')
  

  table[,c(3:5)]<- round(table[,c(3:5)], 4)
  table[,c(6)]<- round(table[,c(6)], 0)
  
  return(table)
}



write.csv(groomit(tabs(pmm_biofilm, 'biofilm'), 'Cooperative trait: biofilm'), 'output/analyses_output/ED_table_2.2.biofilm.csv')

write.csv(groomit(tabs(pmm_QS, 'quorum_sensing'), 'Cooperative trait: quorum sensing'), 'output/analyses_output/ED_table_2.3.quorumSensing.csv')

write.csv(groomit(tabs(pmm_secretion_system, 'secretion_systems'), 'Cooperative trait: secretion systems'), 'output/analyses_output/ED_table_2.4.secretionSystems.csv')

write.csv(groomit(tabs(pmm_siderophores, 'sideorphores'), 'Cooperative trait: siderophores'), 'output/analyses_output/ED_table_2.5.siderophores.csv')

write.csv(groomit(tabs(pmm_ab_degradation, 'ab_degradation'), 'Cooperative trait: antibiotic degradation'), 'output/analyses_output/ED_table_2.6.abDegradation.csv')




# ED TAB 3 ----

# Uses data loaded and functions declared under ED TAB 1 section

mods.R.RA.SPO<- list(ab_degradation = pmm.new.ab_degradation.R.RA.SPO, 
                     siderophores = pmm.new.siderophores.R.RA.SPO,
                     secretion_system_no4 = pmm.new.secretion_system_no4.R.RA.SPO,
                     quorum_sensing = pmm.new.quorum_sensing.R.RA.SPO,
                     biofilm = pmm.new.biofilm.R.RA.SPO,
                     secretome = pmm.new.secretome.R.RA.SPO)


# Table for model including relatedness + abundance + sporulation (ED Table 3)
write.csv(groomit.v2.fullModel(mods.R.RA.SPO$secretome, 'Cooperative trait: secretome', has.gram = TRUE), 'output/analyses_output/ED_table_3.1.secretome.csv')
write.csv(groomit.v2.fullModel(mods.R.RA.SPO$biofilm, 'Cooperative trait: biofilm', has.gram = FALSE), 'output/analyses_output/ED_table_3.2.biofilm.csv')
write.csv(groomit.v2.fullModel(mods.R.RA.SPO$quorum_sensing, 'Cooperative trait: quorum sensing', has.gram = FALSE), 'output/analyses_output/ED_table_3.3.qs.csv')
write.csv(groomit.v2.fullModel(mods.R.RA.SPO$secretion_system_no4, 'Cooperative trait: secretion systems', has.gram = FALSE), 'output/analyses_output/ED_table_3.4.secretionSystems.csv')
write.csv(groomit.v2.fullModel(mods.R.RA.SPO$siderophores, 'Cooperative trait: siderophores', has.gram = FALSE), 'output/analyses_output/ED_table_3.5.siderophores.csv')
write.csv(groomit.v2.fullModel(mods.R.RA.SPO$ab_degradation, 'Cooperative trait: antibiotic degradation', has.gram = FALSE), 'output/analyses_output/ED_table_3.6.abDegradation.csv')

# ED TAB 4 ----

load('/Users/s1687811/Documents/PhD/Research/Hamilton_project/MPGS_December18/output/analyses_output/model_output_relatedness_drivers.RData')


fix<- as.data.frame(summary(m3)$solutions)
fix$pMCMC<- round(fix$pMCMC, 4)


rand<- as.data.frame(summary(m3)$Gcovariances) %>% mutate(pMCMC = c(' ', ' ', ' '), post.mean = posterior.mode(m3$VCV[,c(1, 2, 3)]))
resid<- as.data.frame(summary(m3)$Rcovariances) %>% mutate(pMCMC = ' ', post.mean = posterior.mode(m3$VCV[,c(4)]))

randresid<- as.data.frame(rbind(fix, rand, resid))
row.names(randresid)<- NULL
randresid[,c(1:3)]<- round(randresid[,c(1:3)], 4)

str1<- c('Fixed effects', ' ', ' ', 'Variance structure', ' ', ' ', ' ')
str2<- c('Intercept', 'sporulation score', 'Within host relative abundance', 'species', 'host', 'phylogenetic', 'residual')

tab4<- data.frame(str1, str2, randresid)
tab4$eff.samp<- round(tab4$eff.samp, 0)
colnames(tab4)<- c(' ', ' ', 'Posterior mean/mode', '95% CI lower', '95% CI higher', 'Effective \n sampling', 'pMCMC')

write.csv(x = tab4, file = './output/analyses_output/ED_table_4.csv')


# ED TAB 5 ----

library(metafor)

# a - Simple model, R only ----
datmeta<- read.csv('/Users/s1687811/Documents/PhD/Research/Hamilton_project/MPGS_December18/data/PROCESSED_DATA_TABLES/data_metaanalysis_simple_model.csv', header=TRUE)

datmeta<- datmeta[datmeta$test == 'relatedness only',]

datmeta$trait_plot<- c('Secretome',  'Secretion \n systems', 'Biofilm', 'Quorum \n sensing', 'Siderophores', "Antibiotic \n degradation")
random<-rma(yi=datmeta$effect,sei=datmeta$sd_post) # random effects model fit by REML


value<- c(random$b,
          random$se,
          random$zval,
          random$pval,
          random$ci.lb,
          random$ci.ub,
          summary(random)$fit.stats['ll',2],
          summary(random)$fit.stats['dev',2],
          summary(random)$fit.stats['AIC',2],
          summary(random)$fit.stats['BIC',2],
          summary(random)$fit.stats['AICc',2],
          random$tau2,
          summary(random)$I2,
          summary(random)$H2)


tab3<- data.frame( Statistic = c('estimate', 'se', 'zval', 'pval', 'ci lower', 'ci upper', 'loglik', 'deviance', 'AIC', 'BIC', 'AICc', 'tau^2', 'I^2', 'H^2'),
                   Value = round(value, 3))

write.csv(tab3, 'output/analyses_output/ED_table_5.1.csv')


# b - Complex model, R only ----
datmeta<- read.csv('/Users/s1687811/Documents/PhD/Research/Hamilton_project/MPGS_December18/data/PROCESSED_DATA_TABLES/data_metaanalysis.csv', header=TRUE)

datmeta$trait_plot<- c('Secretome',  'Secretion \n systems', 'Biofilm', 'Quorum \n sensing', 'Siderophores', "Antibiotic \n degradation")
random<-rma(yi=datmeta$effect,sei=datmeta$sd_post) # random effects model fit by REML


value<- c(random$b,
          random$se,
          random$zval,
          random$pval,
          random$ci.lb,
          random$ci.ub,
          summary(random)$fit.stats['ll',2],
          summary(random)$fit.stats['dev',2],
          summary(random)$fit.stats['AIC',2],
          summary(random)$fit.stats['BIC',2],
          summary(random)$fit.stats['AICc',2],
          random$tau2,
          summary(random)$I2,
          summary(random)$H2)


tab3<- data.frame( Statistic = c('estimate', 'se', 'zval', 'pval', 'ci lower', 'ci upper', 'loglik', 'deviance', 'AIC', 'BIC', 'AICc', 'tau^2', 'I^2', 'H^2'),
                   Value = round(value, 3))


write.csv(tab3, 'output/analyses_output/ED_table_5.2.csv')


# c - simple model, R+RA+SPO ----
datmeta<- read.csv('/Users/s1687811/Documents/PhD/Research/Hamilton_project/MPGS_December18/data/PROCESSED_DATA_TABLES/data_metaanalysis_simple_model.csv', header=TRUE)

datmeta<- datmeta[datmeta$test == 'full model' & datmeta$x == 'mean_relatedness' ,]


datmeta$trait_plot<- c('Secretome',  'Secretion \n systems', 'Biofilm', 'Quorum \n sensing', 'Siderophores', "Antibiotic \n degradation")
random<-rma(yi=datmeta$effect,sei=datmeta$sd_post) # random effects model fit by REML


value<- c(random$b,
          random$se,
          random$zval,
          random$pval,
          random$ci.lb,
          random$ci.ub,
          summary(random)$fit.stats['ll',2],
          summary(random)$fit.stats['dev',2],
          summary(random)$fit.stats['AIC',2],
          summary(random)$fit.stats['BIC',2],
          summary(random)$fit.stats['AICc',2],
          random$tau2,
          summary(random)$I2,
          summary(random)$H2)


tab3<- data.frame( Statistic = c('estimate', 'se', 'zval', 'pval', 'ci lower', 'ci upper', 'loglik', 'deviance', 'AIC', 'BIC', 'AICc', 'tau^2', 'I^2', 'H^2'),
                   Value = round(value, 3))


write.csv(tab3, 'output/analyses_output/ED_table_5.3.csv')






