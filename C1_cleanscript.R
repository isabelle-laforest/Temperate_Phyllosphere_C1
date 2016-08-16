##################################################################################################
############################## PHYLLOSPHERE NATURAL TEMPERATE FOREST #############################
################################################################################################## 

## CHAPTER 1
## AUTEUR: Isabelle Laforest-Lapoine

# PACKAGES TO INSTALL

  # install.packages(“picante”, dependencies=TRUE)
  # install.packages("biom", dependencies=TRUE)
  library(picante)
  library(biom)
  library(ggplot2)
  library(Hmisc)
  library(vegan)

# SETTING WORKING DIRECTORY TO IMPORT FILES FROM THE SERVER

  setwd('/data/users/isabelle/Temperate_phyllo_2014/Chimeras/uclust/')

# DATABASES (FULL 45580 or OTU SELECTION 4784)

  otus.biom.uclust.nc<-read_biom("uclust_non_chimeric_otu_table.biom")
  otus.biom.uclust.nc.select<-read_biom("uclust_non_chimeric_otu_table_more_than_20")

##################################################################################################
######################################### ANALYSES FOR H1 ########################################
##################################################################################################  
  
# IMPORT METADATA
  setwd('/data/users/isabelle/Temperate_phyllo_2014/R_analyses')
  samplemetadataH1<-read.delim("Metadata_by_sample_H1_species_location.txt", row.names=1)
  summary(samplemetadataH1)
  dim(samplemetadataH1)
  row.names(samplemetadataH1)

# IMPORT TAXONOMY

  # Note that we have to do some editing “by hand” of the taxonomy assignments file
  # set tab, ; and _ as separators (can do this in a spreadsheet editor)
  otus.taxonomy.forest <- read.delim("tax_blast_non_chimeric_otu_uclust.txt", row.names=1)
  head(otus.taxonomy.forest)
  summary(otus.taxonomy.forest)

# CONVERT BIOM TABLE TO COMMUNITY MATRIX
  
  otus.comm.uclust.nc<-t(as.matrix(biom_data(otus.biom.uclust.nc)))
  otus.comm.uclust.nc<-as.data.frame(otus.comm.uclust.nc)
  dim(otus.comm.uclust.nc)

  otus.comm.uclust.nc.select<-t(as.matrix(biom_data(otus.biom.uclust.nc.select)))
  otus.comm.uclust.nc.select<-as.data.frame(otus.comm.uclust.nc.select)
  dim(otus.comm.uclust.nc.select)

## HOW TO TAKE OUT ALL CHLOROPLASTS FROM OUR DATA

  dim(otus.taxonomy.forest) # 45580 8
  c<-which(otus.taxonomy.forest$CLASS=="Chloroplast")
  otus.taxonomy.forest<-otus.taxonomy.forest[-c,]
  dim(otus.taxonomy.forest) # 45297 8
  dim(otus.comm.uclust.nc) #219 45580
  dim(otus.comm.uclust.nc.select) #219 4784
  otus.comm.uclust.nc<-otus.comm.uclust.nc[,row.names(otus.taxonomy.forest)]
  otus.taxonomy.forest.select<-otus.taxonomy.forest[colnames(otus.comm.uclust.nc.select),]
  #The next one is not necessary because no chloroplast are part of the smaller dataset
  #otus.comm.uclust.nc.select<-otus.comm.uclust.nc.select[,row.names(otus.taxonomy.forest.select)]
  dim(otus.comm.uclust.nc.select) #219 4784
  dim(otus.comm.uclust.nc) # 219 45272

  # Some datanames are wrong

    which(row.names(otus.comm.uclust.nc.select)=="AS9ACSASutton4")
    which(row.names(otus.comm.uclust.nc.select)=="AS10ACSASutton5")
    which(row.names(otus.comm.uclust.nc.select)=="AS11ACSASutton6")
    row.names(otus.comm.uclust.nc.select)[175]<-"AS9ACSASutton1"
    row.names(otus.comm.uclust.nc.select)[114]<-"AS10ACSASutton2"
    row.names(otus.comm.uclust.nc.select)[76]<-"AS11ACSASutton3"

    which(row.names(otus.comm.uclust.nc)=="AS9ACSASutton4")
    which(row.names(otus.comm.uclust.nc)=="AS10ACSASutton5")
    which(row.names(otus.comm.uclust.nc)=="AS11ACSASutton6")
    row.names(otus.comm.uclust.nc)[175]<-"AS9ACSASutton1"
    row.names(otus.comm.uclust.nc)[114]<-"AS10ACSASutton2"
    row.names(otus.comm.uclust.nc)[76]<-"AS11ACSASutton3"

# UNIFRAC
  
  unifrac.ur.unwt <- as.dist(read.delim("unweighted_unifrac_rarefied_otu_table_even4500.txt", row.names=1))
  unifrac.ur.wt <- as.dist(read.delim("weighted_unifrac_rarefied_otu_table_even4500.txt", row.names=1))
  which(row.names(as.matrix(unifrac.ur.unwt))=="AS9ACSASutton4")
  which(row.names(as.matrix(unifrac.ur.unwt))=="AS10ACSASutton5")
  which(row.names(as.matrix(unifrac.ur.unwt))=="AS11ACSASutton6")
  unifrac.ur.unwt<-as.matrix(unifrac.ur.unwt)
  row.names(unifrac.ur.unwt)[175]<-"AS9ACSASutton1"
  row.names(unifrac.ur.unwt)[114]<-"AS10ACSASutton2"
  row.names(unifrac.ur.unwt)[76]<-"AS11ACSASutton3"
  unifrac.ur.unwt<-as.dist(unifrac.ur.unwt)
  which(row.names(as.matrix(unifrac.ur.wt))=="AS9ACSASutton4")
  which(row.names(as.matrix(unifrac.ur.wt))=="AS10ACSASutton5")
  which(row.names(as.matrix(unifrac.ur.wt))=="AS11ACSASutton6")
  unifrac.ur.wt<-as.matrix(unifrac.ur.wt)
  row.names(unifrac.ur.wt)[175]<-"AS9ACSASutton1"
  row.names(unifrac.ur.wt)[114]<-"AS10ACSASutton2"
  row.names(unifrac.ur.wt)[76]<-"AS11ACSASutton3"
  unifrac.ur.wt<-as.dist(unifrac.ur.wt)
  
## MAKE SURE ALL SAMPLES MATCH IN DIFFERENT FILES BEFORE PROCEEDING

  # Verify at this step that there are no NA's in the row.names of comm.sub because of mistmatches (ACSASutton4-5-6)
  comm.sub <- otus.comm.uclust.nc[rownames(samplemetadataH1),]
  dim(comm.sub)
  row.names(comm.sub)
  both <- intersect(rownames(comm.sub), rownames(samplemetadataH1))
  comm.taxo <- otus.taxonomy.forest[colnames(comm.sub),]
  unifrac.wt.sub <- as.dist(as.matrix(unifrac.ur.wt)[both,both])
  unifrac.unwt.sub <- as.dist(as.matrix(unifrac.ur.unwt)[both,both])

  comm.sub.s <- otus.comm.uclust.nc.select[rownames(samplemetadataH1),]
  dim(comm.sub.s)
  row.names(comm.sub.s)

#MAKE CLASS-LEVEL COMMUNITY MATRIX

  comm.taxoClass <- aggregate(t(comm.sub), by=list(class=comm.taxo$CLASS), sum)
  comm.taxoClass <- comm.taxoClass[-1,]
  rownames(comm.taxoClass) <- comm.taxoClass[,1]
  comm.taxoClass <- t(comm.taxoClass[,-1])
  comm.taxoClass <- comm.taxoClass[rownames(comm.sub),]

## RAREFIED TAXONOMY CLASS ANALYSIS
  
  comm.taxoClass.r4500 <- rrarefy(comm.taxoClass, sample=4000)
  comm.taxoClass.r4500.mds <- metaMDS(comm.taxoClass.r4500)
  comm.taxoClass.sub.r4500.agg <- aggregate(decostand(comm.taxoClass.r4500,method="total"), by=list(species=samplemetadataH1$species, site=samplemetadataH1$location), mean)
  write.csv(t(comm.taxoClass.sub.r4500.agg), "comm.taxoClass.r4500.agg.csv")

#MAKE SPECIES-LEVEL COMMUNITY MATRIX

  comm.taxoSpecies <- aggregate(t(comm.sub), by=list(class=comm.taxo$SPECIES), sum)
  comm.taxoSpecies <- comm.taxoSpecies[-1,]
  rownames(comm.taxoSpecies) <- comm.taxoSpecies[,1]
  comm.taxoSpecies <- t(comm.taxoSpecies[,-1])
  comm.taxoSpecies <- comm.taxoSpecies[rownames(comm.sub),]

# GET INFORMATION ABOUT SEQUENCES IN SAMPLE

  hist(apply(comm.sub,1,sum))
  min(apply(comm.sub,1,sum)) #4574
  max(apply(comm.sub,1,sum)) #86276
  
  hist(apply(comm.sub.s,1,sum))
  min(apply(comm.sub.s,1,sum)) #4420
  max(apply(comm.sub.s,1,sum)) #86481

# RAREFACTION OF THE DATA

  comm.sub.rare4500 <- rrarefy(comm.sub, sample=4000)
  comm.sub.s.rare4500 <- rrarefy(comm.sub.s, sample=4000)  

  # THIS IS AN AUTOMATIZED LOOP THAT WILL PRODUCE MANY RAREFACTION TO ASSESS THE ROBUSTNESS OF OUR DATA
  comm.sub.rare4500_2 <- rrarefy(comm.sub, sample=4000)
  test<-list(comm.sub.rare4500,comm.sub.rare4500_2)

  for (i in 1:10){
    test[[i]]<-rrarefy(comm.sub, sample=4000)
  }

  for (i in 11:20){
    test[[i]]<-rrarefy(comm.sub, sample=4000)
  }

  for (i in 21:30){
    test[[i]]<-rrarefy(comm.sub, sample=4000)
  }

  for (i in 31:40){
    test[[i]]<-rrarefy(comm.sub, sample=4000)
  }

  for (i in 41:50){
    test[[i]]<-rrarefy(comm.sub, sample=4000)
  }

  for (i in 51:60){
    test[[i]]<-rrarefy(comm.sub, sample=4000)
  }

  for (i in 61:70){
    test[[i]]<-rrarefy(comm.sub, sample=4000)
  }

  for (i in 71:80){
    test[[i]]<-rrarefy(comm.sub, sample=4000)
  }

  for (i in 81:90){
    test[[i]]<-rrarefy(comm.sub, sample=4000)
  }

  for (i in 91:100){
    test[[i]]<-rrarefy(comm.sub, sample=4000)
  }

# ORDINATION ON THE RAREFIED DATA

  comm.sub.rare4500.mds <- metaMDS(comm.sub.rare4500)
  comm.sub.s.rare4000.mds <- metaMDS(comm.sub.s.rare4500)

  # THIS IS AN AUTOMATIZED LOOP THAT WILL PRODUCE MANY RAREFACTION TO ASSESS THE ROBUSTNESS OF OUR DATA
  comm.sub.rare4500.mds <- metaMDS(test[[1]])
  comm.sub.rare4500.mds_2 <- metaMDS(test[[2]])
  test_2<-list(comm.sub.rare4500.mds,comm.sub.rare4500.mds_2)
  
  for (i in 3:100){
    
    test_2[[i]] <- metaMDS(test[[i]])
    
  }

# ADONIS: TEST ANOVA-LIKE INTERACTIONS AMONG VARIABLES (R2=65%)

  model_test.s<-adonis(comm.sub.rare4500 ~ species*location*moment, data=samplemetadataH1, permutations = 999)
  attributes(model_test.s)
  model_test.s$aov.tab
  
  # THIS IS AN AUTOMATIZED LOOP THAT WILL PRODUCE MANY RAREFACTION TO ASSESS THE ROBUSTNESS OF OUR DATA
  for (i in 1:41){
  model_test<-adonis(test[[i]] ~ species*location*moment, data=samplemetadataH1)
    print(model_test$aov.tab)
  }

##################################################################################################
########################################## PLOTS FOR H1 ##########################################
##################################################################################################

  # Make a plot with all sites and all species (Mess)
  test <- qplot(scores(comm.sub.rare4500.mds)[,1], scores(comm.sub.rare4500.mds)[,2], colour=species, shape=location, group=ind_id, geom="point", xlab="NMDS 1", ylab="NMDS 2", data=samplemetadataH1)
  test + geom_path() + theme_bw() + theme(axis.title=element_text(size="12", color="black"), legend.position="right")

  # Make a plot only for sites
  colnames(samplemetadataH1)[2]<-"Site"
  test <- qplot(scores(comm.sub.s.rare4000.mds)[,1], scores(comm.sub.s.rare4000.mds)[,2], color=Site, shape=Site, geom="point", xlab="NMDS 1", ylab="NMDS 2", data=samplemetadataH1)
  test + stat_ellipse(level=0.68, geom = "polygon", alpha = 1/2, aes(fill=Site)) + theme_bw() + theme(axis.title=element_text(size="12", color="black"), legend.position="right") + scale_color_manual(values = c("dodgerblue", "darkgreen", "darkorange", "red3")) + scale_fill_manual(values = c("dodgerblue", "darkgreen", "darkorange", "red3"))

  test <- qplot(scores(comm.sub.s.rare4000.mds)[,1], scores(comm.sub.s.rare4000.mds)[,2], color=Site, shape=Site, geom="point", xlab="NMDS 1", ylab="NMDS 2", data=samplemetadataH1)
  test + stat_ellipse(level=0.68) + theme_bw() + theme(axis.title=element_text(size="12", color="black"), legend.position="right") + scale_color_manual(values = c("dodgerblue", "darkgreen", "darkorange", "red3")) + scale_fill_manual(values = c("dodgerblue", "darkgreen", "darkorange", "red3"))


  # Make a plot only for species
  colnames(samplemetadataH1)[1]<-"Species"
  test <- qplot(scores(comm.sub.rare4500.mds)[,1], scores(comm.sub.rare4500.mds)[,2], shape=Species, geom="point", xlab="NMDS 1", ylab="NMDS 2", data=samplemetadataH1)
  test + stat_ellipse(level=0.68, geom = "polygon", alpha = 1/2, aes(fill = Species)) + theme_bw() + theme(axis.title=element_text(size="12", color="black"), legend.position="right")

  test <- qplot(scores(comm.sub.rare4500.mds)[,1], scores(comm.sub.rare4500.mds)[,2], shape=Species, geom="point", xlab="NMDS 1", ylab="NMDS 2", data=samplemetadataH1)
  test + stat_ellipse(level=0.68, geom = "polygon", alpha = 1/2, aes(fill = Species)) + theme_bw() + theme(axis.title=element_text(size="12", color="black"), legend.position="right") + scale_fill_manual(values = c("gray90", "black", "gray25", "gray66", "gray85")) + scale_color_manual(values = c("gray0", "gray0", "gray0", "gray0", "gray0")) + scale_size(guide = "none")

  # Make a plot only for months
  colnames(samplemetadataH1)[3]<-"Time"
  test <- qplot(scores(comm.sub.rare4500.mds)[,1], scores(comm.sub.rare4500.mds)[,2], shape=Time, geom="point", xlab="NMDS 1", ylab="NMDS 2", data=samplemetadataH1)
  test + stat_ellipse(level=0.68) + theme_bw() + theme(axis.title=element_text(size="12", color="black"), legend.position="top") + scale_color_manual(values = c("darkorchid", "dodgerblue4", "brown2"))

  test <- qplot(scores(comm.sub.s.rare4000.mds)[,1], scores(comm.sub.s.rare4000.mds)[,2], colour=Time, shape=Time, geom="point", xlab="NMDS 1", ylab="NMDS 2", data=samplemetadataH1)
  test + stat_ellipse(level=0.68) + theme_bw() + theme(axis.title=element_text(size="12", color="black"), legend.position="right") + scale_color_manual(values = c("blue1", "red1", "green3"))

##################################################################################################
######################################## EXTRACT MORE STATS ######################################
################################################################################################## 
  
# NUMBER OF OTUs (45272)

  dim(comm.sub)

# MEAN OTU RICHNESS PER TREE
  x<-matrix(nrow=1, ncol=142)

  for (i in 1:length(rownames(comm.sub))){
    x[1,i]<-length(which(comm.sub[i,]>0))
  }

  mean(x) # 1300.373
  se<-sqrt(var(as.numeric(x))/length(x))
  se#55.77932

# HOW MANY OTUs ARE REALLY IN THE 142 SAMPLES

  z<-decostand(comm.sub, method="pa")
  OTUs<-apply(z,2,sum)
  u<-which(OTUs>0)
  len<-length(u)#32753

##################################################################################################
######################################### CORE MICROBIOME ########################################
##################################################################################################

  comm.taxo$occ<-OTUs/dim(comm.sub)[1]
  sum_seq_per_otus<-apply(comm.sub,2,sum)
  comm.taxo$seq<-sum_seq_per_otus
  tail(comm.taxo[order(comm.taxo$occ),],n=85)
  write.csv(tail(comm.taxo[order(comm.taxo$occ),],n=85),"coremicrobiome.csv")

# How the 20 most common OTUs distribute across samples

  write.csv(comm.sub[,row.names(tail(comm.taxo[order(comm.taxo$occ),],n=20))],"distribution_of_coremicrobiome_per_sample.csv")

# % of OTUs occuring only on one species

# First let's calculate how many OTUs there are in each samples
z<-decostand(comm.sub, method="pa")
order<-apply(z,1,sum)
order_OTUs_eachsample<-order[order(order)]

# Then calculate the percent that are only occuring in one sample
OTU_order<-apply(z,2,sum)
percent_OTU_on_one_sample<-length(which(OTU_order==1))/len
percent_OTU_on_one_sample*100

# How many OTUs are only present in one species (CUSTOM)

beta.comm<-decostand(comm.sub, method="pa")

Acru<-beta.comm[which(samplemetadataH1$species=="ACRU"),]
Acsa<-beta.comm[which(samplemetadataH1$species=="ACSA"),]
Abba<-beta.comm[which(samplemetadataH1$species=="ABBA"),]
Bepa<-beta.comm[which(samplemetadataH1$species=="BEPA"),]
Pisp<-beta.comm[which(samplemetadataH1$species=="Pisp"),]

sumAcru<-decostand(t(as.data.frame(apply(Acru, 2, sum))), method="pa")
sumAcsa<-decostand(t(as.data.frame(apply(Acsa, 2, sum))), method="pa")
sumAbba<-decostand(t(as.data.frame(apply(Abba, 2, sum))), method="pa")
sumBepa<-decostand(t(as.data.frame(apply(Bepa, 2, sum))), method="pa")
sumPisp<-decostand(t(as.data.frame(apply(Pisp, 2, sum))), method="pa")

percent_otu_one_species<-matrix(nrow=5,ncol=45272)
percent_otu_one_species[1,]<-sumAcru
percent_otu_one_species[2,]<-sumAcsa
percent_otu_one_species[3,]<-sumAbba
percent_otu_one_species[4,]<-sumBepa
percent_otu_one_species[5,]<-sumPisp

percent_otu_one_species<-as.data.frame(percent_otu_one_species)
row.names(percent_otu_one_species)<-c("ABBA", "ACRU", "ACSA", "BEPA", "PISP")
colnames(percent_otu_one_species)<-colnames(beta.comm)

sum_percent_otu_one_species<-t(as.data.frame(apply(percent_otu_one_species, 2, sum)))
dim(sum_percent_otu_one_species)

# This calculates the number of OTU that are only on one species
length(which(sum_percent_otu_one_species==6))/length(sum_percent_otu_one_species)*100

##################################################################################################
######################################## COLLECTOR'S CURVE #######################################
##################################################################################################

  x<-specaccum(comm.sub, method='exact', permutations=100)
  x<-specaccum(comm.sub, method='collector', permutations=100)
  plot(x, add=FALSE)
  boxplot(x)

  x<-specaccum(comm.sub.s, method='exact', permutations=100)
  x<-specaccum(comm.sub.s, method='collector', permutations=100)
  plot(x, add=FALSE)
  boxplot(x)

##################################################################################################
################################### SHANNON INDEX APHA DIVERSITY #################################
##################################################################################################
  
  x<-diversity(comm.sub, index="shannon", MARGIN=1, base=exp(1))
  class(x)
  write.csv(x, "shannon_index_allOTUs")
  shannon_test<-read.csv("shannon_diversity_test.csv", header=T, sep=";")
  rownames(shannon_test)<-shannon_test$ID
  test1<-anova(lm(SHANNON~SITE*SPECIES*TIME, data=shannon_test))
  test1
  summary(fm1<-aov(SHANNON~SITE*SPECIES*TIME, data=shannon_test))
  post_oc<-TukeyHSD(fm1, "SPECIES", ordered=TRUE)
  post_oc
   
  plot(post_oc)
  plot(lm(SHANNON~SITE*SPECIES*TIME, data=shannon_test))
  
  post_oc2<-TukeyHSD(fm1, "SITE", ordered=TRUE)
  post_oc2
  plot(post_oc2)
  
  post_oc3<-TukeyHSD(fm1, "TIME", ordered=TRUE)
  post_oc3
  plot(post_oc3)

  post_oc4<-TukeyHSD(fm1, "SITE:SPECIES", ordered=TRUE)
  post_oc4
  plot(post_oc4)

  p<-ggplot(shannon_test, aes(SPECIES, SHANNON))
  p+theme_bw()+geom_boxplot(aes(fill=SPECIES)) + scale_fill_manual(values = c("gray100", "black", "gray25", "gray66", "gray90"))
  p+geom_boxplot()+geom_jitter()
  
  q<-ggplot(shannon_test, aes(SITE, SHANNON) )
  q+geom_boxplot(aes(fill=SITE))  
  q+geom_boxplot()+geom_jitter()
  q+geom_boxplot(aes(fill=SPECIES))
  
  r<-ggplot(shannon_test, aes(TIME, SHANNON) )
  r+geom_boxplot(aes(fill=TIME))  
  r+geom_boxplot()+geom_jitter()
  r+geom_boxplot(aes(fill=SPECIES))
  r+geom_boxplot(aes(fill=SITE))
  
  y<-diversity(comm.sub.s, index="shannon", MARGIN=1, base=exp(1))
  y<-as.data.frame(y)
  colnames(y)<-"SHANNON"
  y$SPECIES<-shannon_test$SPECIES
  y$SITE<-shannon_test$SITE
  y$TIME<-shannon_test$TIME
  test2<-anova(lm(SHANNON~SITE*SPECIES*TIME, data=y))
  plot(lm(SHANNON~SITE*SPECIES*TIME, data=y))

##################################################################################################
############################ REPEATED MESURES MODEL ON ALPHA DIVERSITY ###########################
##################################################################################################
  
  library(lme4)
  m2 <- lmer(SHANNON ~ 1|SUBJECT, shannon_test, REML=FALSE)
  summary(m2)
  m3 <- lmer(SHANNON ~ (1|SUBJECT)+SPECIES, shannon_test, REML=FALSE)
  summary(m3)
  m4 <- lmer(SHANNON ~ (1|SUBJECT)+SPECIES+SITE, shannon_test, REML=FALSE)
  summary(m4)
  m5 <- lmer(SHANNON ~ (1|SUBJECT)+SPECIES+SITE+TIME, shannon_test, REML=FALSE)
  summary(m5)

##################################################################################################
################################ UNIFRAC: PHYLOGENETIC DIVERSITY #################################
##################################################################################################
  
  # Make sure to load from Unifrac_14_Novembre.R
  unifrac.wt.sub
  unifrac.wt.sub<-as.matrix(unifrac.wt.sub)
  unifrac.unwt.sub
  unifrac.r4500.wt.mds <- monoMDS(unifrac.wt.sub)
  unifrac.r4500.unwt.mds <- monoMDS(unifrac.unwt.sub)
  plot(envfit(unifrac.r4500.wt.mds, samplemetadataH1[,1:2], na.rm=TRUE),p.max=0.05)
  ordiplot(unifrac.r4500.unwt.mds)
  plot(envfit(unifrac.r4500.unwt.mds, samplemetadataH1[,7:12], na.rm=TRUE),p.max=0.05)
  
  # Make a plot only for location
  test <- qplot(scores(unifrac.r4500.wt.mds)[,1], scores(unifrac.r4500.wt.mds)[,2], colour=location, geom="point", xlab="MDS 1", ylab="MDS 2", data=samplemetadataH1)
  test + stat_ellipse() + theme_grey() + theme(axis.title=element_text(face="bold.italic", size="12", color="black"), legend.position="right", legend.background=element_rect(fill="grey95")) + scale_color_manual(values = c("red", "darkgreen", "firebrick", "deeppink", "goldenrod", "springgreen"))
  
  test <- qplot(scores(unifrac.r4500.unwt.mds)[,1], scores(unifrac.r4500.unwt.mds)[,2], colour=location, geom="point", xlab="MDS 1", ylab="MDS 2", data=samplemetadataH1)
  test + stat_ellipse() + theme_grey() + theme(axis.title=element_text(face="bold.italic", size="12", color="black"), legend.position="right", legend.background=element_rect(fill="grey95")) + scale_color_manual(values = c("red", "darkgreen", "firebrick", "deeppink", "goldenrod", "springgreen"))
  
  # Make a plot only for species
  test <- qplot(scores(unifrac.r4500.wt.mds)[,1], scores(unifrac.r4500.wt.mds)[,2], colour=species, geom="point", xlab="MDS 1", ylab="MDS 2", data=samplemetadataH1)
  test + stat_ellipse() + theme_grey() + theme(axis.title=element_text(face="bold.italic", size="12", color="black"), legend.position="right", legend.background=element_rect(fill="grey95")) + scale_color_manual(values = c("red", "darkgreen", "firebrick", "deeppink", "goldenrod", "springgreen"))
  
  test <- qplot(scores(unifrac.r4500.unwt.mds)[,1], scores(unifrac.r4500.unwt.mds)[,2], colour=species, geom="point", xlab="MDS 1", ylab="MDS 2", data=samplemetadataH1)
  test + stat_ellipse() + theme_grey() + theme(axis.title=element_text(face="bold.italic", size="12", color="black"), legend.position="right", legend.background=element_rect(fill="grey95")) + scale_color_manual(values = c("red", "darkgreen", "firebrick", "deeppink", "goldenrod", "springgreen"))
  
  # Make a plot without ellipses
  test <- qplot(scores(unifrac.r4500.wt.mds)[,1], scores(unifrac.r4500.wt.mds)[,2], colour=species, shape=location, geom="point", xlab="NMDS 1", ylab="NMDS 2", data=samplemetadataH1)
  test + theme_grey() + theme(axis.title=element_text(face="bold.italic", size="12", color="black"), legend.position="right", legend.background=element_rect(fill="grey95")) + scale_color_manual(values = c("darkorchid", "dodgerblue4", "brown2", "darkgreen", "firebrick"))
  
  test <- qplot(scores(unifrac.r4500.unwt.mds)[,1], scores(unifrac.r4500.unwt.mds)[,2], colour=species, shape=location, geom="point", xlab="NMDS 1", ylab="NMDS 2", data=samplemetadataH1)
  test + theme_grey() + theme(axis.title=element_text(face="bold.italic", size="12", color="black"), legend.position="right", legend.background=element_rect(fill="grey95")) + scale_color_manual(values = c("darkorchid", "dodgerblue4", "brown2", "darkgreen", "firebrick"))

##################################################################################################
##################################### INTRA-MONTH VARIATION ######################################
##################################################################################################

# THIS IS ON ONLY THE SAMPLES THAT HAVE 3 SAMPLES
  three_reps<-matrix(nrow=34*3,ncol=5)
  vec<-attributes(which(summary(samplemetadataH1$ind_id)==3))$names
  rnames<-matrix(nrow=1,ncol=length(vec)*3)
  x=0
  y=1
  for (i in 1:length(vec)){
    t<-which(samplemetadataH1$ind_id==vec[i])
    y<-x+1
    x<-x+3
    as.matrix(samplemetadataH1[t,1:5])->three_reps[y:x,]
    row.names(samplemetadataH1[t,])->rnames[y:x]
  }
  dim(three_reps)
  three_reps<-as.data.frame(three_reps)
  colnames(three_reps)<-colnames(samplemetadataH1[,1:5])
  row.names(three_reps)<-rnames
  comm.three <- otus.comm.uclust.nc[rownames(three_reps),]
  dim(comm.three)
  comm.three.rare4500 <- rrarefy(comm.three, sample=4500)
  comm.three.rare4500.mds <- metaMDS(comm.three.rare4500)

  model_test.three<-adonis(comm.three.rare4500 ~ species*location*moment, data=three_reps)
  model_test.three$aov.tab

  # Make a plot only for months with lines between repetitive measures
  test <- qplot(scores(comm.three.rare4500.mds)[,1], scores(comm.three.rare4500.mds)[,2], colour=species, shape=moment, group=ind_id, geom="point", xlab="NMDS 1", ylab="NMDS 2", data=three_reps)
  test + geom_path() + theme_grey() + theme(axis.title=element_text(size="12", color="black"), legend.position="right") + scale_color_manual(values = c("darkgreen", "orange", "red", "dodgerblue4", "green"))

  test <- qplot(scores(comm.three.rare4500.mds)[,1], scores(comm.three.rare4500.mds)[,2], shape=species, geom="point", xlab="NMDS 1", ylab="NMDS 2", data=three_reps)
  test + geom_path(aes(group=ind_id, colour=species)) + theme_grey() + theme(axis.title=element_text(size="12", color="black"), legend.position="right")

  # THIS FIGURE IS REALLY CONFUSED too many ellipses
  test <- qplot(scores(comm.three.rare4500.mds)[,1], scores(comm.three.rare4500.mds)[,2], shape=species, colour=moment, geom="point", xlab="NMDS 1", ylab="NMDS 2", data=three_reps)
  test + geom_path(aes(group=ind_id)) + stat_ellipse(level=0.68, geom = "polygon", alpha = 1/2, aes(fill = moment)) + theme_grey() + theme(axis.title=element_text(size="12", color="black"), legend.position="right")

  # Make a plot only for sites
  test <- qplot(scores(comm.three.rare4500.mds)[,1], scores(comm.three.rare4500.mds)[,2], colour=location, geom="point", xlab="NMDS 1", ylab="NMDS 2", data=three_reps)
  test + stat_ellipse(level=0.68) + theme_grey() + theme(axis.title=element_text(size="12", color="black"), legend.position="right") + scale_color_manual(values = c("blue", "green", "darkorange", "red3"))

  # Make a plot only for species
  test <- qplot(scores(comm.three.rare4500.mds)[,1], scores(comm.three.rare4500.mds)[,2],color=species, geom="point", xlab="NMDS 1", ylab="NMDS 2", data=three_reps)
  test + stat_ellipse(level=0.68) + theme_grey() + theme(axis.title=element_text(size="12", color="black"), legend.position="right") + scale_color_manual(values=c("darkgreen", "firebrick", "deeppink", "goldenrod", "springgreen"))

  test <- qplot(scores(comm.three.rare4500.mds)[,1], scores(comm.three.rare4500.mds)[,2], shape=species, geom="point", xlab="NMDS 1", ylab="NMDS 2", data=three_reps)
  test + geom_path(aes(group=ind_id, colour=species))+ stat_ellipse(level=0.68, geom = "polygon", alpha = 1/2, aes(fill = species)) + theme_grey() + theme(axis.title=element_text(size="12", color="black"), legend.position="right")

  # Make a plot only for months
  test <- qplot(scores(comm.three.rare4500.mds)[,1], scores(comm.three.rare4500.mds)[,2], colour=moment, geom="point", xlab="NMDS 1", ylab="NMDS 2", data=three_reps)
  test + stat_ellipse(level=0.68) + theme_grey() + theme(axis.title=element_text(size="12", color="black"), legend.position="top") + scale_color_manual(values = c("darkorchid", "dodgerblue4", "brown2"))

  # Plot only BEPA
  bepa<-three_reps[which(three_reps$species=="BEPA"),]
  b<-scores(comm.three.rare4500.mds)[rownames(bepa),]
  test <- qplot(b[,1], b[,2], colour=moment, shape=moment, geom="point", xlab="NMDS 1", ylab="NMDS 2", data=bepa)
  test + geom_path(aes(group=ind_id)) + stat_ellipse(level=0.68,geom = "polygon", alpha = 1/2, aes(fill = moment)) + theme_grey() + theme(axis.title=element_text(size="12", color="black"), legend.position="right")

  # Plot only PISP
  pisp<-three_reps[which(three_reps$species=="Pisp"),]
  pi<-scores(comm.three.rare4500.mds)[rownames(pisp),]
  test <- qplot(pi[,1], pi[,2], colour=species, shape=moment, group=ind_id, geom="point", xlab="NMDS 1", ylab="NMDS 2", data=pisp)
  test + geom_path() + theme_grey() + theme(axis.title=element_text(size="12", color="black"), legend.position="right") + scale_color_manual(values = "darkgreen")

  # Plot only ABBA
  abba<-three_reps[which(three_reps$species=="ABBA"),]
  a<-scores(comm.three.rare4500.mds)[rownames(abba),]
  test <- qplot(a[,1], a[,2], colour=species, shape=moment, group=ind_id, geom="point", xlab="NMDS 1", ylab="NMDS 2", data=abba)
  test + geom_path() + theme_grey() + theme(axis.title=element_text(size="12", color="black"), legend.position="right") + scale_color_manual(values = "forestgreen")

  # Plot only ACSA
  acsa<-three_reps[which(three_reps$species=="ACSA"),]
  ac<-scores(comm.three.rare4500.mds)[rownames(acsa),]
  test <- qplot(ac[,1], ac[,2], colour=species, shape=moment, group=ind_id, geom="point", xlab="NMDS 1", ylab="NMDS 2", data=acsa)
  test + geom_path() + theme_grey() + theme(axis.title=element_text(size="12", color="black"), legend.position="right") + scale_color_manual(values = "red")
  
  # Plot only ACRU
  acru<-three_reps[which(three_reps$species=="ACRU"),]
  acr<-scores(comm.three.rare4500.mds)[rownames(acru),]
  test <- qplot(acr[,1], acr[,2], colour=species, shape=moment, group=ind_id, geom="point", xlab="NMDS 1", ylab="NMDS 2", data=acru)
  test + geom_path() + theme_grey() + theme(axis.title=element_text(face="bold.italic", size="12", color="black"), legend.position="right") + scale_color_manual(values = "darkorange2")

  # THIS IS ON COMPLETE DATA

  # Make a plot only for months with lines between repetitive measures
  test <- qplot(scores(comm.sub.rare4500.mds)[,1], scores(comm.sub.rare4500.mds)[,2], colour=species, shape=moment, group=ind_id, geom="point", xlab="NMDS 1", ylab="NMDS 2", data=samplemetadataH1)
  test + geom_path() + theme_grey() + theme(axis.title=element_text(size="12", color="black"), legend.position="right") + scale_color_manual(values = c("darkgreen", "orange", "red", "dodgerblue4", "green"))

  # Plot only BEPA
  bepa<-samplemetadataH1[which(samplemetadataH1$species=="BEPA"),]
  b<-scores(comm.sub.rare4500.mds)[rownames(bepa),]
  test <- qplot(b[,1], b[,2], colour=species, shape=moment, group=ind_id, geom="point", xlab="NMDS 1", ylab="NMDS 2", data=bepa)
  test + geom_path() + theme_grey() + theme(axis.title=element_text(size="12", color="black"), legend.position="right") + scale_color_manual(values = "dodgerblue4")

  # Plot only PISP
  pisp<-samplemetadataH1[which(samplemetadataH1$species=="Pisp"),]
  pi<-scores(comm.sub.rare4500.mds)[rownames(pisp),]
  test <- qplot(pi[,1], pi[,2], colour=species, shape=moment, group=ind_id, geom="point", xlab="NMDS 1", ylab="NMDS 2", data=pisp)
  test + geom_path() + theme_grey() + theme(axis.title=element_text(size="12", color="black"), legend.position="right") + scale_color_manual(values = "darkgreen")

  # Plot only ABBA
  abba<-samplemetadataH1[which(samplemetadataH1$species=="ABBA"),]
  a<-scores(comm.sub.rare4500.mds)[rownames(abba),]
  test <- qplot(a[,1], a[,2], colour=species, shape=moment, group=ind_id, geom="point", xlab="NMDS 1", ylab="NMDS 2", data=abba)
  test + geom_path() + theme_grey() + theme(axis.title=element_text(size="12", color="black"), legend.position="right") + scale_color_manual(values = "green")

  # Plot only ACSA
  acsa<-samplemetadataH1[which(samplemetadataH1$species=="ACSA"),]
  ac<-scores(comm.sub.rare4500.mds)[rownames(acsa),]
  test <- qplot(ac[,1], ac[,2], colour=species, shape=moment, group=ind_id, geom="point", xlab="NMDS 1", ylab="NMDS 2", data=acsa)
  test + geom_path() + theme_grey() + theme(axis.title=element_text(size="12", color="black"), legend.position="right") + scale_color_manual(values = "red")

  # Plot only ACRU
  acru<-samplemetadataH1[which(samplemetadataH1$species=="ACRU"),]
  acr<-scores(comm.sub.rare4500.mds)[rownames(acru),]
  test <- qplot(acr[,1], acr[,2], colour=species, shape=moment, group=ind_id, geom="point", xlab="NMDS 1", ylab="NMDS 2", data=acru)
  test + geom_path() + theme_grey() + theme(axis.title=element_text(face="bold.italic", size="12", color="black"), legend.position="right") + scale_color_manual(values = "orange")

##################################################################################################
####################################### PERMANOVAS UNIFRAC #######################################
##################################################################################################

  # UNIFRAC

  model_test<-adonis(unifrac.wt.sub ~ species+location+moment+species*location+location*moment, data=samplemetadataH1)
  attributes(model_test)
  model_test$aov.tab

  model_test<-adonis(unifrac.unwt.sub ~ species*location*moment, data=samplemetadataH1)
  attributes(model_test)
  model_test$aov.tab
  
##################################################################################################
######################################## CLIMATIC VARIABLES ######################################
##################################################################################################

  attributes(samplemetadataH1)
  
  data_climate<-matrix(nrow=dim(samplemetadataH1)[1],ncol=2)
  data_climate<-as.data.frame(data_climate)
  colnames(data_climate)<-c("MEAN_TEMP","TOTAL_RAIN")
 rownames(data_climate)<-row.names(samplemetadataH1)
A<-which(samplemetadataH1$moment==levels(samplemetadataH1$moment)[1])
    for(j in 1:length(A)){
      data_climate$MEAN_TEMP[A[j]]<-samplemetadataH1$MEAN_TEMP_JUNE[A[j]]
      data_climate$TOTAL_RAIN[A[j]]<-samplemetadataH1$TOTAL_RAIN_JUNE[A[j]]
    }
A<-which(samplemetadataH1$moment==levels(samplemetadataH1$moment)[2])
for(j in 1:length(A)){
  data_climate$MEAN_TEMP[A[j]]<-samplemetadataH1$MEAN_TEMP_JULY[A[j]]
  data_climate$TOTAL_RAIN[A[j]]<-samplemetadataH1$TOTAL_RAIN_JULY[A[j]]
}
A<-which(samplemetadataH1$moment==levels(samplemetadataH1$moment)[3])
for(j in 1:length(A)){
  data_climate$MEAN_TEMP[A[j]]<-samplemetadataH1$MEAN_TEMP_AUGUST[A[j]]
  data_climate$TOTAL_RAIN[A[j]]<-samplemetadataH1$TOTAL_RAIN_AUGUST[A[j]]
}

  model_test<-adonis(comm.sub.rare4500~MEAN_TEMP+TOTAL_RAIN, permutations = 999, data=samplemetadataH1,strata=species)
  attributes(model_test)
  model_test$aov.tab
  
##################################################################################################
######################################## FUNCTIONAL TRAITS #######################################
##################################################################################################

  samplemetadataH1[,6:12]->traits

  ordiplot(comm.sub.rare4500.mds, type="none", xlab="NMDS Axis 1", ylab="NMDS Axis 2", main="NMDS on Bray-Curtis Distances")
  #rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "gray")
  points(scores(comm.sub.rare4500.mds), pch=19, col="black", cex=1)
  plot(envfit(comm.sub.rare4500.mds, traits[,6:11], na.rm=TRUE),p.max=0.05, col="red", lwd=1.5)
  plot(envfit(comm.sub.rare4500.mds, comm.taxoClass),p.max=0.001, col="blue", lwd=1.5)

  model_test<-adonis(comm.sub.rare4500~wood_density+Nmass+LMA+seed_mass+avg_max_h+shade_tol+drought_tol,permutations = 999, data=samplemetadataH1, strata=location)
  attributes(model_test)
  model_test$aov.tab

##################################################################################################
############################################## THE END ###########################################
##################################################################################################
