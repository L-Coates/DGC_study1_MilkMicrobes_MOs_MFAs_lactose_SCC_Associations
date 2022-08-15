#This script describes the steps taken to determine if any of the microbial
#taxa are predictive of log(SCC)

#Loading the necessary libraries
library(phyloseq)
library(qiime2R)
library(lme4)
library(lmerTest)
library(picante)
library(tidyverse)
library(rmcorr)


#STEP 1:Reading in the metadata table
metadata  <- read.csv("../../../DairyGrandChallenge.Study1.2017-2018/ResearchMaterials.DairyGrandChallenge/Cleaned_For_MetadataFile/metadata.csv")

#keeping only the data that is needed. 
metadata<- metadata[,c(1:15)]
metadata <- metadata[-1,]

#Rename sample 5849-D1 to 5849-D2 because cow 5849 actually missed
#the milk collection on November 25, 2017, but milk was collected from this cow
#on November 26, 2017. The same was true for cow 6233. 
metadata$sample.id <- gsub(pattern="5849-D1", replacement="5849-D2", x=metadata$sample.id)
metadata$sample.id <- gsub(pattern="6233-D1", replacement="6233-D2", x=metadata$sample.id)
metadata[metadata$sample.id=="5849-D2","day"] <- 2
metadata[metadata$sample.id=="6233-D2","day"] <- 2

#STEP 2:Reading in the AgSource data and combining with metadata

#SCC measures from 11/21/17, 11/28/17, 01/02/18, 01/05/18, 03/13/18, and 03/19/18
Nov.21.17 <-readxl::read_xlsx("../../ResearchMaterials.DairyGrandChallenge/Materials.FromKableLab/DairyGrandChallenge_SCCdata/KFK08 112017 112117 - Nov 24.xlsx", col_names=FALSE, col_types=c("date", "text", "skip", "skip", "skip", "skip", "numeric", "skip", "text", "skip", "skip", "skip", "skip"))
Nov.28.17  <-readxl::read_xlsx("../../ResearchMaterials.DairyGrandChallenge/Materials.FromKableLab/DairyGrandChallenge_SCCdata/KFK08 112717 112817 - Dec 1.xlsx", col_names=FALSE, col_types=c("date", "text", "skip", "skip", "skip", "skip", "numeric", "skip", "text", "skip", "skip", "skip", "skip"))
Jan.02.18  <-readxl::read_xlsx("../../ResearchMaterials.DairyGrandChallenge/Materials.FromKableLab/DairyGrandChallenge_SCCdata/KFK08 010218.xlsx", col_names=FALSE, col_types=c("date", "text", "skip", "skip", "skip", "skip", "numeric", "skip", "text", "skip", "skip", "skip", "skip"))
Jan.05.18  <-readxl::read_xlsx("../../ResearchMaterials.DairyGrandChallenge/Materials.FromKableLab/DairyGrandChallenge_SCCdata/KFK08 010518.xlsx", col_names=FALSE, col_types=c("date", "text", "skip", "skip", "skip", "skip", "numeric", "skip", "text", "skip", "skip", "skip", "skip"))
March.13.18  <-readxl::read_xlsx("../../ResearchMaterials.DairyGrandChallenge/Materials.FromKableLab/DairyGrandChallenge_SCCdata/KFK08 031218 031318 - March 16.xlsx", col_names=FALSE, col_types=c("date", "text", "skip", "skip", "skip", "skip", "numeric", "skip", "text", "skip", "skip", "skip", "skip", "skip", "skip", "skip", "skip"))
March.19.18  <-readxl::read_xlsx("../../ResearchMaterials.DairyGrandChallenge/Materials.FromKableLab/DairyGrandChallenge_SCCdata/KFK08 031918 032018 - March 23.xlsx", col_names=FALSE, col_types=c("date", "text", "skip", "skip", "skip", "skip", "numeric", "skip", "text", "skip", "skip", "skip", "skip"))

colnames(Nov.21.17)<- c("date", "cow", "SCCx1000.milking1.beforesample", "milking.number")
colnames(Nov.28.17)<- c("date", "cow", "SCCx1000.milking1.aftersample", "milking.number")
colnames(Jan.02.18)<- c("date", "cow", "SCCx1000.milking1.beforesample", "milking.number")
colnames(Jan.05.18)<- c("date", "cow", "SCCx1000.milking1.aftersample", "milking.number")
colnames(March.13.18)<- c("date", "cow", "SCCx1000.milking1.beforesample", "milking.number")
colnames(March.19.18)<- c("date", "cow", "SCCx1000.milking1.aftersample", "milking.number")

Nov.21.17 <- Nov.21.17[-c(1:9),]
Nov.28.17 <- Nov.28.17[-c(1:9),]
Jan.02.18 <- Jan.02.18[-c(1:9),]
Jan.05.18 <- Jan.05.18[-c(1:9),]
March.13.18 <- March.13.18[-c(1:9),]
March.19.18 <- March.19.18[-c(1:9),]

Nov.21.17$date <-as.Date(Nov.21.17$date)
Nov.21.17 <- Nov.21.17[Nov.21.17$date=="2017-11-21",]
length(unique(Nov.21.17$cow))

Nov.28.17$date <-as.Date(Nov.28.17$date)
Nov.28.17 <- Nov.28.17[Nov.28.17$date=="2017-11-28",]
length(unique(Nov.28.17$cow))

Jan.02.18$date <-as.Date(Jan.02.18$date)
Jan.02.18 <- Jan.02.18[Jan.02.18$date=="2018-01-02",]
length(unique(Jan.02.18$cow))

Jan.05.18$date <-as.Date(Jan.05.18$date)
Jan.05.18 <- Jan.05.18[Jan.05.18$date=="2018-01-05",]
length(unique(Jan.05.18$cow))

March.13.18$date <-as.Date(March.13.18$date)
March.13.18 <- March.13.18[March.13.18$date=="2018-03-13",]
length(unique(March.13.18$cow))

March.19.18$date <-as.Date(March.19.18$date)
March.19.18 <- March.19.18[March.19.18$date=="2018-03-19",]
length(unique(March.19.18$cow))

#Selecting only the M1 measurements (these are from the first milking of the day) from Nov 21st (this is the most complete dataset taken closest to the date of collection of the milk samples used for microbiota analysis.)
Nov.21.17 <- Nov.21.17[Nov.21.17$milking.number=="M1",]
Nov.21.17 <- Nov.21.17[, c("cow", "SCCx1000.milking1.beforesample")]
length(unique(Nov.21.17$cow)) 

Nov.28.17 <- Nov.28.17[Nov.28.17$milking.number=="M1",]
Nov.28.17 <- Nov.28.17[, c("cow", "SCCx1000.milking1.aftersample")]
length(unique(Nov.28.17$cow)) 

Jan.02.18 <- Jan.02.18[Jan.02.18$milking.number=="M1",]
Jan.02.18 <- Jan.02.18[, c("cow", "SCCx1000.milking1.beforesample")]
length(unique(Jan.02.18$cow)) 

Jan.05.18 <- Jan.05.18[Jan.05.18$milking.number=="M1",]
Jan.05.18 <- Jan.05.18[, c("cow", "SCCx1000.milking1.aftersample")]
length(unique(Jan.05.18$cow)) 

March.13.18 <- March.13.18[March.13.18$milking.number=="M1",]
March.13.18 <- March.13.18[, c("cow", "SCCx1000.milking1.beforesample")]
length(unique(March.13.18$cow)) 

March.19.18 <- March.19.18[March.19.18$milking.number=="M1",]
March.19.18 <- March.19.18[, c("cow", "SCCx1000.milking1.aftersample")]
length(unique(March.19.18$cow))


#Matching SCC measurements to metadata measurements by period. 
metadata.pd0 <- metadata[metadata$period=="p_0",]
metadata.pd1 <- metadata[metadata$period=="p_1",]
metadata.pd2 <- metadata[metadata$period=="p_2",]

SCC.pd0 <- merge(Nov.21.17, Nov.28.17, by="cow", all=FALSE)
dim(metadata.pd0)
metadata.pd0 <- merge(metadata.pd0, SCC.pd0, by="cow", all=FALSE)
dim(metadata.pd0)

SCC.pd1 <- merge(Jan.02.18, Jan.05.18, by="cow", all=FALSE)
dim(metadata.pd1)
metadata.pd1 <- merge(metadata.pd1, SCC.pd1, by="cow", all=FALSE)
dim(metadata.pd1)

SCC.pd2 <- merge(March.13.18, March.19.18, by="cow", all=FALSE)
dim(metadata.pd2)
metadata.pd2 <- merge(metadata.pd2, SCC.pd2, by="cow", all=FALSE)
dim(metadata.pd2)

metadata.v2 <- rbind(metadata.pd0, metadata.pd1, metadata.pd2)

#Creating a new column in the metadata file that is the log transformation of SCC
metadata.v2$logSCC.milking1.beforesample <- log10(1000*(metadata.v2$SCCx1000.milking1.beforesample))
metadata.v2$logSCC.milking1.aftersample <- log10(1000*(metadata.v2$SCCx1000.milking1.aftersample))

#STEP 3: Loading in the feature table that's been rarefied and filtered
feature.table <- read_qza("feature-table-rarefied-at-507-samples-for-analysis.qza")
feature.table <- as.data.frame(feature.table$data)

#for the feature table, rename 5849-D1 to 5849-D2 (and sample 6233-D1 to 6233-D2,
#but this sample has already been filtered during the previous steps). 
library(dplyr)
feature.table<-rename(feature.table, "5849-D2"="5849-D1")

#STEP 4: Filter the metadata table down to the samples in the feature table
metadata.v3 <- metadata.v2[which(metadata.v2$sample.id %in% colnames(feature.table)),]
dim(metadata.v3) 
rownames(metadata.v3)=NULL

#STEP 5: Remove cows that don't have at least one sample from each period. 
cow.pd0 <- unique(metadata.v3[metadata.v3$period=="p_0",][["cow"]])
cow.pd1 <- unique(metadata.v3[metadata.v3$period=="p_1",][["cow"]])
cow.pd2 <- unique(metadata.v3[metadata.v3$period=="p_2",][["cow"]])
cows.to.keep <- cow.pd0[which(cow.pd0 %in% cow.pd1)]
cows.to.keep <- cows.to.keep[which(cows.to.keep %in% cow.pd2)]
metadata.v4 <- metadata.v3[which(metadata.v3$cow %in% cows.to.keep),]
dim(metadata.v4) #313 samples
length(unique(metadata.v4$cow)) #60 cows
rownames(metadata.v4)=NULL

#STEP 6: Read in taxonomy
taxonomy <- read_qza("taxonomy-PE-no-singletons-blanks-crmc-decontaminated-no-Bostaurus-mito-chloro-eukary.qza")

#change it to the way phyloseq wants it 
taxonomy2 <- as.data.frame(taxonomy$data) %>% column_to_rownames("Feature.ID") %>% separate(Taxon, c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"), sep = ";") %>% as.matrix()
#remove "Confidence" column
taxonomy2 <- taxonomy2[,-8]

#trim taxonomy to include only ASVs in the ASV table
taxonomy3 <- taxonomy2[c(rownames(taxonomy2)%in%rownames(feature.table)),]

#Replace all empty levels with NA for the sake of consistency
#This prevents inaccurate glomming later
taxonomy3[,grep("Species", colnames(taxonomy3))][grep("s__$", taxonomy3[,grep("Species", colnames(taxonomy3))])] <- NA
taxonomy3[,grep("Genus", colnames(taxonomy3))][grep("g__$", taxonomy3[,grep("Genus", colnames(taxonomy3))])] <- NA
taxonomy3[,grep("Family", colnames(taxonomy3))][grep("f__$", taxonomy3[,grep("Family", colnames(taxonomy3))])] <- NA
taxonomy3[,grep("Order", colnames(taxonomy3))][grep("o__$", taxonomy3[,grep("Order", colnames(taxonomy3))])] <- NA
taxonomy3[,grep("Class", colnames(taxonomy3))][grep("c__$", taxonomy3[,grep("Class", colnames(taxonomy3))])] <- NA
taxonomy3[,grep("Phylum", colnames(taxonomy3))][grep("p__$", taxonomy3[,grep("Phylum", colnames(taxonomy3))])] <- NA

#STEP 7: Read in rooted tree. 
tree <- read_qza("rooted-tree-rarefied507-samplesforanalyses.qza")

#STEP 8: build a phyloseq object
phylobj <- phyloseq(otu_table(feature.table, taxa_are_rows=T),tax_table(taxonomy3), phy_tree(tree$data), sample_data(metadata.v4 %>% as.data.frame() %>% column_to_rownames("sample.id")))

#keep only the taxa that have at least one count
phylobj.filtered <- filter_taxa(phylobj, function(x) sum(x)>0, TRUE)

#STEP 9: Gather taxa at Genus level
phylobj.genus <- tax_glom(phylobj.filtered, 'Genus', NArm=F)
dim(otu_table(phylobj.genus)) #There are 1308 taxa and 313 samples

#Create a list of unique taxa names that can be assigned to the features within the phyloseq object
taxa.names <- tax_table(phylobj.genus)[,6]

#replace empty names with family level identities
taxa.names2 <- tax_table(phylobj.genus)[,5] #make a list of family level designations
taxa.names[c(which(is.na(taxa.names)), grep("g__$", taxa.names))] <-taxa.names2[c(which(is.na(taxa.names)),grep("g__$", taxa.names))]

#replace empty names with order level identities
taxa.names3 <-tax_table(phylobj.genus)[,4] #make a list of order level designations
taxa.names[c(which(is.na(taxa.names2)), grep( "f__$", taxa.names))] <- taxa.names3[c(which(is.na(taxa.names2)), grep("f__$", taxa.names))]

#replace empty names with class level identities
taxa.names4 <- tax_table(phylobj.genus)[,3] #make a list of class level designations
taxa.names[c(which(is.na(taxa.names3)), grep("o__$", taxa.names))] <-taxa.names4[c(which(is.na(taxa.names3)),grep("o__$", taxa.names))]

#replace empty names with phylum level identities
taxa.names5 <- tax_table(phylobj.genus)[,2] #make a list of phylum level designations
taxa.names[c(which(is.na(taxa.names4)),grep("c__$", taxa.names))] <-taxa.names5[c(which(is.na(taxa.names4)),grep("c__$", taxa.names))]

#replace empty names with kingdom level identities
taxa.names6 <- tax_table(phylobj.genus)[,1] #make a list of phylum level designations
taxa.names[c(which(is.na(taxa.names5)),grep("p__$", taxa.names))] <-taxa.names6[c(which(is.na(taxa.names5)),grep("p__$", taxa.names))]

#remove leading white space
taxa.names <- gsub(" ","", taxa.names)

#change the taxa designations to names compatible with R
taxa.names <- make.names(taxa.names, unique = TRUE)

#Make all of the feature names within the phyloseq object correspond to the genus level designation
taxa_names(phylobj.genus) <- taxa.names

#STEP 10: Combine the data at the Family level
phylobj.family <- tax_glom(phylobj.filtered, 'Family', NArm=F)
dim(otu_table(phylobj.family)) #There are 520 family level taxa and 313 samples

#Create a list of unique taxa names that can be assigned to the features within the phyloseq object
taxa.names <- tax_table(phylobj.family)[,5]

#replace empty names with order level identities
taxa.names2 <- tax_table(phylobj.family)[,4] #make a list of order level designations
taxa.names[c(which(is.na(taxa.names)), grep("f__$", taxa.names))] <-taxa.names2[c(which(is.na(taxa.names)),grep("f__$", taxa.names))]

#replace empty names with class level identities
taxa.names3 <-tax_table(phylobj.family)[,3] #make a list of class level designations
taxa.names[c(which(is.na(taxa.names2)), grep("o__$", taxa.names))] <- taxa.names3[c(which(is.na(taxa.names2)), grep("o__$", taxa.names))]

#replace empty names with phylum level identities
taxa.names4 <- tax_table(phylobj.family)[,2] #make a list of phylum level designations
taxa.names[c(which(is.na(taxa.names3)), grep("c__$", taxa.names))] <-taxa.names4[c(which(is.na(taxa.names3)),grep("c__$", taxa.names))]

#replace empty names with kingdom level identities
taxa.names5 <- tax_table(phylobj.family)[,1] #make a list of kingdom level designations
taxa.names[c(which(is.na(taxa.names4)),grep("p__$", taxa.names))] <-taxa.names5[c(which(is.na(taxa.names4)),grep("p__$", taxa.names))]

#remove leading white space
taxa.names <- gsub(" ","", taxa.names)

#change the taxa designations to names compatible with R
taxa.names <- make.names(taxa.names, unique = TRUE)

#Make all of the feature names within the phyloseq object correspond to the family level designation
taxa_names(phylobj.family) <- taxa.names

#STEP 11: Combine all of the data at the Class level
phylobj.class <- tax_glom(phylobj.filtered, 'Class', NArm=F)
dim(otu_table(phylobj.class)) #There are 128 class level taxa and 313 samples

##Create a list of unique taxa names that can be assigned to the features within the phyloseq object
taxa.names <- tax_table(phylobj.class)[,3]

#replace empty names with phylum level identities
taxa.names2 <- tax_table(phylobj.class)[,2] #make a list of phylum level designations
taxa.names[c(which(is.na(taxa.names)), grep("c__$", taxa.names))] <-taxa.names2[c(which(is.na(taxa.names)),grep("c__$", taxa.names))]

#replace empty names with kingdom level identities
taxa.names3 <-tax_table(phylobj.class)[,1] #make a list of kingdom level designations
taxa.names[c(which(is.na(taxa.names2)), grep("p__$", taxa.names))] <- taxa.names3[c(which(is.na(taxa.names2)), grep("p__$", taxa.names))]

#remove leading white space
taxa.names <- gsub(" ","", taxa.names)

#change the taxa designations to names compatible with R
taxa.names <- make.names(taxa.names, unique = TRUE)

#Make all of the feature names within the phyloseq object correspond to the class level designation
taxa_names(phylobj.class) <- taxa.names

#STEP 12: Combine all of the data at the Phylum level
phylobj.phylum <- tax_glom(phylobj.filtered, 'Phylum', NArm=F)
dim(otu_table(phylobj.phylum)) #There are 51 phyla level taxa and 313 samples

##Create a list of unique taxa names that can be assigned to the features within the phyloseq object
taxa.names <- tax_table(phylobj.phylum)[,2]

#replace empty names with kingdom level identities
taxa.names2 <- tax_table(phylobj.phylum)[,1] #make a list of kingdom level designations
taxa.names[c(which(is.na(taxa.names)), grep("p__$", taxa.names))] <-taxa.names2[c(which(is.na(taxa.names)),grep("p__$", taxa.names))]

#remove leading white space
taxa.names <- gsub(" ","", taxa.names)

#change the taxa designations to names compatible with R
taxa.names <- make.names(taxa.names, unique = TRUE)

#Make all of the feature names within the phyloseq object correspond to the phylum level designation
taxa_names(phylobj.phylum) <- taxa.names

#STEP 13: Prepare the phyloseq genus level object and phyloseq family level object
#statistical analyses and visualizations/plots
melted.phylgenus <- psmelt(phylobj.genus)
melted.phylfamily <- psmelt(phylobj.family)
melted.phylclass <- psmelt(phylobj.class)
melted.phylphylum <- psmelt(phylobj.phylum)

#check the data types of each column, and if they aren't what I desire, then change to desired data type.
sapply(melted.phylgenus, class)
sapply(melted.phylfamily, class)
sapply(melted.phylclass, class)
sapply(melted.phylphylum, class)

#STEP 12: Perform linear mixed effects modeling to test if any 
#genera are predictors of log(SCC). 

#Selecting days 11/26, 01/04, and 03/15 for cows; and for cows that don't have
#all samples from these time points remaining after rarefying, then selecting
#earlier time point. In other words: after rarefying, I am choosing the latest 
#time point available for each cow in a period. 
ChooseLatestTimepoint <-function(cow) {
  df <- data.frame()
  for(i in cow) {
    metadata.subsetbycow<-metadata.v4[metadata.v4$cow==i,]
    metadata.subset.period0 <- metadata.subsetbycow[metadata.subsetbycow$period=="p_0",]
    metadata.subset.period1 <- metadata.subsetbycow[metadata.subsetbycow$period=="p_1",]
    metadata.subset.period2 <- metadata.subsetbycow[metadata.subsetbycow$period=="p_2",]
    if(length(metadata.subset.period0$period)>1) {
      metadata.subset.period0<-metadata.subset.period0[order(metadata.subset.period0$day, decreasing=TRUE),]
      metadata.subset.period0.sample <- metadata.subset.period0[match(unique(metadata.subset.period0$period), metadata.subset.period0$period),]
      df <- rbind(df, metadata.subset.period0.sample)
    } else {
      df <- rbind(df, metadata.subset.period0)
    }
    if(length(metadata.subset.period1$period)>1) {
      metadata.subset.period1<-metadata.subset.period1[order(metadata.subset.period1$day, decreasing=TRUE),]
      metadata.subset.period1.sample <- metadata.subset.period1[match(unique(metadata.subset.period1$period), metadata.subset.period1$period),]
      df <- rbind(df, metadata.subset.period1.sample)
    } else {
      df <- rbind(df, metadata.subset.period1)
    }
    if(length(metadata.subset.period2$period)>1) {
      metadata.subset.period2<-metadata.subset.period2[order(metadata.subset.period2$day, decreasing=TRUE),]
      metadata.subset.period2.sample <- metadata.subset.period2[match(unique(metadata.subset.period2$period), metadata.subset.period2$period),]
      df <- rbind(df, metadata.subset.period2.sample)
    } else {
      df <- rbind(df, metadata.subset.period2)
    }
  }
  return(df)
}

metadata.v5 <- ChooseLatestTimepoint(unique(metadata.v4$cow))
dim(metadata.v5) #180 samples
length(unique(metadata.v5$cow)) #60 cows

melted.phylgenus.v2 <- melted.phylgenus[which(melted.phylgenus$Sample %in% metadata.v5$sample.id),]

#Will only test for relationships with the genera that are
#present in 50% of samples and in 90% of cows. 

#looking at the prevalence of the genera
prevalence.per.genus <- melted.phylgenus.v2[melted.phylgenus.v2$Abundance>0, ]
prevalence.per.genus <- data.frame(tapply(prevalence.per.genus$Abundance, prevalence.per.genus$OTU, length))
prevalence.per.genus$OTUnames <- rownames(prevalence.per.genus)
colnames(prevalence.per.genus) <- c("prevalence", "OTU")

#Testing for just genera that are present in at least 50% of the samples
length(unique(melted.phylgenus.v2$Sample))
genera.50percent <-prevalence.per.genus[prevalence.per.genus$prevalence>=(180*0.5),]
genera.50percent <- unique(genera.50percent$OTU)

#And of those genera (present in 50% of samples), selecting the genera that are in 
#at least 90 % of cows (i.e. >= 90 % of cows have that genus for at least one time point)
nonzero.counts.per.cow.per.genus <-melted.phylgenus.v2[which(melted.phylgenus.v2$OTU %in%genera.50percent),]
nonzero.counts.per.cow.per.genus <- nonzero.counts.per.cow.per.genus[nonzero.counts.per.cow.per.genus$Abundance>0, ]
length(unique(melted.phylgenus.v2$cow))

Function.NonZeroCounts.PerCow.PerGenus<-function(genus){
  df <- data.frame()
  for(i in genus){
    OTU.subset <- nonzero.counts.per.cow.per.genus[nonzero.counts.per.cow.per.genus$OTU==i,]
    number.cows<-length(unique(OTU.subset$cow))
    percentage.cows <- (number.cows/60)*100
    OTU.cows <-data.frame(i, percentage.cows)  
    df <- rbind(df, OTU.cows)
  }
  print(df)
}

nonzero.counts.per.cow.per.genus.v2 <- Function.NonZeroCounts.PerCow.PerGenus(unique(nonzero.counts.per.cow.per.genus$OTU))
colnames(nonzero.counts.per.cow.per.genus.v2) <-c("OTU", "percent.of.cows")
genera.50PercentSamples.90PercentCows <- nonzero.counts.per.cow.per.genus.v2[nonzero.counts.per.cow.per.genus.v2$percent.of.cows>=90,][["OTU"]]
length(genera.50PercentSamples.90PercentCows) #40 genera

#Add the full taxonomy to the list of genera tested for relationships with SCC,
#then write to file 
full.taxonomy.genus <- data.frame(tax_table(phylobj.genus))
full.taxonomy.genus$Genus <- rownames(full.taxonomy.genus)
full.taxonomy.genus <- full.taxonomy.genus[which(full.taxonomy.genus$Genus %in% genera.50PercentSamples.90PercentCows), ]
full.taxonomy.genus$collapsed.taxonomy <- paste(full.taxonomy.genus$Kingdom, full.taxonomy.genus$Phylum, full.taxonomy.genus$Class, full.taxonomy.genus$Order, full.taxonomy.genus$Family, full.taxonomy.genus$Genus, sep="")
full.taxonomy.genus <- full.taxonomy.genus[,c(6,8)]
write.csv(full.taxonomy.genus, "Genera_tested_for_relationships_with_SCC.csv")

LME.Microbiota.SCC.genera <- function(genera){
  coefficient.pvalue <- data.frame(messages=NA,coefficient=NA, pvalue=NA, genus=NA)
  for(i in genera){
    OTU.abundance <- melted.phylgenus.v2[melted.phylgenus.v2$OTU==i,]
    lmer.model<-lmer(formula=logSCC.milking1.aftersample~Abundance+(1|cow), data=OTU.abundance)
    allfit.output<-allFit(lmer.model)
    bobyqa.output <- summary(allfit.output[["bobyqa"]])
    summary.subset <- data.frame(messages=paste(bobyqa.output[["optinfo"]][["conv"]][["lme4"]][["messages"]], collapse=" "), coefficient=bobyqa.output[["coefficients"]][2,1], pvalue=bobyqa.output[["coefficients"]][2,5], genus=i)
    coefficient.pvalue <- rbind(coefficient.pvalue, summary.subset)
  }
  print(coefficient.pvalue)
}


genus.SCC.LME.output <- LME.Microbiota.SCC.genera(genera.50PercentSamples.90PercentCows)#before multiple test correction, Staphylococcus and Aerococcus had positive relationships with SCC. 
genus.SCC.LME.output <- genus.SCC.LME.output[-1,]
genus.SCC.LME.output$fdr <- p.adjust(genus.SCC.LME.output$pvalue, method="fdr")

#Creating graphs for Staphylococcus vs. log(SCC) and Aerococcus vs. log(SCC)
staphylococcus <- melted.phylgenus.v2[melted.phylgenus.v2$OTU=="g__Staphylococcus",]
staphylococcus.logSCC.lmer <- allFit(lmer(logSCC.milking1.aftersample~Abundance+(1|cow), data=staphylococcus))
staphylococcus.logSCC.lmer.bobyqa <- staphylococcus.logSCC.lmer[["bobyqa"]]
ggplot(aes(x=Abundance, y=logSCC.milking1.aftersample), data=staphylococcus)+geom_point()+geom_line(aes(y=predict(staphylococcus.logSCC.lmer.bobyqa), group=cow))+xlab("Staphylococcus")+ylab("log(SCC)")+annotate("text", x=300, y= 6.1, label=bquote('1.8 x'~10^-3))
ggsave(filename="staphylococcus.vs.logSCC.jpeg", device=jpeg)

aerococcus <- melted.phylgenus.v2[melted.phylgenus.v2$OTU=="g__Aerococcus",]
aerococcus.logSCC.lmer <- allFit(lmer(logSCC.milking1.aftersample~Abundance+(1|cow), data=aerococcus))
aerococcus.logSCC.lmer.bobyqa <- aerococcus.logSCC.lmer[["bobyqa"]]
ggplot(aes(x=Abundance, y=logSCC.milking1.aftersample), data=aerococcus)+geom_point()+geom_line(aes(y=predict(aerococcus.logSCC.lmer.bobyqa), group=cow))+xlab("Aerococcus")+ylab("log(SCC)")+annotate("text", x=210, y= 5.5, label=bquote('2.3 x'~10^-3))
ggsave(filename="aerococcus.vs.logSCC.jpeg", device=jpeg)

#Running the model with microbial abundance effect varying by cow just to
#confirm the 'direction' of relationship between some genera and logSCC.
LME.Microbiota.SCC.genera.v2 <- function(genera){
    coefficient.pvalue <- data.frame(messages=NA,coefficient=NA, pvalue=NA, genus=NA)
    for(i in genera){
        OTU.abundance <- melted.phylgenus.v2[melted.phylgenus.v2$OTU==i,]
        lmer.model<-lmer(formula=logSCC.milking1.aftersample~Abundance+(Abundance|cow), data=OTU.abundance)
        allfit.output<-allFit(lmer.model)
        bobyqa.output <- summary(allfit.output[["bobyqa"]])
        summary.subset <- data.frame(messages=paste(bobyqa.output[["optinfo"]][["conv"]][["lme4"]][["messages"]], collapse=" "), coefficient=bobyqa.output[["coefficients"]][2,1], pvalue=bobyqa.output[["coefficients"]][2,5], genus=i)
        coefficient.pvalue <- rbind(coefficient.pvalue, summary.subset)
    }
    print(coefficient.pvalue)
}


genus.SCC.LME.output.v2 <- LME.Microbiota.SCC.genera.v2(genera.50PercentSamples.90PercentCows) 
genus.SCC.LME.output.v2 <- genus.SCC.LME.output.v2[-1,]
genus.SCC.LME.output.v2$fdr <- p.adjust(genus.SCC.LME.output.v2$pvalue, method="fdr")


#STEP 13: Perform linear mixed effects modeling to test if any 
#families are predictors of log(SCC).
melted.phylfamily.v2 <- melted.phylfamily[which(melted.phylfamily$Sample %in% metadata.v5$sample.id),]

#looking at the prevalence of the families
prevalence.per.family <- melted.phylfamily.v2[melted.phylfamily.v2$Abundance>0, ]
prevalence.per.family <- data.frame(tapply(prevalence.per.family$Abundance, prevalence.per.family$OTU, length))
prevalence.per.family$OTUnames <- rownames(prevalence.per.family)
colnames(prevalence.per.family) <- c("prevalence", "OTU")

#Testing for just families that are present in at least 50% of the samples
families.50percent <-prevalence.per.family[prevalence.per.family$prevalence>=(180*0.5),]
families.50percent <- unique(families.50percent$OTU)

#And of those families (present in 50% of samples), selecting the families that are in 
#at least 90 % of cows (i.e. >= 90 % of cows have that family for at least one time point)
nonzero.counts.per.cow.per.family <-melted.phylfamily.v2[which(melted.phylfamily.v2$OTU %in%families.50percent),]
nonzero.counts.per.cow.per.family <- nonzero.counts.per.cow.per.family[nonzero.counts.per.cow.per.family$Abundance>0, ]

Function.NonZeroCounts.PerCow.PerFamily<-function(family){
  df <- data.frame()
  for(i in family){
    OTU.subset <- nonzero.counts.per.cow.per.family[nonzero.counts.per.cow.per.family$OTU==i,]
    number.cows<-length(unique(OTU.subset$cow))
    percentage.cows <- (number.cows/60)*100
    OTU.cows <-data.frame(i, percentage.cows)  
    df <- rbind(df, OTU.cows)
  }
  print(df)
}

nonzero.counts.per.cow.per.family.v2 <- Function.NonZeroCounts.PerCow.PerFamily(unique(nonzero.counts.per.cow.per.family$OTU))
colnames(nonzero.counts.per.cow.per.family.v2) <-c("OTU", "percent.of.cows")
families.50PercentSamples.90PercentCows <- nonzero.counts.per.cow.per.family.v2[nonzero.counts.per.cow.per.family.v2$percent.of.cows>=90,][["OTU"]]
length(families.50PercentSamples.90PercentCows) #41 families

LME.Microbiota.SCC.families <- function(family){
    coefficient.pvalue <- data.frame(messages=NA,coefficient=NA, pvalue=NA, family=NA)
    for(i in family){
        OTU.abundance <- melted.phylfamily.v2[melted.phylfamily.v2$OTU==i,]
        lmer.model<-lmer(formula=logSCC.milking1.aftersample~Abundance+(1|cow), data=OTU.abundance)
        allfit.output<-allFit(lmer.model)
        bobyqa.output <- summary(allfit.output[["bobyqa"]])
        summary.subset <- data.frame(messages=paste(bobyqa.output[["optinfo"]][["conv"]][["lme4"]][["messages"]], collapse=" "), coefficient=bobyqa.output[["coefficients"]][2,1], pvalue=bobyqa.output[["coefficients"]][2,5], family=i)
        coefficient.pvalue <- rbind(coefficient.pvalue, summary.subset)
    }
    print(coefficient.pvalue)
}


families.SCC.LME.output <- LME.Microbiota.SCC.families(families.50PercentSamples.90PercentCows)#Before multiple test correction, only Aerococcaceae had a positive relationship with SCC
families.SCC.LME.output <- families.SCC.LME.output[-1,]
families.SCC.LME.output$fdr <- p.adjust(families.SCC.LME.output$pvalue, method="fdr")

LME.Microbiota.SCC.families.v2 <- function(family){
    coefficient.pvalue <- data.frame(messages=NA,coefficient=NA, pvalue=NA, family=NA)
    for(i in family){
        OTU.abundance <- melted.phylfamily.v2[melted.phylfamily.v2$OTU==i,]
        lmer.model<-lmer(formula=logSCC.milking1.aftersample~Abundance+(Abundance|cow), data=OTU.abundance)
        allfit.output<-allFit(lmer.model)
        bobyqa.output <- summary(allfit.output[["bobyqa"]])
        summary.subset <- data.frame(messages=paste(bobyqa.output[["optinfo"]][["conv"]][["lme4"]][["messages"]], collapse=" "), coefficient=bobyqa.output[["coefficients"]][2,1], pvalue=bobyqa.output[["coefficients"]][2,5], family=i)
        coefficient.pvalue <- rbind(coefficient.pvalue, summary.subset)
    }
    print(coefficient.pvalue)
}


families.SCC.LME.output.v2 <- LME.Microbiota.SCC.families.v2(families.50PercentSamples.90PercentCows)#many singular fit warnings and a model convergence failures. 
families.SCC.LME.output.v2 <- families.SCC.LME.output.v2[-1,]
families.SCC.LME.output.v2$fdr <- p.adjust(families.SCC.LME.output.v2$pvalue, method="fdr")


#STEP 14: Perform linear mixed effects modeling to test if any 
#classes are predictors of log(SCC).
melted.phylclass.v2 <- melted.phylclass[which(melted.phylclass$Sample %in% metadata.v5$sample.id),]

#looking at the prevalence of the classes
prevalence.per.class <- melted.phylclass.v2[melted.phylclass.v2$Abundance>0, ]
prevalence.per.class <- data.frame(tapply(prevalence.per.class$Abundance, prevalence.per.class$OTU, length))
prevalence.per.class$OTUnames <- rownames(prevalence.per.class)
colnames(prevalence.per.class) <- c("prevalence", "OTU")

#Testing for just classes that are present in at least 50% of the samples
classes.50percent <-prevalence.per.class[prevalence.per.class$prevalence>=(180*0.5),]
classes.50percent <- unique(classes.50percent$OTU)

#And of those classes (present in 50% of samples), selecting the classes that are in 
#at least 90 % of cows (i.e. >= 90 % of cows have that class for at least one time point)
nonzero.counts.per.cow.per.class <-melted.phylclass.v2[which(melted.phylclass.v2$OTU %in%classes.50percent),]
nonzero.counts.per.cow.per.class <- nonzero.counts.per.cow.per.class[nonzero.counts.per.cow.per.class$Abundance>0, ]

Function.NonZeroCounts.PerCow.PerClass<-function(class){
  df <- data.frame()
  for(i in class){
    OTU.subset <- nonzero.counts.per.cow.per.class[nonzero.counts.per.cow.per.class$OTU==i,]
    number.cows<-length(unique(OTU.subset$cow))
    percentage.cows <- (number.cows/60)*100
    OTU.cows <-data.frame(i, percentage.cows)  
    df <- rbind(df, OTU.cows)
  }
  print(df)
}

nonzero.counts.per.cow.per.class.v2 <- Function.NonZeroCounts.PerCow.PerClass(unique(nonzero.counts.per.cow.per.class$OTU))
colnames(nonzero.counts.per.cow.per.class.v2) <-c("OTU", "percent.of.cows")
classes.50PercentSamples.90PercentCows <- nonzero.counts.per.cow.per.class.v2[nonzero.counts.per.cow.per.class.v2$percent.of.cows>=90,][["OTU"]]

LME.Microbiota.SCC.classes <- function(class){
    coefficient.pvalue <- data.frame(messages=NA,coefficient=NA, pvalue=NA, class=NA)
    for(i in class){
        OTU.abundance <- melted.phylclass.v2[melted.phylclass.v2$OTU==i,]
        lmer.model<-lmer(formula=logSCC.milking1.aftersample~Abundance+(1|cow), data=OTU.abundance)
        allfit.output<-allFit(lmer.model)
        bobyqa.output <- summary(allfit.output[["bobyqa"]])
        summary.subset <- data.frame(messages=paste(bobyqa.output[["optinfo"]][["conv"]][["lme4"]][["messages"]], collapse=" "), coefficient=bobyqa.output[["coefficients"]][2,1], pvalue=bobyqa.output[["coefficients"]][2,5], class=i)
        coefficient.pvalue <- rbind(coefficient.pvalue, summary.subset)
    }
    print(coefficient.pvalue)
}


classes.SCC.LME.output <- LME.Microbiota.SCC.classes(classes.50PercentSamples.90PercentCows) #no significant relationships
classes.SCC.LME.output <- classes.SCC.LME.output[-1,]
classes.SCC.LME.output$fdr <- p.adjust(classes.SCC.LME.output$pvalue, method="fdr")

LME.Microbiota.SCC.classes.v2 <- function(class){
    coefficient.pvalue <- data.frame(messages=NA,coefficient=NA, pvalue=NA, class=NA)
    for(i in class){
        OTU.abundance <- melted.phylclass.v2[melted.phylclass.v2$OTU==i,]
        lmer.model<-lmer(formula=logSCC.milking1.aftersample~Abundance+(Abundance|cow), data=OTU.abundance)
        allfit.output<-allFit(lmer.model)
        bobyqa.output <- summary(allfit.output[["bobyqa"]])
        summary.subset <- data.frame(messages=paste(bobyqa.output[["optinfo"]][["conv"]][["lme4"]][["messages"]], collapse=" "), coefficient=bobyqa.output[["coefficients"]][2,1], pvalue=bobyqa.output[["coefficients"]][2,5], class=i)
        coefficient.pvalue <- rbind(coefficient.pvalue, summary.subset)
    }
    print(coefficient.pvalue)
}


classes.SCC.LME.output.v2 <- LME.Microbiota.SCC.classes.v2(classes.50PercentSamples.90PercentCows) #no significant relationships
classes.SCC.LME.output.v2 <- classes.SCC.LME.output.v2[-1,]
classes.SCC.LME.output.v2$fdr <- p.adjust(classes.SCC.LME.output.v2$pvalue, method="fdr")


#STEP 15: Perform linear mixed effects modeling to test if any 
#phyla are predictors of log(SCC).
melted.phylphyla.v2 <- melted.phylphylum[which(melted.phylphylum$Sample %in% metadata.v5$sample.id),]

#looking at the prevalence of the phyla
prevalence.per.phylum <- melted.phylphyla.v2[melted.phylphyla.v2$Abundance>0, ]
prevalence.per.phylum <- data.frame(tapply(prevalence.per.phylum$Abundance, prevalence.per.phylum$OTU, length))
prevalence.per.phylum$OTUnames <- rownames(prevalence.per.phylum)
colnames(prevalence.per.phylum) <- c("prevalence", "OTU")

#Testing for just phyla that are present in at least 50% of the samples
phyla.50percent <-prevalence.per.phylum[prevalence.per.phylum$prevalence>=(180*0.5),]
phyla.50percent <- unique(phyla.50percent$OTU)

#And of those phyla (present in 50% of samples), selecting the phyla that are in 
#at least 90 % of cows (i.e. >= 90 % of cows have that phylum for at least one time point)
nonzero.counts.per.cow.per.phylum <-melted.phylphyla.v2[which(melted.phylphyla.v2$OTU %in%phyla.50percent),]
nonzero.counts.per.cow.per.phylum <- nonzero.counts.per.cow.per.phylum[nonzero.counts.per.cow.per.phylum$Abundance>0, ]

Function.NonZeroCounts.PerCow.PerPhylum<-function(phylum){
  df <- data.frame()
  for(i in phylum){
    OTU.subset <- nonzero.counts.per.cow.per.phylum[nonzero.counts.per.cow.per.phylum$OTU==i,]
    number.cows<-length(unique(OTU.subset$cow))
    percentage.cows <- (number.cows/60)*100
    OTU.cows <-data.frame(i, percentage.cows)  
    df <- rbind(df, OTU.cows)
  }
  print(df)
}

nonzero.counts.per.cow.per.phylum.v2 <- Function.NonZeroCounts.PerCow.PerPhylum(unique(nonzero.counts.per.cow.per.phylum$OTU))
colnames(nonzero.counts.per.cow.per.phylum.v2) <-c("OTU", "percent.of.cows")
phyla.50PercentSamples.90PercentCows <- nonzero.counts.per.cow.per.phylum.v2[nonzero.counts.per.cow.per.phylum.v2$percent.of.cows>=90,][["OTU"]]

LME.Microbiota.SCC.phyla <- function(phyla){
    coefficient.pvalue <- data.frame(messages=NA,coefficient=NA, pvalue=NA, phyla=NA)
    for(i in phyla){
        OTU.abundance <- melted.phylphyla.v2[melted.phylphyla.v2$OTU==i,]
        lmer.model<-lmer(formula=logSCC.milking1.aftersample~Abundance+(1|cow), data=OTU.abundance)
        allfit.output<-allFit(lmer.model)
        bobyqa.output <- summary(allfit.output[["bobyqa"]])
        summary.subset <- data.frame(messages=paste(bobyqa.output[["optinfo"]][["conv"]][["lme4"]][["messages"]], collapse=" "), coefficient=bobyqa.output[["coefficients"]][2,1], pvalue=bobyqa.output[["coefficients"]][2,5], phyla=i)
        coefficient.pvalue <- rbind(coefficient.pvalue, summary.subset)
    }
    print(coefficient.pvalue)
}


phyla.SCC.LME.output <- LME.Microbiota.SCC.phyla(phyla.50PercentSamples.90PercentCows)#no correlations before multiple test correction

LME.Microbiota.SCC.phyla.v2 <- function(phyla){
    coefficient.pvalue <- data.frame(messages=NA,coefficient=NA, pvalue=NA, phyla=NA)
    for(i in phyla){
        OTU.abundance <- melted.phylphyla.v2[melted.phylphyla.v2$OTU==i,]
        lmer.model<-lmer(formula=logSCC.milking1.aftersample~Abundance+(Abundance|cow), data=OTU.abundance)
        allfit.output<-allFit(lmer.model)
        bobyqa.output <- summary(allfit.output[["bobyqa"]])
        summary.subset <- data.frame(messages=paste(bobyqa.output[["optinfo"]][["conv"]][["lme4"]][["messages"]], collapse=" "), coefficient=bobyqa.output[["coefficients"]][2,1], pvalue=bobyqa.output[["coefficients"]][2,5], phyla=i)
        coefficient.pvalue <- rbind(coefficient.pvalue, summary.subset)
    }
    print(coefficient.pvalue)
}

phyla.SCC.LME.output.v2 <- LME.Microbiota.SCC.phyla.v2(phyla.50PercentSamples.90PercentCows)#numerous singularity warnings. No relationships. 


#STEP 16: Calculate microbial alpha diversity. 

#phyloseq object down to the 180 samples from the 60 cows that have one sample per period. 
phylobj.filtered.v2 <- prune_samples(samples=metadata.v5$sample.id, x=phylobj.filtered)
length(sample_names(phylobj.filtered)) #313 samples
length(sample_names(phylobj.filtered.v2)) #180 samples

#remove taxa that have no counts in the filtered down dataset. 
phylobj.filtered.v3 <- filter_taxa(phylobj.filtered.v2, function(x) sum(x)>0, TRUE)
length(taxa_names(phylobj.filtered.v2)) #10572 taxa
length(taxa_names(phylobj.filtered.v3)) #6951 taxa

#calculate alpha diversity -- shannon and simpson indices
alpha.diversity.microbiota<-estimate_richness(physeq=phylobj.filtered.v3, measures=c("InvSimpson", "Shannon")) #from the phyloseq package

#combine these diversity measures to the metadata.v4 object
colnames(alpha.diversity.microbiota)<- c("microbiota.diversity.shannon", "microbiota.diversity.inversesimpson")
alpha.diversity.microbiota$sample.id <- rownames(alpha.diversity.microbiota)
metadata.v6 <- merge(metadata.v5, alpha.diversity.microbiota, by="sample.id", all=FALSE)
dim(metadata.v6) #180 samples, 21 variables. 

#STEP 17: Run LMM to determine if microbial alpha diversity (shannon, simpson,
#or phylogenetic diversity) has significant effect on log(SCC)
hist(metadata.v6$logSCC.milking1.aftersample)
hist(metadata.v6$microbiota.diversity.inversesimpson) #right skewed
hist(metadata.v6$microbiota.diversity.shannon) #a little left skewed
set.seed(1) 
summary(lmer(logSCC.milking1.aftersample~microbiota.diversity.shannon+(1|cow), data=metadata.v6)) #not significant
summary(lmer(logSCC.milking1.aftersample~microbiota.diversity.shannon+(microbiota.diversity.shannon|cow), data=metadata.v6)) #not significant
metadata.v6$cow <- as.factor(metadata.v6$cow)
rmcorr(participant = cow, measure1=logSCC.milking1.aftersample, measure2=microbiota.diversity.shannon, dataset=metadata.v6) #not significant

#Other functions to calculate alpha diversity from phyloseq object
#include alpha() and diversity() from the microbiome package

#The (Faith's) phylogenetic diversity of the microbiota data
#using the Picante package
#First, I need to extract the list of ASVs from this particular set of samples
#and use this to create a rooted phylogenetic tree with qiime2

ASV.table <- as.data.frame(otu_table(phylobj.filtered.v3))
phylogenetic.diversity.microbiota <- pd(samp=t(ASV.table), tree=phy_tree(phylobj.filtered.v3), include.root=TRUE)
phylogenetic.diversity.microbiota <- as.data.frame(phylogenetic.diversity.microbiota)
phylogenetic.diversity.microbiota$sample.id <- rownames(phylogenetic.diversity.microbiota)
phylogenetic.diversity.microbiota <- phylogenetic.diversity.microbiota[,-2]
colnames(phylogenetic.diversity.microbiota) <- c("microbiota.diversity.PD", "sample.id")
metadata.v7 <- merge(metadata.v6, phylogenetic.diversity.microbiota, by="sample.id", all=FALSE)

#LMM with microbial PD diversity as predictor of log(SCC)
hist(metadata.v7$microbiota.diversity.PD)
set.seed(1)
summary(lmer(logSCC.milking1.aftersample~microbiota.diversity.PD+(1|cow), data=metadata.v7)) #not significant
summary(lmer(logSCC.milking1.aftersample~microbiota.diversity.PD+(microbiota.diversity.PD|cow), data=metadata.v7)) #model failed to converge, not significant
metadata.v7$cow <- as.factor(metadata.v7$cow)
rmcorr(participant = cow, measure1=logSCC.milking1.aftersample, measure2=microbiota.diversity.PD, dataset=metadata.v7) #not significant

#STEP 18: Use GLMM to determine if log(SCC) is a predictor of microbial genera
#first, need to select the sample from each cow that is the earliest in the period. 
ChooseEarliestTimepoint <-function(cow) {
    df <- data.frame()
    for(i in cow) {
        metadata.subsetbycow<-metadata.v4[metadata.v4$cow==i,]
        metadata.subset.period0 <- metadata.subsetbycow[metadata.subsetbycow$period=="p_0",]
        metadata.subset.period1 <- metadata.subsetbycow[metadata.subsetbycow$period=="p_1",]
        metadata.subset.period2 <- metadata.subsetbycow[metadata.subsetbycow$period=="p_2",]
        if(length(metadata.subset.period0$period)>1) {
            metadata.subset.period0<-metadata.subset.period0[order(metadata.subset.period0$day, decreasing=FALSE),]
            metadata.subset.period0.sample <- metadata.subset.period0[match(unique(metadata.subset.period0$period), metadata.subset.period0$period),]
            df <- rbind(df, metadata.subset.period0.sample)
        } else {
            df <- rbind(df, metadata.subset.period0)
        }
        if(length(metadata.subset.period1$period)>1) {
            metadata.subset.period1<-metadata.subset.period1[order(metadata.subset.period1$day, decreasing=FALSE),]
            metadata.subset.period1.sample <- metadata.subset.period1[match(unique(metadata.subset.period1$period), metadata.subset.period1$period),]
            df <- rbind(df, metadata.subset.period1.sample)
        } else {
            df <- rbind(df, metadata.subset.period1)
        }
        if(length(metadata.subset.period2$period)>1) {
            metadata.subset.period2<-metadata.subset.period2[order(metadata.subset.period2$day, decreasing=FALSE),]
            metadata.subset.period2.sample <- metadata.subset.period2[match(unique(metadata.subset.period2$period), metadata.subset.period2$period),]
            df <- rbind(df, metadata.subset.period2.sample)
        } else {
            df <- rbind(df, metadata.subset.period2)
        }
    }
    return(df)
}

metadata.v8 <- ChooseEarliestTimepoint(unique(metadata.v4$cow))
dim(metadata.v8) #180 samples
length(unique(metadata.v8$cow)) #60 cows

#filter genus count table accordingly.
melted.phylgenus.v3 <- melted.phylgenus[which(melted.phylgenus$Sample %in% metadata.v8$sample.id),]

#select genera present in 50% of samples
prevalence.per.genus <- melted.phylgenus.v3[melted.phylgenus.v3$Abundance>0, ]
prevalence.per.genus <- data.frame(tapply(prevalence.per.genus$Abundance, prevalence.per.genus$OTU, length))
prevalence.per.genus$OTUnames <- rownames(prevalence.per.genus)
colnames(prevalence.per.genus) <- c("prevalence", "OTU")

#Testing for just genera that are present in at least 50% of the samples
length(unique(melted.phylgenus.v3$Sample))
genera.50percent <-prevalence.per.genus[prevalence.per.genus$prevalence>=(180*0.5),]
genera.50percent <- unique(genera.50percent$OTU)

#And of those genera (present in 50% of samples), selecting the genera that are in 
#at least 90 % of cows (i.e. >= 90 % of cows have that genus for at least one time point)
nonzero.counts.per.cow.per.genus <-melted.phylgenus.v3[which(melted.phylgenus.v3$OTU %in%genera.50percent),]
nonzero.counts.per.cow.per.genus <- nonzero.counts.per.cow.per.genus[nonzero.counts.per.cow.per.genus$Abundance>0, ]
length(unique(melted.phylgenus.v3$cow))

Function.NonZeroCounts.PerCow.PerGenus<-function(genus){
    df <- data.frame()
    for(i in genus){
        OTU.subset <- nonzero.counts.per.cow.per.genus[nonzero.counts.per.cow.per.genus$OTU==i,]
        number.cows<-length(unique(OTU.subset$cow))
        percentage.cows <- (number.cows/60)*100
        OTU.cows <-data.frame(i, percentage.cows)  
        df <- rbind(df, OTU.cows)
    }
    print(df)
}

nonzero.counts.per.cow.per.genus.v2 <- Function.NonZeroCounts.PerCow.PerGenus(unique(nonzero.counts.per.cow.per.genus$OTU))
colnames(nonzero.counts.per.cow.per.genus.v2) <-c("OTU", "percent.of.cows")
genera.50PercentSamples.90PercentCows <- nonzero.counts.per.cow.per.genus.v2[nonzero.counts.per.cow.per.genus.v2$percent.of.cows>=90,][["OTU"]]
length(genera.50PercentSamples.90PercentCows) #27 genera

#run GLMM
Function.SCCpredictor.GenusResponse.GLMM <-function(microbe) {
  coefficient.pvalue <- data.frame(messages=NA, coefficient=NA, pvalue=NA, taxon=NA)
    for(i in microbe){
      microbe.subset <- melted.phylgenus.v3[melted.phylgenus.v3$OTU==i,]
      glmer.model <- try(glmer(Abundance~logSCC.milking1.beforesample+(1|cow), family=poisson, data=microbe.subset), silent=TRUE)
      model.summary <- try(summary(glmer.model), silent=TRUE)
      glmer.df <- try(data.frame(messages=paste(model.summary[["optinfo"]][["conv"]][["lme4"]][["messages"]], collapse=" "), coefficient=model.summary[["coefficients"]][2,1], pvalue=model.summary[["coefficients"]][2,4], taxon=i), silent=FALSE)
      coefficient.pvalue <- try(rbind(coefficient.pvalue, glmer.df), silent=TRUE)
    }
  print(coefficient.pvalue)
}


results.GLMM.SCC.genera<- Function.SCCpredictor.GenusResponse.GLMM(genera.50PercentSamples.90PercentCows)
results.GLMM.SCC.genera <- results.GLMM.SCC.genera[-1,] #removed first line because it was empty
results.GLMM.SCC.genera$pvalue <- as.numeric(results.GLMM.SCC.genera$pvalue)
results.GLMM.SCC.genera$fdr <- p.adjust(p=results.GLMM.SCC.genera$pvalue, method="fdr")
results.GLMM.SCC.genera$bonferroni <- p.adjust(p=results.GLMM.SCC.genera$pvalue, method="bonferroni")
write.csv(results.GLMM.SCC.genera, file="logSCC.predictor.of.Microbiota.genera.glmm.poisson.csv")

lactobacillus <- melted.phylgenus.v3[melted.phylgenus.v3$OTU=="g__Lactobacillus",]
lactobacillus.logSCC.glmer <- glmer(Abundance~logSCC.milking1.beforesample+(1|cow), data=lactobacillus, family=poisson)
ggplot(aes(x=logSCC.milking1.beforesample, y=Abundance), data=lactobacillus)+geom_point()+geom_line(aes(y=predict(lactobacillus.logSCC.glmer), group=cow))+xlab("log(SCC)")+ylab("Lactobacillus")+annotate("text", x=6.0, y=300,label="-1.3")
ggsave(filename="logSCC.vs.Lactobacillus.jpeg", device="jpeg")

corynebacterium <- melted.phylgenus.v3[melted.phylgenus.v3$OTU=="g__Corynebacterium",]
corynebacterium.logSCC.glmer <- glmer(Abundance~logSCC.milking1.beforesample+(1|cow), data=corynebacterium, family=poisson)
ggplot(aes(x=logSCC.milking1.beforesample, y=Abundance), data=corynebacterium)+geom_point()+geom_line(aes(y=predict(corynebacterium.logSCC.glmer), group=cow))+xlab("log(SCC)")+ylab("Corynebacterium")+annotate("text", x=6.0, y=280, label="1.0")
ggsave("logSCC.vs.Corynebacterium.jpeg", device="jpeg")

#running GLMM with slope varying by cow too, to see if similar results are obtained.  
Function.SCCpredictor.GenusResponse.GLMM.v2 <-function(microbe) {
  coefficient.pvalue <- data.frame(messages=NA, coefficient=NA, pvalue=NA, taxon=NA)
  for(i in microbe){
    microbe.subset <- melted.phylgenus.v3[melted.phylgenus.v3$OTU==i,]
    glmer.model <- try(glmer(Abundance~logSCC.milking1.beforesample+(logSCC.milking1.beforesample|cow), family=poisson, data=microbe.subset), silent=TRUE)
    model.summary <- try(summary(glmer.model), silent=TRUE)
    glmer.df <- try(data.frame(messages=paste(model.summary[["optinfo"]][["conv"]][["lme4"]][["messages"]], collapse=" "), coefficient=model.summary[["coefficients"]][2,1], pvalue=model.summary[["coefficients"]][2,4], taxon=i), silent=FALSE)
    coefficient.pvalue <- try(rbind(coefficient.pvalue, glmer.df), silent=TRUE)
  }
  print(coefficient.pvalue)
}


results.GLMM.SCC.genera.v2<- Function.SCCpredictor.GenusResponse.GLMM.v2(genera.50PercentSamples.90PercentCows)


#STEP 19: Use LMM to determine if log(SCC) is a predictor of microbial alpha diversity.
#phyloseq object down to the 180 samples from the 60 cows that have one sample per period. 
phylobj.filtered.v4 <- prune_samples(samples=metadata.v8$sample.id, x=phylobj.filtered)
length(sample_names(phylobj.filtered)) #313 samples
length(sample_names(phylobj.filtered.v4)) #180 samples

#remove taxa that have no counts in the filtered down dataset. 
phylobj.filtered.v5 <- filter_taxa(phylobj.filtered.v4, function(x) sum(x)>0, TRUE)
length(taxa_names(phylobj.filtered.v4)) #10572 taxa
length(taxa_names(phylobj.filtered.v5)) #7695 taxa

#calculate alpha diversity -- shannon and inverse simpson indices
alpha.diversity.microbiota <-estimate_richness(physeq=phylobj.filtered.v5, measures=c("InvSimpson", "Shannon")) #from the phyloseq package

#combine these diversity measures to the metadata.v8 object
colnames(alpha.diversity.microbiota)<- c("microbiota.diversity.shannon", "microbiota.diversity.inversesimpson")
alpha.diversity.microbiota$sample.id <- rownames(alpha.diversity.microbiota)
metadata.v9 <- merge(metadata.v8, alpha.diversity.microbiota, by="sample.id", all=FALSE)
dim(metadata.v9) #180 samples, 21 variables. 


#The (Faith's) phylogenetic diversity of the microbiota data
#using the Picante package
#First, I need to extract the list of ASVs from this particular set of samples
#and use this to create a rooted phylogenetic tree with qiime2

ASV.table <- as.data.frame(otu_table(phylobj.filtered.v5))
phylogenetic.diversity.microbiota <- pd(samp=t(ASV.table), tree=phy_tree(phylobj.filtered.v5), include.root=TRUE)
phylogenetic.diversity.microbiota <- as.data.frame(phylogenetic.diversity.microbiota)
phylogenetic.diversity.microbiota$sample.id <- rownames(phylogenetic.diversity.microbiota)
phylogenetic.diversity.microbiota <- phylogenetic.diversity.microbiota[,-2]
colnames(phylogenetic.diversity.microbiota) <- c("microbiota.diversity.PD", "sample.id")
dim(metadata.v9) #180 rows, 21 columns
metadata.v9 <- merge(metadata.v9, phylogenetic.diversity.microbiota, by="sample.id", all=FALSE)
dim(metadata.v9) #180 rows, 22 columns

hist(metadata.v9$logSCC.milking1.beforesample) #fairly normal distribution
hist(metadata.v9$microbiota.diversity.shannon) #fairly normal distribution
hist(metadata.v9$microbiota.diversity.inversesimpson) #right skewed.
hist(metadata.v9$microbiota.diversity.PD) #a little right skewed

set.seed(1)
summary(lmer(microbiota.diversity.shannon~logSCC.milking1.beforesample+(1|cow), data=metadata.v9)) #singular fit warning
logSCC.microbiotaShannon<-lmer(microbiota.diversity.shannon~logSCC.milking1.beforesample+(1|cow), data=metadata.v9)
allFit(logSCC.microbiotaShannon) #continue to have singular fit warning
summary(lmer(microbiota.diversity.shannon~logSCC.milking1.beforesample+(logSCC.milking1.beforesample|cow), data=metadata.v9)) #singular fit warning
logSCC.microbiotaShannon.slope<-lmer(microbiota.diversity.shannon~logSCC.milking1.beforesample+(logSCC.milking1.beforesample|cow), data=metadata.v9)
allFit(logSCC.microbiotaShannon.slope) #continue to have singular fit warning
metadata.v9$cow <- as.factor(metadata.v9$cow)
rmcorr(participant=cow, measure1=logSCC.milking1.beforesample, measure2=microbiota.diversity.shannon, dataset=metadata.v9)#not significant

summary(lmer(microbiota.diversity.PD~logSCC.milking1.beforesample+(1|cow), data=metadata.v9)) #not significant
summary(lmer(microbiota.diversity.PD~logSCC.milking1.beforesample+(microbiota.diversity.PD|cow), data=metadata.v9)) #model failed to converge, not significant
rmcorr(participant=cow, measure1=logSCC.milking1.beforesample, measure2=microbiota.diversity.PD, dataset=metadata.v9)#not significant




