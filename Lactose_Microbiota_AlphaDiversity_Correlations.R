#This script describes the steps taken to analyze the relationship between lactose
#and particular microbiota and microbial alpha diversity

#Load needed libraries
library(devtools)
library(tidyverse)
#devtools::install_github("jbisanz/qiime2R")
library(qiime2R)
library(phyloseq)
#install_github("microbiome/microbiome")
library(microbiome)
library(dplyr)
library(plyr)

#Load other needed libraries
library(ggplot2)
library(ggpubr)
library(plotly)
library(psych)
library(cowplot)
library(vegan)
library(picante)
library(lme4)
library(lmerTest)
library(rmcorr)

#STEP 1:Read in metadata file, 
metadata  <- read.csv("../../../DairyGrandChallenge.Study1.2017-2018/ResearchMaterials.DairyGrandChallenge/Cleaned_For_MetadataFile/metadata.csv")

#keeping only the data that is needed. 
metadata <- metadata[,c(1:15)]
metadata <- metadata[-1,]

#STEP 2: Read in the table with number of reads per sample, before rarefaction. 
reads.per.sample <- read.csv("sample-frequency.csv", head=FALSE)
colnames(reads.per.sample) <- c("sample.id", "reads")
metadata.v2 <- merge(metadata, reads.per.sample, by="sample.id", all=FALSE)

#STEP 3: Read in the lactose dataset and clean up to remove empty columns and rows
macronutrient <- readxl::read_xlsx("PMT DVH_Dairy Grand Challenge coded FOSS composition data_101818.xlsx")

#Looking at the first few lines of the file
head(macronutrient, 10) #The intended column headers are on the eighth row of the excel file

#In order to set the correct the column headers set, need to read in the dataset again but without the first row being used as the column headers
macronutrient <- readxl::read_xlsx("PMT DVH_Dairy Grand Challenge coded FOSS composition data_101818.xlsx", col_names=FALSE)

#Set the column headers to the values on the eighth row
colnames(macronutrient)<-macronutrient[8,]
#Remove the first eight rows
macronutrient <- macronutrient[-c(1:8),]

dim(macronutrient) #there are 76 observations and 46 variables in the table
colnames(macronutrient)
head(macronutrient) #The sample IDs are listed under the "Day" columns
#I can tell from the column names that the observations for days 2,3,4,5,6 are listed to the right of day 1 observations instead of below
#I can also see that the first two columns appear to be empty
sum(is.na(macronutrient[,1])) #all entries in the first column are empty/NA
sum(is.na(macronutrient[,2])) #all entries in the second column are empty/NA

#Remove the first two columns
macronutrient <- macronutrient[,-c(1:2)]
colnames(macronutrient) #no longer see those empty column headers

#Remove the "Period 1", "Period 2", "Cow ID", "Treatment 0", "treatment 1", "Treatment 2", "Stolen from", and "Location" columns because these variables will be retrieved from the Kalscheur lab data or at not needed in the compiled dataset
macronutrient <- macronutrient[,-c(1:4, 17, 30, 43, 44)]
head(macronutrient)
colnames(macronutrient)

#Need to subset the by day (i.e. take subset chunks of columns), then add back together under the same column headers so that there are only 
macronutrient.Day1 <- macronutrient[,1:6]
macronutrient.Day2 <- macronutrient[,7:12]
macronutrient.Day3 <- macronutrient[,13:18]
macronutrient.Day4 <- macronutrient[,19:24]
macronutrient.Day5 <- macronutrient[,25:30]
macronutrient.Day6 <- macronutrient[,31:36]

#Change the "Day" column headers in the data subsets to "sample_id"
macronutrient.Day1 <- rename(macronutrient.Day1, sample.id=`Day 1`)
macronutrient.Day2 <- rename(macronutrient.Day2, sample.id=`Day 2`)
macronutrient.Day3 <- rename(macronutrient.Day3, sample.id=`Day 3`)
macronutrient.Day4 <- rename(macronutrient.Day4, sample.id=`Day 4`)
macronutrient.Day5 <- rename(macronutrient.Day5, sample.id=`Day  5`)
macronutrient.Day6 <- rename(macronutrient.Day6, sample.id=`Day 6`)

#All of the column headers should be the same now for the subsets fo the macronutrient data
#Checking to see that they are the same..
colnames(macronutrient.Day1)==colnames(macronutrient.Day2) #All column headers are the same
colnames(macronutrient.Day1)==colnames(macronutrient.Day3) #All column headers are the same
colnames(macronutrient.Day1)==colnames(macronutrient.Day4) #All column headers are the same
colnames(macronutrient.Day1)==colnames(macronutrient.Day5) #All column headers are the same
colnames(macronutrient.Day1)==colnames(macronutrient.Day6) #All column headers are the same

#Now add the subsets together by rows under the same set of 6 column headers
macronutrient <- rbind(macronutrient.Day1, macronutrient.Day2, macronutrient.Day3, macronutrient.Day4, macronutrient.Day5, macronutrient.Day6)
dim(macronutrient) #456 observations and 6 variables (this is what's expected with 76 observations in each subset)

#Check for NA values and duplicate values in the sample_id column
sum(is.na(macronutrient$sample.id)) #the output is zero, so no NA values among the sample_id column
macronutrient.distinct.col1 <- unique(macronutrient$sample.id)
length(macronutrient.distinct.col1) #456 unique observations in sample_id so there are no duplicate sample IDs

#Remove the destination information ("_Wyn") from the sample IDs so that they can be matched and merged to other datasets
macronutrient$sample.id <-gsub(pattern="_Wyn", x=macronutrient$sample.id, replacement="") 
head(macronutrient)

#Add a column specifying the units of measure for the macronutrient data
macronutrient$nut_Units <-"gram/100 gram of milk "

#STEP 4: Prepare lactose and metadata sets for merging. 
#Since the dates are in the sample IDs for macronutrients dataset, 
#I will change them to D1, D2, D3, D4, D5, and D6, so that they match up with 
#the metadata sampe.id's 
macronutrient$sample.id <- gsub(replacement="-D1", pattern="_171125", macronutrient$sample.id)
macronutrient$sample.id <- gsub(replacement="-D2", pattern="_171126", macronutrient$sample.id)
macronutrient$sample.id <- gsub(replacement="-D3", pattern="_180104", macronutrient$sample.id)
macronutrient$sample.id <- gsub(replacement="-D4", pattern="_180105", macronutrient$sample.id)
macronutrient$sample.id <- gsub(replacement="-D5", pattern="_180314", macronutrient$sample.id)
macronutrient$sample.id <- gsub(replacement="-D6", pattern="_180315", macronutrient$sample.id)

#Create a new "sample" column in the lactose and metadata sets that doesn't
#have "-B" designation on the microbiota sample names, so that they can be merged.

macronutrient$sample <- macronutrient$sample.id
metadata.v2$sample <- metadata.v2$sample.id
metadata.v2$sample <- gsub(pattern="-B", replacement="", metadata.v2$sample)

#Merge the macronutrient data with the metadata
metadata.v3 <- merge(metadata.v2, macronutrient, by="sample", all=FALSE)
dim(metadata.v3) #473 samples

#get rid of the "sample", "sample.id.x" and "sample.id.y" column
metadata.v3$sample.id <- metadata.v3$sample.id.x
metadata.v3 <- metadata.v3[,-c(1:2,18)]

#STEP 5: Read in feature table
feature.table <- read_qza("feature-table-rarefied-at-507-samples-for-analysis.qza")
feature.table <- as.data.frame(feature.table$data)

#STEP 6: Filter the metadata samples down to what is the feature table
metadata.v4 <- metadata.v3[which(metadata.v3$sample.id %in% colnames(feature.table)),]
dim(metadata.v4) #334 samples

#STEP 7: Rename sample 5849-D1 to 5849-D2 because cow 5849 actually missed
#the milk collection on November 25, 2017, but milk was collected from this cow
#on November 26, 2017. The same was true for cow 6233..
metadata.v4$sample.id <- gsub(pattern="5849-D1", replacement="5849-D2", x=metadata.v4$sample.id)
metadata.v4$sample.id <- gsub(pattern="6233-D1", replacement="6233-D2", x=metadata.v4$sample.id)
metadata.v4[metadata.v4$sample.id=="5849-D2","day"] <- 2
metadata.v4[metadata.v4$sample.id=="6233-D2","day"] <- 2

feature.table.v2 <- feature.table
library(dplyr)
feature.table.v2 <- rename(feature.table.v2, "5849-D2"="5849-D1")

#STEP 8: Filter the metadata samples to (a) remove lactose outliers (cows #6229 and 5651), 
#(b) remove cows without a sample from a period, and (c), if a cow has more than
#one sample in a period then, select the sample with higher reads (before rarefying) 
#for a cow in a period. 

#(a) remove lactose outlier cows
metadata.v5 <- metadata.v4[metadata.v4$cow!="6229",]
metadata.v5 <- metadata.v5[metadata.v5$cow!="5651",]
dim(metadata.v5) #324 samples

#(b) remove cows that don't have a sample from each period. 
cows.pd0 <- unique(metadata.v5[metadata.v5$period=="p_0",][["cow"]])
cows.pd1 <- unique(metadata.v5[metadata.v5$period=="p_1",][["cow"]])
cows.pd2 <- unique(metadata.v5[metadata.v5$period=="p_2",][["cow"]])
cows.to.keep <- cows.pd0[which(cows.pd0 %in% cows.pd1)]
cows.to.keep.v2 <- cows.to.keep[which(cows.to.keep %in% cows.pd2)]
length(cows.to.keep.v2) #58 cows
metadata.v6 <- metadata.v5[which(metadata.v5$cow %in% cows.to.keep.v2),]
dim(metadata.v6) #303 samples
length(unique(metadata.v6$cow)) #58 cows

#(c)for cows that have more than one sample in a period, select the sample
#that had the higher number of reads before rarefying.

ChooseSampleFromPeriod <-function(cow) {
    NoDuplicates <- data.frame()
    for(i in cow) {
        select.cow <- metadata.v6[metadata.v6$cow==i,]
        period0 <- select.cow[select.cow$period=="p_0",]
        if(length(period0$sample.id)>1) {
            period0<-period0[order(period0$reads, decreasing=TRUE),]
            Higher.Read.Sample <- period0[match(unique(period0$period), period0$period),]
            NoDuplicates <- rbind(NoDuplicates, Higher.Read.Sample)
        } else {
            NoDuplicates <- rbind(NoDuplicates, period0)
        }
        period1 <- select.cow[select.cow$period=="p_1",]
        if(length(period1$sample.id)>1) {
            period1<-period1[order(period1$reads, decreasing=TRUE),]
            Higher.Read.Sample <- period1[match(unique(period1$period), period1$period),]
            NoDuplicates <- rbind(NoDuplicates, Higher.Read.Sample)
        } else {
            NoDuplicates <- rbind(NoDuplicates, period1)
        }
        period2 <- select.cow[select.cow$period=="p_2",]
        if(length(period2$sample.id)>1) {
            period2<-period2[order(period2$reads, decreasing=TRUE),]
            Higher.Read.Sample <- period2[match(unique(period2$period), period2$period),]
            NoDuplicates <- rbind(NoDuplicates, Higher.Read.Sample)
        } else {
            NoDuplicates <- rbind(NoDuplicates, period2)
        }
    }
    return(NoDuplicates)
}

#Run the function "ChooseSampleFromPeriod"
cows <-unique(metadata.v6$cow)
metadata.v7 <- ChooseSampleFromPeriod(cows) 
length(metadata.v7$sample.id) #174 samples
length(unique(metadata.v7$cow)) #58 cows

#STEP 9: Read in taxonomy
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

#STEP 10: Read in rooted tree. 
tree <- read_qza("rooted-tree-rarefied507-samplesforanalyses.qza")

#STEP 11: build a phyloseq object
rownames(metadata.v7)=NULL
phylobj <- phyloseq(otu_table(feature.table.v2, taxa_are_rows=T),tax_table(taxonomy3), phy_tree(tree$data), sample_data(metadata.v7 %>% as.data.frame() %>% column_to_rownames("sample.id")))

#keep only the taxa that have at least one count
phylobj.filtered <- filter_taxa(phylobj, function(x) sum(x)>0, TRUE)

#STEP 12: Calculate alpha diversity for microbiota data..

#The inverse simpson and shannon diversity of the microbiota
alpha.diversity.microbiota<-estimate_richness(physeq=phylobj.filtered, measures=c("InvSimpson", "Shannon")) #from the phyloseq package

#combine these diversity measures to the metadata object
colnames(alpha.diversity.microbiota)<- c("microbiota.diversity.shannon", "microbiota.diversity.inversesimpson")
alpha.diversity.microbiota$sample.id <- rownames(alpha.diversity.microbiota)
metadata.v8 <- merge(metadata.v7, alpha.diversity.microbiota, by="sample.id", all=FALSE)

#Other functions to calculate alpha diversity from phyloseq object
#include alpha() and diversity() from the microbiome package

#The (Faith's) phylogenetic diversity of the microbiota data
#using the Picante package
#First, I need to extract the list of ASVs from this particular set of samples
#and use this to create a rooted phylogenetic tree with qiime2
ASV.table <- as.data.frame(otu_table(phylobj.filtered))
phylogenetic.diversity.microbiota <- pd(samp=t(ASV.table), tree=phy_tree(phylobj.filtered), include.root=TRUE)
phylogenetic.diversity.microbiota <- as.data.frame(phylogenetic.diversity.microbiota)
phylogenetic.diversity.microbiota$sample.id <- rownames(phylogenetic.diversity.microbiota)
phylogenetic.diversity.microbiota <- phylogenetic.diversity.microbiota[,-2]
colnames(phylogenetic.diversity.microbiota) <- c("microbiota.diversity.PD", "sample.id")
metadata.v8 <- merge(metadata.v8, phylogenetic.diversity.microbiota, by="sample.id", all=FALSE)

#Make columns of interest numeric
metadata.v8$Lactose <- as.numeric(metadata.v8$Lactose)
metadata.v8$microbiota.diversity.shannon <- as.numeric(metadata.v8$microbiota.diversity.shannon)
metadata.v8$microbiota.diversity.inversesimpson <- as.numeric(metadata.v8$microbiota.diversity.inversesimpson)
metadata.v8$microbiota.diversity.PD <- as.numeric(metadata.v8$microbiota.diversity.PD)

#STEP 13: Run linear mixed effect model to determine if Lactose is a predictor 
#of microbial alpha diversity
hist(metadata.v8$Lactose) #a little left skewed
hist(metadata.v8$microbiota.diversity.shannon) #a little left skewed
hist(metadata.v8$microbiota.diversity.inversesimpson) #heavily right skewed
hist(metadata.v8$microbiota.diversity.PD) #a little right skewed

summary(lmer(microbiota.diversity.shannon~Lactose+(1|cow), data=metadata.v8)) #singular fit warning, not significant
allFit(lmer(microbiota.diversity.shannon~Lactose+(1|cow), data=metadata.v8)) #singular fit warning, not significant

summary(lmer(microbiota.diversity.PD~Lactose+(1|cow), data=metadata.v8)) #singular fit warning, not significant
(allFit(lmer(microbiota.diversity.PD~Lactose+(1|cow), data=metadata.v8)))[["bobyqa"]]%>% summary  #model failed to converge, not significant

metadata.v8$cow <- as.factor(metadata.v8$cow)
rmcorr(participant=cow, measure1=Lactose, measure2=microbiota.diversity.shannon, data=metadata.v8) #not significant
rmcorr(participant=cow, measure1=Lactose, measure2=microbiota.diversity.PD, data=metadata.v8) #not significant

#STEP 14: Test for correlations between lactose and genus-level taxa
#Combine all of the data at the Genus level
phylobj.genus <- tax_glom(phylobj.filtered, 'Genus', NArm=F)
dim(otu_table(phylobj.genus)) #There are 1162 taxa and 174 samples

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

melted.phylgenus <- psmelt(phylobj.genus)
melted.phylgenus$Lactose <- as.numeric(melted.phylgenus$Lactose)
length(unique(melted.phylgenus$cow)) #58 cows
length(unique(melted.phylgenus$Sample)) #174 samples

#Function to test for correlation between lactose and microbial genus abundance
Lactose.genus.GLMM.function <- function(genus){
  coefficient.pvalue <- data.frame(message=NA,coefficient=NA, pvalue=NA, taxon=NA, BMO=NA)
  for(i in genus){
    OTU.abundance <- melted.phylgenus[melted.phylgenus$OTU==i,]
    output.lactose <- try(summary(glmer(Abundance~Lactose+(1|cow), data=OTU.abundance, family=poisson)), silent=TRUE)
    output.lactose.subset <- try(data.frame(message=paste(output.lactose[["optinfo"]][["conv"]][["lme4"]][["messages"]], collapse=" "), coefficient=output.lactose[["coefficients"]][2,1], pvalue=output.lactose[["coefficients"]][2,4], taxon=i, BMO="Lactose"), silent=TRUE)
    coefficient.pvalue <- try(rbind(coefficient.pvalue, output.lactose.subset), silent=TRUE)
    }
  print(coefficient.pvalue)   
}

#First, looking at the prevalence of the genera
prevalence.per.genus <- melted.phylgenus[melted.phylgenus$Abundance>0, ]
prevalence.per.genus <- data.frame(tapply(prevalence.per.genus$Abundance, prevalence.per.genus$OTU, length))
prevalence.per.genus$OTUnames <- rownames(prevalence.per.genus)
colnames(prevalence.per.genus) <- c("prevalence", "OTU")

#Testing for just genera that are present in at least 50% of the samples
genera.50percent <-prevalence.per.genus[prevalence.per.genus$prevalence>=(174*0.5),]
genera.50percent <- unique(genera.50percent$OTU)

#And of those genera (present in 50% of samples), selecting the genera that are in 
#at least 90 % of cows (i.e. >= 90 % of cows have that genus for at least one time point)
nonzero.counts.per.cow.per.genus <-melted.phylgenus[which(melted.phylgenus$OTU %in%genera.50percent),]
nonzero.counts.per.cow.per.genus <- nonzero.counts.per.cow.per.genus[nonzero.counts.per.cow.per.genus$Abundance>0, ]

Function.NonZeroCounts.PerCow.PerGenus<-function(genus){
  df <- data.frame()
  for(i in genus){
    OTU.subset <- nonzero.counts.per.cow.per.genus[nonzero.counts.per.cow.per.genus$OTU==i,]
    number.cows<-length(unique(OTU.subset$cow))
    percentage.cows <- (number.cows/58)*100
    OTU.cows <-data.frame(i, percentage.cows)  
    df <- rbind(df, OTU.cows)
  }
  print(df)
}

nonzero.counts.per.cow.per.genus.v2 <- Function.NonZeroCounts.PerCow.PerGenus(unique(nonzero.counts.per.cow.per.genus$OTU))
colnames(nonzero.counts.per.cow.per.genus.v2) <-c("OTU", "percent.of.cows")
genera.50PercentSamples.90PercentCows <- nonzero.counts.per.cow.per.genus.v2[nonzero.counts.per.cow.per.genus.v2$percent.of.cows>=90,][["OTU"]]

#Add the full taxonomy to the list of genera tested for relationships with SCC,
#then write to file 
full.taxonomy.genus <- data.frame(tax_table(phylobj.genus))
full.taxonomy.genus$Genus <- rownames(full.taxonomy.genus)
full.taxonomy.genus <- full.taxonomy.genus[which(full.taxonomy.genus$Genus %in% genera.50PercentSamples.90PercentCows), ]
full.taxonomy.genus$collapsed.taxonomy <- paste(full.taxonomy.genus$Kingdom, full.taxonomy.genus$Phylum, full.taxonomy.genus$Class, full.taxonomy.genus$Order, full.taxonomy.genus$Family, full.taxonomy.genus$Genus, sep="")
full.taxonomy.genus <- full.taxonomy.genus[,c(6,8)]
write.csv(full.taxonomy.genus, "Genera_tested_for_relationships_with_Lactose.csv")


results.GLMM.prevalentGenera<-Lactose.genus.GLMM.function(genera.50PercentSamples.90PercentCows)
results.GLMM.prevalentGenera.v2 <- results.GLMM.prevalentGenera[-1,]
results.GLMM.prevalentGenera.v2$pvalue <- as.numeric(results.GLMM.prevalentGenera.v2$pvalue)
#there were a few of the models that failed so trying allFit
UCG10<-melted.phylgenus[melted.phylgenus$OTU=="g__UCG.010",]
(allFit(glmer(Abundance~Lactose+(1|cow), data=UCG10, family=poisson)))[["bobyqa"]]%>%summary #model still failed to converge

Eubacterium<-melted.phylgenus[melted.phylgenus$OTU=="g__.Eubacterium._coprostanoligenes_group",]
(allFit(glmer(Abundance~Lactose+(1|cow), data=Eubacterium, family=poisson)))[["bobyqa"]]%>%summary #model still failed to converge

#removing those two results from the list and performing multiple test correction
results.GLMM.prevalentGenera.v3 <- results.GLMM.prevalentGenera.v2[-(grep(pattern="UCG[.]010|Eubacterium", results.GLMM.prevalentGenera.v2$taxon)),]

results.GLMM.prevalentGenera.v3$fdr <- p.adjust(p=results.GLMM.prevalentGenera.v3$pvalue, method="fdr")
results.GLMM.prevalentGenera.v3$bonferroni <- p.adjust(p=results.GLMM.prevalentGenera.v3$pvalue, method="bonferroni")

write.csv(results.GLMM.prevalentGenera.v3, file="PrevalentGenera.Correlation.Lactose.csv")

#trying GLMM with slope varying by cow just to see if the relationships are similar
Lactose.genus.GLMM.function.v2 <- function(genus){
  coefficient.pvalue <- data.frame(message=NA, coefficient=NA, pvalue=NA, taxon=NA, BMO=NA)
  for(i in genus){
    OTU.abundance <- melted.phylgenus[melted.phylgenus$OTU==i,]
    output.lactose <- try(summary(glmer(Abundance~Lactose+(Lactose|cow), data=OTU.abundance, family=poisson)), silent=TRUE)
    output.lactose.subset <- try(data.frame(message=paste(output.lactose[["optinfo"]][["conv"]][["lme4"]][["messages"]], collapse=" "),coefficient=output.lactose[["coefficients"]][2,1], pvalue=output.lactose[["coefficients"]][2,4], taxon=i, BMO="Lactose"), silent=TRUE)
    coefficient.pvalue <- try(rbind(coefficient.pvalue, output.lactose.subset), silent=TRUE)
  }
  print(coefficient.pvalue)   
}

results.GLMM.prevalentGenera.v4<-Lactose.genus.GLMM.function.v2(genera.50PercentSamples.90PercentCows)
results.GLMM.prevalentGenera.v4 <- results.GLMM.prevalentGenera.v4[-1,]
results.GLMM.prevalentGenera.v4$pvalue <- as.numeric(results.GLMM.prevalentGenera.v4$pvalue)
results.GLMM.prevalentGenera.v4$fdr <- p.adjust(p=results.GLMM.prevalentGenera.v4$pvalue, method="fdr")
results.GLMM.prevalentGenera.v4$bonferroni <- p.adjust(p=results.GLMM.prevalentGenera.v4$pvalue, method="bonferroni")

#STEP 15: Make a barplot of the estimated beta coefficients for lactose
#as a fixed effect on the genera
significant.subset <- results.GLMM.prevalentGenera.v3[results.GLMM.prevalentGenera.v3$bonferroni<=0.01,]
significant.subset$Genus <- significant.subset$taxon
full.taxonomy.genus <- data.frame(tax_table(phylobj.genus))
full.taxonomy.genus$Genus <- rownames(full.taxonomy.genus)
select.taxonomy.genus.lactose <- full.taxonomy.genus[which(full.taxonomy.genus$Genus %in% significant.subset$Genus), ]
significant.subset.v2 <-merge(significant.subset,select.taxonomy.genus.lactose, by="Genus")
significant.subset.v2$Genus <- gsub(pattern="g__", replacement="", significant.subset.v2$Genus)
significant.subset.v2$Genus <- gsub(pattern="f__", replacement="", significant.subset.v2$Genus)
significant.subset.v2$Class <- gsub(pattern="c__", replacement="", significant.subset.v2$Class)

ggplot(aes(x=coefficient, y=Genus), data=significant.subset.v2)+geom_bar(position="stack", stat="identity")+geom_vline(xintercept=0)+xlab("estimated beta coefficient for lactose")+theme(plot.title = element_text(lineheight = 0.5, size=10, hjust=0.45))+facet_grid(vars(Class), scales="free_y", space="free_y")+theme(strip.text.y=element_text(size=8,angle=0.45))
ggsave(filename="Barplot.Genus.Lactose.correlations.jpeg", width=10, height=5)
