#This script contains code used to test for correlations
#between BMO relative abundance and rarefied microbiota counts. 

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
library(lme4)
library(lmerTest)
library(rmcorr)
library(prediction)
library(ggforce)
library(purrr)
library(DHARMa)
library(vegan)
library(picante)

#STEP 1:Reading in the metadata table
metadata  <- read.csv("../../../DairyGrandChallenge.Study1.2017-2018/ResearchMaterials.DairyGrandChallenge/Cleaned_For_MetadataFile/metadata.csv")

#keeping only the data that is needed. 
metadata<- metadata[,c(1:15)]
metadata <- metadata[-1,]

#STEP 2: Read in the table with number of reads per sample, before rarefaction. 
reads.per.sample <- read.csv("sample-frequency.csv", head=FALSE)
colnames(reads.per.sample) <- c("sample.id", "reads")
metadata.v2 <- merge(metadata, reads.per.sample, by="sample.id", all=FALSE)
dim(metadata.v2) #508 samples
length(unique(metadata.v2$cow))#111 cows

#STEP 3: Read in the milk oligosaccharide dataset and edit it. 
BMO <- readxl::read_xlsx("USDA DGC Project 1 - BMO data (Barile Lab) - 11.17.2020 - SDD.xlsx", sheet=1)
BMO <- BMO[,-c(2:7)]
BMO$`Sample ID 1` <- gsub(pattern="_Davis", x=BMO$`Sample ID 1`, replacement="")
BMO$`Sample ID 1` <- gsub(pattern="_171125", x=BMO$`Sample ID 1`, replacement="-D1")
BMO$`Sample ID 1` <- gsub(pattern="_171126", x=BMO$`Sample ID 1`, replacement="-D2")
BMO$`Sample ID 1` <- gsub(pattern="_180104", x=BMO$`Sample ID 1`, replacement="-D3")
BMO$`Sample ID 1` <- gsub(pattern="_180105", x=BMO$`Sample ID 1`, replacement="-D4")
BMO$`Sample ID 1` <- gsub(pattern="_180314", x=BMO$`Sample ID 1`, replacement="-D5")
BMO$`Sample ID 1` <- gsub(pattern="_180315", x=BMO$`Sample ID 1`, replacement="-D6")

BMO<-rename(BMO, sample.id = `Sample ID 1`)

#remove BMO 3_3_1_0_0 in column 11 because it has "NA" values
BMO <- BMO[,-11]

#STEP 4:Merge BMO dataset with metadata
metadata.v2$sample <- metadata.v2$sample.id
metadata.v2$sample <- gsub(pattern="-B", replacement="", metadata.v2$sample)

BMO$sample <- BMO$sample.id
metadata.v3 <- merge(metadata.v2, BMO, by="sample", all=FALSE)
dim(metadata.v3) #368 samples

#can get rid of "sample" column and "sample.id.y" column and rename "sample.id.x"
metadata.v3 <- metadata.v3[,-c(1,18)]
metadata.v3 <- rename(metadata.v3, "sample.id"="sample.id.x")

#STEP 5:Rename sample 5849-D1 to 5849-D2 because cow 5849 actually missed
#the milk collection on November 25, 2017, but milk was collected from this cow
#on November 26, 2017. The same was true for cow 6233. 
metadata.v3$sample.id <- gsub(pattern="5849-D1", replacement="5849-D2", x=metadata.v3$sample.id)
metadata.v3$sample.id <- gsub(pattern="6233-D1", replacement="6233-D2", x=metadata.v3$sample.id)
metadata.v3[metadata.v3$sample.id=="5849-D2","day"] <- 2
metadata.v3[metadata.v3$sample.id=="6233-D2","day"] <- 2

#STEP 6: Loading in the feature table that's been rarefied and filtered
feature.table <- read_qza("feature-table-rarefied-at-507-samples-for-analysis.qza")
feature.table <- as.data.frame(feature.table$data)

#STEP 7: Rename 5849-D1 to 5849-D2
feature.table.v2 <- feature.table
library(dplyr)
feature.table.v2 <- rename(feature.table.v2, "5849-D2"="5849-D1")

#STEP 8: Filter the metadata samples down to what is the feature table
metadata.v4 <- metadata.v3[which(metadata.v3$sample.id %in% colnames(feature.table.v2)),]
dim(metadata.v4) #279 samples


#STEP 9: Filter the metadata samples to (a) remove lactose outliers (cows #6229 and 5651), 
#(b) if a cow has more than one sample in a period then, select the sample with
#higher reads (before rarefying) for a cow in a period. 

#(a) remove lactose outlier cows
metadata.v5 <- metadata.v4[metadata.v4$cow!="6229",]
metadata.v5 <- metadata.v5[metadata.v5$cow!="5651",]
dim(metadata.v5) #279 samples

#(b)for cows that have more than one sample in a period, select the sample
#that had the higher number of reads before rarefying.

ChooseSampleFromPeriod <-function(cow) {
    NoDuplicates <- data.frame()
    for(i in cow) {
        select.cow <- metadata.v5[metadata.v5$cow==i,]
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
cows <-unique(metadata.v5$cow)
metadata.v6 <- ChooseSampleFromPeriod(cows) 
length(metadata.v6$sample.id) #160 samples
length(unique(metadata.v6$cow)) #56 cows

#STEP 10: Read in taxonomy
taxonomy <- read_qza("taxonomy-PE-no-singletons-blanks-crmc-decontaminated-no-Bostaurus-mito-chloro-eukary.qza")

#change it to the way phyloseq wants it 
taxonomy2 <- as.data.frame(taxonomy$data) %>% column_to_rownames("Feature.ID") %>% separate(Taxon, c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"), sep = ";") %>% as.matrix()
#remove "Confidence" column
taxonomy2 <- taxonomy2[,-8]

#trim taxonomy to include only ASVs in the ASV table
taxonomy3 <- taxonomy2[c(rownames(taxonomy2)%in%rownames(feature.table.v2)),]

#Replace all empty levels with NA for the sake of consistency
#This prevents inaccurate glomming later
taxonomy3[,grep("Species", colnames(taxonomy3))][grep("s__$", taxonomy3[,grep("Species", colnames(taxonomy3))])] <- NA
taxonomy3[,grep("Genus", colnames(taxonomy3))][grep("g__$", taxonomy3[,grep("Genus", colnames(taxonomy3))])] <- NA
taxonomy3[,grep("Family", colnames(taxonomy3))][grep("f__$", taxonomy3[,grep("Family", colnames(taxonomy3))])] <- NA
taxonomy3[,grep("Order", colnames(taxonomy3))][grep("o__$", taxonomy3[,grep("Order", colnames(taxonomy3))])] <- NA
taxonomy3[,grep("Class", colnames(taxonomy3))][grep("c__$", taxonomy3[,grep("Class", colnames(taxonomy3))])] <- NA
taxonomy3[,grep("Phylum", colnames(taxonomy3))][grep("p__$", taxonomy3[,grep("Phylum", colnames(taxonomy3))])] <- NA

#STEP 11: Read in rooted tree. 
tree <- read_qza("rooted-tree-rarefied507-samplesforanalyses.qza")

#STEP 12: build a phyloseq object
rownames(metadata.v6)=NULL
phylobj <- phyloseq(otu_table(feature.table.v2, taxa_are_rows=T),tax_table(taxonomy3), phy_tree(tree$data), sample_data(metadata.v6 %>% as.data.frame() %>% column_to_rownames("sample.id")))

#keep only the taxa that have at least one count
phylobj.filtered <- filter_taxa(phylobj, function(x) sum(x)>0, TRUE)

#Note: for now, cows that don't have a sample in every period are still being
#retained in the dataset for beta-diversity analyses because the beta-diversity 
#correlations are being performed by period, separately. In the analyses that
#are not broken up by period though, cows lacking samples from a period will be
#removed from the dataset. 


#STEP 13: Melt the phyloseq object to obtain a dataframe that can be used with vegan
#to generate dissimilarity matrices
melted.phylobj.ASV <- psmelt(phylobj.filtered)

#Change the data frame from long to wide format
melted.phylobj.ASV.wide <- pivot_wider(data=melted.phylobj.ASV, id_cols=Sample, names_from=OTU, values_from=Abundance)


#STEP 14: merge with the metadata to get back the other metadata measures
#including the BMO abundances.
#First, rename the metadata sample.id column name to Sample
metadata.v7 <- rename(metadata.v6, Sample = sample.id)
dim(metadata.v7)
#then merge by "Sample" column name
melted.phylobj.ASV.wide.v2 <- merge(metadata.v7, melted.phylobj.ASV.wide, by="Sample", all=FALSE)
dim(melted.phylobj.ASV.wide.v2) #All 160 retained during merging
length(unique(melted.phylobj.ASV.wide.v2$cow)) #56 cows

#STEP 15: Separate the data frame into periods
period0 <- melted.phylobj.ASV.wide.v2[melted.phylobj.ASV.wide.v2$period=="p_0",]
period1 <- melted.phylobj.ASV.wide.v2[melted.phylobj.ASV.wide.v2$period=="p_1",]
period2 <- melted.phylobj.ASV.wide.v2[melted.phylobj.ASV.wide.v2$period=="p_2",]

#STEP 16:Calculate bray-curtis dissimilarity metrics for BMOs and microbiota 
#in each period. 

period0.BMO <-period0[,17:35] 
period0.ASV <-period0[,36:7442]
period0.BMO <- data.frame(lapply(period0.BMO, as.numeric))
period0.ASV <- data.frame(lapply(period0.ASV, as.numeric))

set.seed(1)
dist.period0.BMO <-vegdist(period0.BMO, method="bray")
dist.period0.ASV <- vegdist(period0.ASV, method="bray")

period1.BMO <-period1[,17:35] 
period1.ASV <-period1[,36:7442]
period1.BMO <- data.frame(lapply(period1.BMO, as.numeric))
period1.ASV <- data.frame(lapply(period1.ASV, as.numeric))

set.seed(1)
dist.period1.BMO <-vegdist(period1.BMO, method="bray")
dist.period1.ASV <- vegdist(period1.ASV, method="bray")

period2.BMO <-period2[,17:35] 
period2.ASV <-period2[,36:7442]
period2.BMO <- data.frame(lapply(period2.BMO, as.numeric))
period2.ASV <- data.frame(lapply(period2.ASV, as.numeric))

set.seed(1)
dist.period2.BMO <-vegdist(period2.BMO, method="bray")
dist.period2.ASV <- vegdist(period2.ASV, method="bray")

#STEP 17: Investigate correlations between milk microbiota beta-diversity 
#and BMO beta-diversity
set.seed(1)
mantel.period0 <- mantel(dist.period0.ASV, dist.period0.BMO, method="spearman", na.rm=FALSE)
mantel.period0 #not significant

set.seed(1)
mantel.period1 <- mantel(dist.period1.ASV, dist.period1.BMO, method="spearman", na.rm=FALSE)
mantel.period1 #not significant

set.seed(1)
mantel.period2 <- mantel(dist.period2.ASV, dist.period2.BMO, method="spearman", na.rm=FALSE)
mantel.period2 #significant

#p-value adjustment
p.adjust(c(mantel.period0[["signif"]],mantel.period1[["signif"]], mantel.period2[["signif"]]),method="fdr")
#not significant after multiple test correction

#STEP 18: Procrustes analysis using the same Bray-Curtis dissimilarity matrices
set.seed(1)
mds.period0.BMO <- monoMDS(dist.period0.BMO, y=cmdscale(dist.period0.BMO))
mds.period0.ASV <- monoMDS(dist.period0.ASV, y=cmdscale(dist.period0.ASV))
protest(X=mds.period0.ASV, Y=mds.period0.BMO, scale=TRUE) #not significant

set.seed(1)
mds.period1.BMO <- monoMDS(dist.period1.BMO, y=cmdscale(dist.period1.BMO))
mds.period1.ASV <- monoMDS(dist.period1.ASV, y=cmdscale(dist.period1.ASV))
protest(X=mds.period1.ASV, Y=mds.period1.BMO, scale=TRUE) #not significant

set.seed(1)
mds.period2.BMO <- monoMDS(dist.period2.BMO, y=cmdscale(dist.period2.BMO))
mds.period2.ASV <- monoMDS(dist.period2.ASV, y=cmdscale(dist.period2.ASV))
protest(X=mds.period2.ASV, Y=mds.period2.BMO, scale=TRUE) #not significant


#STEP 19: Moving on to linear mixed effects modeling and repeated measures 
#correlation analyses, so we now need to remove cows with less than three samples
#(i.e. cows that don't have a sample from each period.)
 
cows.pd0 <- unique(metadata.v7[metadata.v7$period=="p_0",][["cow"]])
cows.pd1 <- unique(metadata.v7[metadata.v7$period=="p_1",][["cow"]])
cows.pd2 <- unique(metadata.v7[metadata.v7$period=="p_2",][["cow"]])
cows.to.keep <- cows.pd0[which(cows.pd0 %in% cows.pd1)]
cows.to.keep.v2 <- cows.to.keep[which(cows.to.keep %in% cows.pd2)]
length(cows.to.keep.v2) #49 cows
metadata.v8 <- metadata.v7[which(metadata.v7$cow %in% cows.to.keep.v2),]
dim(metadata.v8) #147 samples
length(unique(metadata.v8$cow)) #49 cows


#STEP 20:Calculate the inverse simpson and shannon measures of alpha diversity for the
#BMOs in each sample

#First, need to unload the "microbiome" package because it somehow messes up
#the alpha diversity calculation with the "vegan" package. 
unload("microbiome")
metadata.v8$sample.id <- metadata.v8$Sample

alpha.diversity.BMO<- function(sample) {
    metdat <- data.frame()
    for (i in sample) {
        BMO.by.sample <- data.frame(metadata.v8[metadata.v8$sample.id==i,])
        BMO.by.sample2 <- BMO.by.sample[,17:35]
        BMO.by.sample2 <- as.numeric(BMO.by.sample2)
        BMO.by.sample$BMO.diversity.inversesimpson<-diversity(x=BMO.by.sample2, index="invsimpson")
        BMO.by.sample$BMO.diversity.shannon<-diversity(x=BMO.by.sample2, index="shannon")
        metdat <- rbind(metdat, BMO.by.sample) 
    }
    metdat
}

metadata.v9<-alpha.diversity.BMO(metadata.v8$sample.id)

#STEP 21: Build new phyloseq object. 
rownames(metadata.v9) = NULL
#remove the "Sample" column, since the sample.id column is already in the table
metadata.v9 <- metadata.v9[,-1]
phylobj.v2 <- phyloseq(otu_table(feature.table.v2, taxa_are_rows=T),tax_table(taxonomy3), phy_tree(tree$data), sample_data(metadata.v9 %>% as.data.frame() %>% column_to_rownames("sample.id")))
phylobj.v2.filtered <- filter_taxa(phylobj.v2, function(x) sum(x)>0, TRUE)


#STEP 22: Calculate alpha diversity for microbiota data..
#The simpson and shannon diversity of the microbiota
alpha.diversity.microbiota<-estimate_richness(physeq=phylobj.v2.filtered, measures=c("InvSimpson", "Shannon")) #from the phyloseq package

colnames(alpha.diversity.microbiota)<- c("microbiota.diversity.shannon", "microbiota.diversity.inversesimpson")
alpha.diversity.microbiota$sample.id <- rownames(alpha.diversity.microbiota)
metadata.v10 <- merge(metadata.v9, alpha.diversity.microbiota, by="sample.id", all=FALSE)

#Other functions to calculate alpha diversity from phyloseq object
#include alpha() and diversity() from the microbiome package

#The (Faith's) phylogenetic diversity of the microbiota data
#using the Picante package

ASV.table <- as.data.frame(otu_table(phylobj.v2.filtered))
phylogenetic.diversity.microbiota <- pd(samp=t(ASV.table), tree=phy_tree(phylobj.v2.filtered), include.root=TRUE)
phylogenetic.diversity.microbiota <- as.data.frame(phylogenetic.diversity.microbiota)
phylogenetic.diversity.microbiota$sample.id <- rownames(phylogenetic.diversity.microbiota)
phylogenetic.diversity.microbiota <- phylogenetic.diversity.microbiota[,-2]
colnames(phylogenetic.diversity.microbiota) <- c("microbiota.diversity.PD", "sample.id")
metadata.v11 <- merge(metadata.v10, phylogenetic.diversity.microbiota, by="sample.id", all=FALSE)
dim(metadata.v11) #147 samples

#STEP 23: Look at distribution of alpha diversity
hist(metadata.v11$BMO.diversity.shannon) # little left skewed
hist(metadata.v11$BMO.diversity.inversesimpson) #pretty normal
hist(metadata.v11$microbiota.diversity.shannon) #little left skewed
hist(metadata.v11$microbiota.diversity.inversesimpson) # right skewed
hist(metadata.v11$microbiota.diversity.PD) # a little right skewed


#STEP 24: linear mixed effect model
lmer.MicrobeSimpson.BMOsimpson<-lmer(data=metadata.v11, microbiota.diversity.inversesimpson~BMO.diversity.inversesimpson+(1|cow))
allfit.lmer.MicrobeSimpson.BMOsimpson <- allFit(lmer.MicrobeSimpson.BMOsimpson) #singular fit warning
allfit.lmer.MicrobeSimpson.BMOsimpson.summary<-summary(allfit.lmer.MicrobeSimpson.BMOsimpson[["bobyqa"]]) #BMO diversity not predictive of microbiota diversity

lmer.MicrobeShannon.BMOshannon <-lmer(data=metadata.v11, microbiota.diversity.shannon~BMO.diversity.shannon+(1|cow))
allfit.lmer.MicrobeShannon.BMOshannon <- allFit(lmer.MicrobeShannon.BMOshannon) #singular fit warning
allfit.lmer.MicrobeShannon.BMOshannon.summary<-summary(allfit.lmer.MicrobeShannon.BMOshannon[["bobyqa"]]) #BMO diversity not predictive of microbiota diversity

metadata.v11$cow <- as.factor(metadata.v11$cow)
rmcorr(participant=cow, measure1=BMO.diversity.shannon, measure2=microbiota.diversity.shannon, dataset = metadata.v11) #not significant

lmer.MicrobePD.BMOsimpson<-lmer(data=metadata.v11, microbiota.diversity.PD~BMO.diversity.inversesimpson+(1|cow))
allfit.lmer.MicrobePD.BMOsimpson <- allFit(lmer.MicrobePD.BMOsimpson) 
allfit.lmer.MicrobePD.BMOsimpson.summary <-summary(allfit.lmer.MicrobePD.BMOsimpson[["bobyqa"]]) #BMO diversity maybe predictive of BMO diversity
hist(resid(allfit.lmer.MicrobePD.BMOsimpson[["bobyqa"]]))
rmcorr(participant=cow, measure1=BMO.diversity.inversesimpson, measure2=microbiota.diversity.PD, dataset=metadata.v11) #significant relationship

ggplot(data=metadata.v11,aes(x=BMO.diversity.inversesimpson, y=microbiota.diversity.PD, group=cow))+geom_point()+xlab("BMO diversity, inverse simpson")+ylab("microbiota diversity, PD") +geom_line(aes(y=predict(lmer.MicrobePD.BMOsimpson), group=cow))

lmer.MicrobePD.BMOshannon<-lmer(data=metadata.v11, microbiota.diversity.PD~BMO.diversity.shannon+(1|cow))
allfit.lmer.MicrobePD.BMOshannon <- allFit(lmer.MicrobePD.BMOshannon) 
allfit.lmer.MicrobePD.BMOshannon.summary<-summary(allfit.lmer.MicrobePD.BMOshannon[["bobyqa"]]) #BMO diversity not predictive of microbiota diversity
hist(resid(allfit.lmer.MicrobePD.BMOshannon[["bobyqa"]]))
rmcorr(participant=cow, measure1=BMO.diversity.shannon, measure2=microbiota.diversity.PD, dataset=metadata.v11) #significant relationship

ggplot(data=metadata.v11,aes(x=BMO.diversity.shannon, y=microbiota.diversity.PD, group=cow))+geom_point()+xlab("BMO diversity, shannon")+ylab("microbiota diversity, PD") +geom_line(aes(y=predict(lmer.MicrobePD.BMOshannon), group=cow))

#p-value adjustment
p.adjust(c((allfit.lmer.MicrobeSimpson.BMOsimpson.summary[["coefficients"]][2,5]),(allfit.lmer.MicrobeShannon.BMOshannon.summary[["coefficients"]][2,5]),(allfit.lmer.MicrobePD.BMOshannon.summary[["coefficients"]][2,5]), (allfit.lmer.MicrobePD.BMOsimpson.summary[["coefficients"]][2,5])), method="fdr") #not significant after multiple test correction

#STEP 25: Identifying correlations between single BMOs 
#and microbial alpha diversity

#First, change the BMO measures from character to numeric, using a function 
#with a for loop.
Make.BMO.numeric.class <-function(BMO) {
    for(i in BMO) {
        metadata.v11[[i]]<-as.numeric(metadata.v11[[i]])
    } 
    return(metadata.v11)
}

BMO <- colnames(metadata.v11[,17:35])
metadata.v11 <- Make.BMO.numeric.class(BMO)


#Creating a function to evaluate the correlation between each BMO 
#and microbial alpha diversity. 
Function.BMO.MicrobialAlphadiv <- function(diversity, BMO.list){
    coefficient.pvalue <- data.frame(coefficient=NA, pvalue=NA, alpha_diversity_measure=NA, BMO=NA)
    for(i in diversity){
        for(j in BMO.list){
            form <- as.formula(paste0(i,"~",j,"+(1|cow)"))
            lmer.model <- try(lmer(form, data=metadata.v11), silent=TRUE)
            allfit.model <- try(allFit(lmer.model), silent=TRUE)
            model.summary <- try(summary(allfit.model[["bobyqa"]]), silent=TRUE)
            lmer.df <- try(data.frame(messages=paste(model.summary[["optinfo"]][["conv"]][["lme4"]][["messages"]], collapse=" "), coefficient=model.summary[["coefficients"]][2,1], pvalue=model.summary[["coefficients"]][2,5], alpha_diversity_measure=i, BMO=j), silent=FALSE)
            coefficient.pvalue <- try(rbind(coefficient.pvalue, lmer.df), silent=TRUE)
        }
    }
    print(coefficient.pvalue)
}
diversity <- colnames(metadata.v11[,38:40])
BMO.list <- colnames(metadata.v11[,17:35])
results.LME.BMO.AlphaDiv<-Function.BMO.MicrobialAlphadiv(diversity, BMO.list)
results.LME.BMO.AlphaDiv$pvalue <- as.numeric(results.LME.BMO.AlphaDiv$pvalue)
results.LME.BMO.AlphaDiv <- results.LME.BMO.AlphaDiv[!is.na(results.LME.BMO.AlphaDiv$pvalue),]
results.LME.BMO.AlphaDiv$fdr <- p.adjust(p=results.LME.BMO.AlphaDiv$pvalue, method="fdr") #no significant relationships


#STEP 26:Combine all of the data at the Genus level
phylobj.genus <- tax_glom(phylobj.v2.filtered, 'Genus', NArm=F)
dim(otu_table(phylobj.genus)) #There are 1056 taxa and 147 samples

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

#STEP 27: Preparing the phyloseq genus level object for ggplot and statistics
melted.phylgenus <- psmelt(phylobj.genus)

#Change the BMO measurements to numeric variable for count table
melted.phylgenus[,19:37] <- lapply(X=melted.phylgenus[,19:37], FUN=as.numeric)

#Change cow from integer to character
melted.phylgenus$cow <- as.character(melted.phylgenus$cow)

#STEP 28: Cut down to just the genera present in >=50% of samples and >=90% of cows.
prevalence.per.genus <- melted.phylgenus[melted.phylgenus$Abundance>0, ]
prevalence.per.genus <- data.frame(tapply(prevalence.per.genus$Abundance, prevalence.per.genus$OTU, length))
prevalence.per.genus$OTUnames <- rownames(prevalence.per.genus)
colnames(prevalence.per.genus) <- c("prevalence", "OTU")

genera.50percent <-prevalence.per.genus[prevalence.per.genus$prevalence>=(147*0.5),]
genera.50percent <- unique(genera.50percent$OTU)

nonzero.counts.per.cow.per.genus <-melted.phylgenus[which(melted.phylgenus$OTU %in%genera.50percent),]
nonzero.counts.per.cow.per.genus <- nonzero.counts.per.cow.per.genus[nonzero.counts.per.cow.per.genus$Abundance>0, ]

Function.NonZeroCounts.PerCow.PerGenus<-function(genus){
  df <- data.frame()
  for(i in genus){
  OTU.subset <- nonzero.counts.per.cow.per.genus[nonzero.counts.per.cow.per.genus$OTU==i,]
  number.cows<-length(unique(OTU.subset$cow))
  percentage.cows <- (number.cows/49)*100
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
write.csv(full.taxonomy.genus, "Genera_tested_for_relationships_with_BMOs.csv")

#STEP 29: Function that runs a generalized linear mixed effects model with poisson 
#distribution and allowing only intercept to vary by cow
Function.BMOpredictor.GenusResponse.GLME.v1 <-function(microbe, BMO) {
coefficient.pvalue <- data.frame(messages=NA, coefficient=NA, pvalue=NA, taxon=NA, BMO=NA)
    for(i in microbe){
    for(j in BMO){
      microbe.subset <- melted.phylgenus[melted.phylgenus$OTU==i,]
      form<-as.formula(paste0("Abundance~", j, "+(1|cow)"))
      glmer.model <- try(glmer(form, family=poisson, data=microbe.subset), silent=TRUE)
      model.summary <- try(summary(glmer.model), silent=TRUE)
      glmer.df <- try(data.frame(messages=paste(model.summary[["optinfo"]][["conv"]][["lme4"]][["messages"]], collapse=" "), coefficient=model.summary[["coefficients"]][2,1], pvalue=model.summary[["coefficients"]][2,4], taxon=i, BMO=j), silent=FALSE)
      coefficient.pvalue <- try(rbind(coefficient.pvalue, glmer.df), silent=TRUE)
    }
  }
  print(coefficient.pvalue)
}


BMOs <- colnames(melted.phylgenus[,19:37])
results.GLME.BMO.genera.v1<- Function.BMOpredictor.GenusResponse.GLME.v1(microbe=genera.50PercentSamples.90PercentCows, BMO=BMOs)
results.GLME.BMO.genera.v1 <- results.GLME.BMO.genera.v1[-1,] #removed first line because it was empty

#STEP 30: For those (few) models that failed, re-fitting the models 
#with allFit function (spoiler: some of these models still failed after refitting with allFit.) 

Function.BMOpredictor.GenusResponse.GLME.v1.optimizer <-function(microbe, BMO) {
    optimizer.choice <- data.frame(method="", optimizer="bobyqa")
    coefficient.pvalue.bobyqa <- data.frame(messages.bobyqa=NA, coefficient.bobyqa=NA, pvalue.bobyqa=NA, taxon=NA, BMO=NA)
    for(i in microbe){
        for(j in BMO){
            microbe.subset <- melted.phylgenus[melted.phylgenus$OTU==i,]
            form<-as.formula(paste0("Abundance~", j, "+(1|cow)"))
            glmer.model <- try(glmer(form, data=microbe.subset, family=poisson), silent=TRUE)
            glmer.optimization<-try(allFit(glmer.model, meth.tab=optimizer.choice), silent = TRUE)
            
            glmer.bobyqa <- try(summary(glmer.optimization[["bobyqa"]]), silent=TRUE)
            glmer.bobyqa.df <- try(data.frame(messages.bobyqa=paste(glmer.bobyqa[["optinfo"]][["conv"]][["lme4"]][["messages"]], collapse=" "), coefficient.bobyqa=glmer.bobyqa[["coefficients"]][2,1], pvalue.bobyqa=glmer.bobyqa[["coefficients"]][2,4], taxon=i, BMO=j), silent=FALSE)
            coefficient.pvalue.bobyqa <- try(rbind(coefficient.pvalue.bobyqa, glmer.bobyqa.df), silent=TRUE)
        }
    }
    print(coefficient.pvalue.bobyqa)
}

failed.convergence <- results.GLME.BMO.genera.v1
failed.convergence<-failed.convergence[grep(pattern="Model*", x=failed.convergence$messages),]
failed.taxon <- unique(failed.convergence$taxon)
failed.bmo <- unique(failed.convergence$BMO)

results.GLME.BMO.genera.v1.optimizer<- Function.BMOpredictor.GenusResponse.GLME.v1.optimizer(microbe=failed.taxon, BMO=failed.bmo) 
#some of the models still failed

#for the models that were resolved with using the optimizer, putting
#in the coefficient and pvalues for those optimized models, and getting
#rid of the models that failed even after optimization.
failed.convergence.v2 <- merge(failed.convergence, results.GLME.BMO.genera.v1.optimizer, by=c("taxon", "BMO"), all=FALSE)
success.convergence <- failed.convergence.v2[-c(grep(pattern="Model",failed.convergence.v2$messages.bobyqa)),]
success.convergence <- success.convergence[,c(6:8, 1,2)]
colnames(success.convergence) <- c("messages", "coefficient", "pvalue", "taxon", "BMO")

results.GLME.BMO.genera.v1.complete <- results.GLME.BMO.genera.v1[-c(grep(pattern="Model",results.GLME.BMO.genera.v1$messages)),]
results.GLME.BMO.genera.v1.complete <- rbind(results.GLME.BMO.genera.v1.complete, success.convergence)

results.GLME.BMO.genera.v1.complete$pvalue <- as.numeric(results.GLME.BMO.genera.v1.complete$pvalue)
results.GLME.BMO.genera.v1.complete$fdr <- p.adjust(p=results.GLME.BMO.genera.v1.complete$pvalue, method="fdr")
results.GLME.BMO.genera.v1.complete$bonferroni <- p.adjust(p=results.GLME.BMO.genera.v1.complete$pvalue, method="bonferroni")
write.csv(results.GLME.BMO.genera.v1.complete, file="BMOpredictor.GENUSresponse.GlmerPoisson.CowIntercept.csv")


#STEP 31: making plots of the BMOs and genera that are correlated, but just for a few of
#the relationships with the largest estimates for beta coefficients
lactobacillus <- melted.phylgenus[melted.phylgenus$OTU=="g__Lactobacillus",]
model.lactobacillus32000<- glmer(Abundance~X3_2_0_0_0+(1|cow), data=lactobacillus, family=poisson)
summary(model.lactobacillus32000)
ggplot(aes(x=X3_2_0_0_0, y=Abundance), data=lactobacillus)+geom_point()+geom_line(aes(y=predict(model.lactobacillus32000), group=cow))+ylab("Lactobacillus abundance")+xlab("3_2_0_0_0")

enterococcus <- melted.phylgenus[melted.phylgenus$OTU=="g__Enterococcus",]
model.enterococcus.36100 <- glmer(Abundance~X3_6_1_0_0+(1|cow), data=enterococcus, family=poisson)
summary(model.enterococcus.36100)
ggplot(aes(x=X3_6_1_0_0, y=Abundance), data=enterococcus)+geom_point()+geom_line(aes(y=predict(model.enterococcus.36100), group=cow))+ylab("Enterococcus abundance")+xlab("3_6_1_0_0")

paracoccus <- melted.phylgenus[melted.phylgenus$OTU=="g__Paracoccus",]
model.parococcus <- glmer(Abundance~X8_0_0_0_0+(1|cow), data=paracoccus, family=poisson)
summary(model.parococcus)
ggplot(aes(x=X8_0_0_0_0, y=Abundance), data=paracoccus)+geom_point()+geom_line(aes(y=predict(model.parococcus), group=cow))+ylab("Paracoccus abundance")+xlab("8_0_0_0_0")

acinetobacter <- melted.phylgenus[melted.phylgenus$OTU=="g__Acinetobacter",]
model.acinetobacter <- glmer(Abundance~X8_0_0_0_0+(1|cow), data=acinetobacter, family=poisson)
summary(model.acinetobacter)
ggplot(aes(x=X8_0_0_0_0, y=Abundance), data=acinetobacter)+geom_point()+geom_line(aes(y=predict(model.acinetobacter), group=cow))+ylab("Acinetobacter abundance")+xlab("8_0_0_0_0")

#STEP 32: Make stacked bar graphs summarizing the results. 

#Select only tests that resulted in Bonferroni-corrected pvalue <= 0.01. 
significant.subset <- results.GLME.BMO.genera.v1.complete[results.GLME.BMO.genera.v1.complete$bonferroni<=0.01,]
dim(significant.subset) #224 significant results. 

neutral.nonfucosylated.BMOs <- grep(pattern="._._0_0_0", x=unique(significant.subset$BMO), value=TRUE)
neutral.nonfucosylated.BMOs
neutral.fucosylated.BMOs <- grep(pattern="._._[1-9]_0_0", x=unique(significant.subset$BMO), value=TRUE)                                    
neutral.fucosylated.BMOs
sialylated.nonfucosylated.BMOs <- grep(pattern="X3..SL|X6..SL|._._0_[1-9]_[0-9]", x=unique(significant.subset$BMO), value=TRUE)
sialylated.nonfucosylated.BMOs

#For each genus, count the number of relationship with 1) neutral BMOs 2) acidic/sialilated BMOs (and non) 3) fucosylated BMOs 
Function.Summarize.byBMOcategory <- function(genus){
    df <- data.frame(Genus=NA, BMO_category=NA, number=NA)
    for(i in genus){
        taxon.subset <-significant.subset[significant.subset$taxon==i,]
        
        taxon.neutralnonfucosBMO.subset <- taxon.subset[which(taxon.subset$BMO %in% neutral.nonfucosylated.BMOs),]
        taxon.neutralnonfucosBMO.subset.positiveeffect <- taxon.neutralnonfucosBMO.subset[taxon.neutralnonfucosBMO.subset$coefficient>0,]
        taxon.neutralnonfucosBMO.subset.negativeeffect <- taxon.neutralnonfucosBMO.subset[taxon.neutralnonfucosBMO.subset$coefficient<0,]
        df.2a <- data.frame(Genus=i, BMO_category="neutral, non-fucosylated", number=nrow(taxon.neutralnonfucosBMO.subset.positiveeffect))
        df.2b <- data.frame(Genus=i, BMO_category="neutral, non-fucosylated", number=-(nrow(taxon.neutralnonfucosBMO.subset.negativeeffect)))
        df <- rbind(df, df.2a)
        df <- rbind(df, df.2b)
        
        taxon.neutralfucosBMO.subset <- taxon.subset[which(taxon.subset$BMO %in% neutral.fucosylated.BMOs),]
        taxon.neutralfucosBMO.subset.positiveeffect <- taxon.neutralfucosBMO.subset[taxon.neutralfucosBMO.subset$coefficient>0,]
        taxon.neutralfucosBMO.subset.negativeeffect <- taxon.neutralfucosBMO.subset[taxon.neutralfucosBMO.subset$coefficient<0,]
        df.3a <- data.frame(Genus=i, BMO_category="neutral, fucosylated", number=nrow(taxon.neutralfucosBMO.subset.positiveeffect))
        df.3b <- data.frame(Genus=i, BMO_category="neutral, fucosylated", number=-(nrow(taxon.neutralfucosBMO.subset.negativeeffect)))
        df <- rbind(df, df.3a)
        df <- rbind(df, df.3b)
        
        taxon.sialylnonfucosBMO.subset <- taxon.subset[which(taxon.subset$BMO %in% sialylated.nonfucosylated.BMOs),]
        taxon.sialylnonfucosBMO.subset.positiveeffect <- taxon.sialylnonfucosBMO.subset[taxon.sialylnonfucosBMO.subset$coefficient>0,]
        taxon.sialylnonfucosBMO.subset.negativeeffect <- taxon.sialylnonfucosBMO.subset[taxon.sialylnonfucosBMO.subset$coefficient<0,]
        df.4a <- data.frame(Genus=i, BMO_category="sialylated, non-fucosylated", number=nrow(taxon.sialylnonfucosBMO.subset.positiveeffect))
        df.4b <- data.frame(Genus=i, BMO_category="sialylated, non-fucosylated", number=-(nrow(taxon.sialylnonfucosBMO.subset.negativeeffect)))
        df <- rbind(df, df.4a)
        df <- rbind(df, df.4b)
    }
    print(df)
} 

summary.results.byBMO.category <- Function.Summarize.byBMOcategory(unique(significant.subset$taxon))
summary.results.byBMO.category <- summary.results.byBMO.category[!is.na(summary.results.byBMO.category$Genus),]

#Add the Class designations so that genera can be organized by Class
full.taxonomy.genus <- data.frame(tax_table(phylobj.genus))
full.taxonomy.genus$Genus <- rownames(full.taxonomy.genus)
select.taxonomy.genus.BMO <- full.taxonomy.genus[which(full.taxonomy.genus$Genus %in% summary.results.byBMO.category$Genus), ]
summary.results.byBMO.category.v2 <-merge(summary.results.byBMO.category,select.taxonomy.genus.BMO, by="Genus")

#Before graphing, get rid of "g__" before each genus name
summary.results.byBMO.category.v2$Genus <- gsub(pattern="g__", replacement="", summary.results.byBMO.category.v2$Genus)
summary.results.byBMO.category.v2$Genus <- gsub(pattern="f__", replacement="", summary.results.byBMO.category.v2$Genus)
summary.results.byBMO.category.v2$Class <- gsub(pattern="c__", replacement="", summary.results.byBMO.category.v2$Class)

ggplot(aes(x=number, y=Genus, fill=BMO_category), data=summary.results.byBMO.category.v2)+geom_bar(position="stack", stat="identity")+geom_vline(xintercept=0)+labs(fill="MO type")+xlab("number of MOs associated with a genus")+labs(title = "negative associations               positive associations")+theme(plot.title = element_text(lineheight = 0.5, size=10, hjust=0.45))+facet_grid(vars(Class), scales="free_y", space="free_y")+theme(strip.text.y=element_text(size=8,angle=0.45))
ggsave(filename="Barplot.Genus.BMO.correlations.v2.jpeg", width=10, height=5)


#STEP 33: Reading in lactose data and merging with count data to 
#re-run the models with lactose as a fixed effect. 
lactose <- readxl::read_xlsx("PMT DVH_Dairy Grand Challenge coded FOSS composition data_101818.xlsx")

#Looking at the first few lines of the file
head(lactose, 10) #The intended column headers are on the eighth row of the excel file

#In order to set the correct the column headers set, need to read in the dataset again but without the first row being used as the column headers
lactose <- readxl::read_xlsx("PMT DVH_Dairy Grand Challenge coded FOSS composition data_101818.xlsx", col_names=FALSE)

#Set the column headers to the values on the eighth row
colnames(lactose)<-lactose[8,]
#Remove the first eight rows
lactose <- lactose[-c(1:8),]

dim(lactose) #there are 76 observations and 46 variables in the table
colnames(lactose)
head(lactose) #The sample IDs are listed under the "Day" columns
#I can tell from the column names that the observations for days 2,3,4,5,6 are listed to the right of day 1 observations instead of below
#I can also see that the first two columns appear to be empty
sum(is.na(lactose[,1])) #all entries in the first column are empty/NA
sum(is.na(lactose[,2])) #all entries in the second column are empty/NA

#Remove the first two columns
lactose <- lactose[,-c(1:2)]
colnames(lactose) #no longer see those empty column headers

#Remove the "Period 1", "Period 2", "Cow ID", "Treatment 0", "treatment 1", "Treatment 2", "Stolen from", and "Location" columns because these variables will be retrieved from the Kalscheur lab data or at not needed in the compiled dataset
lactose <- lactose[,-c(1:4, 17, 30, 43, 44)]
head(lactose)
colnames(lactose)

#Need to subset the by day (i.e. take subset chunks of columns), then add back together under the same column headers so that there are only 
lactose.Day1 <- lactose[,1:6]
lactose.Day2 <- lactose[,7:12]
lactose.Day3 <- lactose[,13:18]
lactose.Day4 <- lactose[,19:24]
lactose.Day5 <- lactose[,25:30]
lactose.Day6 <- lactose[,31:36]

#Change the "Day" column headers in the data subsets to "sample_id"
lactose.Day1 <- rename(lactose.Day1, sample.id=`Day 1`)
lactose.Day2 <- rename(lactose.Day2, sample.id=`Day 2`)
lactose.Day3 <- rename(lactose.Day3, sample.id=`Day 3`)
lactose.Day4 <- rename(lactose.Day4, sample.id=`Day 4`)
lactose.Day5 <- rename(lactose.Day5, sample.id=`Day  5`)
lactose.Day6 <- rename(lactose.Day6, sample.id=`Day 6`)

#All of the column headers should be the same now for the subsets fo the lactose data
#Checking to see that they are the same..
colnames(lactose.Day1)==colnames(lactose.Day2) #All column headers are the same
colnames(lactose.Day1)==colnames(lactose.Day3) #All column headers are the same
colnames(lactose.Day1)==colnames(lactose.Day4) #All column headers are the same
colnames(lactose.Day1)==colnames(lactose.Day5) #All column headers are the same
colnames(lactose.Day1)==colnames(lactose.Day6) #All column headers are the same

#Now add the subsets together by rows under the same set of 6 column headers
lactose <- rbind(lactose.Day1, lactose.Day2, lactose.Day3, lactose.Day4, lactose.Day5, lactose.Day6)
dim(lactose) #456 observations and 6 variables (this is what's expected with 76 observations in each subset)

#Check for NA values and duplicate values in the sample_id column
sum(is.na(lactose$sample.id)) #the output is zero, so no NA values among the sample_id column
lactose.distinct.col1 <- unique(lactose$sample.id)
length(lactose.distinct.col1) #456 unique observations in sample_id so there are no duplicate sample IDs

#Remove the destination information ("_Wyn") from the sample IDs so that they can be matched and merged to other datasets
lactose$sample.id <-gsub(pattern="_Wyn", x=lactose$sample.id, replacement="") 
head(lactose)

#Replace date with "D1, D2, D3, D4, D5, D6" in sample.id  
lactose$sample.id <- gsub(pattern="_171125", x=lactose$sample.id, replacement="-D1")
lactose$sample.id <- gsub(pattern="_171126", x=lactose$sample.id, replacement="-D2")
lactose$sample.id <- gsub(pattern="_180104", x=lactose$sample.id, replacement="-D3")
lactose$sample.id <- gsub(pattern="_180105", x=lactose$sample.id, replacement="-D4")
lactose$sample.id <- gsub(pattern="_180314", x=lactose$sample.id, replacement="-D5")
lactose$sample.id <- gsub(pattern="_180315", x=lactose$sample.id, replacement="-D6")

#Add a column specifying the units of measure for the lactose data
lactose$nut_Units <-"gram/100 gram of milk "

#Rename sample 5849-D1 to 5849-D2 because cow 5849 actually missed
#the milk collection on November 25, 2017, but milk was collected from this cow
#on November 26, 2017. The same was true for cow 6233. 
lactose$sample.id <- gsub(pattern="5849-D1", replacement="5849-D2", x=lactose$sample.id)
lactose$sample.id <- gsub(pattern="6233-D1", replacement="6233-D2", x=lactose$sample.id)

lactose<-rename(lactose, Sample = sample.id)
#get rid of "-B" from any sample ID's because they won't match with the lactose dataset
melted.phylgenus.v2 <- melted.phylgenus
melted.phylgenus.v2$Sample <- gsub(pattern="-B", replacement="", melted.phylgenus.v2$Sample)
melted.phylgenus.v2 <- merge(melted.phylgenus.v2, lactose, by="Sample", all=FALSE)
dim(melted.phylgenus.v2)
dim(melted.phylgenus)
melted.phylgenus.v2$Lactose <- as.numeric(melted.phylgenus.v2$Lactose)

#STEP 34: Testing for correlations between BMOs and microbiota counts at 
#genus level, AND with addition of lactose as fixed effect

Function.BMO.Lactose.predictors.GenusResponse.GLME <-function(microbe, BMO) {
coefficient.pvalue <- data.frame(messages=NA, coefficient=NA, pvalue=NA, taxon=NA, BMO=NA)
for(i in microbe){
  for(j in BMO){
    microbe.subset <- melted.phylgenus.v2[melted.phylgenus.v2$OTU==i,]
    form<-as.formula(paste0("Abundance~", j,"+Lactose+(1|cow)"))
    glmer.model <- try(glmer(form, family=poisson, data=microbe.subset), silent=TRUE)
    model.summary <- try(summary(glmer.model), silent=TRUE)
    glmer.df <- try(data.frame(messages=paste(model.summary[["optinfo"]][["conv"]][["lme4"]][["messages"]], collapse=" "), coefficient=model.summary[["coefficients"]][2,1], pvalue=model.summary[["coefficients"]][2,4], taxon=i, BMO=j), silent=FALSE)
    coefficient.pvalue <- try(rbind(coefficient.pvalue, glmer.df), silent=TRUE)
  }
}
print(coefficient.pvalue)
}


results.GLME.BMO.lactose.genera.v1<- Function.BMO.Lactose.predictors.GenusResponse.GLME(microbe=genera.50PercentSamples.90PercentCows, BMO=BMOs)
results.GLME.BMO.lactose.genera.v1<- results.GLME.BMO.lactose.genera.v1[-1,] #removed first line because it was empty


#STEP 35: There were some models that failed to converge so for those models, trying 
#allFit with bobyqa optimizer
Function.BMO.Lactose.predictors.GenusResponse.GLME.v1.optimizer <-function(microbe, BMO) {
  optimizer.choice <- data.frame(method="", optimizer="bobyqa")
  coefficient.pvalue.bobyqa <- data.frame(messages.bobyqa=NA, coefficient.bobyqa=NA, pvalue.bobyqa=NA, taxon=NA, BMO=NA)
  for(i in microbe){
    for(j in BMO){
      microbe.subset <- melted.phylgenus.v2[melted.phylgenus.v2$OTU==i,]
      form<-as.formula(paste0("Abundance~", j, "+Lactose+(1|cow)"))
      glmer.model <- try(glmer(form, data=microbe.subset, family=poisson), silent=TRUE)
      glmer.optimization<-try(allFit(glmer.model, meth.tab=optimizer.choice), silent = TRUE)
      
      glmer.bobyqa <- try(summary(glmer.optimization[["bobyqa"]]), silent=TRUE)
      glmer.bobyqa.df <- try(data.frame(messages.bobyqa=paste(glmer.bobyqa[["optinfo"]][["conv"]][["lme4"]][["messages"]], collapse=" "), coefficient.bobyqa=glmer.bobyqa[["coefficients"]][2,1], pvalue.bobyqa=glmer.bobyqa[["coefficients"]][2,4], taxon=i, BMO=j), silent=FALSE)
      coefficient.pvalue.bobyqa <- try(rbind(coefficient.pvalue.bobyqa, glmer.bobyqa.df), silent=TRUE)
    }
  }
  print(coefficient.pvalue.bobyqa)
}

failed.convergence <- results.GLME.BMO.lactose.genera.v1
failed.convergence<-failed.convergence[grep(pattern="Model*", x=failed.convergence$messages),]
failed.taxon <- unique(failed.convergence$taxon)
failed.bmo <- unique(failed.convergence$BMO)

results.GLME.BMO.lactose.genera.v1.optimizer<- Function.BMO.Lactose.predictors.GenusResponse.GLME.v1.optimizer(microbe=failed.taxon, BMO=failed.bmo)
results.GLME.BMO.lactose.genera.v1.optimizer <- results.GLME.BMO.lactose.genera.v1.optimizer[-1,] #remove the first row because it's empty
#all of the models still failed. 

#removing the models that failed, and then performing p-value adjustment for multiple tests.
results.GLME.BMO.lactose.genera.v2  <- results.GLME.BMO.lactose.genera.v1[-c(grep(pattern="Model", results.GLME.BMO.lactose.genera.v1$messages)),]

results.GLME.BMO.lactose.genera.v2$fdr <- p.adjust(results.GLME.BMO.lactose.genera.v2$pvalue, method="fdr")
results.GLME.BMO.lactose.genera.v2$bonferroni <-p.adjust(results.GLME.BMO.lactose.genera.v2$pvalue, method="bonferroni")

significant.subset.v2 <- results.GLME.BMO.lactose.genera.v2[results.GLME.BMO.lactose.genera.v2$bonferroni<0.01,]

comparing.models.with.vs.without.lactose <- merge(significant.subset, significant.subset.v2, by=c("taxon", "BMO"))

write.csv(comparing.models.with.vs.without.lactose, "BMO.vs.Genera.with.vs.without.lactose.fixedeffect.compare.results.csv")
