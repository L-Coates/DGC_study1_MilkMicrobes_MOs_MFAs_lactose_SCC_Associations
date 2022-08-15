#This script demonstrates the steps taken in R to determine relationships between milk fatty acids and milk microbiota from Holstein cows in the Dairy Grand Challenge study 1. 

#First, the libraries needed for the data manipulation, visualization and statistical analyses are loaded
library(dplyr)
library(lubridate)
library(qiime2R)
library(phyloseq)
library(tidyverse)
library(lmerTest)
library(lme4)
library(microbiome)
library(picante)
library(vegan)
library(devtools)
library(ggplot2)
library(ggpubr)
library(plotly)
library(psych)
library(cowplot)
library(lme4)
library(lmerTest)
library(prediction)
library(ggforce)
library(DHARMa)
library(rmcorr)
library(genepop)

#STEP 1:Reading in the metadata table and perform some editing. 
metadata  <- read.csv("../../../DairyGrandChallenge.Study1.2017-2018/ResearchMaterials.DairyGrandChallenge/Cleaned_For_MetadataFile/metadata.csv")

#keeping only the data that is needed. 
metadata<- metadata[,c(1:15)]
metadata <- metadata[-1,]

#STEP 2: Read in the number of reads per sample (before rarefaction) which will 
#be used to select samples from each cow for each period. Merge with metadata.
reads.per.sample <- read.csv("sample-frequency.csv", head=FALSE)
colnames(reads.per.sample) <- c("sample.id", "reads")
metadata.v2 <- merge(metadata, reads.per.sample, by="sample.id", all=FALSE)
dim(metadata.v2) #508 samples
length(unique(metadata.v2$cow))#111 cows

##STEP 3: Reading in the milk fatty acid dataset and edit to prepare
#to merge with other datasets.
fattyacid <- readxl::read_xlsx("Dec 2020 DAPP 1 Fatty acid data.xlsx")
dim(fattyacid)
head(fattyacid)

#The intended column headers are on the second row of the excel file. 
#The "Unique ID" column is what should be used for matching and merging to other datasets. 
#In order to set the correct the column headers set, need to read in the dataset
#again but without the first row being used as the column headers
fattyacid <- readxl::read_xlsx("Dec 2020 DAPP 1 Fatty acid data.xlsx", col_names=FALSE)

#Set the column headers to the values on the second row.
colnames(fattyacid)<-fattyacid[2,]

#Remove the first two rows now.
fattyacid <- fattyacid[-c(1:2),]

#Can also get rid of the "ID", "Cow ID", "Date ID", "Location ID", and 
#"Treatment" columns because they either don't provide unique information, or 
#they are information/metadata that will be gathered from the Kalscheur lab datasets.
fattyacid <- fattyacid[,-c(1, 3:6)]
colnames(fattyacid)


#The Matthew Picklo said that they did some follow up characterization of some 
#of the unknown MFAs in the dataset and said that the following MFAs should be 
#renamed: 1) "Unk 18075 (Putative BCFA 10:0, 8-Me)" should be named "10:1 isomer";
#2) "Unk 21175 (Putative BCFA 12:0, 11-Me)" should be named "12:1 isomer"; 
#3) "14:1 unknown" should be renamed "14:1 isomer"; and 
#4) "Unk_28633 (Putative BCFA 17:0, 15-Me)" should be renamed "17:1 isomer". 
#WARNING: you need to unload tidyverse before running "rename" function, because 
#dplyr is not compatible with tidyverse, at least for this function call. 
fattyacid<-rename(fattyacid, `10:1 isomer`=`Unk 18075 (Putative BCFA 10:0, 8-Me)`)
fattyacid<-rename(fattyacid, `12:1 isomer`=`Unk 21175 (Putative BCFA 12:0, 11-Me)`)
fattyacid<-rename(fattyacid, `14:1 isomer`=`14:1 unknown`)
fattyacid<-rename(fattyacid, `17:1 isomer`=`Unk_28633 (Putative BCFA 17:0, 15-Me)`)

#How many sample IDs are in the "unique ID" column? Are there any NAs in the sample ID columns?
length(fattyacid$`Unique ID`)
sum(is.na(fattyacid$`Unique ID`))
length(grep(pattern="_GF", x=fattyacid$`Unique ID`))

#There are 456 sample IDs out of the 471 entries in Unique ID column. 
#Removing the rows with NA in the Unique ID column
na.rows <- is.na(fattyacid$`Unique ID`)
fattyacid <- fattyacid[!na.rows,]
sum(is.na(fattyacid$`Unique ID`))

#Remove the rows containing repeated column headers, and now all the 456 unique IDs
#are the expected sample ID type.
fattyacid <- fattyacid[fattyacid$`Unique ID`!="Unique ID", ]
dim(fattyacid)
length(grep(pattern="_GF", x=fattyacid$`Unique ID`))

#Making sure that all 456 sample IDs in the "unique ID" column are unique, and 
#not replicated.  
fattyacid.distinct.col1 <-unique(fattyacid$`Unique ID`)
length(fattyacid.distinct.col1)

#Change "Unique ID" to "sample.id" for data merging 
fattyacid <- rename(fattyacid, sample.id=`Unique ID`)

#Removing the destination information ("_GF") from the sample IDs so that they can be matched and merged to other datasets
fattyacid$sample.id <- gsub(pattern="_GF", replacement="", fattyacid$sample.id)

#convert the dates in the sample.id's to D1, D2, D3, D4, D5, and D6 
fattyacid$sample.id <- gsub(pattern="_171125", replacement="-D1", fattyacid$sample.id)
fattyacid$sample.id <- gsub(pattern="_171126", replacement="-D2", fattyacid$sample.id)
fattyacid$sample.id <- gsub(pattern="_180104", replacement="-D3", fattyacid$sample.id)
fattyacid$sample.id <- gsub(pattern="_180105", replacement="-D4", fattyacid$sample.id)
fattyacid$sample.id <- gsub(pattern="_180314", replacement="-D5", fattyacid$sample.id)
fattyacid$sample.id <- gsub(pattern="_180315", replacement="-D6", fattyacid$sample.id)

#Remove the rows with zero values for all fatty acids (i.e. are empty sample place holders) 
fattyacid.v2 <- fattyacid[fattyacid$Total>0,]
dim(fattyacid) #456 samples
dim(fattyacid.v2) #450 samples

#Remove certain fatty acids from the dataset that we don't want to include in analysis
#(e.g. internal standards, contaminants, and fatty acids that are not detected in any samples.)
fattyacid.v3 <- fattyacid.v2
fattyacid.v3<-fattyacid.v3[,-grep(pattern="Internal.Standard", colnames(fattyacid.v3))]
length(grep(pattern="Internal.Standard", colnames(fattyacid.v2)))
length(grep(pattern="Internal.Standard", colnames(fattyacid.v3)))

fattyacid.v3<-fattyacid.v3[,-grep(pattern="contaminant", colnames(fattyacid.v3))]
length(grep(pattern="contaminant", colnames(fattyacid.v2)))
length(grep(pattern="contaminant", colnames(fattyacid.v3)))

fattyacid.v3[,4:81] <- lapply(fattyacid.v3[,4:81], as.numeric)

MFA.sums <- data.frame(colSums(fattyacid.v3[,4:81]))
MFA.sums$MFA.name <- rownames(MFA.sums)
colnames(MFA.sums) <- c("MFAsum", "MFAname")
MFA.zero <- MFA.sums[MFA.sums$MFAsum==0,]
MFA.zero <- MFA.zero$MFAname
fattyacid.v4 <- fattyacid.v3[,-c(which(colnames(fattyacid.v3)%in%MFA.zero))]
dim(fattyacid.v3) #81 columns
dim(fattyacid.v4) #77 columns, so four MFAs were dropped because they weren't present in the samples

#STEP 4: Create new "sample" column with "-B" from duplicate name removed,
#in order to merge metadata with the fatty acid data. 
metadata.v2$sample <- metadata.v2$sample.id
metadata.v2$sample <- gsub(pattern="*-B", replacement="", metadata.v2$sample)
fattyacid.v4$sample <- fattyacid.v4$sample.id

metadata.v3 <- merge(metadata.v2, fattyacid.v4, by="sample", all=FALSE)
dim(metadata.v3) #473 samples

#now get rid of "sample" column, and rename "sample.id.x" to "sample.id"
library(dplyr)
metadata.v3 <- rename(metadata.v3, "sample.id"="sample.id.x")
metadata.v3 <- metadata.v3[,-c(1,18)]

#STEP 5: Rename sample 5849-D1 to 5849-D2 because cow 5849 actually missed
#the milk collection on November 25, 2017, but milk was collected from this cow
#on November 26, 2017. The same was true for cow 6233. 
metadata.v3$sample.id <- gsub(pattern="5849-D1", replacement="5849-D2", x=metadata.v3$sample.id)
metadata.v3$sample.id <- gsub(pattern="6233-D1", replacement="6233-D2", x=metadata.v3$sample.id)
metadata.v3[metadata.v3$sample.id=="5849-D2","day"] <- 2
metadata.v3[metadata.v3$sample.id=="6233-D2","day"] <- 2

#STEP 6: Read in the milk microbiota feature table (which has been rarefied at
#507 reads per sample and filtered to remove antibiotic-treated cows and duplicates). 
feature.table <- read_qza("feature-table-rarefied-at-507-samples-for-analysis.qza")
feature.table <- as.data.frame(feature.table$data)

#STEP 7: Rename 5849-D1 to 5849-D2
feature.table.v2 <- feature.table
library(dplyr)
feature.table.v2 <- rename(feature.table.v2, "5849-D2"="5849-D1")

#STEP 8: Filter the metadata samples down to what is the feature table
metadata.v4 <- metadata.v3[which(metadata.v3$sample.id %in% colnames(feature.table.v2)),]
dim(metadata.v4) #334 samples

#STEP 9: For the cows that have more than one sample from a period remaining in
#the dataset, the sample with the higher number of reads (before rarefying), will
#be kept, and the other removed from the dataset. 
ChooseSampleFromPeriod <-function(cow) {
NoDuplicates <- data.frame()
for(i in cow) {
    select.cow <- metadata.v4[metadata.v4$cow==i,]
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
cows <-unique(metadata.v4$cow)
metadata.v5 <- ChooseSampleFromPeriod(cows) 
length(metadata.v5$sample.id) #191 samples
length(unique(metadata.v5$cow)) #66 cows

#STEP 10: Read in taxonomy
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

#STEP 11: Read in rooted tree. 
tree <- read_qza("rooted-tree-rarefied507-samplesforanalyses.qza")

#STEP 12: build a phyloseq object which will be used for beta-diversity analyses
rownames(metadata.v5)=NULL
phylobj <- phyloseq(otu_table(feature.table.v2, taxa_are_rows=T),tax_table(taxonomy3), phy_tree(tree$data), sample_data(metadata.v5 %>% as.data.frame() %>% column_to_rownames("sample.id")))

#keep only the taxa that have at least one count
phylobj.filtered <- filter_taxa(phylobj, function(x) sum(x)>0, TRUE)

#STEP 13: Calculate beta-diversity of MFAs and microbiota.

#Melt the phyloseq object to obtain a dataframe that can be used with vegan
#to generate dissimilarity matrices
melted.phylobj.ASV <- psmelt(phylobj.filtered)

#Change the data frame from long to wide format
melted.phylobj.ASV.wide <- pivot_wider(data=melted.phylobj.ASV, id_cols=Sample, names_from=OTU, values_from=Abundance)

#merge with the metadata to get back the other metadata measures
#including the MFA abundances.
#First, rename the metadata sample.id column name to Sample
metadata.v6<-rename(metadata.v5, Sample = sample.id)
dim(metadata.v6)

#then merge by "Sample" column name
melted.phylobj.ASV.wide.v2 <- merge(metadata.v6, melted.phylobj.ASV.wide, by="Sample", all=FALSE)
dim(melted.phylobj.ASV.wide.v2) #All 191 retained during merging

#Separate the data frame into periods
period0 <- melted.phylobj.ASV.wide.v2[melted.phylobj.ASV.wide.v2$period=="p_0",]
period1 <- melted.phylobj.ASV.wide.v2[melted.phylobj.ASV.wide.v2$period=="p_1",]
period2 <- melted.phylobj.ASV.wide.v2[melted.phylobj.ASV.wide.v2$period=="p_2",]

#Calculate bray-curtis dissimilarity metrics for MFAs and microbiota 
#in each period. 

period0.MFA <-period0[,19:91] 
period0.ASV <-period0[,93:8525]
period0.MFA <- data.frame(lapply(period0.MFA, as.numeric))
period0.ASV <- data.frame(lapply(period0.ASV, as.numeric))

set.seed(1)
dist.period0.MFA <-vegdist(period0.MFA, method="bray")
dist.period0.ASV <- vegdist(period0.ASV, method="bray")

period1.MFA <-period1[,19:91] 
period1.ASV <-period1[,93:8525]
period1.MFA <- data.frame(lapply(period1.MFA, as.numeric))
period1.ASV <- data.frame(lapply(period1.ASV, as.numeric))

set.seed(1)
dist.period1.MFA <-vegdist(period1.MFA, method="bray")
dist.period1.ASV <- vegdist(period1.ASV, method="bray")

period2.MFA <-period2[,19:91] 
period2.ASV <-period2[,93:8525]
period2.MFA <- data.frame(lapply(period2.MFA, as.numeric))
period2.ASV <- data.frame(lapply(period2.ASV, as.numeric))

set.seed(1)
dist.period2.MFA <-vegdist(period2.MFA, method="bray")
dist.period2.ASV <- vegdist(period2.ASV, method="bray")

#STEP 14: Investigate correlations between milk microbiota beta-diversity 
#and MFA beta-diversity
set.seed(1)
mantel.period0 <- mantel(dist.period0.ASV, dist.period0.MFA, method="spearman", na.rm=FALSE)
mantel.period0 #not significant before multiple test correction. 

set.seed(2)
mantel.period1 <- mantel(dist.period1.ASV, dist.period1.MFA, method="spearman", na.rm=FALSE)
mantel.period1 #not significant

set.seed(3)
mantel.period2 <- mantel(dist.period2.ASV, dist.period2.MFA, method="spearman", na.rm=FALSE)
mantel.period2 #not significant


#STEP 15: Procrustes analysis using the same Bray-Curtis dissimilarity matrices
set.seed(1)
mds.period0.MFA <- monoMDS(dist.period0.MFA, y=cmdscale(dist.period0.MFA))
mds.period0.ASV <- monoMDS(dist.period0.ASV, y=cmdscale(dist.period0.ASV))
protest(X=mds.period0.ASV, Y=mds.period0.MFA, scale=TRUE) #not significant

set.seed(1)
mds.period1.MFA <- monoMDS(dist.period1.MFA, y=cmdscale(dist.period1.MFA))
mds.period1.ASV <- monoMDS(dist.period1.ASV, y=cmdscale(dist.period1.ASV))
protest(X=mds.period1.ASV, Y=mds.period1.MFA, scale=TRUE) #not significant

set.seed(1)
mds.period2.MFA <- monoMDS(dist.period2.MFA, y=cmdscale(dist.period2.MFA))
mds.period2.ASV <- monoMDS(dist.period2.ASV, y=cmdscale(dist.period2.ASV))
protest(X=mds.period2.ASV, Y=mds.period2.MFA, scale=TRUE) #not significant


#STEP 16: Moving on to linear mixed effects modeling and repeated measures 
#correlation analyses, so need to remove cows with less than three samples
#(i.e. cows that don't have a sample from each period.)
cows.pd0 <- unique(metadata.v6[metadata.v6$period=="p_0",][["cow"]])
cows.pd1 <- unique(metadata.v6[metadata.v6$period=="p_1",][["cow"]])
cows.pd2 <- unique(metadata.v6[metadata.v6$period=="p_2",][["cow"]])
cows.to.keep <- cows.pd0[which(cows.pd0 %in% cows.pd1)]
cows.to.keep.v2 <- cows.to.keep[which(cows.to.keep %in% cows.pd2)]
length(cows.to.keep.v2) #60 cows
metadata.v7 <- metadata.v6[which(metadata.v6$cow %in% cows.to.keep.v2),]
dim(metadata.v7) #180 samples
length(unique(metadata.v7$cow)) #60 cows

#Removing MFAs that are not present in any samples in this subset of samples,
#then calculating shannon and simpson MFA alpha diversity in each sample. 
MFA.sums.v2 <- data.frame(colSums(metadata.v7[,19:91]))
MFA.sums.v2$MFA.name <- rownames(MFA.sums.v2)
colnames(MFA.sums.v2) <- c("MFAsum", "MFAname")
MFA.zero.v2 <- MFA.sums.v2[MFA.sums.v2$MFAsum==0,]
MFA.zero.v2 <- MFA.zero.v2$MFAname
metadata.v8 <- metadata.v7[,-c(which(colnames(metadata.v7)%in%MFA.zero.v2))]

#create new phyloseq object. 
metadata.v8 <- rename(metadata.v8, "sample.id"="Sample")
rownames(metadata.v8)=NULL
phylobj.v2 <- phyloseq(otu_table(feature.table.v2, taxa_are_rows=T),tax_table(taxonomy3), phy_tree(tree$data), sample_data(metadata.v8 %>% as.data.frame() %>% column_to_rownames("sample.id")))

#keep only the taxa that have at least one count
phylobj.v2.filtered <- filter_taxa(phylobj.v2, function(x) sum(x)>0, TRUE)


#STEP 17: Looking for a correlation between MFA alpha diversity and microbial alpha diversity.

#calculating alpha diversity
#but first need to unload microbiome package because it somehow interferes
#with vegan package. 
unload("microbiome")

alpha.diversity.MFA<- function(sample) {
    subset <- data.frame()
    for (i in sample) {
        MFA.by.sample <- metadata.v8[metadata.v8$sample.id==i,]
        MFA.by.sample2 <- MFA.by.sample[,19:86]
        MFA.by.sample2 <- as.numeric(MFA.by.sample2)
        MFA.by.sample2[is.na(MFA.by.sample2)] <- 0
        MFA.by.sample$MFA.diversity.inversesimpson<-diversity(x=MFA.by.sample2, index="invsimpson")
        MFA.by.sample$MFA.diversity.shannon<-diversity(x=MFA.by.sample2, index="shannon")
        subset <- rbind(subset, MFA.by.sample) 
    }
    print(subset)
}

metadata.v9 <- alpha.diversity.MFA(metadata.v8$sample.id)

#Calculate microbial phylogenetic alpha diversity
ASV.table <- as.data.frame(otu_table(phylobj.v2.filtered))
phylogenetic.diversity.microbiota <- pd(samp=t(ASV.table), tree=phy_tree(phylobj.v2.filtered), include.root=TRUE)
phylogenetic.diversity.microbiota <- as.data.frame(phylogenetic.diversity.microbiota)
phylogenetic.diversity.microbiota$sample.id <- rownames(phylogenetic.diversity.microbiota)
phylogenetic.diversity.microbiota <- phylogenetic.diversity.microbiota[,-2]
colnames(phylogenetic.diversity.microbiota) <- c("microbiota.diversity.PD", "sample.id")
metadata.v10 <- merge(metadata.v9, phylogenetic.diversity.microbiota, by="sample.id", all=FALSE)

#Calculate the shannon and inverse simpson microbial alpha diversity
alpha.diversity.microbiota <- estimate_richness(physeq=phylobj.v2.filtered, measures=c("InvSimpson", "Shannon")) 
colnames(alpha.diversity.microbiota)<- c("microbiota.diversity.shannon", "microbiota.diversity.inversesimpson")
alpha.diversity.microbiota$sample.id <- rownames(alpha.diversity.microbiota)
metadata.v11 <- merge(metadata.v10, alpha.diversity.microbiota, by="sample.id", all=FALSE)

#Running linear mixed effects modeling to determine if MFA alpha diversity is a 
#predictor (or is correlated with) of microbial alpha diversity. 
hist(metadata.v11$microbiota.diversity.shannon) #somewhat left skewed
hist(metadata.v11$microbiota.diversity.inversesimpson) #heavily right skewed
hist(metadata.v11$microbiota.diversity.PD) #a little right skewed.
hist(metadata.v11$MFA.diversity.shannon) #fairly normal
hist(metadata.v11$MFA.diversity.inversesimpson) #fairly normal


metadata.v11$cow <- as.factor(metadata.v11$cow)
summary(lmer(microbiota.diversity.shannon~MFA.diversity.shannon+(1|cow), data=metadata.v11)) #warning  boundary singular fit, not significant. 
summary(lmer(microbiota.diversity.shannon~MFA.diversity.shannon+(MFA.diversity.shannon|cow), data=metadata.v11)) #warning: model convergence failure,  boundary singular fit.

rmcorr.MFAshannon.microbeshannon <-rmcorr(participant=cow, measure1=MFA.diversity.shannon, measure2=microbiota.diversity.shannon, dataset=metadata.v11)
rmcorr.MFAshannon.microbeshannon #not significant

summary(lmer(microbiota.diversity.PD~MFA.diversity.shannon+(1|cow), data=metadata.v11)) #not significant
summary(lmer(microbiota.diversity.PD~MFA.diversity.shannon+(MFA.diversity.shannon|cow), data=metadata.v11)) #warning  boundary singular fit.

rmcorr.MFAshannon.microbePD<-rmcorr(participant=cow, measure1=MFA.diversity.shannon, measure2=microbiota.diversity.PD, dataset=metadata.v11) #not significant
rmcorr.MFAshannon.microbePD #not significant

summary(lmer(microbiota.diversity.PD~MFA.diversity.inversesimpson+(1|cow), data=metadata.v11)) #not significant
summary(lmer(microbiota.diversity.PD~MFA.diversity.inversesimpson+(MFA.diversity.inversesimpson|cow), data=metadata.v11)) #warning  boundary singular fit.
rmcorr.MFAsimpson.microbePD<-rmcorr(participant=cow, measure1=MFA.diversity.inversesimpson, measure2=microbiota.diversity.PD, dataset = metadata.v11) #not significant
rmcorr.MFAsimpson.microbePD #not significant

#STEP 18: Testing for correlations between single MFAs/total MFAs types and
#microbial alpha diversity, using a function. 

#But first, need to remove MFAs that are present in less than 50 % of the 
#subset of samples because we cannot accurately determine correlation with such 
#a low prevalence MFA. 
Function.FindSparse.MFAs.v2 <- function(MFA) {
    MFA.remove <- data.frame()
    for(i in MFA) {
        if(sum(metadata.v11[[i]]>0)<(180*0.50)){
            MFA.remove <- rbind(MFA.remove, i)
        }
    }
    print(MFA.remove)
}

MFA.list.v2 <- colnames(metadata.v11[,19:87])
MFA.to.remove.v2 <- Function.FindSparse.MFAs.v2(MFA.list.v2)
MFA.to.remove.v2 <- MFA.to.remove.v2[,1]

metadata.v12 <- metadata.v11[,-c(which(colnames(metadata.v11)%in%MFA.to.remove.v2))]

#Remove MFAs that are present in less than 90% of cows
MFA.prevalence.byCow.Function <- function(cow, MFA){
    MFA.keep <- data.frame()
    for(i in cow){
        for(j in MFA){
            cow.subset <-metadata.v12[metadata.v12$cow==i,]
            if(sum(cow.subset[[j]]>0)>0){
                MFA.keep <- rbind(MFA.keep, j)
            }
        }
    }
    print(MFA.keep)
}

cow.list <- unique(metadata$cow)
MFA.list <- colnames(metadata.v12[,19:79])
MFA.occurence <- MFA.prevalence.byCow.Function(cow=cow.list, MFA=MFA.list)
MFA.occurence <-data.frame(table(MFA.occurence))
length(MFA.occurence$MFA.occurence)/sum(MFA.occurence$Freq>=(60*0.9)) #100 percent of the MFAs (that met the prevalence cutoff above) are present in at least 90 %  of the cows. 

#Applying square root transformation of individual MFAs (not total MFAs or MFA types) to resolve the skewed distribution of many of the MFAs.
hist(metadata.v12$Total)
metadata.v13 <- metadata.v12
metadata.v13[,19:79] <- lapply(metadata.v13[,19:79], sqrt)

#running linear mixed effects modeling
Function.MFA.MicrobialShannonDiversity <- function(MFA){
    coefficient.pvalue <- data.frame(message=NA,coefficient=NA, pvalue=NA, MFA=NA)
    for(i in MFA){
        form<-as.formula(paste0("microbiota.diversity.shannon","~`", i, "`+(1|cow)"))
        lmer.output<-try(summary(lmer(form, data=metadata.v13)), silent=TRUE)
        lmer.output2<-try(data.frame(message=paste(lmer.output[["optinfo"]][["conv"]][["lme4"]][["messages"]], collapse=" "),coefficient=lmer.output[["coefficients"]][2,1], pvalue=lmer.output[["coefficients"]][2,5], MFA=i), silent=TRUE)
        coefficient.pvalue <- rbind(coefficient.pvalue, lmer.output2)  
    }
    print(coefficient.pvalue)
}

MFA <- colnames(metadata.v13[,19:80])
MFA.MicrobeShannon.Correlation <- Function.MFA.MicrobialShannonDiversity(MFA)
MFA.MicrobeShannon.Correlation <- MFA.MicrobeShannon.Correlation[-1,] #singular fit warning, not signficant

Function.MFA.MicrobialPDDiversity <- function(MFA){
    coefficient.pvalue <- data.frame(message=NA,coefficient=NA, pvalue=NA, MFA=NA)
    for(i in MFA){
        form<-as.formula(paste0("microbiota.diversity.PD","~`", i, "`+(1|cow)"))
        lmer.output<-try(summary(lmer(form, data=metadata.v13)), silent=TRUE)
        lmer.output2<-try(data.frame(message=paste(lmer.output[["optinfo"]][["conv"]][["lme4"]][["messages"]], collapse=" "), coefficient=lmer.output[["coefficients"]][2,1], pvalue=lmer.output[["coefficients"]][2,5], MFA=i), silent=TRUE)
        coefficient.pvalue <- rbind(coefficient.pvalue, lmer.output2)  
    }
    print(coefficient.pvalue)
}

MFA.MicrobePD.Correlation <- Function.MFA.MicrobialPDDiversity(MFA)
MFA.MicrobePD.Correlation <- MFA.MicrobePD.Correlation[-1,] #not significant


#STEP 19: Make new phyloseq object with the metadata set that has the low
#prevalent MFAs removed and the square root transformed MFAs, then 
#summarize counts by genus. 
phylobj.v3 <- phyloseq(otu_table(feature.table.v2, taxa_are_rows=T),tax_table(taxonomy3), phy_tree(tree$data), sample_data(metadata.v13 %>% as.data.frame() %>% column_to_rownames("sample.id")))

#Keep only the taxa that have at least one count.
phylobj.v3.filtered <- filter_taxa(phylobj.v3, function(x) sum(x)>0, TRUE)

#Combine all of the data at the Genus level

phylobj.genus <- tax_glom(phylobj.v3.filtered, 'Genus', NArm=F)
dim(otu_table(phylobj.genus)) #all 180 samples retained

#Create a list of unique genus names that can be assigned to the features within the phyloseq object
taxa.names <- tax_table(phylobj.genus)[,6]

#Replace empty genus names with family level identities
taxa.names2 <- tax_table(phylobj.genus)[,5] 
taxa.names[c(which(is.na(taxa.names)), grep("g__$", taxa.names))] <-taxa.names2[c(which(is.na(taxa.names)),grep("g__$", taxa.names))]

#Replace empty family names with order level identities
taxa.names3 <-tax_table(phylobj.genus)[,4] 
taxa.names[c(which(is.na(taxa.names2)), grep( "f__$", taxa.names))] <- taxa.names3[c(which(is.na(taxa.names2)), grep("f__$", taxa.names))]

#Replace empty Order names with class level identities
taxa.names4 <- tax_table(phylobj.genus)[,3] 
taxa.names[c(which(is.na(taxa.names3)), grep("o__$", taxa.names))] <-taxa.names4[c(which(is.na(taxa.names3)),grep("o__$", taxa.names))]

#Replace empty Class names with phylum level identities
taxa.names5 <- tax_table(phylobj.genus)[,2] 
taxa.names[c(which(is.na(taxa.names4)),grep("c__$", taxa.names))] <-taxa.names5[c(which(is.na(taxa.names4)),grep("c__$", taxa.names))]

#Replace empty Phylum names with kingdom level identities
taxa.names6 <- tax_table(phylobj.genus)[,1] 
taxa.names[c(which(is.na(taxa.names5)),grep("p__$", taxa.names))] <-taxa.names6[c(which(is.na(taxa.names5)),grep("p__$", taxa.names))]

#Remove leading white space
taxa.names <- gsub(" ","", taxa.names)

#Change the taxa designations to names compatible with R
taxa.names <- make.names(taxa.names, unique = TRUE)

#Make all of the feature names within the phyloseq object correspond to the genus level designation
taxa_names(phylobj.genus) <- taxa.names


## STEP 20: Determine correlations between milk fatty acids and microbial 
#abundance at the genus level.

#Melt genus level phyloseq objects into data frames and change "cow" into factor. 
melted.phylgenus <-psmelt(phylobj.genus)
melted.phylgenus$cow <- as.factor(melted.phylgenus$cow)

#Remove genera that aren't present in at least 50% of samples
prevalence.per.genus <- melted.phylgenus[melted.phylgenus$Abundance>0, ]
prevalence.per.genus <- data.frame(tapply(prevalence.per.genus$Abundance, prevalence.per.genus$OTU, length))
prevalence.per.genus$OTUnames <- rownames(prevalence.per.genus)
colnames(prevalence.per.genus) <- c("prevalence", "OTU")

genera.50percent <-prevalence.per.genus[prevalence.per.genus$prevalence>=(180*0.5),]
genera.50percent <- unique(genera.50percent$OTU) #41 genera

nonzero.counts.per.cow.per.genus <-melted.phylgenus[which(melted.phylgenus$OTU %in%genera.50percent),]
nonzero.counts.per.cow.per.genus <- nonzero.counts.per.cow.per.genus[nonzero.counts.per.cow.per.genus$Abundance>0, ]

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
length(genera.50PercentSamples.90PercentCows) #35 genera

#printing to file the genera tested for relationships with MFAs
#Add the full taxonomy to the list of genera tested for relationships with SCC,
#then write to file 
full.taxonomy.genus <- data.frame(tax_table(phylobj.genus))
full.taxonomy.genus$Genus <- rownames(full.taxonomy.genus)
full.taxonomy.genus <- full.taxonomy.genus[which(full.taxonomy.genus$Genus %in% genera.50PercentSamples.90PercentCows), ]
full.taxonomy.genus$collapsed.taxonomy <- paste(full.taxonomy.genus$Kingdom, full.taxonomy.genus$Phylum, full.taxonomy.genus$Class, full.taxonomy.genus$Order, full.taxonomy.genus$Family, full.taxonomy.genus$Genus, sep="")
full.taxonomy.genus <- full.taxonomy.genus[,c(6,8)]
write.csv(full.taxonomy.genus, "Genera_tested_for_relationships_with_MFAs.csv")

#Functions are created that performs generalized linear mixed effects modeling 
#on each genus and fatty acid combination and adjusts the p-value for multiple comparisons.  
Function.MFApredictor.GenusResponse.GLME.v1 <-function(microbe, MFA) {
    coefficient.pvalue <- data.frame(messages=NA, coefficient=NA, pvalue=NA, taxon=NA, MFA=NA)
    for(i in microbe){
        for(j in MFA){
            microbe.subset <- melted.phylgenus[melted.phylgenus$OTU==i,]
            form<-as.formula(paste0("Abundance~`", j, "`+(1|cow)"))
            glmer.model <- try(glmer(form, family=poisson, data=microbe.subset), silent=TRUE)
            model.summary <- try(summary(glmer.model), silent=TRUE)
            glmer.df <- try(data.frame(messages=paste(model.summary[["optinfo"]][["conv"]][["lme4"]][["messages"]], collapse=" "), coefficient=model.summary[["coefficients"]][2,1], pvalue=model.summary[["coefficients"]][2,4], taxon=i, MFA=j), silent=FALSE)
            coefficient.pvalue <- try(rbind(coefficient.pvalue, glmer.df), silent=TRUE)
        }
    }
    print(coefficient.pvalue)
}

MFA.list.v3 <- colnames(melted.phylgenus[,21:81])
results.MFApredictor.GenusResponse.v1 <-Function.MFApredictor.GenusResponse.GLME.v1(microbe=genera.50PercentSamples.90PercentCows, MFA=MFA.list.v3)
results.MFApredictor.GenusResponse.v1 <- results.MFApredictor.GenusResponse.v1[-1,]

#There were some models that failed to converge so for those models, trying 
#allFit with bobyqa optimizer
Function.MFApredictor.GenusResponse.GLME.v1.optimizer <-function(microbe, MFA) {
    optimizer.choice <- data.frame(method="", optimizer="bobyqa")
    coefficient.pvalue <- data.frame(messages=NA, coefficient=NA, pvalue=NA, taxon=NA, MFA=NA)
    for(i in microbe){
        for(j in MFA){
            microbe.subset <- melted.phylgenus[melted.phylgenus$OTU==i,]
            form<-as.formula(paste0("Abundance~`", j, "`+(1|cow)"))
            glmer.model <- try(glmer(form, data=microbe.subset, family=poisson), silent=TRUE)
            glmer.optimization<-try(allFit(glmer.model, meth.tab=optimizer.choice), silent = TRUE)
            
            glmer.bobyqa <- try(summary(glmer.optimization[["bobyqa"]]), silent=TRUE)
            glmer.bobyqa.df <- try(data.frame(messages=paste(glmer.bobyqa[["optinfo"]][["conv"]][["lme4"]][["messages"]], collapse=" "), coefficient=glmer.bobyqa[["coefficients"]][2,1], pvalue=glmer.bobyqa[["coefficients"]][2,4], taxon=i, MFA=j), silent=FALSE)
            coefficient.pvalue <- try(rbind(coefficient.pvalue, glmer.bobyqa.df), silent=TRUE)
        }
    }
    print(coefficient.pvalue)
}

failed.convergence <- results.MFApredictor.GenusResponse.v1
failed.convergence<-failed.convergence[grep(pattern="Model*", x=failed.convergence$messages),]
failed.taxon <- unique(failed.convergence$taxon)
failed.mfa <- unique(failed.convergence$MFA)

results.GLME.MFA.genera.v1.optimizer<- Function.MFApredictor.GenusResponse.GLME.v1.optimizer(microbe=failed.taxon, MFA=failed.mfa)
results.GLME.MFA.genera.v1.optimizer <- results.GLME.MFA.genera.v1.optimizer[-1,] #remove the first row because it's empty

#some of the models still failed with the optimizer, so removing those
results.GLME.MFA.genera.v1.optimizer.v2 <- results.GLME.MFA.genera.v1.optimizer[-(grep(pattern="Model", results.GLME.MFA.genera.v1.optimizer$messages)),]

#all failed with bobyqa optimizer. Removing those models from the dataset, and 
#calculating fdr
results.MFApredictor.GenusResponse.v1.b <- results.MFApredictor.GenusResponse.v1[-grep(pattern="Model*", x=results.MFApredictor.GenusResponse.v1$messages),]
results.MFApredictor.GenusResponse.v1.c <-rbind(results.MFApredictor.GenusResponse.v1.b, results.GLME.MFA.genera.v1.optimizer.v2)

results.MFApredictor.GenusResponse.v1.c$fdr <- p.adjust(p=results.MFApredictor.GenusResponse.v1.c$pvalue, method="fdr")
results.MFApredictor.GenusResponse.v1.c$bonferroni <- p.adjust(p=results.MFApredictor.GenusResponse.v1.c$pvalue, method="bonferroni")
write.csv(results.MFApredictor.GenusResponse.v1.c, file="MFApredictor.GENUSresponse.GlmerPoisson.CowIntercept.csv")

#STEP 21: Summarize results (significant relationships) by category of MFA
significant.subset <- results.MFApredictor.GenusResponse.v1.c[results.MFApredictor.GenusResponse.v1.c$bonferroni<=0.01,]
dim(significant.subset) #925 significant results. 


#For each genus, count the number of relationship with 1) SFA - SCFA 
#2) SFA – non-SCFA and non-BCFA, 3) SFA - BCFA, 4) MUFA, and 5) PUFA. 
length(unique(significant.subset$MFA)) #61 total MFAs

SFA.SCFA <- unique(grep(pattern="X4\\.0|X6\\.0",significant.subset$MFA, value=TRUE))
SFA.SCFA
length(SFA.SCFA) #2 total

SFA.BCFA <- unique(grep(pattern="Me+",significant.subset$MFA, value=TRUE))
SFA.BCFA
length(SFA.BCFA) #7 total

SFA.EvenChain <- unique(grep(pattern="X8\\.0$|X[1-2][0,2,4,6,8]\\.0$",significant.subset$MFA, value=TRUE))
SFA.EvenChain
length(SFA.EvenChain) #9 total

SFA.OddChain <- unique(grep(pattern="X[7,9]\\.0$|X[1-2][1,3,5,7,9]\\.0$",significant.subset$MFA, value=TRUE))
SFA.OddChain
length(SFA.OddChain) #8 total

MUFA <- unique(grep(pattern="X[0-9]\\.1|X[0-2][0-9]\\.1", significant.subset$MFA, value=TRUE))
MUFA
length(MUFA) #17 total

PUFA <- unique(grep(pattern="X[0-9]\\.[2-9]|X[0-2][0-9]\\.[2-9]", significant.subset$MFA, value=TRUE))
PUFA
length(PUFA) #18 total

sum(length(SFA.SCFA)+length(SFA.BCFA)+length(SFA.EvenChain)+length(SFA.OddChain)+length(MUFA)+length(PUFA)) #61 total
#any MFA duplicated in the lists/categories?
length(unique(c(SFA.SCFA, SFA.BCFA, SFA.EvenChain, SFA.OddChain, MUFA, PUFA))) #61

Function.Summarize.byMFAcategory <- function(genus){
    df <- data.frame(Genus=NA, MFA_category=NA, number=NA)
    for(i in genus){
        taxon.subset <-significant.subset[significant.subset$taxon==i,]
        
        taxon.SFA.SCFA.subset <- taxon.subset[which(taxon.subset$MFA %in% SFA.SCFA),]
        taxon.SFA.SCFA.subset.positiveeffect <- taxon.SFA.SCFA.subset[taxon.SFA.SCFA.subset$coefficient>0,]
        taxon.SFA.SCFA.subset.negativeeffect <- taxon.SFA.SCFA.subset[taxon.SFA.SCFA.subset$coefficient<0,]
        df.1a <- data.frame(Genus=i, MFA_category="SFA: SCFA", number=nrow(taxon.SFA.SCFA.subset.positiveeffect))
        df.1b <- data.frame(Genus=i, MFA_category="SFA: SCFA", number=-(nrow(taxon.SFA.SCFA.subset.negativeeffect)))
        df <- rbind(df, df.1a)
        df <- rbind(df, df.1b)
        
        taxon.SFA.BCFA.subset <- taxon.subset[which(taxon.subset$MFA %in% SFA.BCFA),]
        taxon.SFA.BCFA.subset.positiveeffect <- taxon.SFA.BCFA.subset[taxon.SFA.BCFA.subset$coefficient>0,]
        taxon.SFA.BCFA.subset.negativeeffect <- taxon.SFA.BCFA.subset[taxon.SFA.BCFA.subset$coefficient<0,]
        df.2a <- data.frame(Genus=i, MFA_category="SFA: BCFA", number=nrow(taxon.SFA.BCFA.subset.positiveeffect))
        df.2b <- data.frame(Genus=i, MFA_category="SFA: BCFA", number=-(nrow(taxon.SFA.BCFA.subset.negativeeffect)))
        df <- rbind(df, df.2a)
        df <- rbind(df, df.2b)
        
        taxon.SFA.EvenChain <- taxon.subset[which(taxon.subset$MFA %in% SFA.EvenChain),]
        taxon.SFA.EvenChain.subset.positiveeffect <- taxon.SFA.EvenChain[taxon.SFA.EvenChain$coefficient>0,]
        taxon.SFA.EvenChain.subset.negativeeffect <- taxon.SFA.EvenChain[taxon.SFA.EvenChain$coefficient<0,]
        df.3a <- data.frame(Genus=i, MFA_category="SFA: even chain", number=nrow(taxon.SFA.EvenChain.subset.positiveeffect))
        df.3b <- data.frame(Genus=i, MFA_category="SFA: even chain", number=-(nrow(taxon.SFA.EvenChain.subset.negativeeffect)))
        df <- rbind(df, df.3a)
        df <- rbind(df, df.3b)
        
        taxon.SFA.OddChain <- taxon.subset[which(taxon.subset$MFA %in% SFA.OddChain),]
        taxon.SFA.OddChain.subset.positiveeffect <- taxon.SFA.OddChain[taxon.SFA.OddChain$coefficient>0,]
        taxon.SFA.OddChain.subset.negativeeffect <- taxon.SFA.OddChain[taxon.SFA.OddChain$coefficient<0,]
        df.4a <- data.frame(Genus=i, MFA_category="SFA: odd chain", number=nrow(taxon.SFA.OddChain.subset.positiveeffect))
        df.4b <- data.frame(Genus=i, MFA_category="SFA: odd chain", number=-(nrow(taxon.SFA.OddChain.subset.negativeeffect)))
        df <- rbind(df, df.4a)
        df <- rbind(df, df.4b)
       
        taxon.MUFA <- taxon.subset[which(taxon.subset$MFA %in% MUFA),]
        taxon.MUFA.positiveeffect <- taxon.MUFA[taxon.MUFA$coefficient>0,]
        taxon.MUFA.negativeeffect <- taxon.MUFA[taxon.MUFA$coefficient<0,]
        df.5a <- data.frame(Genus=i, MFA_category="MUFA", number=nrow(taxon.MUFA.positiveeffect))
        df.5b <- data.frame(Genus=i, MFA_category="MUFA", number=-(nrow(taxon.MUFA.negativeeffect)))
        df <- rbind(df, df.5a)
        df <- rbind(df, df.5b) 
        
        taxon.PUFA <- taxon.subset[which(taxon.subset$MFA %in% PUFA),]
        taxon.PUFA.positiveeffect <- taxon.PUFA[taxon.PUFA$coefficient>0,]
        taxon.PUFA.negativeeffect <- taxon.PUFA[taxon.PUFA$coefficient<0,]
        df.6a <- data.frame(Genus=i, MFA_category="PUFA", number=nrow(taxon.PUFA.positiveeffect))
        df.6b <- data.frame(Genus=i, MFA_category="PUFA", number=-(nrow(taxon.PUFA.negativeeffect)))
        df <- rbind(df, df.6a)
        df <- rbind(df, df.6b)
    }
    print(df)
} 

summary.results.byMFA.category <- Function.Summarize.byMFAcategory(unique(significant.subset$taxon))
summary.results.byMFA.category <- summary.results.byMFA.category[!is.na(summary.results.byMFA.category$Genus),]
dim(summary.results.byMFA.category) #396 rows

#add the name of the Class for each genus
full.taxonomy.genus <- data.frame(tax_table(phylobj.genus))
full.taxonomy.genus$Genus <- rownames(full.taxonomy.genus)
select.taxonomy.genus.MFA <- full.taxonomy.genus[which(full.taxonomy.genus$Genus %in% summary.results.byMFA.category$Genus), ]
summary.results.byMFA.category.v2 <-merge(summary.results.byMFA.category,select.taxonomy.genus.MFA, by="Genus")

#Before graphing, get rid of "g__" before each genus name, "f__" before each family,
#and "c__" before the class name.
summary.results.byMFA.category.v2$Genus <- gsub(pattern="g__", replacement="", summary.results.byMFA.category.v2$Genus)
summary.results.byMFA.category.v2$Genus <- gsub(pattern="f__", replacement="", summary.results.byMFA.category.v2$Genus)
summary.results.byMFA.category.v2$Class <- gsub(pattern="c__", replacement="", summary.results.byMFA.category.v2$Class)

ggplot(aes(x=number, y=Genus,fill=MFA_category), data=summary.results.byMFA.category.v2)+geom_bar(position="stack", stat="identity")+geom_vline(xintercept=0)+xlab("number of MFAs associated with a genus")+labs(fill="fatty acid type")+labs(title = "negative associations               positive associations")+theme(plot.title = element_text(lineheight = 0.5, size=10, hjust=0.45))+facet_grid(vars(Class), scales="free_y", space="free_y")+theme(strip.text.y=element_text(size=8,angle=0.45))
ggsave(filename="Barplot.Genus.MFA.correlations.jpeg", width=10, height=5)


#STEP 22: graphing some of the statistically significant relationships from 
#models that had large estimates for beta coefficient generate any warnings.

acinetobacter <- melted.phylgenus[melted.phylgenus$OTU=="g__Acinetobacter",]
model.acinetobacter <- glmer(Abundance~X16.1.E.mixture+(1|cow),family=poisson, data=acinetobacter)
ggplot(aes(x=X16.1.E.mixture, y=Abundance), data=acinetobacter)+geom_point()+geom_line(aes(y=predict(model.acinetobacter), group=cow))+ylab("Acinetobacter")+xlab("square root 16:1 E mixture")

ggsave(filename="MFA16.1.E.vs.Acinetobacter.jpeg")

hist(residuals(model.acinetobacter)) #pretty normal distribution
median(residuals(model.acinetobacter)) #kind of close to zero 

enterococcus <- melted.phylgenus[melted.phylgenus$OTU=="g__Enterococcus",]
model.enterococcus <- glmer(Abundance~X21.0+(1|cow),family=poisson, data=enterococcus)
ggplot(aes(x=X21.0, y=Abundance), data=enterococcus)+geom_point()+geom_line(aes(y=predict(model.enterococcus), group=cow))+ylab("Enterococcus")+xlab("square root 21:0")

ggsave(filename="MFA21.0.vs.Enterococcus.jpeg")

hist(residuals(model.enterococcus)) 
median(residuals(model.enterococcus)) #kind of close to zero 


muribaculaceae <- melted.phylgenus[melted.phylgenus$OTU=="g__Muribaculaceae",]
model.muribaculaceae <- glmer(Abundance~X20.3n.3+(1|cow),family=poisson, data=muribaculaceae)
ggplot(aes(x=X20.3n.3, y=Abundance), data=muribaculaceae)+geom_point()+geom_line(aes(y=predict(model.muribaculaceae), group=cow))+ylab("Muribaculaceae")+xlab("square root 20:3,n-3")

ggsave(filename="MFA20.3.n3.vs.Muribaculaceae.jpeg")

hist(residuals(model.muribaculaceae)) 
median(residuals(model.muribaculaceae)) #close to zero 

#STEP 23: Since some of these genera seemed to have all negative or all positive
#associations with the individual milk fatty acid types, we think those same 
#genera may have a similar association with the total fatty acid measure.

Function.TotalMFApredictor.GenusResponse.GLME <-function(microbe) {
    coefficient.pvalue <- data.frame(messages=NA, coefficient=NA, pvalue=NA, taxon=NA)
    for(i in microbe){
            microbe.subset <- melted.phylgenus[melted.phylgenus$OTU==i,]
            glmer.model <- try(glmer(Abundance~Total+(1|cow), family=poisson, data=microbe.subset), silent=TRUE)
            model.summary <- try(summary(glmer.model), silent=TRUE)
            glmer.df <- try(data.frame(messages=paste(model.summary[["optinfo"]][["conv"]][["lme4"]][["messages"]], collapse=" "), coefficient=model.summary[["coefficients"]][2,1], pvalue=model.summary[["coefficients"]][2,4], taxon=i), silent=FALSE)
            coefficient.pvalue <- try(rbind(coefficient.pvalue, glmer.df), silent=TRUE)
    }
    print(coefficient.pvalue)
}

significant.genera <- unique(summary.results.byMFA.category.v2$Genus)
results.TotalMFApredictor.GenusResponse <-Function.TotalMFApredictor.GenusResponse.GLME(microbe=significant.genera)
results.TotalMFApredictor.GenusResponse <- results.TotalMFApredictor.GenusResponse[-1,] #one model failed
#re-run the failed model with an optimizer
Function.TotalMFApredictor.GenusResponse.GLME.optimizer <-function(microbe) {
    optimizer.choice <- data.frame(method="", optimizer="bobyqa")
    coefficient.pvalue <- data.frame(messages=NA, coefficient=NA, pvalue=NA, taxon=NA)
    for(i in microbe){
            microbe.subset <- melted.phylgenus[melted.phylgenus$OTU==i,]
            form<-as.formula(paste0("Abundance~Total+(1|cow)"))
            glmer.model <- try(glmer(form, data=microbe.subset, family=poisson), silent=TRUE)
            glmer.optimization<-try(allFit(glmer.model, meth.tab=optimizer.choice), silent = TRUE)
            
            glmer.bobyqa <- try(summary(glmer.optimization[["bobyqa"]]), silent=TRUE)
            glmer.bobyqa.df <- try(data.frame(messages=paste(glmer.bobyqa[["optinfo"]][["conv"]][["lme4"]][["messages"]], collapse=" "), coefficient=glmer.bobyqa[["coefficients"]][2,1], pvalue=glmer.bobyqa[["coefficients"]][2,4], taxon=i), silent=FALSE)
            coefficient.pvalue <- try(rbind(coefficient.pvalue, glmer.bobyqa.df), silent=TRUE)
    }
    print(coefficient.pvalue)
}

model.TotalMFA.UCG10.omptimzed<-Function.TotalMFApredictor.GenusResponse.GLME.optimizer(microbe="g__UCG.010")
#still failed, so remove from list and adjust pvalue for multiple tests
results.TotalMFApredictor.GenusResponse.v2 <- results.TotalMFApredictor.GenusResponse[-c(grep(pattern="Model",results.TotalMFApredictor.GenusResponse$messages)),]
results.TotalMFApredictor.GenusResponse.v2$fdr <- p.adjust(results.TotalMFApredictor.GenusResponse.v2$pvalue, method="fdr")
results.TotalMFApredictor.GenusResponse.v2$bonferroni <- p.adjust(results.TotalMFApredictor.GenusResponse.v2$pvalue, method="bonferroni")

#STEP 24: make graph for the associations between Total MFAs and genera
results.TotalMFApredictor.GenusResponse.v2 <- rename(results.TotalMFApredictor.GenusResponse.v2, "Genus"="taxon")
results.TotalMFApredictor.GenusResponse.v3 <-merge(results.TotalMFApredictor.GenusResponse.v2,select.taxonomy.genus.MFA, by="Genus")

#Before graphing, get rid of "g__" before each genus name, "f__" before each family,
#and "c__" before the class name.
results.TotalMFApredictor.GenusResponse.v3$Genus <- gsub(pattern="g__", replacement="", results.TotalMFApredictor.GenusResponse.v3$Genus)
results.TotalMFApredictor.GenusResponse.v3$Genus <- gsub(pattern="f__", replacement="", results.TotalMFApredictor.GenusResponse.v3$Genus)
results.TotalMFApredictor.GenusResponse.v3$Class <- gsub(pattern="c__", replacement="", results.TotalMFApredictor.GenusResponse.v3$Class)

ggplot(aes(x=coefficient, y=Genus), data=results.TotalMFApredictor.GenusResponse.v3)+geom_bar(stat="identity")+geom_vline(xintercept=0)+xlab("beta coefficient estimate for [total MFA]")+facet_grid(vars(Class), scales="free_y", space="free_y")+theme(strip.text.y=element_text(size=8,angle=0.45))
ggsave(filename="Barplot.TotalMFA.betacoefficient.Genus.jpeg", width=10, height=5)
