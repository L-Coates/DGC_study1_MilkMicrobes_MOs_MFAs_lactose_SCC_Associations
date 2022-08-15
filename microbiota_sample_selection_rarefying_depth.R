#This scripts describes the steps involved in selecting samples and
#determining rarefaction depth (15th percentile)

#Load libraries
library(devtools)
library(tidyverse)
#devtools::install_github("jbisanz/qiime2R")
library(qiime2R)
library(phyloseq)
#install_github("microbiome/microbiome")
library(microbiome)
library(dplyr)
library(plyr)
library(ggplot2)

#Import list of reads per sample (which is a file downloaded from
#the feature table summary) from the feature table that has had: 
#(1) contaminants removed (2) Bos taurus mitochondrion ASV removed
#(3) mitochondria, chloroplast, and eukarya ASVs removed. 
read.totals <- read.csv("sample-frequency.csv", head=FALSE)
colnames(read.totals) <- c("sample.id", "reads")

#STEP 1: Remove commercial raw milk
read.totals.v2 <- read.totals[-c(grep(pattern="CR", x=read.totals$sample.id)),]

#STEP 2: Remove samples with fewer reads than blanks
blanks <- read.totals[c(grep(pattern="BL",x=read.totals$sample.id)),]
summary(blanks$reads) #max number of reads among blanks is 49. 
read.totals.v3 <- read.totals.v2[read.totals.v2$reads>49,]

#STEP 3: Remove antibiotic-treated cows
antibiotic.cows <- c("5020", "5212", "5297", "5405", "5697", "6076", "6215", "6232", "6233", "6241")
read.totals.v3$cow <- gsub(pattern="-D[123456]|-B", replacement="", x=read.totals.v3$sample.id)
length(unique(read.totals.v3$cow)) #76 cows
read.totals.v4 <- read.totals.v3[-c(which(read.totals.v3$cow %in% antibiotic.cows)),]
length(unique(read.totals.v4$cow)) #66 cows

# STEP 4: Select duplicate with higher number of reads
read.totals.v5 <- read.totals.v4
read.totals.v5$sample <- read.totals.v5$sample.id
read.totals.v5$sample <- gsub(pattern="-B", replacement="", x=read.totals.v5$sample)

ChooseDuplicate <-function(sample) {
    NoDuplicates <- data.frame()
    for(i in sample) {
        SubsetBySample<-read.totals.v5[read.totals.v5$sample==i,]
        if(length(SubsetBySample$sample)>1) {
            SubsetBySample<-SubsetBySample[order(SubsetBySample$reads, decreasing=TRUE),]
            OneSample <- SubsetBySample[match(unique(SubsetBySample$sample), SubsetBySample$sample),]
            NoDuplicates <- rbind(NoDuplicates, OneSample)
        } else {
            NoDuplicates <- rbind(NoDuplicates, SubsetBySample)
        }
    }
    return(NoDuplicates)
}

#Run the function "ChooseDuplicate"
sample <- unique(read.totals.v5$sample)
read.totals.v5 <- ChooseDuplicate(sample) 
dim(read.totals.v5) #393 samples

#STEP 5: Remove samples below the 15th percentile
#Calculate 15th percentile
quantile(read.totals.v5$reads, 0.15) #506.4 reads is the 15th percentile

#Remove samples that have fewer than 506.4 reads
read.totals.v6 <- read.totals.v5[read.totals.v5$reads>=507,]
dim(read.totals.v6) #334 samples

#Note: I am still keeping cows that may be missing a sample from a period because
#I would still like to include those cows in the Beta-diversity analyses which I 
#am performing within each period, separately. 

#STEP 6: write list of sample ids to text file to use in qiime to filter the 
#feature table, the representative sequences, and to build a phylogenetic tree
sample.list <- data.frame(read.totals.v6$sample.id)
colnames(sample.list) <- "sample-id"
rownames(sample.list)=NULL
write.table(sample.list, quote=FALSE, row.names=FALSE, file="samples.for.microbiota.analyses.second.manuscript.txt")
