#This script includes the steps taken to determine relationships between 
#milk somatic cell counts (SCC) and milk fatty acids. 

## Load libraries that are needed.
library(dplyr)
library(lme4)
library(lmerTest)
library(MASS)
library(vegan)
library(ggplot2)
library(rmcorr)

## Step 1: Read in the milk fatty acid (MFA) and somatic cell count (SCC) datasets and match them. 
fattyacid <- readxl::read_xlsx("Dec 2020 DAPP 1 Fatty acid data.xlsx", col_names=FALSE)

### Set the column headers to the values on the second row.
colnames(fattyacid)<-fattyacid[2,]

### Remove the first two rows now.
fattyacid <- fattyacid[-c(1:2),]

### Can also get rid of the "ID" and "Location ID" columns since they aren't needed.
fattyacid <- fattyacid[,-c(1,5)]
colnames(fattyacid)

### The Matthew Picklo said that they did some follow up characterization of some of the unknown MFAs in the dataset and said that the following MFAs should be renamed: 1) "Unk 18075 (Putative BCFA 10:0, 8-Me)" should be named "10:1 isomer"; 2) "Unk 21175 (Putative BCFA 12:0, 11-Me)" should be named "12:1 isomer"; 3) "14:1 unknown" should be renamed "14:1 isomer"; and 4) "Unk_28633 (Putative BCFA 17:0, 15-Me)" should be renamed "17:1 isomer". 
fattyacid<-rename(fattyacid, `10:1 isomer`=`Unk 18075 (Putative BCFA 10:0, 8-Me)`)
fattyacid<-rename(fattyacid, `12:1 isomer`=`Unk 21175 (Putative BCFA 12:0, 11-Me)`)
fattyacid<-rename(fattyacid, `14:1 isomer`=`14:1 unknown`)
fattyacid<-rename(fattyacid, `17:1 isomer`=`Unk_28633 (Putative BCFA 17:0, 15-Me)`)

### How many sample IDs are in the "unique ID" column? Are there any NAs in the sample ID columns?
length(fattyacid$`Unique ID`)
sum(is.na(fattyacid$`Unique ID`))
length(grep(pattern="_GF", x=fattyacid$`Unique ID`))

### There are 456 sample IDs out of the 471 entries in Unique ID column. Removing the rows with NA in the Unique ID column
na.rows <- is.na(fattyacid$`Unique ID`)
fattyacid <- fattyacid[!na.rows,]
sum(is.na(fattyacid$`Unique ID`))

### Remove the rows containing repeated column headers, and now all the 456 unique IDs are the expected sample ID type.
fattyacid <- fattyacid[fattyacid$`Unique ID`!="Unique ID", ]
dim(fattyacid)
length(grep(pattern="_GF", x=fattyacid$`Unique ID`)) 

### Making sure that all 456 sample IDs in the "unique ID" column are unique, and not replicated.  
fattyacid.distinct.col1 <-unique(fattyacid$`Unique ID`)
length(fattyacid.distinct.col1)

### Change "Unique ID" to "sample.id" for data merging 
fattyacid <- rename(fattyacid, sample.id=`Unique ID`)

### Removing the destination information ("_GF") from the sample IDs since that bit of information isn't needed in this analysis. 
fattyacid$sample.id <- gsub(pattern="_GF", replacement="", fattyacid$sample.id)

### Removing fatty acids that are "internal standards" or "contaminants"
fattyacid <- fattyacid[,-c(grep("*contaminant|*Internal Standard", colnames(fattyacid)))]

### Removing any samples that have zeros for all fatty acids, or in other words, have zero for "Total" milk fatty acids. 
dim(fattyacid)
fattyacid$Total <- as.numeric(fattyacid$Total)
fattyacid <- fattyacid[fattyacid$Total!=0,]
dim(fattyacid)

### Cows 5849 and 6233 were missed on the first collection timepoint (11/25/2017)
#But milk samples were collected from them on 11/26/2017, so the date in the
#sample id is incorrect and will be changed to "171126"
fattyacid$sample.id <- gsub(pattern="5849_171125", replacement="5849_171126", x=fattyacid$sample.id)
fattyacid$sample.id <- gsub(pattern="6233_171125", replacement="6233_171126", x=fattyacid$sample.id)
fattyacid[fattyacid$sample.id=="5849_171126","Date ID"] = "171126"
fattyacid[fattyacid$sample.id=="6233_171126","Date ID"] = "171126"

### Reading in the AgSource data for the SCC measures from 11/21/17, 11/28/17, 01/02/18, 01/05/18, 03/13/18, and 03/19/18 (which are the same agsource SCC measures that were used for the analysis of SCC & microbiota data.)
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


### Selecting only the M1 measurements (these are from the first milking of the day) from Nov 21st (this is the most complete dataset taken closest to the date of collection of the milk samples used for microbiota analysis.)
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


### Matching SCC measurements to MFA measurements by period. 
fattyacid<-rename(fattyacid,cow=`Cow ID`)
fattyacid.pd0 <- rbind(fattyacid[fattyacid$`Date ID`=="171125",], fattyacid[fattyacid$`Date ID`=="171126",])
fattyacid.pd1 <- rbind(fattyacid[fattyacid$`Date ID`=="180104",], fattyacid[fattyacid$`Date ID`=="180105",])
fattyacid.pd2 <- rbind(fattyacid[fattyacid$`Date ID`=="180314",], fattyacid[fattyacid$`Date ID`=="180315",])

SCC.pd0 <- merge(Nov.21.17, Nov.28.17, by="cow", all=FALSE)
fattyacid.pd0 <- merge(fattyacid.pd0, SCC.pd0, by="cow", all=FALSE)

SCC.pd1 <- merge(Jan.02.18, Jan.05.18, by="cow", all=FALSE)
fattyacid.pd1 <- merge(fattyacid.pd1, SCC.pd1, by="cow", all=FALSE)

SCC.pd2 <- merge(March.13.18, March.19.18, by="cow", all=FALSE)
fattyacid.pd2 <- merge(fattyacid.pd2, SCC.pd2, by="cow", all=FALSE)

fattyacid.v2 <- rbind(fattyacid.pd0, fattyacid.pd1, fattyacid.pd2)

## Step 2: Determine if there is a correlation between MFA alpha diversity and SCC.
fattyacid.v2[,7:84] <- lapply(fattyacid.v2[,7:84], as.numeric)

### Change the column names of the MFAs so that they can be identified in making a formula for lmer(). 
colnames(fattyacid.v2)<-gsub(pattern="-", replacement=".", colnames(fattyacid.v2))
colnames(fattyacid.v2)<-gsub(pattern=":", replacement=".", colnames(fattyacid.v2))
colnames(fattyacid.v2)<-gsub(pattern=" ", replacement=".", colnames(fattyacid.v2))
colnames(fattyacid.v2)<-gsub(pattern="  ", replacement=".", colnames(fattyacid.v2))
colnames(fattyacid.v2)<-gsub(pattern=",", replacement=".", colnames(fattyacid.v2))
colnames(fattyacid.v2)<-gsub(pattern="/", replacement=".", colnames(fattyacid.v2))
colnames(fattyacid.v2)<-gsub(pattern="[()]", replacement=".", colnames(fattyacid.v2))

colnames(fattyacid.v2)[7:83]<-paste0("X", colnames(fattyacid.v2[7:83]))


### Combining MFA abundances by MFA type/class (e.g. saturated fatty acid, monounsaturated fatty acid, polyunsaturated fatty acid, omega 3 fatty acids, omega 6 fatty acids). 
SFA <- grep(pattern="^..\\.0|^...\\.0|BCFA",colnames(fattyacid.v2[,7:83]), value=TRUE)

MUFA <- grep(pattern="^..\\.1|^...\\.1", colnames(fattyacid.v2[,7:83]), value=TRUE)

PUFA <- grep(pattern="^..\\.[2-9]|^...\\.[2-9]", colnames(fattyacid.v2[,7:83]), value=TRUE)

omega3 <- grep(pattern="n.3", colnames(fattyacid.v2[,7:83]), value=TRUE)
omega6 <- grep(pattern="n.6", colnames(fattyacid.v2[,7:83]), value=TRUE)
BCFA <- grep(pattern="Me", colnames(fattyacid.v2[,7:83]), value=TRUE)

sum(length(SFA)+length(MUFA)+length(PUFA))
table(SFA%in%MUFA)
table(MUFA%in%PUFA)
table(SFA%in%PUFA)

fattyacid.v2$SFA<-rowSums(fattyacid.v2[,c(SFA)])
fattyacid.v2$MUFA <- rowSums(fattyacid.v2[,c(MUFA)])
fattyacid.v2$PUFA <- rowSums(fattyacid.v2[,c(PUFA)])
fattyacid.v2$omega3 <- rowSums(fattyacid.v2[,c(omega3)])
fattyacid.v2$omega6 <- rowSums(fattyacid.v2[,c(omega6)])
fattyacid.v2$BCFA <- rowSums(fattyacid.v2[,c(BCFA)])


### Log transforming SCC in the metadata object
fattyacid.v2$LogSCC.milking1.beforesample <- log10(1000*(fattyacid.v2$SCCx1000.milking1.beforesample))
fattyacid.v2$LogSCC.milking1.aftersample <- log10(1000*(fattyacid.v2$SCCx1000.milking1.aftersample))

###Remove the antibiotic-treated cows
antibiotic.cows <- c("5020", "5212","5297", "5405", "5697", "6076", "6215", "6232", "6233", "6241")
length(unique(fattyacid.v2$cow)) #76 cows
fattyacid.v2 <- fattyacid.v2[-c(which(fattyacid.v2$cow %in% antibiotic.cows)),]
length(unique(fattyacid.v2$cow)) #66 cows

### Separating the fatty acid data into two datasets -- one with the 
#first sampling date in the period, and the other with the second sampling 
#date in the period, since the SCC measure is the same for both samples from the
#same period.
unique(fattyacid.v2$Date.ID)
first.sample.dates <- c("171125", "180104", "180314")
second.sample.dates <- c("171126", "180105", "180315")
fattyacid.v2.firstsample <- fattyacid.v2[which(fattyacid.v2$Date.ID%in%first.sample.dates),]
fattyacid.v2.secondsample <- fattyacid.v2[which(fattyacid.v2$Date.ID%in%second.sample.dates),]

### Retain cows in the dataset that have a sample for each period. 
cow.counts.firstsample<-data.frame(table(fattyacid.v2.firstsample$cow))
toRemove <- cow.counts.firstsample[cow.counts.firstsample$Freq<3,]
toRemove$Var1 <- as.character(toRemove$Var1)
toRemove <- toRemove$Var1
length(unique(fattyacid.v2.firstsample$cow))
fattyacid.v2.firstsample.v2 <- fattyacid.v2.firstsample[!(fattyacid.v2.firstsample$cow%in%toRemove),]
length(unique(fattyacid.v2.firstsample.v2$cow)) #one was dropped.

cow.counts.secondsample<-data.frame(table(fattyacid.v2.secondsample$cow))
toRemove <- cow.counts.secondsample[cow.counts.secondsample$Freq<3,]
toRemove$Var1 <- as.character(toRemove$Var1)
toRemove <- toRemove$Var1
length(unique(fattyacid.v2.secondsample$cow)) 
fattyacid.v2.secondsample.v2 <- fattyacid.v2.secondsample[!(fattyacid.v2.secondsample$cow%in%toRemove),]
length(unique(fattyacid.v2.secondsample.v2$cow)) #none were dropped

### Remove MFAs that are absent in all samples
fattyacid.totals.firstsample <- data.frame(colSums(fattyacid.v2.firstsample.v2[,7:84]))
fattyacid.totals.firstsample$fattyacid <-rownames(fattyacid.totals.firstsample)
colnames(fattyacid.totals.firstsample) <- c("column.total", "fattyacid")
MFAs.remove.from.firstsample <- fattyacid.totals.firstsample[fattyacid.totals.firstsample$column.total==0,]
MFAs.remove.from.firstsample <- MFAs.remove.from.firstsample$fattyacid #
length(MFAs.remove.from.firstsample) #8 MFAs to remove
dim(fattyacid.v2.firstsample.v2) #195 rows, 94 columns
fattyacid.v2.firstsample.v2 <- fattyacid.v2.firstsample.v2[,-c(which(colnames(fattyacid.v2.firstsample.v2) %in% MFAs.remove.from.firstsample))]
dim(fattyacid.v2.firstsample.v2) #195 rows, 86 columns

fattyacid.totals.secondsample <- data.frame(colSums(fattyacid.v2.secondsample.v2[,7:84]))
fattyacid.totals.secondsample$fattyacid <-rownames(fattyacid.totals.secondsample)
colnames(fattyacid.totals.secondsample) <- c("column.total", "fattyacid")
MFAs.remove.from.secondsample <- fattyacid.totals.secondsample[fattyacid.totals.secondsample$column.total==0,]
MFAs.remove.from.secondsample <- MFAs.remove.from.secondsample$fattyacid #
length(MFAs.remove.from.secondsample) #9 MFAs to remove
dim(fattyacid.v2.secondsample.v2) #198 rows, 94 columns
fattyacid.v2.secondsample.v2 <- fattyacid.v2.secondsample.v2[,-c(which(colnames(fattyacid.v2.secondsample.v2) %in% MFAs.remove.from.secondsample))]
dim(fattyacid.v2.secondsample.v2) #195 rows, 85 columns

### Calculating the shannon and inverse simpson diversity of the MFAs in each sample. 
fattyacid.v2.firstsample.v3 <- fattyacid.v2.firstsample.v2
fattyacid.v2.firstsample.v3$MFA.diversity.inversesimpson<-diversity(x=fattyacid.v2.firstsample.v3[,7:75], index="invsimpson")
fattyacid.v2.firstsample.v3$MFA.diversity.shannon<-diversity(x=fattyacid.v2.firstsample.v3[,7:75], index="shannon")

fattyacid.v2.secondsample.v3 <- fattyacid.v2.secondsample.v2
fattyacid.v2.secondsample.v3$MFA.diversity.inversesimpson<-diversity(x=fattyacid.v2.secondsample.v3[,7:74], index="invsimpson")
fattyacid.v2.secondsample.v3$MFA.diversity.shannon<-diversity(x=fattyacid.v2.secondsample.v3[,7:74], index="shannon")

### Running linear mixed effects modeling 
hist(fattyacid.v2.firstsample.v3$MFA.diversity.inversesimpson)#normal
hist(fattyacid.v2.firstsample.v3$MFA.diversity.shannon)#normal
hist(fattyacid.v2.firstsample.v3$LogSCC.milking1.beforesample) #slightly right skewed
summary(lmer(MFA.diversity.inversesimpson~LogSCC.milking1.beforesample+(1|cow), data=fattyacid.v2.firstsample.v3)) #not significant
summary(lmer(MFA.diversity.inversesimpson~LogSCC.milking1.beforesample+(LogSCC.milking1.beforesample|cow), data=fattyacid.v2.firstsample.v3)) #singularity warning
summary(lmer(MFA.diversity.shannon~LogSCC.milking1.beforesample+(1|cow), data=fattyacid.v2.firstsample.v3)) #not significant
summary(lmer(MFA.diversity.shannon~LogSCC.milking1.beforesample+(LogSCC.milking1.beforesample|cow), data=fattyacid.v2.firstsample.v3)) #warning -- model failed to converge

hist(fattyacid.v2.secondsample.v3$MFA.diversity.inversesimpson)#normal
hist(fattyacid.v2.secondsample.v3$MFA.diversity.shannon)#normal
hist(fattyacid.v2.secondsample.v3$LogSCC.milking1.aftersample) #slightly right skewed

summary(lmer(LogSCC.milking1.aftersample~MFA.diversity.inversesimpson+(1|cow), data=fattyacid.v2.secondsample.v3)) #not significant
summary(lmer(LogSCC.milking1.aftersample~MFA.diversity.inversesimpson+(MFA.diversity.inversesimpson|cow), data=fattyacid.v2.secondsample.v3)) #warning-- model failed to converge

summary(lmer(LogSCC.milking1.aftersample~MFA.diversity.shannon+(1|cow), data=fattyacid.v2.secondsample.v3)) #not significant
summary(lmer(LogSCC.milking1.aftersample~MFA.diversity.shannon+(MFA.diversity.shannon|cow), data=fattyacid.v2.secondsample.v3)) #not significant

### Running repeated measures correlation
fattyacid.v2.firstsample.v3$cow <- as.factor(fattyacid.v2.firstsample.v3$cow)
rmcorr(participant=cow, measure1=LogSCC.milking1.beforesample, measure2=MFA.diversity.shannon, dataset=fattyacid.v2.firstsample.v3) #not significant
rmcorr(participant=cow, measure1=LogSCC.milking1.beforesample, measure2=MFA.diversity.inversesimpson, dataset=fattyacid.v2.firstsample.v3) #not significant

fattyacid.v2.secondsample.v3$cow <- as.factor(fattyacid.v2.secondsample.v3$cow)
rmcorr(participant=cow, measure1=LogSCC.milking1.aftersample, measure2=MFA.diversity.shannon, dataset=fattyacid.v2.secondsample.v3) #not significant
rmcorr(participant=cow, measure1=LogSCC.milking1.aftersample, measure2=MFA.diversity.inversesimpson, dataset=fattyacid.v2.secondsample.v3) #not significant

### Step 3: Determine if there is a correlation between any of the individual MFAs (that are present in least 75 % of samples) and log(SCC)
# Remove the MFAs that are present in fewer than 50 % of the samples. 
Function.FindSparse.MFAs.InFirstSamples <- function(MFA) {
    MFA.remove <- data.frame()
    for(i in MFA) {
        if(sum(fattyacid.v2.firstsample.v3[[i]]>0)<((length(fattyacid.v2.firstsample.v3$sample.id))*0.50)){
            MFA.remove <- rbind(MFA.remove, i)
        }
    }
    print(MFA.remove)
}

MFA.list <- colnames(fattyacid.v2.firstsample.v3[,7:75])
MFA.to.remove <- Function.FindSparse.MFAs.InFirstSamples(MFA.list)
MFA.to.remove <- MFA.to.remove[,1]

fattyacid.v2.firstsample.v4 <- fattyacid.v2.firstsample.v3[,-c(which(colnames(fattyacid.v2.firstsample.v3)%in%MFA.to.remove))]

Function.FindSparse.MFAs.InSecondSamples <- function(MFA) {
    MFA.remove <- data.frame()
    for(i in MFA) {
        if(sum(fattyacid.v2.secondsample.v3[[i]]>0)<((length(fattyacid.v2.secondsample.v3$sample.id))*0.5)){
            MFA.remove <- rbind(MFA.remove, i)
        }
    }
    print(MFA.remove)
}

MFA.list <- colnames(fattyacid.v2.secondsample.v3[,7:74])
MFA.to.remove <- Function.FindSparse.MFAs.InSecondSamples(MFA.list)
MFA.to.remove <- MFA.to.remove[,1]
fattyacid.v2.secondsample.v4 <- fattyacid.v2.secondsample.v3[,-c(which(colnames(fattyacid.v2.secondsample.v3)%in%MFA.to.remove))]


# Remove MFAs that are present in less than 90 % of cows
Function.NonZero.PerCow.PerMFA.FirstSamples<-function(MFAs){
    df <- data.frame()
    for(i in MFAs){
        MFA.subset <- fattyacid.v2.firstsample.v4[,c("cow",i)]
        MFA.subset.v2 <- MFA.subset[MFA.subset[[i]]>0,]
        number.cows<-length(unique(MFA.subset.v2$cow))
        percentage.cows <- (number.cows/66)*100
        MFAs.cows <-data.frame(i, percentage.cows)  
        df <- rbind(df, MFAs.cows)
    }
    print(df)
}

MFAs <- colnames(fattyacid.v2.firstsample.v4[,7:67])
nonzero.per.cow.per.MFA.FirstSample <- Function.NonZero.PerCow.PerMFA.FirstSamples(MFAs)
colnames(nonzero.per.cow.per.MFA.FirstSample) <-c("MFA", "percent.of.cows")
MFAs.50PercentSamples.90PercentCows.FirstSamples <- nonzero.per.cow.per.MFA.FirstSample[nonzero.per.cow.per.MFA.FirstSample$percent.of.cows>=90,][["MFA"]]

#Repeating for the second sampling dates from each period. 
Function.NonZero.PerCow.PerMFA.SecondSamples<-function(MFAs){
  df <- data.frame()
  for(i in MFAs){
    MFA.subset <- fattyacid.v2.secondsample.v4[,c("cow",i)]
    MFA.subset.v2 <- MFA.subset[MFA.subset[[i]]>0,]
    number.cows<-length(unique(MFA.subset.v2$cow))
    percentage.cows <- (number.cows/65)*100
    MFAs.cows <-data.frame(i, percentage.cows)  
    df <- rbind(df, MFAs.cows)
  }
  print(df)
}

MFAs <- colnames(fattyacid.v2.secondsample.v4[,7:67])
nonzero.per.cow.per.MFA.SecondSample <- Function.NonZero.PerCow.PerMFA.SecondSamples(MFAs)
colnames(nonzero.per.cow.per.MFA.SecondSample) <-c("MFA", "percent.of.cows")
MFAs.50PercentSamples.90PercentCows.SecondSamples <- nonzero.per.cow.per.MFA.SecondSample[nonzero.per.cow.per.MFA.SecondSample$percent.of.cows>=90,][["MFA"]]

### Running a function to test for a correlation between each MFA (response) and log(SCC) (predictor). 
#And since several MFAs have skewed distribution, the MFA measures are also 
#being transformed by taking the square root and it is the square root measure 
#that will be used in the linear mixed effect model.  
fattyacid.v2.firstsample.v4[,c(7:68,71:76)] <- sqrt(fattyacid.v2.firstsample.v4[,c(7:68,71:76)])

Function.MFA.SCC.FirstSample <- function(MFA){
    coefficient.pvalue <- data.frame(coefficient=NA, pvalue=NA, MFA=NA)
    for(i in MFA){
        form<-as.formula(paste0(i,"~","LogSCC.milking1.beforesample", "+(1|cow)"))
        lmer.output<-try(summary(lmer(form, data=fattyacid.v2.firstsample.v4)), silent=TRUE)
        lmer.output2<-try(data.frame(coefficient=lmer.output[["coefficients"]][2,1], pvalue=lmer.output[["coefficients"]][2,5], MFA=i), silent=TRUE)
        coefficient.pvalue <- rbind(coefficient.pvalue, lmer.output2)  
    }
    print(coefficient.pvalue)
}

MFA <- colnames(fattyacid.v2.firstsample.v4[,c(7:68,71:76)])
MFA.SCC.correlations.FirstSample <- Function.MFA.SCC.FirstSample(MFA)
MFA.SCC.correlations.FirstSample <- MFA.SCC.correlations.FirstSample[-1,]
MFA.SCC.correlations.FirstSample$FDR <- p.adjust(p=MFA.SCC.correlations.FirstSample$pvalue, method="fdr") #no significant relationships

fattyacid.v2.secondsample.v4[,c(7:68, 71:76)] <- sqrt(fattyacid.v2.secondsample.v4[,c(7:68,71:76)])

Function.MFA.SCC.SecondSample <- function(MFA){
    coefficient.pvalue <- data.frame(coefficient=NA, pvalue=NA, MFA=NA)
    for(i in MFA){
        form<-as.formula(paste0("LogSCC.milking1.aftersample","~", i, "+(1|cow)"))
        lmer.output<-try(summary(lmer(form, data=fattyacid.v2.secondsample.v4)), silent=TRUE)
        lmer.output2<-try(data.frame(coefficient=lmer.output[["coefficients"]][2,1], pvalue=lmer.output[["coefficients"]][2,5], MFA=i), silent=TRUE)
        coefficient.pvalue <- rbind(coefficient.pvalue, lmer.output2)  
    }
    print(coefficient.pvalue)
}

MFA <- colnames(fattyacid.v2.secondsample.v4[,c(7:68,71:76)])
MFA.SCC.correlations.SecondSample <- Function.MFA.SCC.SecondSample(MFA) #no significant relationsihps

#None of the (square root transformed) MFAs were correlated with logSCC when 
#running LMM. 

#Trying repeated measures correlation (RMCORR)
Function.RMCORR.MFA.SCC.FirstSample <- function(MFA){
  coefficient.pvalue <- data.frame(coefficient=NA, pvalue=NA, MFA=NA)
  for(i in MFA){
    rmcorr.output<-try(rmcorr(measure1=i, measure2=LogSCC.milking1.beforesample, participant=cow, dataset=fattyacid.v2.firstsample.v4), silent=TRUE)
    rmcorr.output2<-try(data.frame(coefficient=rmcorr.output[["r"]], pvalue=rmcorr.output[["p"]], MFA=i), silent=TRUE)
    coefficient.pvalue <- rbind(coefficient.pvalue, rmcorr.output2)  
  }
  print(coefficient.pvalue)
}

MFA <- colnames(fattyacid.v2.firstsample.v4[,c(7:68,71:76)])
RMCORR.MFA.SCC.correlations.FirstSample <- Function.RMCORR.MFA.SCC.FirstSample(MFA)
RMCORR.MFA.SCC.correlations.FirstSample <- RMCORR.MFA.SCC.correlations.FirstSample[-1,]
RMCORR.MFA.SCC.correlations.FirstSample$FDR <- p.adjust(p=RMCORR.MFA.SCC.correlations.FirstSample$pvalue, method="fdr") #no significant relationships

Function.RMCORR.MFA.SCC.SecondSample <- function(MFA){
  coefficient.pvalue <- data.frame(coefficient=NA, pvalue=NA, MFA=NA)
  for(i in MFA){
    rmcorr.output<-try(rmcorr(measure1=i, measure2=LogSCC.milking1.aftersample, participant=cow, dataset=fattyacid.v2.secondsample.v4), silent=TRUE)
    rmcorr.output2<-try(data.frame(coefficient=rmcorr.output[["r"]], pvalue=rmcorr.output[["p"]], MFA=i), silent=TRUE)
    coefficient.pvalue <- rbind(coefficient.pvalue, rmcorr.output2)  
  }
  print(coefficient.pvalue)
}

MFA <- colnames(fattyacid.v2.secondsample.v4[,c(7:68,71:76)])
RMCORR.MFA.SCC.correlations.SecondSample <- Function.RMCORR.MFA.SCC.SecondSample(MFA)
RMCORR.MFA.SCC.correlations.SecondSample <- RMCORR.MFA.SCC.correlations.SecondSample[-1,]
RMCORR.MFA.SCC.correlations.SecondSample$FDR <- p.adjust(p=RMCORR.MFA.SCC.correlations.SecondSample$pvalue, method="fdr") #no significant relationships

