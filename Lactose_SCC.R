#This script describes the steps taken to determine if there is a relationship
#between lactose and SCC. 

#Load libraries
library(dplyr)
library(lme4)
library(lmerTest)
library(ggplot2)
library(lattice)
library(DHARMa)
library(rmcorr)

## STEP 1: read in the data files and merge them. 

#Reading in the AgSource data for the SCC measures from 11/21/17, 11/28/17, 01/02/18, 01/05/18, 03/13/18, and 03/19/18 (which are the same agsource SCC measures that were used for the analysis of SCC & microbiota data.)
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

#Reading in Peggy Tomasula's macronutrient dataset that includes lactose measures
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

#Need to subset the by day (i.e. take subset chunks of columns), then add back together under the same column headers so that there are only 
macronutrient.Day1 <- macronutrient[,c(3,5,8)]
macronutrient.Day2 <- macronutrient[,c(3,11,14)]
macronutrient.Day3 <- macronutrient[,c(3,18,21)]
macronutrient.Day4 <- macronutrient[,c(3,24,27)]
macronutrient.Day5 <- macronutrient[,c(3,31,34)]
macronutrient.Day6 <- macronutrient[,c(3,37,40)]

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
macronutrient <- rename(macronutrient, cow=`Cow ID`)

#Check for NA values and duplicate values in the sample_id column
sum(is.na(macronutrient$sample.id)) #the output is zero, so no NA values among the sample_id column
length(unique(macronutrient$sample.id)) #456 unique observations in sample_id so there are no duplicate sample IDs

#Remove the destination information ("_Wyn") from the sample IDs so that they can be matched and merged to other datasets
macronutrient$sample.id <-gsub(pattern="_Wyn", x=macronutrient$sample.id, replacement="") 
head(macronutrient)

#Add a column specifying the units of measure for the lactose measurements
macronutrient$lactose_Units <-"gram/100 gram of milk "

#Remove samples with "NA" for lactose measure
table(is.na(macronutrient$Lactose)) #six samples do not have lactose measure
macronutrient <- macronutrient[!is.na(macronutrient$Lactose),]
table(is.na(macronutrient$Lactose)) #all were removed. 

#Cows 5849 and 6233 were missed on the first collection timepoint (11/25/2017)
#But milk samples were collected from them on 11/26/2017, so the date in the
#sample id is incorrect and will be changed to "171126"
macronutrient$sample.id <- gsub(pattern="5849_171125", replacement="5849_171126", x=macronutrient$sample.id)
macronutrient$sample.id <- gsub(pattern="6233_171125", replacement="6233_171126", x=macronutrient$sample.id)

#Remove the cows that the Barile lab determined are outliers for lactose
macronutrient.v2 <- macronutrient[macronutrient$cow!="6229",]
macronutrient.v2 <- macronutrient.v2[macronutrient.v2$cow!="5651",]

#Matching SCC measurements by period 
SCC.pd0 <- merge(Nov.21.17, Nov.28.17, by="cow", all=FALSE)
length(unique(SCC.pd0$cow)) #76
SCC.pd1 <- merge(Jan.02.18, Jan.05.18, by="cow", all=FALSE)
length(unique(SCC.pd1$cow)) #76
SCC.pd2 <- merge(March.13.18, March.19.18, by="cow", all=FALSE)
length(unique(SCC.pd2$cow)) #74

#Separating lactose measures by period
macronutrient.v2$period<-macronutrient.v2$sample.id
macronutrient.v2$period <-gsub(pattern="^...._171125|^...._171126", replacement="0", macronutrient.v2$period)
macronutrient.v2$period <-gsub(pattern="^...._180104|^...._180105", replacement="1", macronutrient.v2$period)
macronutrient.v2$period <-gsub(pattern="^...._180314|^...._180315", replacement="2", macronutrient.v2$period)
table(macronutrient.v2$period)

macronutrient.v2.pd0 <- macronutrient.v2[macronutrient.v2$period=="0",]
macronutrient.v2.pd1 <- macronutrient.v2[macronutrient.v2$period=="1",]
macronutrient.v2.pd2 <- macronutrient.v2[macronutrient.v2$period=="2",]

#Merge SCC and lactose data by cow
macronutrient.v2.pd0 <- merge(macronutrient.v2.pd0, SCC.pd0, by="cow", all=FALSE)
macronutrient.v2.pd1 <- merge(macronutrient.v2.pd1, SCC.pd1, by="cow", all=FALSE)
macronutrient.v2.pd2 <- merge(macronutrient.v2.pd2, SCC.pd2, by="cow", all=FALSE)
macronutrient.v3 <- rbind(macronutrient.v2.pd0, macronutrient.v2.pd1, macronutrient.v2.pd2)

#Remove antibiotic-treated cows
length(unique(macronutrient.v3$cow)) #74 cows
antibiotic.cows <- c("5020", "5212", "5297", "5405", "5697", "6076", "6215", "6232", "6233", "6241")
macronutrient.v3 <- macronutrient.v3[-c(which(macronutrient.v3$cow %in% antibiotic.cows)),]
length(unique(macronutrient.v3$cow)) #64 cows

#Since samples from the same period have the same SCC measures, going to break
#up the dataset into samples from first versus second collection day of a period. 
macronutrient.v3.firstcollection <- macronutrient.v3[grep(pattern="^...._171125|^...._180104|^...._180314", macronutrient.v3$sample.id),]
macronutrient.v3.secondcollection <- macronutrient.v3[grep(pattern="^...._171126|^...._180105|^...._180315", macronutrient.v3$sample.id),]

#Choose only cows that have three samples per period
macronutrient.v3.firstcollection.v2 <- macronutrient.v3.firstcollection[!is.na(macronutrient.v3.firstcollection$Lactose),]
macronutrient.v3.firstcollection.v2 <- macronutrient.v3.firstcollection[!is.na(macronutrient.v3.firstcollection$SCCx1000.milking1.aftersample),]
samples.percow.firstcollection<-data.frame(table(macronutrient.v3.firstcollection.v2$cow))
cows.keep.firstcollection <- samples.percow.firstcollection[samples.percow.firstcollection$Freq==3,]["Var1"]
macronutrient.v3.firstcollection.v3 <- macronutrient.v3.firstcollection.v2[c(which(macronutrient.v3.firstcollection.v2$cow%in%cows.keep.firstcollection$Var1)),]

macronutrient.v3.secondcollection.v2 <- macronutrient.v3.secondcollection[!is.na(macronutrient.v3.secondcollection$Lactose),]
macronutrient.v3.secondcollection.v2 <- macronutrient.v3.secondcollection.v2[!is.na(macronutrient.v3.secondcollection.v2$SCCx1000.milking1.aftersample),]
samples.percow.secondcollection<-data.frame(table(macronutrient.v3.secondcollection.v2$cow))
cows.keep.secondcollection <- samples.percow.secondcollection[samples.percow.secondcollection$Freq==3,]["Var1"]
cows.keep.secondcollection <- as.character(cows.keep.secondcollection$Var1)
macronutrient.v3.secondcollection.v3 <- macronutrient.v3.secondcollection.v2[c(which(macronutrient.v3.secondcollection.v2$cow%in%cows.keep.secondcollection)),]

## STEP 2: Determine if lactose is a significant predictor of SCC, and if SCC is
#a significant predictor of lactose.

#Look to see if the distribution of lactose is heavily skewed. 
macronutrient.v3.firstcollection.v3$Lactose<-as.numeric(macronutrient.v3.firstcollection.v3$Lactose)
hist(macronutrient.v3.firstcollection.v3$Lactose) #Slightly skewed but not extremely
macronutrient.v3.secondcollection.v3$Lactose<-as.numeric(macronutrient.v3.secondcollection.v3$Lactose)
hist(macronutrient.v3.secondcollection.v3$Lactose) #maybe slightly left skewed. 

#SCC is in thousands, so multiply by 1000 then look at distribution
hist(1000*macronutrient.v3.firstcollection.v3$SCCx1000.milking1.beforesample) #heavily skewed
hist(1000*macronutrient.v3.firstcollection.v3$SCCx1000.milking1.aftersample) #heavily skewed
hist(1000*macronutrient.v3.secondcollection.v3$SCCx1000.milking1.beforesample) #heavily skewed
hist(1000*macronutrient.v3.secondcollection.v3$SCCx1000.milking1.aftersample) #heavily skewed

#Performing log transformation of SCC data
macronutrient.v3.firstcollection.v3$LogSCC.milking1.beforesample<-log10(1000*macronutrient.v3.firstcollection.v3$SCCx1000.milking1.beforesample)
macronutrient.v3.firstcollection.v3$LogSCC.milking1.aftersample<-log10(1000*macronutrient.v3.firstcollection.v3$SCCx1000.milking1.aftersample)
macronutrient.v3.secondcollection.v3$LogSCC.milking1.beforesample<-log10(1000*macronutrient.v3.secondcollection.v3$SCCx1000.milking1.beforesample)
macronutrient.v3.secondcollection.v3$LogSCC.milking1.aftersample<-log10(1000*macronutrient.v3.secondcollection.v3$SCCx1000.milking1.aftersample)

#Performing linear mixed effect modeling
summary(lmer(Lactose~LogSCC.milking1.beforesample+(1|cow), macronutrient.v3.firstcollection.v3)) #significant, negative correlation
summary(lmer(Lactose~LogSCC.milking1.beforesample+(LogSCC.milking1.beforesample|cow), macronutrient.v3.firstcollection.v3)) #singular fit warning

summary(lmer(LogSCC.milking1.aftersample~Lactose+(1|cow), macronutrient.v3.secondcollection.v3)) #significant negative relationship
summary(lmer(LogSCC.milking1.aftersample~Lactose+(Lactose|cow), macronutrient.v3.secondcollection.v3)) #model failure warning


#Graph Lactose vs. Log(SCC)
results.firstcollection.Intercept <- lmer(Lactose~LogSCC.milking1.beforesample+(1|cow), data=macronutrient.v3.firstcollection.v3)
ggplot(aes(x=LogSCC.milking1.beforesample, y=Lactose), data=macronutrient.v3.firstcollection.v3)+geom_point()+xlab("log(SCC)")+ylab("lactose (g/100g milk))")+geom_line(aes(y=predict(results.firstcollection.Intercept), group=cow))+annotate("text",x=6.2, y=5.2, label="-0.096")
ggsave(filename="LogSCC.vs.Lactose.firstcollectiondates.jpeg", device="jpeg")

results.secondcollection.Intercept <- lmer(LogSCC.milking1.aftersample~Lactose+(1|cow), data=macronutrient.v3.secondcollection.v3)
ggplot(aes(x=Lactose, y=LogSCC.milking1.aftersample), data=macronutrient.v3.secondcollection.v3)+geom_point()+ylab("log(SCC)")+xlab("lactose (g/100g milk)")+geom_line(aes(y=predict(results.secondcollection.Intercept), group=cow))+annotate("text", x=5.2, y=6.4, label="-0.79")
ggsave(filename="Lactose.vs.LogSCC.secondcollectiondates.jpeg")

#determine if the model is a good fit for lactose and SCC
#Does the model meet the assumptions?

#Plotting the residuals against fitted values to visually determine if
#there is a non-linear relationship between the predictors and the response
#variables and if there is unequal variances (and can be used to identify outliers)
plot(fitted(results.firstcollection.Intercept), resid(results.firstcollection.Intercept)) 

#The points showed a consistent and fairly even spread around the residual=0 line
#so this suggests that the assumptions of linearity and equal variances has not been violated.
#Looking for normality of residuals with a histogram
hist(resid(results.firstcollection.Intercept)) #residuals look like they have fairly normal distribution

#Looking at QQ plot to confirm that the observations come from a normal distribution
qqmath(results.firstcollection.Intercept) #observations mostly fall on the line but not near the ends


#Also try running repeated measure correlation with plots
macronutrient.v3.firstcollection.v3$cow <- as.factor(macronutrient.v3.firstcollection.v3$cow)
macronutrient.v3.secondcollection.v3$cow <- as.factor(macronutrient.v3.secondcollection.v3$cow)
rmcor.firstcollectionLactose.beforesampleSCC <- rmcorr(participant=cow, measure1=LogSCC.milking1.beforesample, measure2=Lactose, dataset=macronutrient.v3.firstcollection.v3) #significant negative correlation
rmcor.secondcollectionLactose.aftersampleSCC <- rmcorr(participant=cow, measure1=LogSCC.milking1.aftersample, measure2=Lactose, dataset=macronutrient.v3.secondcollection.v3) #significant negative correlation
plot(rmcor.firstcollectionLactose.beforesampleSCC)
plot(rmcor.secondcollectionLactose.aftersampleSCC)
