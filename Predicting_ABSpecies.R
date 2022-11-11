#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

## need to load models from "RFmodel_ABSpecies.RData"
load("RFmodel_ABSpecies.RData")
########################################################################################

data_in <- read.csv(paste(args[1], "analook.csv", sep="/")) #call data into R - use the csv filename also as part of the output filenames
#summary(data_in)
#head(data_in)
#names(data_in)

# data <- data_in[c(5,7:17)] ##create dataframe of only parameters for the model - don't include st, Qual, Pmc, or TBC
# The python script will do this trimming for us

#summary(data)
names(data)
#table(data_in$Date,data_in$Site)

######################################################################
### Set 4
library(caret) ##load package
library(randomForest)


rownames(models4)

dataM <- predict(models4, newdata=data_in, "prob") ##use models to predict bat calls
#testPred <- predict(models4, newdata=data) 
#lapply(testPred, function (x) x[1:5]) 

write.table(dataM,file=paste(args[1], "M1_analook.csv", sep="/"),sep=",",row.names=F) ##export certainty score csv file

######################################################################
## Aggregrating pulses
library(plyr) ##load package 

M1 <- read.csv(paste(args[1], "M1_analook.csv", sep="/")) ##read back in certainty score csv file
head(M1)
M1$Filename <- data_in$Filename ##add columns to M1
M1$Site <- data_in$Site ##add columns to M1
M1$Date <- data_in$Date ##add columns to M1

write.table(M1,file=paste(args[1], "M1_analook.csv", sep="/"),sep=",",row.names=F) ##export more complete certainty score csv file

names(M1)
nrow(M1)
head(M1)
#
#
### run through decision tree for bat call predictions - if below certainty score identify the call as "unknown"
#M1$cs60 <- as.factor(ifelse(M1$rf4.noise>=0.60,"noise",
#                            ifelse(M1$rf4.EPFU>=0.60,"EPFU",
#                                   ifelse(M1$rf4.LABO>=0.60,"LABO",
#                                          ifelse(M1$rf4.LACI>=0.60,"LACI",
#                                                 ifelse(M1$rf4.LANO>=0.60,"LANO",
#                                                        ifelse(M1$rf4.MYCA>=0.60,"MYCA",
#                                                        ifelse(M1$rf4.MYCI>=0.60,"MYCI",
#                                                               ifelse(M1$rf4.MYEV>=0.60,"MYEV",
#                                                                      ifelse(M1$rf4.MYLU>=0.60,"MYLU",
#                                                                            ifelse(M1$rf4.MYSE>=0.60,"MYSE",
#                                                                                    ifelse(M1$rf4.MYVO>=0.60,"MYVO",
#                                                                                           ifelse(M1$rf4.EPFU>=0.30 & M1$rf4.EPFU<0.60 & M1$rf4.LANO>=0.30 & M1$rf4.LANO<0.60 ,"EPFU-LANO",
#                                                                                                  ifelse(M1$rf4.LABO>=0.30 & M1$rf4.LABO<0.60 & M1$rf4.MYLU>=0.30 & M1$rf4.MYLU<0.60 ,"LABO-MYLU",
#                                                                                                         ifelse(M1$rf4.MYEV>=0.30 & M1$rf4.MYEV<0.60 & M1$rf4.MYSE>=0.30 & M1$rf4.MYSE<0.60 ,"MYEV-MYSE",
#                                                                                                                ifelse(M1$rf4.MYCI>=0.30 & M1$rf4.MYCI<0.60 | M1$rf4.MYLU>=0.30 & M1$rf4.MYLU<0.60 | M1$rf4.MYVO>=0.30 & M1$rf4.MYVO<0.60 ,"MYOTIS40K","unknown"))))))))))))))))
#
#
#
##summary(M1)
##levels(M1$cs60)
#
#########################################
####rf4
#
###create the function to aggregate the pulses
#rf4.fun <- function(x){
#  tbl <- table(x$cs60)
#  nmax <- sum(tbl == max(tbl))
#  if (nmax == 1)
#    x$rf4.freq <- rep(names(tbl)[which.max(tbl)],nrow(x))
#  else
#    x$rf4.freq <- "unknown"
#  x
#}
#
##nrow(M1)
#M1.rf4 <- ddply(M1,.(Filename),.fun=rf4.fun) ##apply the function
#
##summary(M1.rf4)
##nrow(M1.rf4)
#
#M1.rf4$Frf4.freq <- as.factor(M1.rf4$rf4.freq) ##convert the column to factor variable
##levels(M1.rf4$Frf4.freq)
#
#
#########################################
##head(data_in)
#Filename <-data_in$Filename ##creat Filename vector
#
#Fileuniq <- unique(Filename) ##strip out duplicates
#Fileuniq <-as.data.frame(Fileuniq) ##convert vector into data frame
##head(Fileuniq)
##nrow(Fileuniq)
#
#SumData <- Fileuniq ##change the name of the dataframe
##names(SumData)
##head(SumData)
#
#SumData$rf4 <- M1.rf4$Frf4.freq[match(SumData$Fileuniq, M1.rf4$Filename)] ##match bat identificaiton to each unique Filename
#SumData$rf4 <- as.factor(SumData$rf4) ##convert to factor variable
#
#SumData$Site <- M1$Site[match(SumData$Fileuniq, M1$Filename)] ##add in columns
#SumData$Date <- M1$Date[match(SumData$Fileuniq, M1$Filename)] ##add in columns
#
#summary(SumData)
##levels(SumData$rf4)
#
#table(SumData$rf4, SumData$Site)
##table(SumData$Date, SumData$Site)
##summary(SumData$Site)
#
#write.table(SumData,file=paste(args[1], "sumdata.csv", sep="/"),sep=",",row.names=F) ##export consolidated data as csv file
