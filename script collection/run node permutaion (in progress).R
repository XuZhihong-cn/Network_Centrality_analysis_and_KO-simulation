
library(dplyr)
library(glmmTMB)

IndividualData <- read.csv("raw_data/IndData.csv",  sep = ";")
header2 <- as.character(read.csv("raw_data/IndData.csv", header = FALSE, fileEncoding="UTF-8-BOM", sep = ";")[1, ])
IBD <- IndividualData
colnames(IBD) <- header2

# Load the data
SampleData <- read.csv("raw_data/SampleData.csv",  sep = ";")
header1 <- as.character(read.csv("raw_data/SampleData.csv", header = FALSE, fileEncoding="UTF-8-BOM", sep = ";")[1, ])
SBD <- SampleData
colnames(SBD) <- header1


# "Clean" the data, transform variables with their correct identity (e.g. "character" into "numeric")
# set data as numeric or factors. 
SBD[, c("D", "S", "EV", "A", "EL")] <- apply(SBD[, c("D", "S", "EV", "A", "EL")], 2, function(x)as.numeric(scale(x)))
SBD$EPG <- round(SBD$EPG)
SBD$CD<-as.factor(SBD$CD)
SBD$ID<-as.factor(SBD$ID)
SBD$SP<-as.factor(SBD$Species)
SBDaf<-droplevels(subset(droplevels(subset(SBD, AG=="adult")),SEX=="female"))
SBDjuv<-droplevels(subset(SBD, AG=="juvenile"))



randomall <- function(i){
  # create a temporary dataframe
  s <- sample(1:52)
  IBD1 <- IBD[,c("ID","D","S","EV")]
  IBD1$D <- IBD$D[s]
  IBD1$S <- IBD$S[s]
  IBD1$EV <- IBD$EV[s]
  SBD1<- SBD[,c(-4:-6)]
  SBD1<-left_join(SBD1,IBD1,by="ID")
  SBD1[, c("D", "S", "EV")] <- apply(SBD1[, c("D", "S", "EV")], 2, function(x)as.numeric(scale(x)))
  # run the glmms  
  try(tempD <- glmmTMB(EPG ~ D+A+SEX+EL+(1|CD)+(1|ID)+(1|SP), ziformula = ~ D+A+SEX+EL+(1|CD)+(1|ID)+(1|SP), data = SBD1, family = 'nbinom2',doFit=TRUE))
  try(tempS <- glmmTMB(EPG ~ S+A+SEX+EL+(1|CD)+(1|ID)+(1|SP), ziformula = ~ S+A+SEX+EL+(1|CD)+(1|ID)+(1|SP), data = SBD1, family ='nbinom2',doFit=TRUE))
  try(tempEV <- glmmTMB(EPG ~ EV+A+SEX+EL+(1|CD)+(1|ID)+(1|SP), ziformula = ~ EV+A+SEX+EL+(1|CD)+(1|ID)+(1|SP), data = SBD1 ,family = 'nbinom2',doFit=TRUE))
  
  # extract confidence intervals and estimates   
  try(temp <- c(confint(tempEV)[2,1],tempEV$fit$par[2],confint(tempEV)[2,2],confint(tempS)[2,1],tempS$fit$par[2],confint(tempS)[2,2],confint(tempD)[2,1],tempD$fit$par[2],confint(tempD)[2,2]))
  try(temp)
  
} 