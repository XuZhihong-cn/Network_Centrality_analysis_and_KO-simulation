# this script is for running node permutations.
# the nodes in the social network is re-arranged randomly, then the glmm recalculated and the results are compared with origin (see ms text)
# here we need the individual data "IndData" as extra data + the previous script "run GLMMs"

library(dplyr)
library(glmmTMB)
library(foreach)
library(doParallel)
library(doSNOW)


# Load the data
IndividualData <- read.csv("raw_data/IndData.csv",  sep = ";")
header2 <- as.character(read.csv("raw_data/IndData.csv", header = FALSE, fileEncoding="UTF-8-BOM", sep = ";")[1, ])
IBD <- IndividualData
colnames(IBD) <- header2
IBDaf<-droplevels(subset(droplevels(subset(IBD, AG=="adult")),SEX=="female"))
IBDjuv<-droplevels(subset(IBD, AG=="juvenile"))


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


## node permutation functions   ##
permall <- function(i){
  # create the re-arranged data frame, in which the nodes and the individual they represented are re-arranged.
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


permaf <- function(i){
  # create the re-arranged data frame, in which the nodes and the individual they represented are re-arranged.
  s <- sample(1:20)
  IBDaf1 <- IBDaf[,c("ID","D","S","EV")]
  IBDaf1$D <- IBDaf$D[s]
  IBDaf1$S <- IBDaf$S[s]
  IBDaf1$EV <- IBDaf$EV[s]
  SBDaf1<- SBDaf[,c(-4:-6)]
  SBDaf1<-left_join(SBDaf1,IBDaf1,by="ID")
  SBDaf1[, c("D", "S", "EV")] <- apply(SBDaf1[, c("D", "S", "EV")], 2, function(x)as.numeric(scale(x)))
  # run the glmms  
  try(tempD <- glmmTMB(EPG ~ D+A+SEX+EL+(1|CD)+(1|ID)+(1|SP), ziformula = ~ D+A+SEX+EL+(1|CD)+(1|ID)+(1|SP), data = SBDaf1, family = 'nbinom2',doFit=TRUE))
  try(tempS <- glmmTMB(EPG ~ S+A+SEX+EL+(1|CD)+(1|ID)+(1|SP), ziformula = ~ S+A+SEX+EL+(1|CD)+(1|ID)+(1|SP), data = SBDaf1, family ='nbinom2',doFit=TRUE))
  try(tempEV <- glmmTMB(EPG ~ EV+A+SEX+EL+(1|CD)+(1|ID)+(1|SP), ziformula = ~ EV+A+SEX+EL+(1|CD)+(1|ID)+(1|SP), data = SBDaf1 ,family = 'nbinom2',doFit=TRUE))
  
  # extract confidence intervals and estimates   
  try(temp <- c(confint(tempEV)[2,1],tempEV$fit$par[2],confint(tempEV)[2,2],confint(tempS)[2,1],tempS$fit$par[2],confint(tempS)[2,2],confint(tempD)[2,1],tempD$fit$par[2],confint(tempD)[2,2]))
  try(temp)
  
} 


permjuv <- function(i){
  # create the re-arranged data frame, in which the nodes and the individual they represented are re-arranged.
  s <- sample(1:24)
  IBDjuv1 <- IBDjuv[,c("ID","D","S","EV")]
  IBDjuv1$D <- IBDjuv$D[s]
  IBDjuv1$S <- IBDjuv$S[s]
  IBDjuv1$EV <- IBDjuv$EV[s]
  SBDjuv1<- SBDjuv[,c(-4:-6)]
  SBDjuv1<-left_join(SBDjuv1,IBDjuv1,by="ID")
  SBDjuv1[, c("D", "S", "EV")] <- apply(SBDjuv1[, c("D", "S", "EV")], 2, function(x)as.numeric(scale(x)))
  # run the glmms  
  try(tempD <- glmmTMB(EPG ~ D+A+SEX+EL+(1|CD)+(1|ID)+(1|SP), ziformula = ~ D+A+SEX+EL+(1|CD)+(1|ID)+(1|SP), data = SBDjuv1, family = 'nbinom2',doFit=TRUE))
  try(tempS <- glmmTMB(EPG ~ S+A+SEX+EL+(1|CD)+(1|ID)+(1|SP), ziformula = ~ S+A+SEX+EL+(1|CD)+(1|ID)+(1|SP), data = SBDjuv1, family ='nbinom2',doFit=TRUE))
  try(tempEV <- glmmTMB(EPG ~ EV+A+SEX+EL+(1|CD)+(1|ID)+(1|SP), ziformula = ~ EV+A+SEX+EL+(1|CD)+(1|ID)+(1|SP), data = SBDjuv1 ,family = 'nbinom2',doFit=TRUE))
  
  # extract confidence intervals and estimates   
  try(temp <- c(confint(tempEV)[2,1],tempEV$fit$par[2],confint(tempEV)[2,2],confint(tempS)[2,1],tempS$fit$par[2],confint(tempS)[2,2],confint(tempD)[2,1],tempD$fit$par[2],confint(tempD)[2,2]))
  try(temp)
  
} 



# Start the node permutation procedure
T1 <- Sys.time()
perm <- 3 

# Set up progress bar
pb <- txtProgressBar(max = perm, style = 3)
progress <- function(n)  setTxtProgressBar(pb, n)
opts <- list(progress = progress)

# Split cores to work in parallel
core <- detectCores() 
cl<- makeCluster(core-1)   
registerDoParallel(cl)
registerDoSNOW(cl)

# Run the permutation
PMres <- foreach(i = 1:perm,.packages = c('glmmTMB','dplyr'), .options.snow = opts) %dopar% permall (i)  
PMresaf <- foreach(i = 1:perm,.packages = c('glmmTMB','dplyr'), .options.snow = opts) %dopar% permaf (i)  
PMresjuv <- foreach(i = 1:perm,.packages = c('glmmTMB','dplyr'), .options.snow = opts) %dopar% permjuv (i)  

# Save the results
PMres <- as.data.frame(do.call(rbind, PMres))
PMresaf <- as.data.frame(do.call(rbind, PMresaf))
PMresjuv <- as.data.frame(do.call(rbind, PMresjuv))

write.csv(na.omit(PMres),file = "PMres.csv")
write.csv(na.omit(PMresaf),file = "PMresaf.csv")
write.csv(na.omit(PMresjuv),file = "PMresjuv.csv")

# End of process
stopImplicitCluster()
stopCluster(cl) 
T2 <- Sys.time()
T2-T1