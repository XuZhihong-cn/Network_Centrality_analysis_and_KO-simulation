# This script is for running the knock-out simulations
# scripts needed before: "run GLMMS" and "run randomisations"

# steps to take: 
# 1. modify the simulation function so as to set up the number of individuals to remove from the network, here we removed either everybody BUT adult females, everybody BUT juveniles, or 5% (3 individuals), 10% (6 individuals), 25% (13 individuals), and 50% (24 individuals) of the group, 
# 2. run the simulation, 
# 3. save the results for plotting

library(glmmTMB)
library(car)
library(dotwhisker)
library(igraph)
library(foreach)
library(doParallel)
library(doSNOW)

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


IndividualData <- read.csv("raw_data/IndData.csv",  sep = ";")
header2 <- as.character(read.csv("raw_data/IndData.csv", header = FALSE, fileEncoding="UTF-8-BOM", sep = ";")[1, ])
IBD <- IndividualData
colnames(IBD) <- header2


NetData <- read.csv("raw_Data/NetData.csv",  sep = ";")
# removed the first row of individual names
NetData <- NetData[,-1]
# turn the dataframe into a matrix object and name the col/row as individual ID. 
NetData <- as.matrix(NetData)
is.matrix(NetData)
rownames(NetData)<- colnames(NetData)


##  the simulation function   ##
simulation <- function(i)  {
  
  # create a list of names to exclude from the list of all individuals of interest
  # change the size to change the amount of individuals to remove
  # here is an example with 3 individuals removed from the network
  sam = sample(IBD$ID, size = 3)
  # create another list of IDs to delete individuals in the network measure dataframe. Here is the example with individuals without fecal samples collected in the study
  NS <- c('takana','mushi','neji','uso')
  # delete corresponding rows in the original dataframe 
  temp <- droplevels(SBD[-(which(SBD$ID %in% sam)),])    
  
  # delete corresponding individuals in the network matrix (you might have to reload it, see script "calculate network measures")
  NetData1 = NetData[-(which(rownames(NetData) %in% sam)),-(which(colnames(NetData) %in% sam))]
  
  # create network
  NetTemp <- graph.adjacency(NetData1, mode = "undirected", weighted = TRUE, diag = FALSE)   
  
  # recalculate eigenvector and save it into the original dataframe
  D <- degree(NetTemp)
  S <- strength(NetTemp)
  EV <-eigen_centrality(NetTemp)$vector
  
  D <- data.frame(ID = names(D), D = D)
  D <- droplevels(D[-(which(D$ID %in% NS)),])
  S <- data.frame(ID = names(S), S = S)
  S <- droplevels(S[-(which(S$ID %in% NS)),])  
  EV <- data.frame(ID = names(EV), EV = EV)
  EV <- droplevels(EV[-(which(EV$ID %in% NS)),])
  
  # put the new centrality measures back in original data frame
  temp$EV <- scale(EV$EV[as.numeric(temp$ID)])
  temp$D <- scale(D$D[as.numeric(temp$ID)])
  temp$S <- scale(S$S[as.numeric(temp$ID)])
  
  # run the models
  try(tempEV <- glmmTMB(EPG ~ EV+A+SEX+EL+(1|CD)+(1|ID)+(1|SP), ziformula = ~ EV+A+SEX+EL+(1|CD)+(1|ID)+(1|SP), data = temp,family = 'nbinom2',doFit=TRUE))
  try(tempS <- glmmTMB(EPG ~ S+A+SEX+EL+(1|CD)+(1|ID)+(1|SP), ziformula = ~ S+A+SEX+EL+(1|CD)+(1|ID)+(1|SP), data = temp, family ='nbinom2',doFit=TRUE))
  try(tempD <- glmmTMB(EPG ~ D+A+SEX+EL+(1|CD)+(1|ID)+(1|SP), ziformula = ~ D+A+SEX+EL+(1|CD)+(1|ID)+(1|SP), data = temp, family = 'nbinom2',doFit=TRUE))
  
  try(temp <- c(confint(tempEV)[2,1],tempEV$fit$par[2],confint(tempEV)[2,2],confint(tempS)[2,1],tempS$fit$par[2],confint(tempS)[2,2],confint(tempD)[2,1],tempD$fit$par[2],confint(tempD)[2,2]))
  
  try(temp)
  
}


# start the simulation procedure
T1 <- Sys.time()
perm <- 10

# Set up progress bar
pb <- txtProgressBar(max = perm, style = 3)
progress <- function(n)  setTxtProgressBar(pb, n)
opts <- list(progress = progress)

# Split cores to work in parallel if troubles
core <- detectCores()   
cl<- makeCluster(core-1)
registerDoParallel(cl)
registerDoSNOW(cl)

# Run the knock-out simulations - repeat this for each removal so at the end there is 4 "KOres", KOres5p, KOres10p, KOres25p and KOres50p
KOres <- foreach(i = 1:perm,.packages = c('glmmTMB','igraph','readxl'), .options.snow = opts) %dopar% simulation (i)  

# Save the results after each simulation run so there will be 4 datasets, KOres5p, KOres10p, KOres25p, KOres50p
KOres <- as.data.frame(do.call(rbind, KOres))
write.csv(na.omit(KOres5p),file = "KOres.csv")

# End of process
stopImplicitCluster()
stopCluster(cl) 
T2 <- Sys.time()
T2-T1


