# this script is for running randomizations 
# the social network is rewired according to a range of probabilities, then network measures are re-calculated, included again in the glmm and results updated and compared (see ms text)
# here we need the individual data "IndData" as extra data + the previous script "run GLMMs"

library(foreach)
library(doParallel)
library(doSNOW)
library(igraph)
library(glmmTMB)
# Load the data
IBD <- read.csv("raw_data/IndData.csv",  sep = ";") # be mindful of separator type
SBD <- read.csv("raw_data/SampleData.csv",  sep = ";")  # be mindful of separator type

# use these to reset the colnames if the first letter became garbled
header1 <- as.character(read.csv("raw_data/SampleData.csv", header = FALSE, fileEncoding="UTF-8-BOM", sep = ";")[1, ])
header2 <- as.character(read.csv("raw_data/IndData.csv", header = FALSE, fileEncoding="UTF-8-BOM", sep = ";")[1, ])
colnames(SBD) <- header1
colnames(IBD) <- header2

# "Clean" the data, transform variables with their correct identity (e.g. "character" into "numeric")
# set data as numeric or factors. 
SBD[, c("D", "S", "EV", "A", "EL")] <- apply(SBD[, c("D", "S", "EV", "A", "EL")], 2, function(x)as.numeric(scale(x)))
SBD$EPG <- round(SBD$EPG)
SBD$CD<-as.factor(SBD$CD)
SBD$ID<-as.factor(SBD$ID)
SBD$SP<-as.factor(SBD$Species)
SBDaf<-droplevels(subset(droplevels(subset(SBD, AG=="adult")),SEX=="female"))
SBDjuv<-droplevels(subset(SBD, AG=="juvenile"))

# create complete network
NetData <- read.csv("raw_Data/NetData.csv",  sep = ";") # be mindful of separator types
# removed the first row of individual names
NetData <- NetData[,-1]
# turn the dataframe into a matrix object and name the col/row as individual ID. 
NetData <- as.matrix(NetData)
is.matrix(NetData)
rownames(NetData)<- colnames(NetData)
netall <- graph.adjacency(NetData, mode = "undirected", weighted = TRUE, diag = FALSE)

# create adult female network
IBD_AF <- droplevels(subset(droplevels(subset(IBD, AG=="adult")),SEX=="female"))
#subset the proximity matrix to have an adult female only matrix
sam <- IBD_AF$ID
NetDataAF<- NetData[(which(rownames(NetData) %in% sam)),(which(colnames(NetData) %in% sam))]
netAF <- graph.adjacency(NetDataAF, mode = "undirected", weighted = TRUE, diag = FALSE)

# create juvenile network
IBD_juv <-droplevels(subset(IBD, AG=="juvenile"))
#subset the proximity matrix to have an juveniles only matrix
sam <- IBD_juv$ID
NetDataJUV<- NetData[(which(rownames(NetData) %in% sam)),(which(colnames(NetData) %in% sam))]
netJUV <- graph.adjacency(NetDataJUV, mode = "undirected", weighted = TRUE, diag = FALSE)

# make a list of individuals that do not have fecal sample
NS <- c('takana','mushi','neji','uso')


## randomisation function   ##

randomall <- function(i){
# create a temporary dataframe
  tempdata <- data.frame(EPG = SBD$EPG, ID = SBD$ID, A = SBD$A, SEX = SBD$SEX, EL = SBD$EL,EV=rep(0,length(SBD$ID)),D=rep(0,length(SBD$ID)),S=rep(0,length(SBD$ID)),CD=SBD$CD,SP=SBD$SP,NO=SBD$NO)
# create the randomized network  
  tempnet <- rewire(netall, with = each_edge(prob = runif(1,0,1), loops = FALSE, multiple = TRUE))
# calculate the new network measures
  D<-degree(tempnet, mode="all")
  S<-strength(tempnet,mode="all")
  EV<-eigen_centrality(tempnet, directed=FALSE)$vector
# put the new network measures in the temporary dataframe
  EV1 <- data.frame(ID = IBD$ID, EV = EV)
  D1 <- data.frame(ID = IBD$ID, D = D)
  S1 <- data.frame(ID = IBD$ID, S = S)
# remove "no sample" individuals  
  EV1 <- droplevels(EV1[-(which(EV1$ID %in% NS)),])
  S1 <- droplevels(S1[-(which(S1$ID %in% NS)),])
  D1 <- droplevels(D1[-(which(D1$ID %in% NS)),])

  tempdata$EV <- scale(EV1$EV[as.numeric(tempdata$ID)])
  tempdata$S <- scale(S1$S[as.numeric(tempdata$ID)])
  tempdata$D <- scale(D1$D[as.numeric(tempdata$ID)])
# run the glmms  
  try(tempD <- glmmTMB(EPG ~ D+A+SEX+EL+SP+(1|CD)+(1|ID)+(1|NO), ziformula = ~ D+A+SEX+EL+SP+(1|CD)+(1|ID)+(1|NO), data = tempdata, family = 'nbinom2',doFit=TRUE))
  try(tempS <- glmmTMB(EPG ~ S+A+SEX+EL+SP+(1|CD)+(1|ID)+(1|NO), ziformula = ~ S+A+SEX+EL+SP+(1|CD)+(1|ID)+(1|NO), data = tempdata, family ='nbinom2',doFit=TRUE))
  try(tempEV <- glmmTMB(EPG ~ EV+A+SEX+EL+SP+(1|CD)+(1|ID)+(1|NO), ziformula = ~ EV+A+SEX+EL+SP+(1|CD)+(1|ID)+(1|NO), data = tempdata ,family = 'nbinom2',doFit=TRUE))
  
# extract confidence intervals and estimates   
  try(temp <- c(confint(tempEV)[2,1],tempEV$fit$par[2],confint(tempEV)[2,2],confint(tempS)[2,1],tempS$fit$par[2],confint(tempS)[2,2],confint(tempD)[2,1],tempD$fit$par[2],confint(tempD)[2,2]))
  try(temp)

} 
######################work until here
################################################
#################################################
randomaf <- function(i){
  
  tempdataaf <- data.frame(EPG = SBDaf$EPG, ID = SBDaf$ID, A = SBDaf$A, EL = SBDaf$EL, EV = rep(0, length(SBDaf$ID)), D = rep(0, length(SBDaf$ID)), S = rep(0, length(SBDaf$ID)), CD = SBDaf$CD, SP = SBDaf$SP, NO = SBDaf$NO)
  
  tempnetaf <- rewire(netAF, with = each_edge(prob = runif(1,0,1), loops = FALSE, multiple = TRUE))
  
  D<-degree(tempnetaf, mode="all")
  S<-strength(tempnetaf,mode="all")
  EV<-eigen_centrality(tempnetaf, directed=FALSE)$vector
  
  Daf <- data.frame(ID = rownames(NetDataAF), D = D)
  Saf <- data.frame(ID = rownames(NetDataAF), S = S)
  EVaf <- data.frame(ID = rownames(NetDataAF), EV = EV)
  
  tempdataaf$EV <- scale(EVaf$EV[as.numeric(tempdataaf$ID)])
  tempdataaf$S <- scale(Saf$S[as.numeric(tempdataaf$ID)])
  tempdataaf$D <- scale(Daf$D[as.numeric(tempdataaf$ID)])
  
  try(tempafD <- glmmTMB(EPG ~ D+A+EL+SP+(1|CD)+(1|ID)+(1|NO), ziformula = ~ D+A+EL+SP+(1|CD)+(1|ID)+(1|NO), data = tempdataaf, family = 'nbinom2',doFit=TRUE))
  try(tempafS <- glmmTMB(EPG ~ S+A+EL+SP+(1|CD)+(1|ID)+(1|NO), ziformula = ~ S+A+EL+SP+(1|CD)+(1|ID)+(1|NO), data = tempdataaf, family ='nbinom2',doFit=TRUE))
  try(tempafEV <- glmmTMB(EPG ~ EV+A+EL+SP+(1|CD)+(1|ID)+(1|NO), ziformula = ~ EV+A+EL+SP+(1|CD)+(1|ID)+(1|NO), data = tempdataaf ,family = 'nbinom2',doFit=TRUE))
  
  try(tempaf <- c(confint(tempafEV)[2,1],tempafEV$fit$par[2],confint(tempafEV)[2,2],confint(tempafS)[2,1],tempafS$fit$par[2],confint(tempafS)[2,2],confint(tempafD)[2,1],tempafD$fit$par[2],confint(tempafD)[2,2]))
  
  try(tempaf)
} 

randomjuv <- function(i){
 
  tempdatajuv <- data.frame(EPG = SBDjuv$EPG, ID = SBDjuv$ID, A = SBDjuv$A, SEX = SBDjuv$SEX, EL = SBDjuv$EL, EV = rep(0, length(SBDjuv$ID)), D = rep(0, length(SBDjuv$ID)), S = rep(0, length(SBDjuv$ID)), CD = SBDjuv$CD, SP = SBDjuv$SP, NO = SBDjuv$NO)
  
  tempnetjuv <- rewire(netJUV, with = each_edge(prob = runif(1,0,1), loops = FALSE, multiple = TRUE))
  
  D<-degree(tempnetjuv, mode="all")
  S<-strength(tempnetjuv,mode="all")
  EV<-eigen_centrality(tempnetjuv, directed=FALSE)$vector
 
  Djuv <- data.frame(ID = rownames(NetDataJUV), D = D)
  Sjuv <- data.frame(ID = rownames(NetDataJUV), S = S)
  EVjuv <- data.frame(ID = rownames(NetDataJUV), EV = EV)
  
  Djuv <- droplevels(Djuv[-(which(Djuv$ID %in% NS)),])
  Sjuv <- droplevels(Sjuv[-(which(Sjuv$ID %in% NS)),])
  EVjuv <- droplevels(EVjuv[-(which(EVjuv$ID %in% NS)),])
  
  tempdatajuv$D <- scale(Djuv$D[as.numeric(tempdatajuv$ID)])
  tempdatajuv$S <- scale(Sjuv$S[as.numeric(tempdatajuv$ID)])
  tempdatajuv$EV <- scale(EVjuv$EV[as.numeric(tempdatajuv$ID)])
  
  try(tempjuvD <- glmmTMB(EPG ~ D+A+SEX+EL+SP+(1|CD)+(1|ID)+(1|NO), ziformula = ~ D+A+SEX+EL+SP+(1|CD)+(1|ID)+(1|NO), data = tempdatajuv, family = 'nbinom2',doFit=TRUE))
  try(tempjuvS <- glmmTMB(EPG ~ S+A+SEX+EL+SP+(1|CD)+(1|ID)+(1|NO), ziformula = ~ S+A+SEX+EL+SP+(1|CD)+(1|ID)+(1|NO), data = tempdatajuv, family ='nbinom2',doFit=TRUE))
  try(tempjuvEV <- glmmTMB(EPG ~ EV+A+SEX+EL+SP+(1|CD)+(1|ID)+(1|NO), ziformula = ~ EV+A+SEX+EL+SP+(1|CD)+(1|ID)+(1|NO), data = tempdatajuv ,family = 'nbinom2',doFit=TRUE))
  
  try(tempjuv <- c(confint(tempjuvEV)[2,1],tempjuvEV$fit$par[2],confint(tempjuvEV)[2,2],confint(tempjuvS)[2,1],tempjuvS$fit$par[2],confint(tempjuvS)[2,2],confint(tempjuvD)[2,1],tempjuvD$fit$par[2],confint(tempjuvD)[2,2]))
  
  try(tempjuv)
} 


# Start the randomisation procedure
T1 <- Sys.time()
perm <- 3 

# Set up progress bar
pb <- txtProgressBar(max = perm, style = 3)
progress <- function(n)  setTxtProgressBar(pb, n)
opts <- list(progress = progress)
 
# Split cores to work in parallel if troubles
core <- detectCores() 
cl<- makeCluster(core-1)   
registerDoParallel(cl)
registerDoSNOW(cl)

# Run the randomisation
RDres <- foreach(i = 1:perm,.packages = c('glmmTMB','igraph'), .options.snow = opts) %dopar% randomall (i)  
RDresaf <- foreach(i = 1:perm,.packages = c('glmmTMB','igraph'), .options.snow = opts) %dopar% randomaf (i)  
RDresjuv <- foreach(i = 1:perm,.packages = c('glmmTMB','igraph'), .options.snow = opts) %dopar% randomjuv (i)  

# Save the results
RDres <- as.data.frame(do.call(rbind, RDres))
RDresaf <- as.data.frame(do.call(rbind, RDresaf))
RDresjuv <- as.data.frame(do.call(rbind, RDresjuv))

write.csv(na.omit(RDres),file = "RDres.csv")
write.csv(na.omit(RDresaf),file = "RDresaf.csv")
write.csv(na.omit(RDresjuv),file = "RDresjuv.csv")

# End of process
stopImplicitCluster()
stopCluster(cl) 
T2 <- Sys.time()
T2-T1
