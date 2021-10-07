# This script is for running models

library(fitdistrplus)
library(glmmTMB)
library(performance)
library(car)
library(effects)
library(multcomp)
library(MuMIn)
library(DHARMa)
library(ggplot2)
library(igraph)

# Load the data 

SampleData <- read.csv("raw_data/SampleData.csv",  sep = ";") # be mindful of separator type
IndividualData <- read.csv("raw_data/IndData.csv",  sep = ";") # be mindful of separator type
NetData <- read.csv("raw_data/NetData.csv",  sep = ";") # be mindful of separator type
# note for variable names: NO = sample ID,	ID = individual ID,	CD = sample collection day,	D	= degree, S	= strength, EV	= eignevector, A = age,	AG = age category,	SEX	= sex, EL = Elo-rating,	FOA = number of focal scans,	FSA	= number of focal sample, EPG = egg count per gram of feces,	SP =  worm species

# Make a list of individuals that do not have fecal sample (useful for later)
NS <- c('takana','mushi','neji','uso')

# "Clean" the data if necessary, delete first explanatory row and transform variables with their correct identity (e.g. "character" into "numeric")
NetData <- NetData[,-1]
NetData <- as.matrix(NetData)
is.matrix(NetData)
rownames(NetData)<- colnames(NetData)
IBD <- IndividualData
SBD <- SampleData

# use these to reset the colnames if the first letter became garbled
header1 <- as.character(read.csv("raw_data/SampleData.csv", header = FALSE, fileEncoding="UTF-8-BOM", sep = ";")[1, ])
header2 <- as.character(read.csv("raw_data/IndData.csv", header = FALSE, fileEncoding="UTF-8-BOM", sep = ";")[1, ])
colnames(SBD) <- header1
colnames(IBD) <- header2
SBD[, c("D", "S", "EV", "A", "EL")] <- apply(SBD[, c("D", "S", "EV", "A", "EL")], 2, function(x)as.numeric(scale(x)))
SBD$EPG <- round(SBD$EPG)
SBD$CD<-as.factor(SBD$CD)
SBD$ID<-as.factor(SBD$ID)
SBD$SP<-as.factor(SBD$Species)

# EPG is the response variable of interest. It is count data and contains a lot of 0s. 
hist(SBD$EPG)
sum(SBD$EPG==0)

# check the distribution of the response variable
# plot the possible distributions
descdist(SBD$EPG, discrete = TRUE)
# most likely distribution fit = negative binomial
# because we have zero-inflated data, we will also use a model fit that takes it into account


## WHOLE-GROUP ##

# run the models with 3 centrality measures separately but the EPG from the 3 parasite species together, variables are scaled as well for easier modelling and comparability (see ms text)

# with degree
glmmD <- glmmTMB(EPG ~ D+A+SEX+EL+(SP)+(1|CD)+(1|ID)+(1|NO), ziformula = ~ D+A+SEX+EL+(SP)+(1|CD)+(1|ID)+(1|NO), data = SBD, family = 'nbinom2',doFit=TRUE)
# with strength
glmmS <- glmmTMB(EPG ~ S+A+SEX+EL+(SP)+(1|CD)+(1|ID)+(1|NO), ziformula = ~ S+A+SEX+EL+(SP)+(1|CD)+(1|ID)+(1|NO), data = SBD, family ='nbinom2',doFit=TRUE)
# with eignevector
glmmEV <- glmmTMB(EPG ~ EV+A+SEX+EL+(SP)+(1|CD)+(1|ID)+(1|NO), ziformula = ~ EV+A+SEX+EL+(SP)+(1|CD)+(1|ID)+(1|NO), data = SBD ,family = 'nbinom2',doFit=TRUE)


# Check overdispersion
testDispersion(glmmD)
testDispersion(glmmS)
testDispersion(glmmEV)

# Check residual distribution
sglmmD <- simulateResiduals(fittedModel = glmmD, plot = T)
sglmmS <- simulateResiduals(fittedModel = glmmS, plot = T)
sglmmEV <- simulateResiduals(fittedModel = glmmEV, plot = T)

sum(residuals(sglmmD)<0.5)/sum(residuals(sglmmD)<10)
sum(residuals(sglmmS)<0.5)/sum(residuals(sglmmS)<10)
sum(residuals(sglmmEV)<0.5)/sum(residuals(sglmmEV)<10)

# Check for multicollinearity
check_collinearity(glmmD)
check_collinearity(glmmS)
check_collinearity(glmmEV)
plot(check_collinearity(glmmEV))

# Run further diagnostics (e.g. Cook's distance and leverage) with an adaptation of the "influence.me" package for glmmTMB (time consuming to run!)

source(system.file("other_methods","influence_mixed.R", package="glmmTMB"))

glmmDinf <- influence_mixed(glmmD,groups = "ID")
glmmSinf <- influence_mixed(glmmS,groups = "ID")
glmmEVinf <- influence_mixed(glmmEV,groups = "ID")

infIndexPlot(glmmDinf)
infIndexPlot(glmmSinf)
infIndexPlot(glmmEVinf)

infD <- as.data.frame(glmmDinf[["fixed.effects[-ID]"]])
infS <- as.data.frame(glmmSinf[["fixed.effects[-ID]"]])
infEV <- as.data.frame(glmmEVinf[["fixed.effects[-ID]"]])

infD <- transform(infD,
                  id=rownames(infD),
                  cooks=cooks.distance(glmmDinf))
infD$ord <- rank(infD$cooks)

infS <- transform(infS,
                  id=rownames(infS),
                  cooks=cooks.distance(glmmSinf))
infS$ord <- rank(infS$cooks)

infEV <- transform(infEV,
                   id=rownames(infEV),
                   cooks=cooks.distance(glmmEVinf))
infEV$ord <- rank(infEV$cooks)


if (require(reshape2)) {
  infD_long <- melt(infD, id.vars=c("ord","id"))
  gg_inflD <- (ggplot(infD_long,aes(ord,value))
               + geom_point()
               + facet_wrap(~variable, scale="free_y")
               + scale_x_reverse(expand=expansion(mult=0.15))
               + scale_y_continuous(expand=expansion(mult=0.15))
               + geom_text(data=subset(infD_long,ord>24),
                           aes(label=id),vjust=-1.05)
  )
  print(gg_inflD)
}
if (require(reshape2)) {
  infS_long <- melt(infS, id.vars=c("ord","id"))
  gg_inflS <- (ggplot(infS_long,aes(ord,value))
               + geom_point()
               + facet_wrap(~variable, scale="free_y")
               + scale_x_reverse(expand=expansion(mult=0.15))
               + scale_y_continuous(expand=expansion(mult=0.15))
               + geom_text(data=subset(infS_long,ord>24),
                           aes(label=id),vjust=-1.05)
  )
  print(gg_inflS)
}
if (require(reshape2)) {
  infEV_long <- melt(infEV, id.vars=c("ord","id"))
  gg_inflEV <- (ggplot(infEV_long,aes(ord,value))
                + geom_point()
                + facet_wrap(~variable, scale="free_y")
                + scale_x_reverse(expand=expansion(mult=0.15))
                + scale_y_continuous(expand=expansion(mult=0.15))
                + geom_text(data=subset(infEV_long,ord>24),
                            aes(label=id),vjust=-1.05)
  )
  print(gg_inflEV)
}

# Run reduced null models with all variables except centrality measures
glmmnull <- glmmTMB(EPG ~ A+SEX+EL+(SP)+(1|CD)+(1|ID)+(1|NO), ziformula = ~ A+SEX+EL+(SP)+(1|CD)+(1|ID)+(1|NO), data = SBD ,family = 'nbinom2',doFit=TRUE)

# Check the full models against the reduced null model with Likelihood ratio tests
anova(glmmD, glmmnull, test="Chisq")
anova(glmmS, glmmnull, test="Chisq")
anova(glmmEV, glmmnull, test="Chisq")

# Summary of the models
summary(glmmD)
summary(glmmS)
summary(glmmEV)

# Plot the results
plot(allEffects(glmmD))
plot(allEffects(glmmS))
plot(allEffects(glmmEV))

# averaging the models to have a overall estimate of shared factors

glmm <- model.avg(glmmD,glmmS,glmmEV)
summary(glmm)


## ADULT FEMALES ##

# subset the data set to select adult females only 
SampleData <- read.csv("raw_data/SampleData.csv",  sep = ";") # be mindful of separator type
SBD <- SampleData

# use these to reset the colnames if the first letter became garbled
header1 <- as.character(read.csv("raw_data/SampleData.csv", header = FALSE, fileEncoding="UTF-8-BOM", sep = ";")[1, ])
colnames(SBD) <- header1

SBDaf<-droplevels(subset(droplevels(subset(SBD, AG=="adult")),SEX=="female"))
head(SBDaf)
IBD_AF <- droplevels(subset(droplevels(subset(IBD, AG=="adult")),SEX=="female"))
SBDaf$EPG <- round(SBDaf$EPG)
SBDaf$CD<-as.factor(SBDaf$CD)
SBDaf$ID<-as.factor(SBDaf$ID)
SBDaf$SP<-as.factor(SBDaf$Species)
#subset the proximity matrix to have an adult female only matrix
sam <- IBD_AF[,1]
NetDataAF<- NetData[(which(rownames(NetData) %in% sam)),(which(colnames(NetData) %in% sam))]

# create a network object
netAF <- graph.adjacency(NetDataAF, mode = "undirected", weighted = TRUE, diag = FALSE)

#calculate network measures
#degree
degAF<-degree(netAF, mode="all")
#strength
strengthAF<-strength(netAF,mode="all")
#eigenvector
evcentAF<-evcent(netAF)$vector

# Put the new centrality measures in a dataframe
degAF <- data.frame(ID = as.factor(names(degAF)), D = degAF)
strengthAF <- data.frame(ID = as.factor(names(strengthAF)), S = strengthAF)
evcentAF <- data.frame(ID = as.factor(names(evcentAF)), EV = evcentAF)



# Check if IDs match between network measure dataframe and analysis dataframe 
if(!identical(levels(evcentAF$ID), levels(droplevels(SBDaf$ID)))) stop("something wrong with factor levels in adult-only")

# If yes, then put the new centrality measures back into the analysis dataframe 

SBDaf$EVaf <- scale(evcentAF$EV[as.numeric(SBDaf$ID)])
SBDaf$Saf <- scale(strengthAF$S[as.numeric(SBDaf$ID)])
SBDaf$Daf <- scale(degAF$D[as.numeric(SBDaf$ID)])
SBDaf[, c( "A", "EL")] <- apply(SBDaf[, c( "A", "EL")], 2, function(x)as.numeric(scale(x)))

# Run the GLMMs
# with degree
glmmAF_D <- glmmTMB(EPG ~ Daf+A+EL+(SP)+(1|CD)+(1|ID)+(1|NO), ziformula = ~ Daf+A+EL+(SP)+(1|CD)+(1|ID)+(1|NO), data = SBDaf, family = 'nbinom2',doFit=TRUE)
# with strength
glmmAF_S <- glmmTMB(EPG ~ Saf+A+EL+(SP)+(1|CD)+(1|ID)+(1|NO), ziformula = ~ Saf+A+EL+(SP)+(1|CD)+(1|ID)+(1|NO), data = SBDaf, family ='nbinom2',doFit=TRUE)
# with eigenvector
glmmAF_EV <- glmmTMB(EPG ~ EVaf+A+EL+(SP)+(1|CD)+(1|ID)+(1|NO), ziformula = ~ EVaf+A+EL+(SP)+(1|CD)+(1|ID)+(1|NO), data = SBDaf ,family = 'nbinom2',doFit=TRUE)

# Run all the checks and diagnostics as before (not repeated here, copy and paste code from above)

# Run reduced null models with all variables except centrality measures
glmmAF_null <- glmmTMB(EPG ~ A+EL+(SP)+(1|CD)+(1|ID)+(1|NO), ziformula = ~ A+EL+(SP)+(1|CD)+(1|ID)+(1|NO), data = SBDaf ,family = 'nbinom2',doFit=TRUE)

# Check the full models against the reduced null model with Likelihood ratio tests
anova(glmmAF_D, glmmAF_null, test="Chisq")
anova(glmmAF_S, glmmAF_null, test="Chisq")
anova(glmmAF_EV, glmmAF_null, test="Chisq")

# Summary of the models
summary(glmmAF_D)
summary(glmmAF_S)
summary(glmmAF_EV)

# Plot the results
plot(allEffects(glmmAF_D))
plot(allEffects(glmmAF_S))
plot(allEffects(glmmAF_EV))

# averaging the models to have a overall estimate of shared factors

glmmAF <- model.avg(glmmAF_D,glmmAF_S,glmmAF_EV)
summary(glmmAF)

## JUVENILES ##

# subset the data set to select juveniles only 

SampleData <- read.csv("raw_data/SampleData.csv",  sep = ";") # be mindful of separator type
SBD <- SampleData

# use these to reset the colnames if the first letter became garbled
header1 <- as.character(read.csv("raw_data/SampleData.csv", header = FALSE, fileEncoding="UTF-8-BOM", sep = ";")[1, ])
colnames(SBD) <- header1

SBDjuv<-droplevels(subset(SBD, AG=="juvenile"))
head(SBDjuv)
IBD_juv <-droplevels(subset(IBD, AG=="juvenile"))
SBDjuv$EPG <- round(SBDjuv$EPG)
SBDjuv$CD<-as.factor(SBDjuv$CD)
SBDjuv$ID<-as.factor(SBDjuv$ID)
SBDjuv$SP<-as.factor(SBDjuv$Species)

#subset the proximity matrix to have an juveniles only matrix
sam <- IBD_juv[,1]
NetDataJUV<- NetData[(which(rownames(NetData) %in% sam)),(which(colnames(NetData) %in% sam))]

# create a network object
netJUV <- graph.adjacency(NetDataJUV, mode = "undirected", weighted = TRUE, diag = FALSE)

#calculate network measures
#degree
degJUV<-degree(netJUV, mode="all")
#strength
strengthJUV<-strength(netJUV,mode="all")
#eigenvector
evcentJUV<-evcent(netJUV)$vector

# Put the new centrality measures in a dataframe
degJUV <- data.frame(ID = as.factor(names(degJUV)), D = degJUV)
degJUV <- droplevels(degJUV[-(which(degJUV$ID %in% NS)),])
strengthJUV <- data.frame(ID = as.factor(names(strengthJUV)), S = strengthJUV)
strengthJUV <- droplevels(strengthJUV[-(which(strengthJUV$ID %in% NS)),])
evcentJUV <- data.frame(ID = as.factor(names(evcentJUV)), EV = evcentJUV)
evcentJUV <- droplevels(evcentJUV[-(which(evcentJUV$ID %in% NS)),])

# Check if IDs match between network measure dataframe and analysis dataframe 
if(!identical(levels(evcentJUV$ID), levels(droplevels(SBDjuv$ID)))) stop("something wrong with factor levels in adult-only")

# If yes, then put the new centrality measures back into the analysis dataframe 
SBDjuv$Djuv <- scale(degJUV$D[as.numeric(SBDjuv$ID)])
SBDjuv$Sjuv <- scale(strengthJUV$S[as.numeric(SBDjuv$ID)])
SBDjuv$EVjuv <- scale(evcentJUV$EV[as.numeric(SBDjuv$ID)])
SBDjuv[, c( "A", "EL")] <- apply(SBDjuv[, c( "A", "EL")], 2, function(x)as.numeric(scale(x)))

# Run the GLMMs
# with degree
glmmJUV_D <- glmmTMB(EPG ~ Djuv+A+SEX+EL+SP+(1|CD)+(1|ID)+(1|NO), ziformula = ~ Djuv+A+SEX+EL+SP+(1|CD)+(1|ID)+(1|NO), data = SBDjuv, family = 'nbinom2',doFit=TRUE)
# with strength
glmmJUV_S <- glmmTMB(EPG ~ Sjuv+A+SEX+EL+SP+(1|CD)+(1|ID)+(1|NO), ziformula = ~ Sjuv+A+SEX+EL+SP+(1|CD)+(1|ID)+(1|NO), data = SBDjuv, family ='nbinom2',doFit=TRUE)
# with eigenvector
glmmJUV_EV <- glmmTMB(EPG ~ EVjuv+A+SEX+EL+SP+(1|CD)+(1|ID)+(1|NO), ziformula = ~ EVjuv+A+SEX+EL+SP+(1|CD)+(1|ID)+(1|NO), data = SBDjuv ,family = 'nbinom2',doFit=TRUE)

# Run all the checks and diagnostics as before (not repeated here, copy and paste code from above)

# Run reduced null models with all variables except centrality measures
glmmJUV_null <- glmmTMB(EPG ~ A+SEX+EL+SP+(1|CD)+(1|ID)+(1|NO), ziformula = ~ A+SEX+EL+SP+(1|CD)+(1|ID)+(1|NO), data = SBDjuv ,family = 'nbinom2',doFit=TRUE)

# Check the full models against the reduced null model with Likelihood ratio tests
anova(glmmJUV_D, glmmJUV_null, test="Chisq")
anova(glmmJUV_S, glmmJUV_null, test="Chisq")
anova(glmmJUV_EV, glmmJUV_null, test="Chisq")

# Summary of the models
summary(glmmJUV_D)
summary(glmmJUV_S)
summary(glmmJUV_EV)

# Plot the results
plot(allEffects(glmmJUV_D))
plot(allEffects(glmmJUV_S))
plot(allEffects(glmmJUV_EV))

# averaging the models to have a overall estimate of shared factors

glmmJUV <- model.avg(glmmJUV_D,glmmJUV_S,glmmJUV_EV)
summary(glmmJUV)