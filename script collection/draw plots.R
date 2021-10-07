# This script details the plotting used in the paper

## NETWORK PLOT (partial credits @ Christof Neumann)

library(igraph)

# Function for ploting multiple networks together for comparison (credits @ Christof Neumann)
plotfoo <- function(pdata = pdata,
                    graphobject = net,
                    plist = plist,
                    maxradius = 10,
                    innerrad = 5,
                    lim = 12,
                    nind = nrow(pdata),
                    cexfac = 3,
                    cols = hcl.colors(length(plist) + 1, "zissou1")) {
# Assign the location of each node in order to get a network in perfect circle
    pdata$xloc <- maxradius * cos(seq(2/nind * pi, 2*pi, length.out = nind))
  pdata$yloc <- maxradius * sin(seq(2/nind * pi, 2*pi, length.out = nind))
# Set the size of nodes
      pdata$cexval <- sqrt(pdata$eigenrank / pi * cexfac)
# Plot the most out ring of nodes 
   plot(0, 0, type="n", 
         xlim = c(-lim, lim), ylim = c(-5*lim/3, 5*lim/3),
         asp = 1, ann = FALSE, axes = FALSE)
# Draw edges
  edata <- as_long_data_frame(graphobject)
  apply(edata[, 1:2], 1, function(x)segments(x0 = pdata$xloc[x[1]], y0 = pdata$yloc[x[1]], x1 = pdata$xloc[x[2]], y1 = pdata$yloc[x[2]], lwd = 0.5))
  
  points(pdata$xloc, pdata$yloc, cex = pdata$cexval, bg = cols[1], pch = 21)
# Do the same thing for nodes from inner ring networks  
  i=1
  for (i in seq_len(length(plist))) {
    r <- seq(maxradius, innerrad, length.out = length(plist) + 1)[-1][i]
    pd <- pdata
    rvals <- seq(2/nind * pi, 2*pi, length.out = nind)
    pd$xloc <- r * cos(rvals)
    pd$yloc <- r * sin(rvals)
    
    idsubset <- plist[[i]]$ids
    pd <- pd[pd$ids %in% idsubset, ]
    pd[plist[[i]]$ids, "eigenrank"] <- plist[[i]]$eigenrank
    pd$cexval <- sqrt(pd$eigenrank / pi * cexfac)
    
    points(pd$xloc, pd$yloc, pch = 21, bg = cols[i + 1], cex = pd$cexval)
  }
}

# Load the data (if necessary again)
# "Clean" the data, delete first explanatory row and transform variables with their correct identity (e.g. "character" into "numeric")

IndividualData <- read.csv("raw_data/IndData.csv",  sep = ";") # be mindful of separator type
NetData <- read.csv("raw_data/NetData.csv",  sep = ";") # be mindful of separator type

NS <- c('takana','mushi','neji','uso')

NetData <- NetData[,-1]
NetData <- as.matrix(NetData)
is.matrix(NetData)
rownames(NetData)<- colnames(NetData)
IBD <- IndividualData
header2 <- as.character(read.csv("raw_data/IndData.csv", header = FALSE, fileEncoding="UTF-8-BOM", sep = ";")[1, ])
colnames(IBD) <- header2


# Subset the proximity matrix to have an adult female only matrix
IBD_AF <- droplevels(subset(droplevels(subset(IBD, AG=="adult")),SEX=="female"))
sam <- IBD_AF$ID
NetDataAF<- NetData[(which(rownames(NetData) %in% sam)),(which(colnames(NetData) %in% sam))]

# Subset the proximity matrix to have an adult female only matrix
IBD_juv <-droplevels(subset(IBD, AG=="juvenile"))
sam <- IBD_juv$ID
NetDataJUV<- NetData[(which(rownames(NetData) %in% sam)),(which(colnames(NetData) %in% sam))]

# Create a list of names to exclude from the list of all individuals of interest
sam1 = sample(IBD$ID,3)
sam2 = sample(IBD$ID,5)
sam3 = sample(IBD$ID,13)
sam4 = sample(IBD$ID,26)

# Exclude the individuals from the matrix
NetData5 = NetData[-(which(rownames(NetData) %in% sam1)),-(which(colnames(NetData) %in% sam1))]
NetData10 = NetData[-(which(rownames(NetData) %in% sam2)),-(which(colnames(NetData) %in% sam2))]
NetData25 = NetData[-(which(rownames(NetData) %in% sam3)),-(which(colnames(NetData) %in% sam3))]
NetData50 = NetData[-(which(rownames(NetData) %in% sam4)),-(which(colnames(NetData) %in% sam4))]

# Create the network
NetALL <- igraph::graph.adjacency(NetData,mode= "undirected",weighted = TRUE, diag = FALSE)
NetAF <- igraph::graph.adjacency(NetDataAF,mode= "undirected",weighted = TRUE, diag = FALSE)
NetJUV <- igraph::graph.adjacency(NetDataJUV,mode= "undirected",weighted = TRUE, diag = FALSE)
Net5 <- igraph::graph.adjacency(NetData5,mode= "undirected",weighted = TRUE, diag = FALSE)
Net10 <- igraph::graph.adjacency(NetData10,mode= "undirected",weighted = TRUE, diag = FALSE)
Net25 <- igraph::graph.adjacency(NetData25,mode= "undirected",weighted = TRUE, diag = FALSE)
Net50 <- igraph::graph.adjacency(NetData50,mode= "undirected",weighted = TRUE, diag = FALSE)

# Create data frame for the most out ring network
NetALLR <- data.frame(ids = colnames(NetData), eigencent = eigen_centrality(NetALL)$vector)
NetALLR$eigenrank <- rank(NetALLR$eigencent )

# Exclude ranks from individuals that are excluded
nub5<-NetALLR[which(NetALLR[,1] %in% sam1),3]
nub10<-NetALLR[which(NetALLR[,1] %in% sam2),3]
nub25<-NetALLR[which(NetALLR[,1] %in% sam3),3]
nub50<-NetALLR[which(NetALLR[,1] %in% sam4),3]
nubjuv<-NetALLR[which(NetALLR[,1] %in% colnames(NetDataJUV)),3]
nubaf<-NetALLR[which(NetALLR[,1] %in% colnames(NetDataAF)),3]

rk5 <-c(1:52)[-nub5] 
rk10 <-c(1:52)[-nub10] 
rk25 <-c(1:52)[-nub25] 
rk50 <-c(1:52)[-nub50] 
rkjuv <-c(1:52)[-nubjuv] 
rkaf <-c(1:52)[-nubaf] 

# Create data frame for inner ring networks
NetAFR <- data.frame(ids = colnames(NetDataAF),eigenrank = rank(eigen_centrality(NetAF)$vector))
NetJUVR <- data.frame(ids = colnames(NetDataJUV),eigenrank = rank(eigen_centrality(NetJUV)$vector))
Net5R <- data.frame(ids = colnames(NetData5),eigenrank = rank(eigen_centrality(Net5)$vector))
Net10R <- data.frame(ids = colnames(NetData10),eigenrank = rank(eigen_centrality(Net10)$vector))
Net25R <- data.frame(ids = colnames(NetData25),eigenrank = rank(eigen_centrality(Net25)$vector))
Net50R <- data.frame(ids = colnames(NetData50),eigenrank = rank(eigen_centrality(Net50)$vector))

# Re-rank the inner ring networks with excluded rank list, so that it's comparable to the original network.
NetAFR$eigenrank<-rkaf[NetAFR$eigenrank]
NetJUVR$eigenrank<-rkjuv[NetJUVR$eigenrank]
Net5R$eigenrank<-rk5[Net5R$eigenrank]
Net10R$eigenrank<-rk10[Net10R$eigenrank]
Net25R$eigenrank<-rk25[Net25R$eigenrank]
Net50R$eigenrank<-rk50[Net50R$eigenrank]

# Create the inner ring netlist
NetList <- list(Net5R,Net10R,Net25R,Net50R,NetJUVR,NetAFR)

# Plot the multiple network.    
# "pdata" for the original network data, "graphobject" for the original network object (which edges come from), "plist" for the network list, "maxradius" for the outest ring radius, "innerrad" for the innest ring radius, "cexfac" for the node size, "lim" for the ploting area size.
plot<- plotfoo(pdata = NetALLR, graphobject = NetALL, plist = NetList, maxradius = 50, innerrad = 22, cexfac = 2, lim=30)

## GLMM RESULTS PLOTS (credits @ Christof Neumann)
# scripts needed: "run GLMMs"
# each variable gets its own plot

library(effects)
require(MuMIn)
require(glmmTMB)

# if needed load the dataset and the network matrix again 

SBD <- read.csv("raw_data/SampleData.csv",  sep = ";") # be mindful of separator type
header1 <- as.character(read.csv("raw_data/SampleData.csv", header = FALSE, fileEncoding="UTF-8-BOM", sep = ";")[1, ])
colnames(SBD) <- header1
SBD[, c("D", "S", "EV", "A", "EL")] <- apply(SBD[, c("D", "S", "EV", "A", "EL")], 2, function(x)as.numeric(scale(x)))
SBD$EPG <- round(SBD$EPG)
SBD$CD<-as.factor(SBD$CD)
SBD$ID<-as.factor(SBD$ID)
SBD$SP<-as.factor(SBD$Species)

# if needed run the models again
glmmD <- glmmTMB(EPG ~ D+A+SEX+EL+(SP)+(1|CD)+(1|ID)+(1|NO), ziformula = ~ D+A+SEX+EL+(SP)+(1|CD)+(1|ID)+(1|NO), data = SBD, family = 'nbinom2',doFit=TRUE)
glmmS <- glmmTMB(EPG ~ S+A+SEX+EL+(SP)+(1|CD)+(1|ID)+(1|NO), ziformula = ~ S+A+SEX+EL+(SP)+(1|CD)+(1|ID)+(1|NO), data = SBD, family ='nbinom2',doFit=TRUE)
glmmEV <- glmmTMB(EPG ~ EV+A+SEX+EL+(SP)+(1|CD)+(1|ID)+(1|NO), ziformula = ~ EV+A+SEX+EL+(SP)+(1|CD)+(1|ID)+(1|NO), data = SBD ,family = 'nbinom2',doFit=TRUE)


# Create the dataframe for Effect
pdataD <- data.frame(Effect("D", glmmD, xlevels = list(D = seq(min(SBD$D), max(SBD$D), length.out = 51))))
pdataS <- data.frame(Effect("S", glmmS, xlevels = list(S = seq(min(SBD$S), max(SBD$S), length.out = 51))))
pdataEV <- data.frame(Effect("EV", glmmEV, xlevels = list(EV = seq(min(SBD$EV), max(SBD$EV), length.out = 51))))
pdataA1 <- data.frame(Effect("A", glmmD, xlevels = list(A = seq(min(SBD$A), max(SBD$A), length.out = 51))))
pdataA2 <- data.frame(Effect("A", glmmS, xlevels = list(A = seq(min(SBD$A), max(SBD$A), length.out = 51))))
pdataA3 <- data.frame(Effect("A", glmmEV, xlevels = list(A = seq(min(SBD$A), max(SBD$A), length.out = 51))))
pdataEL1 <- data.frame(Effect("EL", glmmD, xlevels = list(EL = seq(min(SBD$EL), max(SBD$EL), length.out = 51))))
pdataEL2 <- data.frame(Effect("EL", glmmS, xlevels = list(EL = seq(min(SBD$EL), max(SBD$EL), length.out = 51))))
pdataEL3 <- data.frame(Effect("EL", glmmEV, xlevels = list(EL = seq(min(SBD$EL), max(SBD$EL), length.out = 51))))

# Turn sex into binary numeric variable 1/0, edit text later
pdataSEX1 <- data.frame(Effect("SEX", glmmD,xlevels = list(SBD$SEX)))
pdataSEX11 <- data.frame("SEX" = numeric(2),"fit" = numeric (2),"se" = numeric (2),"lower" = numeric (2),"upper" = numeric (2))
pdataSEX11[,2:5] <- pdataSEX1[,2:5]
pdataSEX11[2,1] <- 1
pdataSEX2 <- data.frame(Effect("SEX", glmmS, xlevels = list(SBD$SEX)))
pdataSEX21 <- data.frame("SEX" = numeric(2),"fit" = numeric (2),"se" = numeric (2),"lower" = numeric (2),"upper" = numeric (2))
pdataSEX21[,2:5] <- pdataSEX2[,2:5]
pdataSEX21[2,1] <- 1
pdataSEX3 <- data.frame(Effect("SEX", glmmEV, xlevels = list(SBD$SEX)))
pdataSEX31 <- data.frame("SEX" = numeric(2),"fit" = numeric (2),"se" = numeric (2),"lower" = numeric (2),"upper" = numeric (2))
pdataSEX31[,2:5] <- pdataSEX3[,2:5]
pdataSEX31[2,1] <- 1


require(ggplot2)
require(gridExtra)

plotD <- ggplot()+
  geom_point(aes(SBD$D, SBD$EPG), pch =16, color = "black", size =5,alpha=0.6)+
  geom_polygon(aes(x = c(pdataD$D, rev(pdataD$D)), y = c(pdataD$lower, rev(pdataD$upper))), fill =  hcl.colors(5, "zissou1", alpha = 0.4)[3])+
  geom_line(aes(pdataD$D, pdataD$fit), lwd = 1.5, color = hcl.colors(5, "zissou1", alpha = 1)[4])+
  theme_bw()+
  ylim(c(0,25000))+
  xlim(c(-3,2))+
  theme(panel.grid.major.x = element_blank(),axis.line = element_line(color = "black",size=1.5),axis.ticks.length = unit(.5,"cm"),axis.ticks = element_line(color = "black",size = 1.5))+
  theme(panel.grid = element_blank(),axis.title = element_blank(),axis.text.y = element_text(size=25),axis.text.x = element_text(size=25))
     

plotS <- ggplot()+
  geom_point(aes(SBD$S, SBD$EPG), pch =16, color = "black", size =5,alpha=0.6)+
  geom_polygon(aes(x = c(pdataS$S, rev(pdataS$S)), y = c(pdataS$lower, rev(pdataS$upper))), fill =  hcl.colors(5, "zissou1", alpha = 0.4)[3])+
  geom_line(aes(pdataS$S, pdataS$fit), lwd = 1.5, color = hcl.colors(5, "zissou1", alpha = 1)[4])+
  theme_bw()+
  ylim(c(0,25000))+
  xlim(c(-3,3))+
  theme(panel.grid.major.x = element_blank(),axis.line = element_line(color = "black",size=1.5),axis.ticks.length = unit(.5,"cm"),axis.ticks = element_line(color = "black",size = 1.5))+
  theme(panel.grid = element_blank(),axis.title = element_blank(),axis.text.y = element_text(size=25),axis.text.x = element_text(size=25))



plotEV <- ggplot()+
  geom_point(aes(SBD$EV, SBD$EPG), pch =16, color = "black", size =5,alpha=0.6)+
  geom_polygon(aes(x = c(pdataEV$EV, rev(pdataEV$EV)), y = c(pdataEV$lower, rev(pdataEV$upper))), fill =  hcl.colors(5, "zissou1", alpha = 0.4)[3])+
  geom_line(aes(pdataEV$EV, pdataEV$fit), lwd = 1.5, color = hcl.colors(5, "zissou1", alpha = 1)[4])+
  theme_bw()+
  ylim(c(0,25000))+
  xlim(c(-1.5,2.5))+
  theme(panel.grid.major.x = element_blank(),axis.line = element_line(color = "black",size=1.5),axis.ticks.length = unit(.5,"cm"),axis.ticks = element_line(color = "black",size = 1.5))+
  theme(panel.grid = element_blank(),axis.title = element_blank(),axis.text.y = element_text(size=25),axis.text.x = element_text(size=25))


# The different models are plotted together
plotA <- ggplot()+
  geom_point(aes(SBD$A, SBD$EPG), pch = 16, color = "black", size =5,alpha=0.6)+
  geom_polygon(aes(x= c(pdataA1$A,rev(pdataA1$A)),y = c(pdataA1$lower,rev(pdataA1$upper))), fill =  hcl.colors(5, "zissou1", alpha = 0.7)[3])+
  geom_polygon(aes(x= c(pdataA2$A,rev(pdataA2$A)),y = c(pdataA2$lower,rev(pdataA2$upper))), fill =  hcl.colors(5, "zissou1", alpha = 0.5)[3])+
  geom_polygon(aes(x= c(pdataA3$A,rev(pdataA3$A)),y = c(pdataA3$lower,rev(pdataA3$upper))), fill =  hcl.colors(5, "zissou1", alpha = 0.3)[3])+
  geom_line(aes(pdataA1$A, pdataA1$fit), lwd = 2,col = hcl.colors(5, "zissou1", alpha = 1)[1])+
  geom_line(aes(pdataA2$A, pdataA2$fit), lwd = 1.5,col = hcl.colors(5, "zissou1", alpha = 1)[2])+
  geom_line(aes(pdataA3$A, pdataA3$fit), lwd = 1,col = hcl.colors(5, "zissou1", alpha = 1)[4])+
  theme_bw()+
  ylim(c(0,25000))+
  xlim(c(-1.5,3))+
  theme(panel.grid.major.x = element_blank(),axis.line = element_line(color = "black",size=1.5),axis.ticks.length = unit(.5,"cm"),axis.ticks = element_line(color = "black",size = 1.5))+
  theme(panel.grid = element_blank(),axis.title = element_blank(),axis.text.y = element_text(size=25),axis.text.x = element_text(size=25))



plotEL <- ggplot()+
  geom_point(aes(SBD$EL, SBD$EPG),pch = 16, color = "black", size =5,alpha=0.6)+
  geom_polygon(aes(x= c(pdataEL1$EL,rev(pdataEL1$EL)),y = c(pdataEL1$lower,rev(pdataEL1$upper))), fill =  hcl.colors(5, "zissou1", alpha = 0.7)[3])+
  geom_polygon(aes(x= c(pdataEL2$EL,rev(pdataEL2$EL)),y = c(pdataEL2$lower,rev(pdataEL2$upper))), fill =  hcl.colors(5, "zissou1", alpha = 0.5)[3])+
  geom_polygon(aes(x= c(pdataEL3$EL,rev(pdataEL3$EL)),y = c(pdataEL3$lower,rev(pdataEL3$upper))), fill =  hcl.colors(5, "zissou1", alpha = 0.3)[3])+
  geom_line(aes(pdataEL1$EL, pdataEL1$fit), lwd = 2,col = hcl.colors(5, "zissou1", alpha = 1)[1])+
  geom_line(aes(pdataEL2$EL, pdataEL2$fit), lwd = 1.5,col = hcl.colors(5, "zissou1", alpha = 1)[2])+
  geom_line(aes(pdataEL3$EL, pdataEL3$fit), lwd = 1,col = hcl.colors(5, "zissou1", alpha = 1)[4])+
  theme_bw()+
  ylim(c(0,25000))+
  xlim(c(-2,3.5))+
  scale_x_continuous(breaks=c(-2,-1,0,1,2,3))+
  theme(panel.grid.major.x = element_blank(),axis.line = element_line(color = "black",size=1.5),axis.ticks.length = unit(.5,"cm"),axis.ticks = element_line(color = "black",size = 1.5))+
  theme(panel.grid = element_blank(),axis.title = element_blank(),axis.text.y = element_text(size=25),axis.text.x = element_text(size=25))


# Turn sex into binary numeric variable
SBD1 <- SBD
SBD1$SEX[which(SBD1$SEX=="female")]<-0
SBD1$SEX[which(SBD1$SEX=="male")]<-1

plotsex <- ggplot()+
  geom_point(aes(as.numeric(SBD1$SEX), SBD1$EPG), pch = 16, color = "black", size =5,alpha=0.6)+
  geom_polygon(aes(x =  c(-0.1,-0.1,0.1,0.1), y = c(pdataSEX11$lower[1],pdataSEX11$lower[2],pdataSEX11$lower[2],pdataSEX11$lower[1])), fill = hcl.colors(5, "zissou1", alpha = 0.7)[3])+
  geom_polygon(aes(x =  c(0.9,0.9,1.1,1.1), y = c(pdataSEX11$upper[1],pdataSEX11$upper[2],pdataSEX11$upper[2],pdataSEX11$upper[1])), fill = hcl.colors(5, "zissou1", alpha = 0.7)[3])+
  geom_polygon(aes(x =  c(-0.2,-0.2,0.2,0.2), y = c(pdataSEX21$lower[1],pdataSEX21$lower[2],pdataSEX21$lower[2],pdataSEX21$lower[1])), fill = hcl.colors(5, "zissou1", alpha = 0.5)[3])+
  geom_polygon(aes(x =  c(0.8,0.8,1.2,1.2), y = c(pdataSEX21$upper[1],pdataSEX21$upper[2],pdataSEX21$upper[2],pdataSEX21$upper[1])), fill = hcl.colors(5, "zissou1", alpha = 0.5)[3])+
  geom_polygon(aes(x =  c(-0.3,-0.3,0.3,0.3), y = c(pdataSEX31$lower[1],pdataSEX31$lower[2],pdataSEX31$lower[2],pdataSEX31$lower[1])), fill = hcl.colors(5, "zissou1", alpha = 0.3)[3])+
  geom_polygon(aes(x =  c(0.7,0.7,1.3,1.3), y = c(pdataSEX31$upper[1],pdataSEX31$upper[2],pdataSEX31$upper[2],pdataSEX31$upper[1])), fill = hcl.colors(5, "zissou1", alpha = 0.3)[3])+
  geom_point(aes(0, pdataSEX11$fit[1]),pch = 16,size = 5,col = hcl.colors(5, "zissou1", alpha = .7)[1])+
  geom_point(aes(1, pdataSEX11$fit[2]),pch = 16,size = 5,col = hcl.colors(5, "zissou1", alpha = .7)[1])+
  geom_point(aes(0, pdataSEX21$fit[1]),pch = 16,size = 5,col = hcl.colors(5, "zissou1", alpha = .7)[2])+
  geom_point(aes(1, pdataSEX21$fit[2]),pch = 16,size = 5,col = hcl.colors(5, "zissou1", alpha = .7)[2])+
  geom_point(aes(0, pdataSEX31$fit[1]),pch = 16,size = 5,col = hcl.colors(5, "zissou1", alpha = .7)[4])+
  geom_point(aes(1, pdataSEX31$fit[2]),pch = 16,size = 5,col = hcl.colors(5, "zissou1", alpha = .7)[4])+
    theme_bw()+
  ylim(c(0,25000))+
  scale_x_continuous(breaks=c(0,1),limits=c(-0.3,1.3))+
  theme(panel.grid.major.x = element_blank(),axis.line = element_line(color = "black",size=1.5),axis.ticks.length = unit(.5,"cm"),axis.ticks = element_line(color = "black",size = 1.5))+
  theme(panel.grid = element_blank(),axis.title = element_blank(),axis.text.y = element_text(size=25),axis.text.x = element_text(size=25))


ggsave('Degree.png', plot = plotD, device = 'png',path = "C:/Users/Xu Zhihong/Desktop/testing area/Parts of the paper 1/review round 1/Regression plot", dpi = 300,width = 5.29, height = 4.30,units= "in")
ggsave('Strength.png', plot = plotS, device = 'png',path = "C:/Users/Xu Zhihong/Desktop/testing area/Parts of the paper 1/review round 1/Regression plot", dpi = 300,width = 5.29, height = 4.30,units= "in")
ggsave('Eigenvector.png', plot = plotEV, device = 'png',path = "C:/Users/Xu Zhihong/Desktop/testing area/Parts of the paper 1/review round 1/Regression plot", dpi = 300,width = 5.29, height = 4.30,units= "in")
ggsave('Elo-rating.png', plot = plotEL, device = 'png',path = "C:/Users/Xu Zhihong/Desktop/testing area/Parts of the paper 1/review round 1/Regression plot", dpi = 300,width = 5.29, height = 4.30,units= "in")
ggsave('Age.png', plot = plotA, device = 'png',path = "C:/Users/Xu Zhihong/Desktop/testing area/Parts of the paper 1/review round 1/Regression plot", dpi = 300,width = 5.29, height = 4.30,units= "in")
ggsave('Sex.png', plot = plotsex, device = 'png',path = "C:/Users/Xu Zhihong/Desktop/testing area/Parts of the paper 1/review round 1/Regression plot", dpi = 300,width = 5.29, height = 4.30,units= "in")

# if you want only the bands instead (rather than the polygon)
# points(pdata$EV, pdata$lower, type = "l", lty = 2, lwd = 3, col = hcl.colors(2, "zissou1")[2])
# points(pdata$EV, pdata$upper, type = "l", lty = 2, lwd = 3, col = hcl.colors(2, "zissou1")[2])

grid.arrange(plotD, plotS, plotEV, nrow = 1)

## RANDOMISATION AND KO SIMULATION PLOTS
# This part is based on the result output of "run randomisations.R" or "run KO simulations.R". Their output are in same format so that they can be both applied on this script.
# Here we provide one example based on "RDres_Example.csv"
require(dplyr)
# Load the data if necessary
RDres <- read.csv("raw_data/RDres_Example.csv",  sep = ";") # be mindful of separator type
header3 <- as.character(read.csv("raw_data/RDres_Example.csv", header = FALSE, fileEncoding="UTF-8-BOM", sep = ";")[1, ])
colnames(RDres) <- header3
arrange(RDres,EV)
# Create data frame for plotting models considering strength factors. Other factors can be plotted separately following the same logic.

# CP for Coefficient Plot
CPdata_RD_S <- data.frame(Min = RDres$MinS,Estimate = RDres$S,Max = RDres$MaxS,V = "S")
CPdata_RD_S <- arrange(CPdata_RD_S,Min)
CPdata_RD_S$NO <- c(1:15)
CPdata_RD_S$P <- "0.025"
CPdata_RD_S$E <- CPdata_RD_S$Min
CPdata_RD_S[16:30,] <- CPdata_RD_S[1:15,]
CPdata_RD_S[16:30,6] <- "0.5"
CPdata_RD_S[16:30,7] <- CPdata_RD_S$Estimate[1:15]
CPdata_RD_S[31:45,] <- CPdata_RD_S[1:15,]
CPdata_RD_S[31:45,6] <- "0.975"
CPdata_RD_S[31:45,7] <- CPdata_RD_S$Max[1:15]
CPdata_RD_S$tag = "0"
CPdata_RD_S$tag[which(CPdata_RD_S$Max < 0)] = "1"
CPdata_RD_S$tag[which(CPdata_RD_S$Min > 0)] = "1"


Plot <- ggplot(CPdata_RD_S)  +  
# First plot a line showing the coefficient of each model
    geom_line(aes(x=E,y=as.factor(NO)),size=1,data = subset(CPdata_RD_S, P %in% c('0.025','0.975')),color='grey')+     
# Then plot the estimates of each model
    geom_point(data = subset(CPdata_RD_S,P %in% '0.5'),aes(x=E,y=as.factor(NO),color = as.factor(tag)),size =2) +
  scale_color_manual("",values = c("1" = 'black', "0" = grey(.8)))+
  geom_vline(xintercept=0,linetype=1,color = 'black')+
# Here the numbers are estimates and confidence intervals of original model  
  geom_vline(xintercept=0.09821095,linetype=2,color='black',lwd=2)+
  geom_vline(xintercept=0.47406758,linetype=2,color='black',lwd=2)+
  geom_vline(xintercept=0.2861392670,linetype=1,color='black',lwd=2)+
  theme(legend.background=element_blank(),panel.background = element_rect(fill = 'white'))+
  coord_cartesian(xlim=c(-0.6,0.6))+
  theme_void()+
  theme(legend.position = "none")+
  theme(panel.grid.major.x = element_blank(),axis.line.x = element_line(color = "black"))+
  theme(axis.title = element_blank(),axis.text.y = element_blank(),axis.ticks.y = element_blank(),axis.text.x = element_text(size=20))
Plot
