# This script details the plotting used in the paper


###################  Eigenvector change across partial network

library(igraph)


IndividualData <- read.csv("raw_data/IndData.csv",  sep = ";") # be mindful of separator type
NetData <- read.csv("raw_data/NetData.csv",  sep = ";") # be mindful of separator type

NS <- c('takana','mushi','neji','uso')

NetData <- NetData[,-1]
NetData <- as.matrix(NetData)
is.matrix(NetData)
rownames(NetData)<- colnames(NetData)
IBD <- IndividualData
header2 <- as.character(unlist(read.csv("raw_data/IndData.csv", header = FALSE, fileEncoding="UTF-8-BOM", sep = ";")[1, ]))
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
NetALLR <- data.frame(ids = colnames(NetData), EV = scale(eigen_centrality(NetALL)$vector), EV5=NA,EV10=NA,EV25=NA,EV50=NA,EVaf=NA,EVjuv=NA)

sam1x<-NetALLR$ids[-which(NetALLR$ids %in% sam1)]
sam2x<-NetALLR$ids[-which(NetALLR$ids %in% sam2)]
sam3x<-NetALLR$ids[-which(NetALLR$ids %in% sam3)]
sam4x<-NetALLR$ids[-which(NetALLR$ids %in% sam4)]

NetALLR$EV5[which(NetALLR$ids %in% sam1x)]<-scale(eigen_centrality(Net5)$vector)
NetALLR$EV10[which(NetALLR$ids %in% sam2x)]<-scale(eigen_centrality(Net10)$vector)
NetALLR$EV25[which(NetALLR$ids %in% sam3x)]<-scale(eigen_centrality(Net25)$vector)
NetALLR$EV50[which(NetALLR$ids %in% sam4x)]<-scale(eigen_centrality(Net50)$vector)
NetALLR$EVaf[which(NetALLR$ids %in% colnames(NetDataAF))]<-scale(eigen_centrality(NetAF)$vector)
NetALLR$EVjuv[which(NetALLR$ids %in% colnames(NetDataJUV))]<-scale(eigen_centrality(NetJUV)$vector)


require(ggplot2)
Example <- NetALLR[c(3,9,11,15,32,33,47,49),]
Example <-as.data.frame(t(Example)[-1,])
Example$Net <- NA
Example$Net <-  factor(c("100%","95%","90%","75%","50%","AF","JUV"))
levels(Example$Net) <- c("100%","95%","90%","75%","50%","AF","JUV")
Example[,1][which(Example[,1] %in% NA)] <- -1.2
Example[,2][which(Example[,2] %in% NA)] <- -1.2
Example[,3][which(Example[,3] %in% NA)] <- -1.2
Example[,4][which(Example[,4] %in% NA)] <- -1.2
Example[,5][which(Example[,5] %in% NA)] <- -1.2
Example[,6][which(Example[,6] %in% NA)] <- -1.2
Example[,7][which(Example[,7] %in% NA)] <- -1.2
Example[,8][which(Example[,8] %in% NA)] <- -1.2


plot1 <-ggplot(Example)+
  geom_point(aes(x=Net,y=as.numeric(kei),shape=ifelse(as.numeric(kei) %in% -1.2,"1","2")),size=6)+
  scale_shape_manual(values = c(13,16))+
  theme_void()+
  theme(panel.grid.major.x = element_blank(),axis.line = element_line(color = "black",size=1.5),axis.ticks.length = unit(.5,"cm"),axis.ticks = element_line(color = "black",size = 1.5))+
  theme(legend.position ="none" ,panel.grid = element_blank(),axis.title = element_blank(),axis.text.y = element_text(size=25),axis.text.x = element_text(size=15))+
  scale_y_continuous(breaks=c(-1,-0,0,1,2),limits = c(-1.2,2.4))



plot2 <-ggplot(Example)+
  geom_point(aes(x=Net,y=as.numeric(minku),shape=ifelse(as.numeric(minku) %in% -1.2,"1","2")),size=6)+
  scale_shape_manual(values = c(13,16))+
  theme_void()+
  theme(panel.grid.major.x = element_blank(),axis.line = element_line(color = "black",size=1.5),axis.ticks.length = unit(.5,"cm"),axis.ticks = element_line(color = "black",size = 1.5))+
  theme(legend.position ="none" ,panel.grid = element_blank(),axis.title = element_blank(),axis.text.y = element_text(size=25),axis.text.x = element_text(size=15))+
  scale_y_continuous(breaks=c(-1,-0,0,1,2),limits = c(-1.2,2.4))


plot3 <-ggplot(Example)+
  geom_point(aes(x=Net,y=as.numeric(muku),shape=ifelse(as.numeric(muku) %in% -1.2,"1","2")),size=6)+
  scale_shape_manual(values = c(13,16))+
  theme_void()+
  theme(panel.grid.major.x = element_blank(),axis.line = element_line(color = "black",size=1.5),axis.ticks.length = unit(.5,"cm"),axis.ticks = element_line(color = "black",size = 1.5))+
  theme(legend.position ="none" ,panel.grid = element_blank(),axis.title = element_blank(),axis.text.y = element_text(size=25),axis.text.x = element_text(size=15))+
  scale_y_continuous(breaks=c(-1,-0,0,1,2),limits = c(-1.2,2.4))


plot4 <-ggplot(Example)+
  geom_point(aes(x=Net,y=as.numeric(okura),shape=ifelse(as.numeric(okura) %in% -1.2,"1","2")),size=6)+
  scale_shape_manual(values = c(13,16))+
  theme_void()+
  theme(panel.grid.major.x = element_blank(),axis.line = element_line(color = "black",size=1.5),axis.ticks.length = unit(.5,"cm"),axis.ticks = element_line(color = "black",size = 1.5))+
  theme(legend.position ="none" ,panel.grid = element_blank(),axis.title = element_blank(),axis.text.y = element_text(size=25),axis.text.x = element_text(size=15))+
  scale_y_continuous(breaks=c(-1,-0,0,1,2),limits = c(-1.2,2.4))



plot5 <-ggplot(Example)+
  geom_point(aes(x=Net,y=as.numeric(botan),shape=ifelse(as.numeric(botan) %in% -1.2,"1","2")),size=6)+
  scale_shape_manual(values = c(13,16))+
  theme_void()+
  theme(panel.grid.major.x = element_blank(),axis.line = element_line(color = "black",size=1.5),axis.ticks.length = unit(.5,"cm"),axis.ticks = element_line(color = "black",size = 1.5))+
  theme(legend.position ="none" ,panel.grid = element_blank(),axis.title = element_blank(),axis.text.y = element_text(size=25),axis.text.x = element_text(size=15))+
  scale_y_continuous(breaks=c(-1,-0,0,1,2),limits = c(-1.2,2.4))



plot6 <-ggplot(Example)+
  geom_point(aes(x=Net,y=as.numeric(hado),shape=ifelse(as.numeric(hado) %in% -1.2,"1","2")),size=6)+
  scale_shape_manual(values = c(13,16))+
  theme_void()+
  theme(panel.grid.major.x = element_blank(),axis.line = element_line(color = "black",size=1.5),axis.ticks.length = unit(.5,"cm"),axis.ticks = element_line(color = "black",size = 1.5))+
  theme(legend.position ="none" ,panel.grid = element_blank(),axis.title = element_blank(),axis.text.y = element_text(size=25),axis.text.x = element_text(size=15))+
  scale_y_continuous(breaks=c(-1,-0,0,1,2),limits = c(-1.2,2.4))



plot7 <-ggplot(Example)+
  geom_point(aes(x=Net,y=as.numeric(takana),shape=ifelse(as.numeric(takana) %in% -1.2,"1","2")),size=6)+
  scale_shape_manual(values = c(13,16))+
  theme_void()+
  theme(panel.grid.major.x = element_blank(),axis.line = element_line(color = "black",size=1.5),axis.ticks.length = unit(.5,"cm"),axis.ticks = element_line(color = "black",size = 1.5))+
  theme(legend.position ="none" ,panel.grid = element_blank(),axis.title = element_blank(),axis.text.y = element_text(size=25),axis.text.x = element_text(size=15))+
  scale_y_continuous(breaks=c(-1,-0,0,1,2),limits = c(-1.2,2.4))



plot8 <-ggplot(Example)+
  geom_point(aes(x=Net,y=as.numeric(yomogi),shape=ifelse(as.numeric(yomogi) %in% -1.2,"1","2")),size=6)+
  scale_shape_manual(values = c(13,16))+
  theme_void()+
  theme(panel.grid.major.x = element_blank(),axis.line = element_line(color = "black",size=1.5),axis.ticks.length = unit(.5,"cm"),axis.ticks = element_line(color = "black",size = 1.5))+
  theme(legend.position ="none" ,panel.grid = element_blank(),axis.title = element_blank(),axis.text.y = element_text(size=25),axis.text.x = element_text(size=15))+
  scale_y_continuous(breaks=c(-1,-0,0,1,2),limits = c(-1.2,2.4))



ggsave('Net1.png', plot = plot1, device = 'png',path = "C:/Users/Xu Zhihong/Desktop/testing area/Parts of the paper 1/review round 2/network plot", dpi = 300,width = 4.32, height = 4.30,units= "in")
ggsave('Net2.png', plot = plot2, device = 'png',path = "C:/Users/Xu Zhihong/Desktop/testing area/Parts of the paper 1/review round 2/network plot", dpi = 300,width = 4.32, height = 4.30,units= "in")
ggsave('Net3.png', plot = plot3, device = 'png',path = "C:/Users/Xu Zhihong/Desktop/testing area/Parts of the paper 1/review round 2/network plot", dpi = 300,width = 4.32, height = 4.30,units= "in")
ggsave('Net4.png', plot = plot4, device = 'png',path = "C:/Users/Xu Zhihong/Desktop/testing area/Parts of the paper 1/review round 2/network plot", dpi = 300,width = 4.32, height = 4.30,units= "in")
ggsave('Net5.png', plot = plot5, device = 'png',path = "C:/Users/Xu Zhihong/Desktop/testing area/Parts of the paper 1/review round 2/network plot", dpi = 300,width = 4.32, height = 4.30,units= "in")
ggsave('Net6.png', plot = plot6, device = 'png',path = "C:/Users/Xu Zhihong/Desktop/testing area/Parts of the paper 1/review round 2/network plot", dpi = 300,width = 4.32, height = 4.30,units= "in")
ggsave('Net7.png', plot = plot7, device = 'png',path = "C:/Users/Xu Zhihong/Desktop/testing area/Parts of the paper 1/review round 2/network plot", dpi = 300,width = 4.32, height = 4.30,units= "in")
ggsave('Net8.png', plot = plot8, device = 'png',path = "C:/Users/Xu Zhihong/Desktop/testing area/Parts of the paper 1/review round 2/network plot", dpi = 300,width = 4.32, height = 4.30,units= "in")


##
## GLMM RESULTS PLOTS (credits @ Christof Neumann)
# scripts needed: "run GLMMs"
# each variable gets its own plot
library(effects)
require(MuMIn)
require(glmmTMB)

# if needed load the dataset and the network matrix again 

SBD <- read.csv("raw_data/SampleData.csv",  sep = ";") # be mindful of separator type
header1 <- as.character(unlist(read.csv("raw_data/SampleData.csv", header = FALSE, fileEncoding="UTF-8-BOM", sep = ";")[1, ]))
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
  geom_point(aes(SBD$D, SBD$EPG), pch =16, color = "#6F6F6F", size =3,alpha=0.6)+
  geom_polygon(aes(x = c(pdataD$D, rev(pdataD$D)), y = c(pdataD$lower, rev(pdataD$upper))), fill =  "#91CBDC", alpha = 0.7)+
  geom_line(aes(pdataD$D, pdataD$fit), lwd = 1.5, color ="#C0FAFF")+
  theme_bw()+
  ylim(c(0,24000))+
  xlim(c(-3,2))+
  theme(panel.grid.major.x = element_blank(),axis.line = element_line(color = "black",size=1.5),axis.ticks.length = unit(.5,"cm"),axis.ticks = element_line(color = "black",size = 1.5))+
  theme(panel.grid = element_blank(),axis.title = element_blank(),axis.text.y = element_text(size=25),axis.text.x = element_text(size=25))


plotS <- ggplot()+
  geom_point(aes(SBD$S, SBD$EPG), pch =16, color = "#6F6F6F", size =3,alpha=0.6)+
  geom_polygon(aes(x = c(pdataS$S, rev(pdataS$S)), y = c(pdataS$lower, rev(pdataS$upper))), fill ='#B87A9B',alpha=0.7)+
  geom_line(aes(pdataS$S, pdataS$fit), lwd = 1.5, color = "#F1AFD1")+
  theme_bw()+
  ylim(c(0,24000))+
  xlim(c(-3,3))+
  theme(panel.grid.major.x = element_blank(),axis.line = element_line(color = "black",size=1.5),axis.ticks.length = unit(.5,"cm"),axis.ticks = element_line(color = "black",size = 1.5))+
  theme(panel.grid = element_blank(),axis.title = element_blank(),axis.text.y = element_text(size=25),axis.text.x = element_text(size=25))



plotEV <- ggplot()+
  geom_point(aes(SBD$EV, SBD$EPG), pch =16, color = "#6F6F6F", size =3,alpha=0.6)+
  geom_polygon(aes(x = c(pdataEV$EV, rev(pdataEV$EV)), y = c(pdataEV$lower, rev(pdataEV$upper))), fill = '#FFE373',alpha=0.7)+
  geom_line(aes(pdataEV$EV, pdataEV$fit), lwd = 1.5, col = "#FFFD00")+
  theme_bw()+
  ylim(c(0,24000))+
  xlim(c(-1.5,2.5))+
  theme(panel.grid.major.x = element_blank(),axis.line = element_line(color = "black",size=1.5),axis.ticks.length = unit(.5,"cm"),axis.ticks = element_line(color = "black",size = 1.5))+
  theme(panel.grid = element_blank(),axis.title = element_blank(),axis.text.y = element_text(size=25),axis.text.x = element_text(size=25))


# The different models are plotted together
plotA <- ggplot()+
  geom_point(aes(SBD$A, SBD$EPG), pch = 16, color = "#6F6F6F", size =3,alpha=0.6)+
  geom_polygon(aes(x= c(pdataA1$A,rev(pdataA1$A)),y = c(pdataA1$lower,rev(pdataA1$upper))), fill = '#91CBDC',alpha=0.7)+
  geom_polygon(aes(x= c(pdataA2$A,rev(pdataA2$A)),y = c(pdataA2$lower,rev(pdataA2$upper))), fill =  '#B87A9B',alpha=0.7)+
  geom_polygon(aes(x= c(pdataA3$A,rev(pdataA3$A)),y = c(pdataA3$lower,rev(pdataA3$upper))), fill =  '#FFE373',alpha=0.7)+
  geom_line(aes(pdataA1$A, pdataA1$fit), lwd = 1.2,col = "#C0FAFF")+
  geom_line(aes(pdataA2$A, pdataA2$fit), lwd = 1.2,col = "#F1AFD1")+
  geom_line(aes(pdataA3$A, pdataA3$fit), lwd = 1.2,col = "#FFFD00")+
  theme_bw()+
  ylim(c(0,24000))+
  xlim(c(-1.5,3))+
  theme(panel.grid.major.x = element_blank(),axis.line = element_line(color = "black",size=1.5),axis.ticks.length = unit(.5,"cm"),axis.ticks = element_line(color = "black",size = 1.5))+
  theme(panel.grid = element_blank(),axis.title = element_blank(),axis.text.y = element_text(size=25),axis.text.x = element_text(size=25))



plotEL <- ggplot()+
  geom_point(aes(SBD$EL, SBD$EPG),pch = 16, color = "#6F6F6F", size =3,alpha=0.6)+
  geom_polygon(aes(x= c(pdataEL1$EL,rev(pdataEL1$EL)),y = c(pdataEL1$lower,rev(pdataEL1$upper))), fill =  '#91CBDC',alpha=0.7)+
  geom_polygon(aes(x= c(pdataEL2$EL,rev(pdataEL2$EL)),y = c(pdataEL2$lower,rev(pdataEL2$upper))), fill =  '#B87A9B',alpha=0.7)+
  geom_polygon(aes(x= c(pdataEL3$EL,rev(pdataEL3$EL)),y = c(pdataEL3$lower,rev(pdataEL3$upper))), fill =  '#FFE373',alpha=0.7)+
  geom_line(aes(pdataEL1$EL, pdataEL1$fit), lwd = 1.2,col = "#C0FAFF")+
  geom_line(aes(pdataEL2$EL, pdataEL2$fit), lwd = 1.2,col = "#F1AFD1")+
  geom_line(aes(pdataEL3$EL, pdataEL3$fit), lwd = 1.2,col = "#FFFD00")+
  theme_bw()+
  ylim(c(0,24000))+
  xlim(c(-2,3.5))+
  scale_x_continuous(breaks=c(-2,-1,0,1,2,3))+
  theme(panel.grid.major.x = element_blank(),axis.line = element_line(color = "black",size=1.5),axis.ticks.length = unit(.5,"cm"),axis.ticks = element_line(color = "black",size = 1.5))+
  theme(panel.grid = element_blank(),axis.title = element_blank(),axis.text.y = element_text(size=25),axis.text.x = element_text(size=25))


# Turn sex into binary numeric variable
SBD1 <- SBD

pdataSEX<- data.frame('SEX'=c('female','male','female','male','female','male'),"Model1"=c('D','D','S','S','EV','EV'),"Model2"=c('D','D','S','S','EV','EV'),"upper"=c(pdataSEX1$upper,pdataSEX2$upper,pdataSEX3$upper),"lower"=c(pdataSEX1$lower,pdataSEX2$lower,pdataSEX3$lower),"fit"=c(pdataSEX1$fit,pdataSEX2$fit,pdataSEX3$fit))
require(ggnewscale)

plotsex <- ggplot(SBD1)+
  ylim(c(0,24000))+
  geom_violin(aes(x=SEX, y=EPG), fill = "#6F6F6F",color = "#6F6F6F",adjust=2)+
  scale_color_manual("Model1",values=c('#91CBDC', '#B87A9B','#FFE373'))+
  geom_errorbar(data=pdataSEX,aes(x=SEX,ymin=lower,ymax=upper,col=Model1),width=0.2,lwd=1.5,position =position_dodge(width=0.4))+
  new_scale_color()+
  scale_color_manual("Model2",values=c("#C0FAFF", "#F1AFD1", "#FFFD00"))+
  geom_point(data=pdataSEX,aes(x=SEX,y=fit,col=Model2),shape=19,size=3,position=position_dodge(width=0.4))+
  theme_bw()+
  theme(panel.grid.major.x = element_blank(),axis.line = element_line(color = "black",size=1.5),axis.ticks.length = unit(.5,"cm"),axis.ticks = element_line(color = "black",size = 1.5))+
  theme(legend.position ="none",panel.grid = element_blank(),axis.title = element_blank(),axis.text.y = element_text(size=25),axis.text.x = element_text(size=25))

ggsave('Degree.png', plot = plotD2,device = 'png',path = , dpi = 300,width = 5.29, height = 5,units= "in")
ggsave('Strength.png', plot = plotS2, device = 'png', dpi = 300,width = 5.29, height = 5,units= "in")
ggsave('Eigenvector.png', plot = plotEV2, device = 'png', dpi = 300,width = 5.29, height = 5,units= "in")
ggsave('Elo-rating.png', plot = plotEL2, device = 'png', dpi = 300,width = 5.29, height = 5,units= "in")
ggsave('Age.png', plot = plotA2, device = 'png', dpi = 300,width = 5.29, height = 5,units= "in")
ggsave('Sex.png', plot = plotsex2, device = 'png', dpi = 300,width = 5.29, height = 5,units= "in")



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
