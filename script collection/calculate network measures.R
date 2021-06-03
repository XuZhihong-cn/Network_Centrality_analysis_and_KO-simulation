# this script aims at calculating network measures for each individual


library(igraph)

#read the data
NetData <- read.csv("raw_data/NetData.csv",sep=";") # be mindful of separator type


#turn the dataframe into a matrix object and name the col/row as individual ID
NetData <- NetData[,-1]
NetData <- as.matrix(NetData)
is.matrix(NetData)
rownames(NetData)<- colnames(NetData)

#turn data matrix into network object
netall <- graph.adjacency(NetData, mode = "undirected", weighted = TRUE, diag = FALSE)

#calculate network measures
#degree
degall<-degree(netall, mode="all")
#strength
strengthall<-strength(netall,mode="all")
#eigenvector
evcentall<-evcent(netall)$vector

# see also the "run GLMMs" script to build and calculate adult-female or juvenile only network measures


