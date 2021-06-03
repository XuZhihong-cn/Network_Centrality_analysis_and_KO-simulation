# this script aims at calculating dominance ranks (Elo-rating) for each individual

library(EloRating)
AggData <- read.csv("raw_data/AggData.csv",sep=";") # be mindful of separator type
# Check that the sequence of interactions for missing data and such
seqcheck(winner=AggData$winner, loser=AggData$loser, Date = AggData$Date)
# Calculate Elo-Ratings
res <- elo.seq(winner=AggData$winner, loser=AggData$loser, Date =AggData$Date,runcheck=TRUE)
# Extract the Elo-rating of each individual at the last day of data collection
extract_elo(res, "2017-06-21")
