
source("mood.R") # library("Rccp"); library("parallel") # options$mc.cores

## system options ##
options("mc.cores"=2) #the number of cores to be used

######################
# Example dataset
library(fda)
data(growth)

data = growth$hgtf

data <- t(read.csv("../vars100000.csv", header=F))

cur.wd <- getwd()
moodIdx.classic <- computeMoodIndices.classic(data)
setwd("src/")
#moodIdx.classic <- computeMoodIndices.rcpp(data)
moodIdx.rcpp <- computeMoodIndices.rcpp(data)
setwd(cur.wd)

boxplot(abs(moodIdx.rcpp - moodIdx.classic))
(err <- sqrt(colSums(abs(moodIdx.rcpp - moodIdx.classic) ^ 2)))
# > x.e-13 y.e-11 z.e-13
