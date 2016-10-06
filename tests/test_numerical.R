
source("mood.R") # library("Rccp"); library("parallel") # options$mc.cores

## system options ##
options("mc.cores"=2) #the number of cores to be used

######################
# Example dataset
library(fda)
data(growth)

data = growth$hgtf

cur.wd <- getwd()
moodIdx.classic <- computeMoodIndices.classic(data)
setwd("src/")
moodIdx.rcpp <- computeMoodIndices.rcpp(data)
setwd(cur.wd)

boxplot(moodIdx.rcpp - moodIdx.classic)
err <- sqrt(colSums(abs(moodIdx.rcpp - moodIdx.classic) ^ 2))
# > x.e-15 y.e-13 z.e-15
