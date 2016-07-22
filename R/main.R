
### Example CODE ###

source("mood.R") # library("Rccp"); library("parallel") # options$mc.cores

## system options ##
options("mc.cores"=2) #the number of cores to be used

######################
# Example dataset
library(fda)
data(growth)

data = growth$hgtf

# REMARK: sometimes data has to be trasposed, since the function expects
#   a matrix whose columns are the 'users' (and variables as rows)
#

# plot data as FDA
domain = growth$age
m = length(domain)
# Obtaining curves at a domain at equidistant points
domain.equi = seq(from=1,to=18,length.out=m)
X = t(growth$hgtf)
X.hat <- t(apply(X, 1, function(x){
  x.sfun <- splinefun(domain, x, method="natural")
  x.sfun(domain.equi, deriv = 0)
}))
# Since I do not need X anymore:
X = X.hat
main.plot = "Growth - girls"
xlab.plot = "years"
ylab.plot = "heights"
plot(NULL, xlim=c(domain.equi[1], domain.equi[m]), ylim = c(min(X), max(X))
     , type = "n", xlab = xlab.plot, ylab = ylab.plot, main = main.plot)

addLines <- function(X, domain.equi, color="black", lwd=1){
  trash <- apply(X, 1, function(x){
    lines(domain.equi, x, col=color, lwd=lwd)
  })
}

addLines(X, domain.equi)
######################


# force NaN's (e.g., missing data); still works!
# data[4,7] <- NaN

# Run algorithm: default options
outliers <- getOutliers(data) # JUST THE ONLY FUNCTION YOU HAVE TO CALL
# outliers
outliers$shape
outliers$amplitude
outliers$magnitude

# done!

######################
# Continue Example
# plot outliers
addLines(t(data[,outliers$shape]), domain.equi, color="yellow", lwd=2)
addLines(t(data[,outliers$amplitude]), domain.equi, color="red", lwd=2)
addLines(t(data[,outliers$magnitude]), domain.equi, color="blue", lwd=2)
