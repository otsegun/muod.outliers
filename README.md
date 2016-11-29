# MOOD
The R implementation of MOOD (outliers)

## Example usage with `mtcars`

```r
#install.packages("/path/to/mood.outiers.git/", repos=NULL, type="source")
devtools::install_github("luisfo/mood.outliers")
library(mood)

data(mtcars)
data <- as.matrix(mtcars)

options("mc.cores"=4) # Let's use 4 cores

# This is the only function to call
## it expects a numeric matrix as first parameter
outliers <- getOutliers(t(data), slope=1)

ylim <- range(data)
xlim <- c(1, ncol(data))

# plot
plot(NULL, type='n', ylim=ylim, xlim=xlim)
## plot all
apply(data, 1, lines)
## plot shape outliers in red
apply(data[outliers$shape,], 1, lines, col="red", lwd=2)
## plot amplitude outliers in blue
apply(data[outliers$amplitude,], 1, lines, col="blue", lwd=2)
## plot magnitude outliers in green
apply(data[outliers$magnitude,], 1, lines, col="green", lwd=2)
```

### Alternative use with a local parallel threads cluster
```r
library(mood)
library(parallel)

# configure
options("mc.cores"=4) # Let's use 4 cores
# before loading data, instantiate cluster
cl <- makeCluster(getOption("mc.cores", 1))

# load data
data(mtcars)
data <- as.matrix(mtcars)

# This is the only function to call
## it expects a numeric matrix as first parameter
outliers <- getOutliers(t(data), slope=1, parClus=cl)

# stop cluster
stopCluster(cl)

ylim <- range(data)
xlim <- c(1, ncol(data))

# plot
plot(NULL, type='n', ylim=ylim, xlim=xlim)
## plot all
apply(data, 1, lines)
## plot shape outliers in red
apply(data[outliers$shape,], 1, lines, col="red", lwd=2)
## plot amplitude outliers in blue
apply(data[outliers$amplitude,], 1, lines, col="blue", lwd=2)
## plot magnitude outliers in green
apply(data[outliers$magnitude,], 1, lines, col="green", lwd=2)
```
