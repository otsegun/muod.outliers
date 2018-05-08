# MUOD
The R implementation of MUOD (Massive Unsupervised Outliers Detection)

## How To Cite Us
Azcorra, Arturo, Luis F. Chiroque, Rubén Cuevas, Antonio Fernández Anta, Henry Laniado, Rosa Elvira Lillo, Juan Romo, and Carlo Sguera. "Unsupervised Scalable Statistical Method for Identifying Influential Users in Online Social Networks." Scientific Reports (2018).

[doi:10.1038/s41598-018-24874-2](https://www.nature.com/articles/s41598-018-24874-2)

***

## Example usage with `mtcars`

```r
#install.packages("/path/to/muod.outiers.git/", repos=NULL, type="source")
devtools::install_github("luisfo/muod.outliers")
library(muod)

data(mtcars)
data <- as.matrix(mtcars)

options("mc.cores"=4) # Let's use 4 cores

# This is the only function to call
## it expects a numeric matrix as first parameter
outliers <- getOutliers(t(data), slope=1)
```

### Alternative use with a local parallel threads cluster
```r
library(muod)
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
```

### Then, you can color outiliers
```r
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

## Bugs and Feedback

We will appreaciate you report us bugs or give us your feedback.
Luis F. Chiroque's email: 
[luisfo89@gmail.com](mailto:luisfo89@gmail.com) 

