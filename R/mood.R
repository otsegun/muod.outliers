
#
# Title: MOOD algorithm (Metric Outlier Online Detection) %title might change
# Author: Luis F. Chiroque
# Affiliation: IMDEA Networks Institute
# e-mail: lf.chiroque@imdea.org
#
#   Algorithm originally designed by The Statistics Department \
# at Universidad Carlos III de Madrid
#
# version: 0.9
#

library("parallel") # mclapply(); options("mc.cores")
library("Rcpp") #sourceCpp()

#compiling the sources; this may take a while
#sourceCpp("../src/cor_cov_blockwise.cpp")

#
# auxiliar functions for indices cuttoff method
#

# from hot-spots paper
# source: http://www.nature.com/articles/srep05276
find.tangent.X.intercept.orig <- function(x, y, newx=x[which.max(diff(y))+1], plot=F){
  spl <- smooth.spline(y ~ x)
  #newx <- x[which.max(diff(y))+1] # which.max(x)
  pred0 <- predict(spl, x=newx, deriv=0)
  pred1 <- predict(spl, x=newx, deriv=1)

  #slope correction
  pred1$y <- max(pred1$y, 1.5)

  y.intercept <- pred0$y - (pred1$y*newx)
  x.intercept <- -y.intercept/pred1$y

  if( plot ) {
    plot(x, y, type="l", ylim=c(0,max(y)))
    abline(h=0, col=8)
    lines(spl, col=2) # spline
    points(pred0, col=2, pch=19) # point to predict tangent 
    lines(x, y.intercept + pred1$y*x, col=3) # tangent (1st deriv. of spline at newx)
    points(x.intercept, 0, col=3, pch=19) # x intercept
  }
  x.intercept
}

find.tangent.X.intercept <- function(x, y, which.x=which.max(diff(y))+1, plot=F){
  spl <- smooth.spline(y ~ x)
  #which.x <- which.max(diff(y))+1 # which.max(x)
  newx.0 <- x[which.min(diff(diff(y)))]
  newx.1 <- mean(x[c(which.x-1, which.x)], na.rm=TRUE)
  pred0 <- predict(spl, x=newx.0, deriv=0)
  pred1 <- predict(spl, x=newx.1, deriv=1)

  #slope correction
  pred1$y <- max(pred1$y, 1.5)

  y.intercept <- pred0$y - (pred1$y*newx.0)
  x.intercept <- -y.intercept/pred1$y

  if( plot ) {
    plot(x, y, type="l", ylim=c(0,max(y)))
    abline(h=0, col=8)
    lines(spl, col=2) # spline
    points(pred0, col=2, pch=19) # point to predict tangent 
    lines(x, y.intercept + pred1$y*x, col=3) # tangent (1st deriv. of spline at newx)
    points(x.intercept, 0, col=3, pch=19) # x intercept
  }
  x.intercept
}

compute_tangent_cutoff <- function(metric, plot=F){
  cpoint <- find.tangent.X.intercept(seq(0,1,length=length(metric))
                                   , metric
                                   #, newx = newx
                                   , plot = plot)
  ceiling(cpoint * length(metric))
}


compute_exp_elbow.last <- function(metric, base=2, exp=2
                              , exp.limit=max(3, ceiling(log(sqrt(length(metric))))+1)
                              , slope=1, h.from=1-1/base^2, h.to=1, rSide=T
                              , plot=F, debug=F) {
  # 1st derivative
  h <- 1 / base^exp
  step <- ceiling(length(metric) * h)
  #if( metric[1] > metric[length(metric)] ) {
  deriv.1.seq <- (step + 1):length(metric) / length(metric)
  deriv.1 <- diff(metric, lag=step) / h
  #} else {
  #  deriv.1.seq <- 0:(length(metric) - step-1) / length(metric)
  #  deriv.1 <- diff(metric, lag=step) * length(metric)
  #}
  middle <- (h.to - h.from) / 2 + h.from
  if ( plot ) {
    #plot(deriv.1.seq, deriv.1
    #     , type="l", xlim=c(h.from * .9, min(1, h.to * 1.1)))
    plot(deriv.1.seq, deriv.1, main=paste0("h=1/", base, "^", exp)
         , type="l", log="", xlim=c(h.from * .999, min(1, h.to * 1.001)))
    abline(h=slope, lty=2)
    abline(v=h.from, lty=3)
    abline(v=h.to, lty=3)
    abline(v=middle, lty=4)
  }
  if ( exp < exp.limit && slope > min(deriv.1[which(deriv.1.seq >= h.from)]) && slope < max(deriv.1[which(deriv.1.seq <= h.to)]) ) {
    # iteraitve step
    if ( slope > min(deriv.1[which(deriv.1.seq >= (h.to - h.from) / 2 + h.from)]) ) {
      h.from <- h.to - (h.to - h.from) / 2 # at the right
      new.rSide <- T
    }else{
      h.to <- h.from + (h.to - h.from) / 2 # at the left
      new.rSide <- F
    }
    if( rSide==T && new.rSide==F ){
      return( middle )
    }
    solution <- compute_exp_elbow.last(metric, exp=exp+1, exp.limit = exp.limit, slope=slope, h.from=h.from, h.to=h.to, plot=plot, rSide=rSide)
    if ( solution == Inf ) {
      solution <- max(deriv.1.seq[which(deriv.1 <= slope)])
      if( debug ) print(paste0("h=1/", base, "^", exp))
    }
  }else{
    solution <- Inf
  }
  return ( solution )
}

minmax.consecutive.seq <- function(seq) {
  if( length(seq) == 0 ) return(seq)
  s <- cumsum(c(F, diff(seq)!=1))
  lst.consec <- lapply(which(table(s)>1)-1
                       , function(id) {seq[s %in% id]} )
  if( length(lst.consec) == 0 ) return(c())
  ifelse(max(seq) %in% lst.consec[[length(lst.consec)]]
         , min(lst.consec[[length(lst.consec)]])
         , max(seq))
}

# it computes the derivative of a given curve as:
#   f'(x) = lim_{h -> 0} ( f(x) - f(x-h*) ) / h
# h*: discrete
# h: continuous
# Assumptions: x \in [0, 1]
#
compute_exp_elbow <- function(metric, base=2, exp=1 #max(1, ceiling(log(length(metric)/log(10))))
                              , exp.limit=max(2, floor(log(length(metric))))
                              , slope=1, enhanced=T, plot=F, debug=F) {
  # 1st derivative
  #h.discrete <- base ^ exp
  h.discrete <- ceiling(length(metric) / base ^ exp )
  h <- h.discrete / length(metric)
  #deriv1.seq <- 1:(length(metric) - h)
  deriv1.seq <- (h.discrete+1):length(metric)
  deriv1 <- diff(metric, lag=h.discrete) / h
  derivThrs <- slope # / (2*h) / 2
  #   if( norm ){
  #     derivThrs <- derivThrs / h
  #   }
  #if( debug ) print(derivThrs)
  offset <- floor(length(deriv1) / 2)
  subset <- offset:length(deriv1)
  if( enhanced ){
    cutoff <- min(minmax.consecutive.seq(which(deriv1[subset] >= derivThrs))+offset-1+h.discrete-1, Inf)
  }else{
    cutoff <- min(which(deriv1[subset] >= derivThrs)+offset-1+h.discrete-1, Inf)
  }
  if ( plot ) {
    plot(deriv1.seq, deriv1, main=paste0("h = ", base, "^", exp, " cutoff=", cutoff)
         , type="l", log=""
         , xlim=c(1, length(metric)), ylim=c(min(deriv1, derivThrs), max(deriv1, derivThrs)))
    abline(h=derivThrs, lty=2)
  }
  if( is.infinite(cutoff) ){
    if( exp < exp.limit ){
      cutoff <- compute_exp_elbow(metric, base, exp+1, exp.limit, slope, enhanced, plot, debug)
      #       if( ! is.infinite(cutoff.deep) ){
      #         cutoff <- cutoff.deep
      #       }
    }
  }
  cutoff
}

###### elbow/ROC method
#http://stackoverflow.com/questions/2018178/finding-the-best-trade-off-point-on-a-curve/2022348#2022348
getFurtherPoint <- function(curve)
{
  # get coordinates of all the points
  nPoints <- length(curve)
  allCoord <- cbind(1:nPoints, curve)
  
  # pull out first point
  firstPoint <- allCoord[1,]
  
  # get vector between first and last point - this is the line
  lineVec <- allCoord[nPoints,] - firstPoint
  
  # normalize the line vector
  lineVecN <- lineVec / sqrt(sum(lineVec^2))
  
  # find the distance from each point to the line:
  # vector between all points and first point
  vecFromFirst <- t(t(allCoord) - firstPoint)
  
  # To calculate the distance to the line, we split vecFromFirst into two 
  # components, one that is parallel to the line and one that is perpendicular 
  # Then, we take the norm of the part that is perpendicular to the line and 
  # get the distance.
  # We find the vector parallel to the line by projecting vecFromFirst onto 
  # the line. The perpendicular vector is vecFromFirst - vecFromFirstParallel
  # We project vecFromFirst by taking the scalar product of the vector with 
  # the unit vector that points in the direction of the line (this gives us 
  # the length of the projection of vecFromFirst onto the line). If we 
  # multiply the scalar product by the unit vector, we have vecFromFirstParallel
  scalarProduct <- rowSums(t(t(vecFromFirst) * lineVecN))
  vecFromFirstParallel <- scalarProduct %*% t(lineVecN)
  vecToLine <- vecFromFirst - vecFromFirstParallel
  
  # distance to line is the norm of vecToLine
  distToLine <- sqrt(rowSums(vecToLine^2))
  which.max(distToLine)
}

###
# Compute indices cutoff to get outliers
###
getOutlierCutoff <- function(curve, method=c("tangent", "deriv", "deriv-enh", "deriv.old", "roc")
                              , slope=slope, curveName="", plot=FALSE){
  method <- match.arg(method)
  
  curve.seq <- seq(1 / length(curve), 1, by = 1 / length(curve))
  
  if( method == "tangent" ){
    which_best <- compute_tangent_cutoff(curve, plot = plot)
  } else if( method == "roc" ){ #compute the furthest point between the extreme points
    which_best <- getFurtherPoint(curve)
  }else if( method == "deriv" ){ # compute elbow point using 1st derivative
    #curve.norm2 <- curve / sqrt(sum(curve^2)) # norm-2
    curve.norm2 <- curve / max(curve) # norm infty
    targetPoint <- compute_exp_elbow(curve.norm2, slope=slope, enhanced=F, plot=plot) # normalized curve as param
    which_best <- min(targetPoint, length(curve)+1)
  }else if( method == "deriv-enh" ){ # compute elbow point using 1st derivative
    #curve.norm2 <- curve / sqrt(sum(curve^2)) # norm-2
    curve.norm2 <- curve / max(curve) # norm infty
    targetPoint <- compute_exp_elbow(curve.norm2, slope=slope, enhanced=T, plot=plot) # normalized curve as param
    which_best <- min(targetPoint, length(curve)+1)
  }else if( method == "deriv.old" ){
    curve.norm2 <- curve / max(curve) # norm infty
    elbow.point <- compute_exp_elbow.last(curve.norm2, slope=slope, plot=plot) # normalized curve as param
    metric.seq <- seq(1 / length(curve), 1, by = 1 / length(curve))
    which_best <- min(which(metric.seq > elbow.point), length(curve)+1)
  }
  elbow.point <- curve.seq[which_best]
  cutoff <- curve[which_best]
  if( plot ) {
    plot(curve.seq, curve, type="l", log=""
         , main=paste0("elbow cut @", (1 - elbow.point) * 100
                       , "% [", curveName, ifelse(exists("threshold"), paste0("; >=",threshold), ""), "]"))
    points(elbow.point, cutoff, col="red", pch=4)
  }
  cutoff
}


###
# Mood algorithm implementations
###

# computes the 3 indices for a range of columns ($i)
meanCorLSM.new <- function(i, mtx, means, vars, sds){
  covMtx <- cov(mtx, mtx[,i], use="pairwise.complete.obs")
  tmp <- covMtx / vars
  preOutl <- cbind("shape"=colMeans(covMtx / sds, na.rm=T)
                   , "magnitude"=colMeans(tmp * means, na.rm=T)
                   , "amplitude"=colMeans(tmp, na.rm=T))
  preOutl[,1] <- preOutl[,1] / sds[i] # shape; ~cor()
  preOutl[,2] <- means[i] - preOutl[,2] # beta0; b = y - ax
  preOutl
}

# computes the 3 indices for a range of columns ($i)
meanCorLSM <- function(i, mtx, means, vars, sds){
  preOutl <- apply(cov(mtx, mtx[,i], use="pairwise.complete.obs"), 2
              , function(col) {
                  tmp <- col / vars
                  c("shape"=mean(col / sds, na.rm=T) #PRECOMPUTED
                    , "magnitude"=mean(tmp * means, na.rm=T) # PRECOMPUTED beta0; b, from y=ax+b
                    , "amplitude"=mean(tmp, na.rm=T)) #beta1; a, from y=ax+b
                }
              )
  preOutl[1,] <- preOutl[1,] / sds[i] # shape; ~cor()
  preOutl[2,] <- means[i] - preOutl[2,] # beta0; b = y - ax
  t(preOutl)
}

# computes the LS mean coefficients for a range of columns ($i) [magnitude, amplitude]
# Linear Algebra method
# source: http://stats.stackexchange.com/questions/22718/what-is-the-difference-between-linear-regression-on-y-with-x-and-x-with-y
meanLSM <- function(i, mtx, means, vars){
  beta.1 <- apply(cov(mtx[,i], mtx, use="pairwise.complete.obs"), 1, function(col) col / vars)
  beta.0 <- means[i] - colMeans(apply(beta.1, 2, function(col) col * means))
  cbind(beta.0, colMeans(beta.1))
}

# computes the column mean of the correlation matrix for a range of columns ($i) [shape]
meanCor <- function(i, mtx){
  #y <- mtx[,i]
  #apply(cor(mtx, y, use="pairwise.complete.obs"), 2, mean, na.rm=T)
  #apply(cor(mtx, mtx[,i], use="pairwise.complete.obs"), 2, mean, na.rm=T)
  colMeans(cor(mtx, mtx[,i], use="pairwise.complete.obs"), na.rm=T)
}


# computes the LS mean coefficients for a range of columns ($i) [magnitude, amplitude]
# lm() method
meanLSM.old <- function(i, mtx){
  #y <- mtx[,i]
  matrix(
          #apply(apply(mtx, 2, function(x) lm(y~x)$coeff), 1, mean, na.rm=T),
          #apply(apply(mtx, 2, function(x) lm(mtx[,i]~x)$coeff), 1, mean, na.rm=T),
          rowMeans(apply(mtx, 2, function(x) lm(mtx[,i]~x)$coeff), na.rm=T),
          nrow = length(i), ncol = 2, byrow = T
        )
}

# computes the 3 indices for a range of columns ($i)
meanCorLSM.rcpp <- function(i, mtx2, means, vars, sds){
  preOutl <- c_cov(mtx2, means, vars, sds, i[1], i[length(i)])
  preOutl[,1] <- preOutl[,1] / sds[i]
  preOutl[,2] <- means[i] - preOutl[,2]
  preOutl
}


#
# First parallel version
# Features:
#   -Scalable, parallel
#
computeMoodIndices.old <- function(data, n=ncol(data), nGroups=n / getOption("mc.cores", 1)){
  # compute the column splits/partition for parallel processing
  splits <- split(1:n, sort(rank(1:n) %% max(1, as.integer(n/nGroups))))
  # auxiliar vars
  data.means <- colMeans(data, na.rm=T)
  data.vars <- apply(data, 2, var, na.rm=T)
  # compute the outlier values
  shape <- do.call(c, mclapply(splits, meanCor, data))
  pMean <- do.call(rbind, mclapply(splits, meanLSM, data, data.means, data.vars))
  # set the result column names
  colnames(pMean) <- c("magnitude", "amplitude") # (b, a) from  y = ax + b
  
  data.frame(shape, pMean)
}

#
# Parallel full version
# Features:
#   -Scalable, parallel
#
computeMoodIndices <- function(data, n=ncol(data), nGroups=n / getOption("mc.cores", 1)){
  # compute the column splits/partition for parallel processing
  splits <- split(1:n, sort(rank(1:n) %% max(1, as.integer(n/nGroups))))
  # compute auxiliar support data
  data.means <- colMeans(data, na.rm=T)
  data.vars <- apply(data, 2, var, na.rm=T)
  data.sds <- apply(data, 2, sd, na.rm=T)
  # compute Outliers
  vectors <- do.call(rbind, mclapply(splits, meanCorLSM, data, data.means, data.vars, data.sds))
  # alternative version using a different routine 'meanCorLSM.new'
  #vectors <- do.call(rbind, mclapply(splits, meanCorLSM.new, data, data.means, data.vars, data.sds))

  vectors <- data.frame(vectors)
  colnames(vectors) <- c("shape", "magnitude", "amplitude")
  vectors
}

#
# Classical version using correlation & linear regression (LSM)
# Features:
#   -Scalable, parallel
#
computeMoodIndices.classic <- function(data, n=ncol(data), nGroups=n / getOption("mc.cores", 1)){
  # compute the column splits/partition for parallel processing
  splits <- split(1:n, sort(rank(1:n) %% max(1, as.integer(n/nGroups))))
  # compute the outlier values
  shape <- do.call(c, mclapply(splits, meanCor, data))
  pMean <- do.call(rbind, mclapply(splits, meanLSM.old, data))
  # set the result column names
  colnames(pMean) <- c("magnitude", "amplitude") # (b, a) from  y = ax + b

  data.frame(shape, pMean)
}

#
# Rcpp full version
# features:
#   -Scalable, parallel
#   -Low memory consumption (=> more scalable)
#
computeMoodIndices.rcpp <- function(data, n=ncol(data), nGroups=n / getOption("mc.cores", 1)){
  # compute the column splits/partition for parallel processing
  splits <- split(1:n, sort(rank(1:n) %% max(1, as.integer(n/nGroups))))
  # compute auxiliar support data
  data.means <- colMeans(data, na.rm=T)
  data.vars <- apply(data, 2, var, na.rm=T)
  data.sds <- apply(data, 2, sd, na.rm=T)
  data2 <- t(t(data) - data.means) #pre computed mean-distance data
  # compute indices 
  #sourceCpp("cor_cov_blockwise.cpp") #in case
  vectors <- do.call(rbind, mclapply(splits, meanCorLSM.rcpp
                                     , data2, data.means, data.vars, data.sds))
  vectors <- data.frame(vectors)
  colnames(vectors) <- c("shape", "magnitude", "amplitude")
  vectors
}

# clean NA's
sanitize_data <- function(data){
  t(apply(data, 1, function(x) {
    x[which(is.na(x))] <- mean(x, na.rm=T)
    x
  }))
}

###
# Main function to obtain outliers using the Mood algorithm
###
getOutliers <- function(data, n=ncol(data), nGroups=n / getOption("mc.cores", 1)
                          , method=c("rcpp", "classic", "normal", "old")
                          , outl.method=c("tangent", "deriv", "deriv-enh", "deriv.old", "roc")
                          #, slope=2, benchmark=c(1, 0, 1), plotCutoff=F) {
                          , slope=2, benchmark=c(1, 0, 1)) {
  method <- match.arg(method)
  outl.method <- match.arg(outl.method)
  outl.type <- c("shape", "amplitude", "magnitude")
  
  data <- sanitize_data(data)
  
  # get index values
  if( method == "rcpp" ){
    pre.indices <- computeMoodIndices.rcpp(data, n, nGroups)
  }else if( method == "classic" ){
    pre.indices <- computeMoodIndices.classic(data, n, nGroups)
  }else if( method == "normal" ){
    pre.indices <- computeMoodIndices(data, n, nGroups)
  }else if( method == "old" ){
    pre.indices <- computeMoodsIndices.old(data, n, nGroups)
  }
  # apply benchmark
  indices <- abs(as.matrix(pre.indices - matrix(benchmark, nrow(pre.indices), 3, byrow = T)))
  # compute outliers cutoff
  sapply(outl.type, function(outl.name, plot=F) {
    # sort the curve
    metric <- sort(indices[,outl.name])
    # compute elbow point
    cutoff <- getOutlierCutoff(metric, method=outl.method
                               , slope=slope
                               , curveName=outl.name
                               , plot=plot)
    # filter outliers
    outl <- which(indices[,outl.name] > cutoff)
    outl[order(indices[outl,outl.name], decreasing=T)]
  }, plot=T, simplify=F, USE.NAMES=T)
}

