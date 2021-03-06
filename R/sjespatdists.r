## General code for  computing spatial distributions.
library(splancs)
library(sjevor)
library(sjedmin)
library(sjedrp)

## Originally in count.opp.R, now duplicated.
##source("~/mosaics/beta_rgc/count.opp.R")

## Some useful constants
density.um.mm2 <- 1e6                   #change density from um^2 to mm^2

######################################################################
## new functions.

sje.vorarea <- function(pts, w, breaks, need.v=FALSE) {
  ## MAX.AREA is the max area.
  ## BREAKS are the histogram breaks to use.
  ## NEED.V is TRUE if we want to use Voronoi information for
  ## further processing.
  v <- vorcr(pts[,1], pts[,2], w[1], w[2], w[3], w[4])
  
  ## Keep just the non-negative areas.
  a <- v$info[,4]
  validareas <- a[a >= 0.0]

  vd.max <- breaks[ length(breaks)]
  ## Use recyling to set an upper limit for the validareas.
  validareas <- pmin(validareas, vd.max)
  if (max(validareas) > vd.max) {
    warning(paste("maximum area is", max(validareas)))
  }
  
  fs <- hist(validareas, breaks=breaks, plot=FALSE)
  ahy <- cumsum(fs$counts)/ sum(fs$counts)

  if (!need.v)
    v=NULL
  res <- list(x=fs$mids, y=ahy, v=v)
  res
}

sje.dellens <- function(pts, w, v=NULL, ds.breaks) {
  ## If V is already calculated by vorarea, we can reuse that else we
  ## tesselate using pts and w.
  ## 
  if (is.null(v)) {
    v <- vorcr(pts[,1], pts[,2], w[1], w[2], w[3], w[4])
  }
  
  ds.max <- ds.breaks[length(ds.breaks)]
  dellens.acc <- vorcr.dellens(v, v$delacc)
  dellens <- pmin(dellens.acc, ds.max)

  ## Now bin the segments into a histogram.
  ds <- hist(dellens, breaks=ds.breaks, plot=FALSE)
  dellencdf <- cumsum(ds$counts)/ sum(ds$counts)

  ## By naming 
  res <- list (x=ds$mids, y=dellencdf)

}

sje.internalangles <- function(pts, w, v=NULL) {
  ## Compute the internal angles of Voronoi polygon.  If V is already
  ## known, we can re-use that.
  if (is.null(v)) {
    v <- vorcr(pts[,1], pts[,2], w[1], w[2], w[3], w[4])
  }
  
  ia <- ianglesplot(v$iangles,show=F)
}


sjespatdists <- function (pts, w, note, plot=F, param=NULL) {
  ## General analysis routine.
  ## pts[npts,2] is a 2-d array of size NPTS * 2.
  ## w is the bounding box window (xmin, xmax, ymin, ymax)
  ## the boundary of the rectangular region where the points lie.
  ## Assume pts already stores the data points, restricted to just inside the
  ## analysis region.
  ## New version
  ## PARAM is a list of relevant parameters that we want to compute, together
  ## with their arguments (e.g. distance params)
  ## Functions (extra params needed in brackets)
  ## g (steps)
  ## f (steps)
  ## l (steps)
  ## vd (vd.breaks)
  ## ds (ds.breaks)

  ## steps - vector of distances.
  ## vd, ds - like steps, but for Voronoi and Delaunay information.
  

  if( length(w) != 4)
    stop(paste("w (", paste(w, collapse=' '), ") should be of length 4."))
  else {
    xmin=w[1]; xmax=w[2]; ymin=w[3]; ymax=w[4]; 
  }
  
  ht <- ymax - ymin
  wid <- xmax - xmin

  steps <- param$steps
  
  datapoly <- spoints( c(xmin,ymin, xmin, ymax,   xmax, ymax,  xmax,ymin))
  ## Check to see that all points are within the polygon.
  outside <- pip(pts, datapoly, out=T, bound=TRUE)
  if (dim(outside)[1]>0) {
    print("some points outside polygon\n")
    browser()
  }

  null.xylist <-  list(x=NULL, y=NULL)
  npts <- length(pts[,1])
  grid.pts <- gridpts(datapoly, npts) #Diggle+Gratton(84) suggest using ~npts.

  if (!is.null(param$distribs$g))
    g <- list(x=steps, y=Ghat(pts,steps))
  else
    g <- null.xylist

  if (!is.null(param$distribs$f))
    f <- list(x=steps, y=Fhat(pts,grid.pts,steps))
  else
    f <- null.xylist

  if (!is.null(param$distribs$l)) {
    l <- list(x=steps, y= sqrt(khat(pts, datapoly, steps)/pi))
  } else {
    l <- null.xylist
  }

  ##polymap(datapoly); pointmap(pts, add=T)#show the points.

  ## Voronoi areas.
  if (!is.null(param$distribs$vd)) {
    vd <- sje.vorarea(pts, w, param$vd.breaks, need.v=TRUE)
  } else{
    vd <- list(x=NULL, y=NULL, v=NULL)
  }

  if (!is.null(param$distribs$ri)) {
    ri <- calc.ri(pts, w)
  } else{
    ri <- NULL
  }


  ## Central angles:
  if (!is.null(param$distribs$ia)) {
    ia <- sje.internalangles(pts, w, vd$v)
  } else {
    ia <- null.xylist
  }

  ## Delaunay segment lengths. xxx
  if (!is.null(param$distribs$ds)) {
    ds <- sje.dellens(pts, w, vd$v, param$ds.breaks)
  } else {
    ds <- null.xylist
  }

  ## Get NN distances from Voronoi information.
  ##nn <- a$v$info[,3];   nn <- nn[nn>=0.0]

  ## Before returning results, remove Voronoi plot, don't think we
  ## neeed to return that.
  vd$v <- NULL
  
  res <- list(
              note = note,
              pts = pts,
              w = w,
              npts = npts,
              g = g, f = f, l = l,
              vd = vd,
              ia = ia,
              ds = ds,
              ri = ri,
              param=param
              )

  class(res) <- "sjespatdists"

  res

}

plot.sjespatdists <- function(res) {
  ##  Plot the spatial distributions, class function.
  distribs <- res$param$distribs

  par(mfrow=c(2,3), bty='n')

  if (!is.null(distribs$g))
    plot(res$g, type='l', col='red', ylab='g')
  
  if (!is.null(distribs$f))
    plot(res$f, type='l', col='red', ylab='f')

  if (!is.null(distribs$l))
    plot(res$l, type='l', col='red', ylab='l')

  if (!is.null(distribs$vd))
    plot(res$vd, type='l', col='red', ylab='VD')

  if (!is.null(distribs$ia))
    plot(res$ia, type='l', col='red', ylab='IA')

  if (!is.null(distribs$ds))
    plot(res$ds, type='l', col='red', ylab='DS')

}


new.dist.arr <- function( one.dist, nsims) {
  ## NSIMS is the number of simulations planned.
  ## Fill up the arrays and then go for it!
  
  distribs <- one.dist$param$distribs
  null.distribs <- sapply(distribs, is.null)
  if ( any(null.distribs) ){
    distribs <- distribs[ - which(null.distribs)]   #remove NULL elements
  }
      
  arrs <- list()
  for (name in names(distribs)) {

    this.distrib <- one.dist[[name]]    #is a list with (x,y) components.
    mat <- matrix(0, length(this.distrib$y), nrow=nsims+1)
    mat[1,] <- this.distrib$y
    colnames(mat) <- this.distrib$x
    attr(mat, "name") <- name
    ## Set the class for the matrix.
    class(mat) <- "spat.array"
    if(name == "opp")
      class(mat) <- "spat.opp"
    if(name == "ri3")
      class(mat) <- "spat.ri3"

    if(name == "ri")
      class(mat) <- "spat.ri"

    new <- list( mat); names(new) <- name
    arrs <- c(arrs, new)
  }

  list(
       get = function() {
         arrs
       },
       set.row = function(sim.dist, row) {
         ## set one simulation row.
         for (name in names(arrs)) {
           this.distrib <- sim.dist[[name]]
           arrs[[name]][row,] <<- this.distrib$y
         }
       }
       ,
       ranks = function() {
         ## Compute the ranking of each matrix.
         sapply(arrs, ranking)
       }
       ,
       plot = function() {
         ## Plot each matrix.
         sapply(arrs, plot)
       }
       ,
       clear.sim = function() {
         ## Remove the simulations from each matrix.
         ## Not necessary usually, but tidy.
         for (name in names(arrs)) {
           arrs[[name]][-1,] <<- 0
         }
       }
       )
}

plot.spat.array <- function(arr, r=NULL, ylab, ...) {
  ## Plot a spatial distribution (K, F, G).
  ## Real data is shown in red; simulation envelope in black.
  if (is.null(r)) {
    r <- colnames(arr)
  }
  if (missing (ylab)) {
    ylab <- attributes(arr)$name
    if (is.null(ylab))
      ylab <- deparse(substitute(arr))
  }
  plot(r, arr[1,], col='red', type='l', bty='n',
       main=ranking(arr), 
       ylab=ylab, ...)
  ##title(main= paste("p =", ranking(arr)), line=-1)
  lines(r, apply(arr[-1,], 2, min), lty=1)
  lines(r, apply(arr[-1,], 2, max), lty=1)
  lines(r, apply(arr[-1,], 2, mean), lty=1)

  ## The following section is new.
  if (substring(ylab,1,1) == "l") {
    ## this is an L function.
    max.r <- max(as.numeric(r))
    segments(0, 0, max.r, max.r, lwd=1, lty=2)
  }
}


plot.spat.opp <- function(arr,cex=0.5, real.col='black', ...) {
  ## Plot the fraction of opposites
  stripchart(list(arr[-1,1], arr[-1,2], arr[-1,3], arr[-1,4]),
             method="jitter", pch=19, vertical=TRUE,
             ylim=c(min(arr), 1), 
             ##group.names=c(expression(1^{st}), "2", "3", "all"),
             group.names=rep("", 4),
             ylab="fraction opposite type", cex=cex, ...)

  axis(1, at=1:4, labels=c(expression(1^{st}),
                    expression(2^{nd}),
                    expression(3^{rd}),
                    "all"))
  dx <- 0.3; i <- 1:4
  segments(i-dx, arr[1,], i+dx, arr[1,], lwd=0.6, col=real.col)
  median.sim <- apply(arr[-1,], 2, median)
  segments(i-dx, median.sim, i+dx, median.sim, lwd=0.6, lty=2)

}

calc.ri <- function(pts, w) {
  ## Return the regularity index of PTS with bounding box W.  This is
  ## potentially wasteful in that Voronoi information is computed
  ## twice (once here, and once by sjespatdists).
  list(x=1, y=vorcr2(pts, w)$cr)
}

## Bivariate regularity indexes.
calc.ri3 <- function(on, of, w) {
  ## Calculate the regularity index of the ON, OFF and ON+OFF population.
  xmin <- w[1]; xmax <- w[2];
  ymin <- w[3]; ymax <- w[4];
  v.n <- vorcr(on[,1], on[,2], xmin, xmax, ymin, ymax)
  v.f <- vorcr(of[,1], of[,2], xmin, xmax, ymin, ymax)
  b <- rbind(on, of)
  v.0 <- vorcr( b[,1],  b[,2], xmin, xmax, ymin, ymax)

  y <- c(v.n$cr, v.f$cr, v.0$cr)
  x <- c("ON", "OFF", "ON+OFF")
  res <- list(x=x, y=y)
}

plot.spat.ri3 <- function(ri3, cex=0.5, ylim=range(ri3),
                          ylab='regularity index',
                          real.col = 'red', ...) {
  ## Plot the regularity indexes.
  res <- list(on=ri3[-1,1], of=ri3[-1,2],on.off=ri3[-1,3])
  stripchart(res, vert=T, pch=19, method="jitter",
             cex=cex,
             ylim=ylim,
             ##group.names=c("ON", "OFF", "ON+OFF"),
             group.names=c("1", "2", "1+2"),
             main="",
             ylab=ylab)
  
  median.sim <- apply(ri3[-1,], 2, median)
  i <- 1:3; dx <- 0.3;
  segments(i-dx, ri3[1,], i+dx, ri3[1,], lwd=0.6, col=real.col)
  segments(i-dx, median.sim, i+dx, median.sim, lwd=0.6, lty=2)
  ##legend(x=1, y=3.5, lty=c(1,2),
  ##       legend=c("experimental RI", "median RI of sims"))
  
}

plot.spat.ri <- function(ri, cex=0.5, ylim=range(ri),
                          real.col = 'red', ...) {
  ## Plot the regularity index for univariate pattern.
  res <- ri[-1,1]
  stripchart(res, vert=T, pch=19, method="jitter",
             cex=cex,
             ylim=ylim,
             main="",
             ylab="regularity index")
  
  median.sim <- median(res)
  i <- 1:3; dx <- 0.3;
  segments(i-dx, ri[1,], i+dx, ri[1,], lwd=0.6, col=real.col)
  segments(i-dx, median.sim, i+dx, median.sim, lwd=0.6, lty=2)
  
}

ranking <- function(arr) {
  ## Evaluate fit of row 1 (the real data) with remaining rows (simulations)
  ## using equations 6.3--6.5 from (Diggle, 2002, page 89).
  n.s <- nrow(arr)
  u <- rep(0, n.s)
  for (i in 1:n.s) {
    ave.i <- apply( arr[-i,], 2, sum) / (n.s - 1)
    u[i] <- sum((ave.i - arr[i,])^2)
    ##u[i] <- max( abs( ave.i - arr[i,]) )
  }
  signif( (rank(u)[1])/n.s, 5)
}

sjespatdists.biv <- function (pts1, pts2, w, note, plot=F, param=NULL) {
  ## GENERAL ANALYSIS ROUTINE, BIVARIATE VERSION.
  ## pts[npts,2] is a 2-d array of size NPTS * 2.
  ## w is the bounding box window (xmin, xmax, ymin, ymax)
  ## the boundary of the rectangular region where the points lie.
  ## Assume pts already stores the data points, restricted to just inside the
  ## analysis region.
  ## New version
  ## PARAM is a list of relevant parameters that we want to compute, together
  ## with their arguments (e.g. distance params)

  if( length(w) != 4)
    stop(paste("w (", paste(w, collapse=' '), ") should be of length 4."))
  else {
    xmin=w[1]; xmax=w[2]; ymin=w[3]; ymax=w[4]; 
  }
  
  ht <- ymax - ymin
  wid <- xmax - xmin


  ## This routine was written for analysing bivariate patterns where both
  ## sets were of roughtly the same size.  (e.g. on and off RGCs occur in
  ## roughly equal numbers). However, sometimes we have different sizes
  ## in the two types of cell, (e.g. horizontal cells).

  
  either <- function(v2, v1) {
    ## If v2 is not null, return v2, else return v1.
    ## Simple helper function for allowing type 1 and type 2 axes to differ.
    if (is.null(v2))
      v1
    else
      v2
  }

  
  steps   <- param$steps
  steps2  <- either(param$steps2,steps)
  steps12 <- either(param$steps12,steps)
  steps0  <- either(param$steps0,steps)

  pts0 <- rbind(pts1, pts2)
  
  datapoly <- spoints( c(xmin,ymin, xmin, ymax,   xmax, ymax,  xmax,ymin))
  ## Check to see that all points are within the polygon.
  outside <- pip(pts0, datapoly, out=T, bound=TRUE)
  if (dim(outside)[1]>0) {
    print("some points outside polygon\n")
    browser()
  }

  
  null.xylist <-  list(x=NULL, y=NULL)

  g0 <- g1 <- g2 <- null.xylist
  if (!is.null(param$distribs$g0))
    g0 <- list(x=steps, y=Ghat(pts0,steps))
  
  if (!is.null(param$distribs$g1))
    g1 <- list(x=steps, y=Ghat(pts1,steps))
  
  if (!is.null(param$distribs$g2))
    g2 <- list(x=steps, y=Ghat(pts2,steps))


  f0 <- f1 <- f2 <- null.xylist
  if (!is.null(param$distribs$f0))
    f0 <- list(x=steps2, y=Fhat(pts0, gridpts(datapoly, dim(pts0)[1]),steps2))

  if (!is.null(param$distribs$f1))
    f1 <- list(x=steps2, y=Fhat(pts1, gridpts(datapoly, dim(pts1)[1]),steps2))
  
  if (!is.null(param$distribs$f2))
    f2 <- list(x=steps2, y=Fhat(pts2, gridpts(datapoly, dim(pts2)[1]),steps2))

  l0 <- l1 <- l2 <- l12 <- null.xylist
  if (!is.null(param$distribs$l0))
    l0 <- list(x=steps0, y= sqrt(khat(pts0, datapoly, steps0)/pi))
    
  if (!is.null(param$distribs$l1))
    l1 <- list(x=steps, y= sqrt(khat(pts1, datapoly, steps)/pi))

  if (!is.null(param$distribs$l2))
    l2 <- list(x=steps2, y= sqrt(khat(pts2, datapoly, steps2)/pi))

  if (!is.null(param$distribs$l12))
    l12 <- list(x=steps12, y=sqrt(k12hat(pts1, pts2, datapoly, steps12)/pi))

  ## Voronoi areas.
  vd0 <- vd1 <- vd2 <- null.xylist
  if (!is.null(param$distribs$vd0))
    vd0 <- sje.vorarea(pts0, w, param$vd0.breaks, need.v=TRUE)

  if (!is.null(param$distribs$vd1))
    vd1 <- sje.vorarea(pts1, w, param$vd1.breaks, need.v=TRUE)

  if (!is.null(param$distribs$vd2))
    vd2 <- sje.vorarea(pts2, w, either(param$vd2.breaks,param$vd1.breaks),
                                       need.v=TRUE)

                                 
  ## Central angles:
  ia0 <- ia1 <- ia2 <- null.xylist
  if (!is.null(param$distribs$ia0))
    ia0 <- sje.internalangles(pts0, w, vd0$v)

  if (!is.null(param$distribs$ia1)) 
    ia1 <- sje.internalangles(pts1, w, vd1$v)

  if (!is.null(param$distribs$ia2))
    ia2 <- sje.internalangles(pts2, w, vd2$v)


  ## Delaunay segment lengths.
  ds0 <- ds1 <- ds2 <- null.xylist
  if (!is.null(param$distribs$ds0))
    ds0 <- sje.dellens(pts0, w, vd0$v, param$ds0.breaks)

  if (!is.null(param$distribs$ds1))
    ds1 <- sje.dellens(pts1, w, vd1$v, param$ds1.breaks)

  if (!is.null(param$distribs$ds2))
    ds2 <- sje.dellens(pts2, w, vd2$v,
                       either(param$ds2.breaks,param$ds1.breaks))

  ## Before returning results, remove Voronoi plot, don't think we
  ## neeed to return that.
  vd0$v <- NULL;  vd1$v <- NULL; vd2$v <- NULL


  ## Could reuse v0 here.  Or should pass in w.
  opp <- NULL
  if (!is.null(param$distribs$opp)) {
    opp <- count.opp(pts1, pts2, w)$res
  }

  ri3 <- NULL
  if (!is.null(param$distribs$ri3)) {
    ri3 <- calc.ri3(pts1, pts2, w)
  }
  
  res <- list(
              note = note,
              pts0 = pts0, pts1=pts1, pts2=pts2,
              w = w,
              g0=g0, g1=g1, g2=g2,
              f0=f0, f1=f1, f2=f2,
              l0=l0, l1=l1, l2=l2, l12=l12,
              vd0=vd0, vd1=vd1, vd2=vd2,
              ia0=ia0, ia1=ia1, ia2=ia2,
              ds0=ds0, ds1=ds1, ds2=ds2,
              opp=opp, ri3=ri3,
              param=param
              )

  class(res) <- "sjespatdistsbiv"

  res
}

plot.sjespatdistsbiv <- function(x) {

  distribs <- names(x$param$distribs)
  print(distribs)
  for (distrib in distribs) {
    if (!is.null(x$param$distribs[[distrib]]))
      plot(x[[distrib]], xlab='', ylab=distrib, type='l', col='red')
  }
}

count.opp <- function(a, b, w=NULL, which.labels=0, plot=FALSE,
                      debug=FALSE) {
  ## A, B are matrices of the cell points
  ## w is the bounding box of the points.
  
  labels <- c( rep(1, dim(a)[1]), rep(2, dim(b)[1]) )

  ## combined matrix of all points.
  c0 <- rbind(a, b)
  c0.n <- dim(c0)[1]

  if (is.null(w)) {
    w <- real(4)
    w[1:2] <- range(c0[,1]); w[3:4] <- range(c0[,2])
  }
  
  v0 <- vorcr( c0[,1], c0[,2], w[1], w[2], w[3], w[4])
              
  if (plot) {
    plot(v0, type='n')
    text(c0, cex=0.8, col=ifelse(labels==1, "green", "orangered"))
  }
    
  ## Find the neighs, and convert -1 to NA
  neighs <- v0$neighs
  neighs[which(neighs == -1, arr.ind=T)] <- NA

  ## Also remove any neighbour information for cells that we want to ignore.
  if (which.labels > 0) {
    neighs[-which(labels==which.labels),] <- NA
  }
  
  ## convert to 1/2 for type.
  opps.label <- t(apply(neighs, 1, function(x) { labels[x]}))
  opps.same <- t(apply(cbind(labels, opps.label), 1,
                       function(x) {x[1] == x[-1]}))
  counts.opp <- apply(opps.same, 2, function(c) { length(which(c==FALSE))})
  
  valid.neighs <- apply(opps.same, 2, function(c) { length(which(!is.na(c)))})
  
  opps.percent <- counts.opp / valid.neighs
  names(opps.percent) <- paste("nn", 1:length(opps.percent), sep='')
  
  ## For each of nth neighbours, print % of cells that were opposite class.
  ## e.g. the 2nd element shows the % of 2nd nearest neighbours that
  ## are opposite type.
  if (debug)
    print(opps.percent)

  ## Now check each row, to see how many per cell are opposite type.
  opps.cell <-  apply(opps.same, 1, function(c) { length(which(c==FALSE))})

  cell.opps.info <- cbind(1:length(labels),
                          v0$numneighs, 
                        opps.cell,
                          opps.cell / v0$numneighs
                          )

  ## remove rows that we don't want.
  cell.opps.info <- cell.opps.info[ - which(opps.cell == 0),]
  colnames(cell.opps.info) <- c("id", "num neighs", "num opp", "percent opp")
  ## should not just take mean since this distribution is not normal.
  summary(cell.opps.info)

  ## although better to compute median, it produces more stratified results
  ## than the mean.
  mean.nn <- mean(cell.opps.info[,4])
  res <- list(x=c(1,2,3,0), y=c(opps.percent[1:3], mean.nn))
  list(res=res, other=list(
                  v0=v0,
                  cell.opps.info = cell.opps.info,
                  opps.percent=opps.percent,
                  nvalid=dim(cell.opps.info)[1]))
}

######################################################################
## older functions.

sjespatdists.old <- function (hpts, xmin, xmax, ymin, ymax, note, plot=F) {
  ## General analysis routine.
  ## hpts[npts,2] is a 2-d array of size NPTS * 2. xmin,ymin, xmax, ymax set
  ## the boundary of the rectangular region where the points lie.
  ## Assume hpts already stores the data points, restricted to just inside the
  ## analysis region.
  ## This used to be called goboth().
  ht <- ymax - ymin
  wid <- xmax - xmin
  
  datapoly <- spoints( c(xmin,ymin, xmin, ymax,   xmax, ymax,  xmax,ymin))
  ## Check to see that all points are within the polygon.
  outside <- pip(hpts, datapoly, out=T, bound=TRUE)
  if (dim(outside)[1]>0) {
    print("some points outside polygon\n")
    browser()
  }
  
  npts <- length(hpts[,1])
  grid.pts <- gridpts(datapoly, npts) #Diggle+Gratton(84) suggest using ~npts.
  
  g <- Ghat(hpts,sje.steps)
  f <- Fhat(hpts,grid.pts,sje.steps)
  k <- khat(hpts, datapoly, sje.steps);
  l <- sqrt(k/pi)

  ##polymap(datapoly); pointmap(hpts, add=T)#show the points.

  ## Voronoi analysis now...
  v <- vorcr(hpts[,1], hpts[,2], xmin,xmax, ymin,ymax)
  
  ## Keep just the non-negative areas.
  a <- v$info[,4]
  validareas <- a[a >= 0.0]

  ## Use recyling to set an upper limit for the validareas.
  validareas <- pmin(validareas, sje.areamax)
  if (max(validareas) > sje.areamax) {
    warning(paste("maximum area is", max(validareas)))
  }
  
  fs <- hist(validareas, breaks=sje.vorareabins, plot=FALSE)
  ahx <- fs$mids; ahy <- cumsum(fs$counts)/ sum(fs$counts)

  ## Central angles:
  ia <- ianglesplot(v$iangles,show=F)

  ## Delaunay segment lengths.
  dellens.acc <- vorcr.dellens(v, v$delacc)
  dellens <- pmin(dellens.acc, sje.dellenmax)
  
  ##cat(paste("dellens number of elements: ", length(dellens), "\n"))
  ds <- hist(dellens, breaks=sje.dellenbins, plot=FALSE )
  dellencdf <- cumsum(ds$counts)/ sum(ds$counts)

  ## Get NN distances from Voronoi information.
  nn <- v$info[,3]
  nn <- nn[nn>=0.0]

  res <- list(
              note = note,
              ##retnum = retnum,
              pts = hpts,
              xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax,
              ht = ht, wid = wid,
              npts = npts,
              g = g, f = f,
              l = l,
              v = v,
              ahx = ahx, ahy = ahy, ia = ia,
              areas=validareas,
              iay = ia$y,                      #return separately so we can
                                        #easily test p-value.
              dellencdf = dellencdf,
              dellens = dellens,
              nn=nn
              )

  if (plot) 
    plotsumresp(s=NULL, res, title=NULL, ps=NULL,show.means=F)

  res

}

## Ranking and summary functions

sje.ranking <- function (distribs, method=1) {
  ## Sun 07 Jan 2001:
  ## distribs is a 2-d array; each row is a different distribution.
  ## Compare row 1 with the other rows.
  ## Method: choice of function to compare two distributions:
  ## 1: max absolute distance.
  ## 2: integrated distance.
  ##
  ## Simple test:
  ## raw <- matrix ( c( 1,2,1,4,  5,6,3,7,  7,2,4,8), nrow=3, ncol=4,byrow=T)
  ## s <- ranking(raw)


  rankabsdistint <- function (x, y) {
    ## Compute the integrated distance between the two distributions, x and y.
    ## Helper function.
    if (length(x) != length(y))
      stop("rankabsdistint: arguments of different lengths")
    sum( abs(x-y))
  }

  rankabsdistmax <- function (x, y) {
    ## Compute  maximum abs distance between two distributions, x and y.
    if (length(x) != length(y))
      stop("rankabsdistint: arguments of different lengths")
    max( abs(x-y))
  }
  
  nres <- dim(distribs)[1]              #number of rows = # of distributions.
  score <- double(nres)                 #will store all the scores.
  for (i in 1:nres) {
    ## Take the mean of everything except row i.
    copydistribs <- distribs
    copydistribs[i,] <- 0               #zero row i so it doesn't affect sum.
    sums <- apply(copydistribs, 2, sum)
    means <- sums / (nres-1)

    ## try couple of different measures for ranking.
    if (method==1) {
      score[i] <- rankabsdistmax( distribs[i,], means)
    } else {
      score[i] <- rankabsdistint( distribs[i,], means)
    }
  }
  ## now rank the scores and compute p-val etc.
  r <- rank(score)[1]
  p <- r / nres;
  list (p = p, r =r, score=score)
}

sje.findallranks <- function(exptl, sims, ...)
{
  ## Find the ranking probability of the exptl versus several
  ## simulations for each of the distributions examined.  List
  ## returned giving the goodness of fit probability of the exptl
  ## mosaic within the simulations.
  allfields <- c("g", "f", "l", "ahy", "iay", "dellencdf")
  nfields <- length(allfields)
  nreps <- length(sims)
  ps <- double(length=nfields)
  for (f in 1:nfields) {
    field <- allfields[f]
    exptldata <- exptl[[field]]
    distribs <- matrix( nrow = 1 + nreps, ncol = length(exptldata))
    
    distribs[1,] <- exptldata
    for (i in 1:nreps) {distribs[i+1,] <- sims[[i]][[field]]}
    s <- sje.ranking(distribs, ...)
    ps[f] <- s$p
  }
  ps
}

sje.makesum <- function(results)  {
  ## Make the mean summary of the results from each field.

  ## Helper function
  myenv <- function(a) {
    ## find upper,mean,lower bounds of the array A for drawing envelope.
    u <- apply(a, 2, max);                #upper
    l <- apply(a, 2, min);                #lower
    m <- apply(a, 2, mean);               #mean
    list(m = m, u = u, l = l)
  }

  allg <-  t(sapply(results, function(x) x$g))
  allf <-  t(sapply(results, function(x) x$f))
  alll <-  t(sapply(results, function(x) x$l))
  alla <-  t(sapply(results, function(x) x$ahy))
  alli <-  t(sapply(results, function(x) x$ia$y))
  alldl<-  t(sapply(results, function(x) x$dellencdf))  

  genv <- myenv(allg);
  fenv <- myenv(allf);
  lenv <- myenv(alll);
  aenv <- myenv(alla);
  ienv <- myenv(alli);
  dlenv <- myenv(alldl);

  list(genv = genv, fenv = fenv, lenv = lenv, aenv=aenv, ienv = ienv,
       dlenv = dlenv)
}

plotsumresp <- function (s=NULL, r, title=NULL, ps=NULL,
                         show.means=F)
{
  ## Plot the summary (s) of several mosaics versus one result (r).
  ## Do the significance testing as well.
  ## If S is not given, summary is not drawn.
  ## ps stores the p values for each plot.  If PS is null, no
  ## probabilities are plotted.
  ## If SHOW.MEANS is true, the mean of the summary (as well as the
  ## envelope) is drawn.

  par(mfrow=c(2,3))

  show.envelope <- !is.null(s)
  show.p <- !is.null(ps)

  ##  G function.
  plot(sje.steps, r$g,
       col="red", type="l",
       xlab = expression(paste("distance (", mu,"m)")),
       ylab = "G: cumulative probability of NND",axes=F)
  axis(1, labels=T, sje.steps.tics)
  axis(2, labels=T, c(seq(0,1,0.2)))

  if (show.envelope) {
    lines(sje.steps, s$genv$u);lines(sje.steps, s$genv$l)
    if (show.means) lines(sje.steps, s$genv$m)
  }
  title(paste("m", r$note,
              ##"e",  if (show.envelope) deparse(substitute(s)), title))
              "e",  title))
  if (show.p) title(sub=paste("p ", ps[1]))


  ## F function.
  plot(sje.steps, r$f, xlab = expression(paste("distance (", mu,"m)")),
       ylab = "F: cumulative probability of distance to grid points",
       col="red", type="l",
       axes=F)
  if (show.envelope) {
    lines(sje.steps, s$fenv$u);lines(sje.steps, s$fenv$l)
    if (show.means) lines(sje.steps, s$fenv$m)
  }
  if (show.p) title(sub=paste("p ", ps[2]))
  axis(1, labels=T, sje.steps.tics)
  axis(2, labels=T, c(seq(0,1,0.2)))


  ## L function.
  plot(sje.steps, r$l, xlab = expression(paste("distance (", mu,"m)")),
       col="red", type="l",
       ylab = "L function", ylim=c(0, sje.steps.max), axes=F);
  if (show.envelope) {
    lines(sje.steps, s$lenv$u);lines(sje.steps, s$lenv$l)
    if (show.means) lines(sje.steps, s$lenv$m)
  }
  maxstep <- sje.steps[length(sje.steps)];
  lines( c(0, maxstep), c(0, maxstep))
  
  if (show.p) title(sub=paste("p ", ps[3]))
  axis(1, labels=T, sje.steps.tics)
  axis(2, labels=T, sje.steps.tics)


  ## Voronoi area
  plot(sje.vorareax,
       col="red", type="l",
       r$ahy,xlab = expression(paste("Voronoi area (",mu,m^2,")")),
       ylab = "cumulative probability", axes=F);
  if (show.envelope) {
    lines(sje.vorareax, s$aenv$u); lines(sje.vorareax, s$aenv$l);
    if (show.means) lines(sje.vorareax, s$aenv$m)
  }
  if (show.p) title(sub=paste("p ", ps[4]))
  axis(1, labels=T, sje.vorarea.tics)
  axis(2, labels=T, c(seq(0,1,0.2)))
  

  ## Internal angles.
  plot(sje.ianglesx, r$ia$y, ylim=c(0,1),
       col="red", type="l",
       xlab = expression(paste("internal angle (",
      degree,")")), ylab="cumulative probability", axes=F);
  axis(1, labels=T, c(seq(0,180,30)))
  axis(2, labels=T, c(seq(0,1,0.2)))
  if (show.envelope) {
    lines(sje.ianglesx, s$ienv$u); lines(sje.ianglesx, s$ienv$l)
    if (show.means) lines(sje.ianglesx, s$ienv$m)
  }
  if (show.p) title(sub=paste("p ", ps[5]))


  ## Delauanay segment length.
  plot(sje.dellenx, r$dellencdf, ylim=c(0,1),
       col="red", type="l",
       xlab = expression(paste("Delaunay segment length (",mu, m, ")")),
       ylab="cumulative probability", axes=F)
  if (show.envelope) {
    lines(sje.dellenx, s$dlenv$u); lines(sje.dellenx, s$dlenv$l);
    if (show.means) lines(sje.dellenx, s$dlenv$m)
  }
  axis(2, labels=T, c(seq(0,1,0.2)))
  axis(1, labels=T, c(seq(0,sje.dellenmax,length=5)))
  if (show.p) title(sub=paste("p ", ps[6]))
}


## * Printing out the exptl measures into a data.frame
summary.exptl <- function(res, file) {
  ## Helper functions.
  num.or.na <- function(x) { if (is.null(x)) NA else x}
  minvec <- function(vec) { if (is.null(vec)) NA else min(vec[which(vec>0)])}
  maxvec <- function(vec) { if (is.null(vec)) NA else max(vec[which(vec>0)])}

  onpts <- sapply(res, function(x) num.or.na(x$onpts))
  owid  <- sapply(res, function(x) num.or.na(x$owid))
  oht   <- sapply(res, function(x) num.or.na(x$oht))
  oden  <- signif(sapply(res, function(x) num.or.na(x$oden)), digits=5)
  npts  <- sapply(res, function(x) num.or.na(x$npts))

  den  <- signif(sapply(res, function(x) num.or.na(x$den)), digits=5)
  cr    <- signif(sapply(res, function(x) num.or.na(x$v$cr)), 3)
  pf    <- sapply(res, function(x) num.or.na(x$pf))
  effrad <- sapply(res, function(x) num.or.na(x$effrad))
  
  minnnd <- signif(sapply(res, function(x) minvec(x$v$info[,3])),3)
  maxnnd <- signif(sapply(res, function(x) maxvec(x$v$info[,3])),3)
  
  minarea <- signif(sapply(res, function(x) minvec(x$v$info[,4])),3)
  maxarea <- signif(sapply(res, function(x) maxvec(x$v$info[,4])),3)

  mindell <- signif(sapply(res, function(x) minvec(x$dellens)),3)
  maxdell <- signif(sapply(res, function(x) maxvec(x$dellens)),3)
  
  mh <- data.frame(onpts=onpts, owid=owid, oht=oht, oden=oden,
                   npts=npts, den=den, cr=cr,
                   minnnd=minnnd, maxnnd=maxnnd,
                   minarea=minarea, maxarea=maxarea,
                   mindell=mindell, maxdell = maxdell,
                   pf=pf, effrad=effrad,
                   row.names=sapply(res, function(x) x$note))
  ## todo - rownames may not work with Hor cell data if
  ## x$note not defined ...
  write.table(mh, file, sep='\t', quote=F, col.names=NA)
  ##write.table(mh, file, sep='\t', quote=F)
  mh
}
