## General code for  computing spatial distributions.
library(splancs)
library(sjevor)
library(sjedmin)
library(sjedrp)

## Some useful constants
density.um.mm2 <- 1e6                   #change density from um^2 to mm^2

sjespatdists <- function (hpts, xmin, xmax, ymin, ymax, note, plot=F) {
  ## General analysis routine.
  ## hpts[npts,2] is a 2-d array of size NPTS * 2. xmin,ymin, xmax, ymax set
  ## the boundary of the rectangular region where the points lie.
  ## Assume hpts already stores the data points, restricted to just inside the
  ## analysis region.
  ## This used to be called goboth().
  ht <- ymax - ymin
  wid <- xmax - xmin
  
  datapoly <- spoints( c(xmin,ymin, xmin, ymax,   xmax, ymax,  xmax,ymin))
  
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
  
  fs <- hist(validareas, breaks=sje.vorareabins, plot=FALSE, freq=FALSE )
  ahx <- fs$mids; ahy <- cumsum(fs$counts)/ sum(fs$counts)

  ## Central angles:
  ia <- ianglesplot(v$iangles,show=F)

  ## Delaunay segment lengths.
  dellens.acc <- vorcr.dellens(v, v$delacc)
  dellens <- pmin(dellens.acc, sje.dellenmax)
  
  ##cat(paste("dellens number of elements: ", length(dellens), "\n"))
  ds <- hist(dellens, breaks=sje.dellenbins, plot=FALSE, freq=FALSE )
  dellencdf <- cumsum(ds$counts)/ sum(ds$counts)


  res <- list(
       note = note,
       ##retnum = retnum,
       pts = hpts,
       xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax,
       ht = ht, wid = wid,
       ##todo - do we need to return these next two distribs?
       ##sje.steps.max = sje.steps.max,
       ##steps = steps,
       npts = npts,
       g = g, f = f,
       l = l,
       v = v,
       ahx = ahx, ahy = ahy, ia = ia,
       iay = ia$y,                      #return separately so we can
                                        #easily test p-value.
       dellencdf = dellencdf,
       dellens = dellens
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
  

  ## Internal angles
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
  axis(1, labels=T, c(seq(0,sje.dellenmax,25)))
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
