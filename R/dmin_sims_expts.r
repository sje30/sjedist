## Typical run for the dmin analysis.
## Can run in batch on unix with  R BATCH dmin_sims_expts.r
source("sjespatdists.r")

sizes.gs <- function() {
  ## Set the sizes for the various distributions.
  ## To analyse other data sets, these variables probably should be adapted,
  ## best to take a copy of this function and call it...

  ## Distances at which G, F, L are evaluated.
  sje.steps.max <<- 45 #was 60 um.
  sje.steps <<- seq(0, sje.steps.max, by=1)  #of length sje.steps.max+1.
  sje.steps.tics <<- seq(0, sje.steps.max, 15) #tic marks on graph.

  ## Voronoi areas.
  sje.areamax <<- 2000
  sje.vorareabins <<- seq(0,sje.areamax, length=50) #was 20
  sje.vorareax <<- hist(0, breaks=sje.vorareabins, plot=F, freq=F)$mids
  sje.vorarea.tics <<- seq(0, sje.areamax, 500)

  ## Just generate the x-points for the interior angles.
  sje.ianglesx <<- ianglesplot(0, show=FALSE)$x

  ## Delaunay segment lengths
  sje.dellenmax <<- 75                         
  ## evaulate distance every um.
  sje.dellenbins <<- seq(0, sje.dellenmax, length=sje.dellenmax+1)
  sje.dellenx <<- hist(0,  breaks=sje.dellenbins, plot=FALSE, freq=FALSE )$mids

  ## DRP -- not used here.
  sje.drp.nbins <<- 20;  sje.drp.r <<- 10
}

sizes.gs()                              #set the sizes for the analysis.

## This is the location of an experimental data file.
exptl.file <- "gs_e8.1.dat"

## Read in the experimental data point, and save as a (N,2) matrix
## where N is the number of points.
exptl.pts <- as.matrix(read.table(exptl.file, header=F,sep="\t"))
plot(exptl.pts, asp=1, main=exptl.file)

## get data max and min values.
xmin <- floor(min(exptl.pts[,1]));   ymin <- floor(min(exptl.pts[,2]))
xmax <- ceiling(max(exptl.pts[,1])); ymax <- ceiling(max(exptl.pts[,2]))

note <- "gs field e8.1"
expt.res <- sjespatdists(exptl.pts, xmin, xmax, ymin, ymax, note, plot=T)


######################################################################
## Try a simulation
## p1,p2 are values of dmin mean and s.d.
p1 <- 8; p2 <- 1.5
nreps <- 99
simresults <- list()
for (i in 1:nreps) {
  npts <- expt.res$npts
  dsim <- dminl(npts=npts, wid=expt.res$wid,
                ht = expt.res$ht, dmin=p1, dminsd=p2)
  dsimpts <- cbind(dsim$x, dsim$y)
  simresults[i] <- list(sjespatdists(dsimpts,
                                     0, expt.res$wid,
                                     0, expt.res$ht,
                                     dsim$note, plot=F))
}

ps <- sje.findallranks(expt.res, simresults); ps
simsum <- sje.makesum(simresults)
plotsumresp(simsum, expt.res, simresults[[1]]$note, ps=ps,show.means=T)
