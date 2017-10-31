
rm(list=ls())
# Install and library the following R packages
install.packages("TDA")
install.packages("gdistance")
install.packages("raster")
install.packages("vegan")
install.packages("shape")
install.packages('R.matlab')
biocLite('destiny')

library("TDA")
library("gdistance")
library("raster")
library("vegan")
library("shape")
library('R.matlab')
library("destiny")


# Read the data files
path<-('～/densitypath')
pathname<-file.path(path,'testdata.mat')
treedata<-readMat(pathname)
X1<-as.matrix(treedata$testdata)
XX<-X1[,1:2]

# Set the parameters of density clustering
k<-50
h<-0.3

# Set the path of graphic output
SI_fig_dir <- "～/densitypath"


# the function of DensityPath 
DensityPath<-function(XX,k,h,SI_fig_dir){
  densitypath<-list()
  
  # step1:Reduce the dimensionality of scRNAseq data
  if(ncol(XX)!=2){
    XX<-DiffusionMap(XX,k=2)
  }
  #######################################################################
  # step2:Estimate the density function (landscape) and the level sets.
  
  # call the clusterTree function in the TDA package, then get the density clusters idKDE
  TreeKDE <- clusterTree(XX, k = k, h = h, density = "kde",
                         printProgress = FALSE)
  densityKDE<-TreeKDE$density
  
  
  idKDE<-setdiff(TreeKDE$id,TreeKDE$parent)
  numKDEleaves<-length(idKDE)
  
  ########################################################################
  
  # stpe3:Select high density clusters as the representative cell states
  l<-1
  KDEdensitypeaks<-matrix(1,numKDEleaves,2)
  for (i in idKDE){
    clusterkde<-TreeKDE$DataPoints[[i]]
    KDEdensitypeaks[l,]<-XX[clusterkde[which.max(densityKDE[clusterkde])],]
    l<-l+1
  }
  
  #########################################################################
  # step4:Construct the cell state-transition path
  
  # divide the grid
  xmin <-min(XX)
  xmax<-max(XX)
  ymin<-min(XX)
  ymax<-max(XX)
  
  
  numcell<-nrow(XX)
  
  if (numcell<=10000)
  {
    numx<-101
    numy<-101
  }
  
  if(numcell>10000)
  {
    numx<-501
    numy<-501
  }
  
  
  Xlim <- c(xmin,xmax);  Ylim <- c(ymin, ymax);  
  by <- (xmax-xmin)/(numx-1)
  
  
  Xseq <- seq(Xlim[1], Xlim[2], by = by)
  Yseq <- seq(Ylim[2],Ylim[1], by = -by)
  Grid <- expand.grid(Xseq, Yseq)
  
  # calculate the density on the grid,and form a density surface
  KDE <- kde(X = XX, Grid = Grid, h = h)
  
  r <- raster(nrows=numy, ncols=numx, xmn=xmin, xmx=xmax, ymn=ymin, ymx=ymax,crs="+proj=utm +units=m")
  r[] <- KDE
  T <- transition(r, function(x) mean(x), 8)
  T <- geoCorrection(T)
  C <-KDEdensitypeaks
  D <-KDEdensitypeaks
  
  # use the costDistance function in the gdistance package to calculate the geodesic distance between any two density peaks of the mesh surface
  dis<- costDistance(T, C, D)
  
  # calculate the minimun spanning tree by calling the spantree function in the vegan package using the geodesic distance matrix.
  spanningtree <- spantree(dis)
  
  #######################################################################
  # Output of the graph
  setwd(SI_fig_dir)
  pdf("DensityPath.pdf")
  par(mfrow = c(2,2))
  par(mai = c(0.7,0.4,0.4,0.4),oma=c(0.7,0.7,0.7,0.7))
  
  # 2-dimensional projection
  plot(XX, pch = 19, cex = 0.6, main = "2D mapping of single cell points",xlab="",ylab="")
  
  # 3-dimensional density surface
  Xseq <- seq(Xlim[1], Xlim[2], by = by)
  Yseq <- seq(Ylim[1],Ylim[2], by = by)
  Grid <- expand.grid(Xseq, Yseq)
  zh<- matrix(KDE,ncol=length(Yseq),nrow=length(Xseq))
  
  KDE <- kde(X = XX, Grid = Grid, h = h)
  op <- par(bg = "white")
  res<-persp(Xseq, Yseq, zh, theta = 345, phi = 30,
             expand = 0.5, col = "white",d =1,
             r=90,
             ltheta = 90,
             shade = 0, 
             ticktype = "detailed",
             #xlab = "X", ylab = "Y", zlab = "Sinc( r )" ,
             box = TRUE,
             border=NA,main="Density landscape",
             axes=F
  )
  
  # discrete density clusters
  plot(1:5,1:5,xlim=c(xmin,xmax),ylim=c(ymin,ymax),type = "n",main = "Density clusters",xlab="",ylab="")
  c<-1
  for (i in idKDE){
    points(matrix(XX[TreeKDE[["DataPoints"]][[i]], ], ncol = 2), col = c,pch = 19, cex = 0.6)
    c<-c+1
  }
  
  # density path
  p<-matrix(1,2,2)
  plot(r,xlim=c(xmin,xmax),ylim=c(ymin,ymax),main="Density path",xlab="",ylab="")
  KDEidspantree<-spanningtree$kid
  numKDEedge<-numKDEleaves-1
  KDEspantree<-matrix(1,numKDEedge,2)
  KDEspantree[,1]<-as.matrix(2:numKDEleaves,numKDEleaves,1)
  KDEspantree[,2]<-t(KDEidspantree)
  for (i in 1:numKDEedge){
    p1<-KDEspantree[i,1]
    p2<-KDEspantree[i,2]
    p[1,]<-KDEdensitypeaks[p1,]
    p[2,]<-KDEdensitypeaks[p2,]
    #lines(p, col="red", lwd=2)
    p1top2 <- shortestPath(T, p[1,], p[2,], output="SpatialLines")
    lines(p1top2, col="red", lwd=2)
  }
  for (i in 1:numKDEleaves){
    points(KDEdensitypeaks[i,1],KDEdensitypeaks[i,2],pch = 19, cex = 1,col="blue")
  }
  
  minspantreepath<-list()
  for (i in 1:numKDEedge){
    p1<-KDEspantree[i,1]
    p2<-KDEspantree[i,2]
    p[1,]<-KDEdensitypeaks[p1,]
    p[2,]<-KDEdensitypeaks[p2,]
    #lines(p, col="red", lwd=2)
    p1top2 <- shortestPath(T, p[1,], p[2,], output="SpatialLines")
    adjpath<-as.matrix(((((p1top2@lines)[[1]])@Lines)[[1]])@coords)
    nr<-nrow(adjpath)
    a1<-min(adjpath[,1])-by
    a2<-max(adjpath[,1])+by
    b1<-min(adjpath[,2])-by
    b2<-max(adjpath[,2])+by
    Xseq1 <- seq(a1, a2, by = by)
    Yseq1 <- seq(b2, b1, by = -by)
    Grid1 <- expand.grid(Xseq1, Yseq1)
    KDE1 <- kde(X = XX, Grid = Grid1, h = h)
    pathpoints<-matrix(0,nrow(adjpath),3)
    pathpoints[,1:2]<-adjpath
    for(m in 1:nr){
      l<-which(abs(adjpath[m,1]-Grid1[,1])+abs(adjpath[m,2]-Grid1[,2])==min(abs(adjpath[m,1]-Grid1[,1])+abs(adjpath[m,2]-Grid1[,2])))
      pathpoints[m,3]<-KDE1[l]
    }
    minspantreepath<-c(minspantreepath,list(pathpoints))
  }
  dev.off()
  
  #######################################################################                
  # Return value                
  # return the density path, which contains the estimated density of the sample points(densityKDE),
  #               the two-dimensional coordinates of the density peaks(KDEdensitypeaks), 
  #               the geodesic distance between the density peaks(dis), 
  #               and the path of the minimum spanning tree on the three-dimensional density surface(minspantreepath).
  densitypath<-c(densitypath,densityKDE=list(densityKDE))
  densitypath<-c(densitypath,KDEdensitypeaks=list(KDEdensitypeaks))
  densitypath<-c(densitypath,peaksdistance=list(dis))
  densitypath<-c(densitypath,minspantreepath=list(minspantreepath))
  
  return(densitypath)
}

####################################################################### 
# Examples of calling function
datadensitypath<-DensityPath(XX,k,h,SI_fig_dir)

