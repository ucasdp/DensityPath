
rm(list=ls())
# Install and library the following R packages
install.packages("TDA")
install.packages("gdistance")
install.packages("raster")
install.packages("vegan")
install.packages("shape")
source("https://bioconductor.org/biocLite.R")
biocLite('destiny')

library("TDA")
library("gdistance")
library("raster")
library("vegan")
library("shape")
library("destiny")


# Read the data files
X1<-read.table('～/densitypath/testdata.csv',sep=',',header=FALSE)
XX<-as.matrix(X1[,1:2])

# Set the parameters of density clustering
numcells<-nrow(XX)
if (numcells<1000){
  k<-10
}
if (numcells>=1000 & numcells<10000){
  k<-50
}
if (numcells>10000){
  k<-100
}

h<-mean(dist(XX))/10

# Set the path of graphic output
SI_fig_dir <- "～/densitypath"

#Select 62 points in XX for mapping
X<-matrix(0,62,2)
for(i in 1:62){
  X[i,]=XX[40*(i-1)+1,]
}

#################################################################################
# the function of DensityPath 
DensityPath<-function(XX,k,h,SI_fig_dir,mappingpoints){
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
  
  maxfKDE<-rep(1,numKDEleaves)
  minfKDE<-rep(1,numKDEleaves)
  maxfKDE2<-rep(1,numKDEleaves)
  minfKDE2<-rep(1,numKDEleaves)
  
  l<-1
  for (i in idKDE) {
    maxfKDE[l]<-TreeKDE[["lambdaTop"]][[i]]
    minfKDE[l]<-TreeKDE[["lambdaBottom"]][[i]]
    maxfKDE2[l]<-TreeKDE[["kappaTop"]][[i]]
    minfKDE2[l]<-TreeKDE[["kappaBottom"]][[i]]
    l<-l+1
  }
  
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
  

  
  ######################################################################
  #step5:Map the single cells onto the cell state-transition path.
  
  if (missing(mappingpoints)==FALSE){
    numadjpath=length(minspantreepath)
    eachadjpath<-rep(1,(numadjpath+1))
    for(i in 1:numadjpath){
      eachadjpath[i+1]<-nrow(minspantreepath[[i]])
    }
    n<-sum(eachadjpath)-1
    globaladjpath<-matrix(0,n,3)
    for(i in 1:numadjpath){
      n1<-sum(eachadjpath[1:i])
      n2<-sum(eachadjpath[1:(i+1)])-1
      globaladjpath[n1:n2,]<-minspantreepath[[i]]
    }
    globalpoints<-globaladjpath[,1:2]
    
    #draw the density path and mapping path
    pdf("MappingPath.pdf")
    par(mfrow = c(1,2))
    par(mai = c(0.7,0.4,0.4,0.4),oma=c(9,0.7,9,0.7))
    plot(r,axes=FALSE,main="Mapping Path",xlab="",ylab="")
    #plot(r,xlim=c(xmin,xmax),ylim=c(ymin,ymax),main="Mapping Path",xlab="",ylab="")
    for (i in 1:numKDEedge){
      p1<-KDEspantree[i,1]
      p2<-KDEspantree[i,2]
      p[1,]<-KDEdensitypeaks[p1,]
      p[2,]<-KDEdensitypeaks[p2,]
      p1top2 <- shortestPath(T, p[1,], p[2,], output="SpatialLines")
      lines(p1top2, col="red", lwd=2)
    }
    
    for (i in 1:numKDEleaves){
      points(KDEdensitypeaks[i,1],KDEdensitypeaks[i,2],pch = 19, cex = 1,col="blue")
    }
    
    minadjpaths<-list()
    for(i in (1:(nrow(mappingpoints)))){
      mappingpoint<-mappingpoints[i,]
      #Look for the nearest point on the density Path that is closest to mappingpoint
      dis<- costDistance(T, mappingpoint, globalpoints)
      correpoint<-globalpoints[min((which(dis==min(dis)))),]
      p3top4 <- shortestPath(T, mappingpoint,correpoint, output="SpatialLines")
      
      
      lines(p3top4, col="cyan", lwd=2)
      points(mappingpoint[1],mappingpoint[2],pch = 19, cex = 1,col="cyan")
      points(correpoint[1],correpoint[2],pch = 19, cex = 1,col="cyan")
      minadjpath<-as.matrix(((((p3top4@lines)[[1]])@Lines)[[1]])@coords)
      minadjpath2<-matrix(0,nrow(minadjpath)+2,2)
      minadjpath2[1,]<-mappingpoint
      minadjpath2[2:(nrow(minadjpath)+1),]<-minadjpath
      minadjpath2[(nrow(minadjpath)+2),]<-correpoint
      minadjpaths<-c(minadjpaths,list(minadjpath2))
      }
    
    
    Xseq <- seq(Xlim[1], Xlim[2], by = by)
    Yseq <- seq(Ylim[1],Ylim[2], by = by)
    Grid <- expand.grid(Xseq, Yseq)
    zh<- matrix(KDE,ncol=length(Yseq),nrow=length(Xseq))
    op <- par(bg = "white")
    res<-persp(Xseq, Yseq, zh, theta = 165, phi = 60,
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
    
    for(i in 1:numadjpath){
      pathpoints<-minspantreepath[[i]]
      lines(trans3d(pathpoints[,1], pathpoints[,2], pathpoints[,3], pmat = res), col = "red", lwd=3)
    }
    
    for (i in 1:numKDEleaves){
      points(trans3d(KDEdensitypeaks[i,1],KDEdensitypeaks[i,2],maxfKDE[i], pmat = res),pch = 19, cex = 1,col="blue")
    }
    
    
    for(i in 1:length(minadjpaths)){
      minadjpath<-minadjpaths[[i]]
      nr<-nrow(minadjpath)
      a1<-min(minadjpath[,1])-by
      a2<-max(minadjpath[,1])+by
      b1<-min(minadjpath[,2])-by
      b2<-max(minadjpath[,2])+by
      Xseq1 <- seq(a1, a2, by = by)
      Yseq1 <- seq(b2, b1, by = -by)
      Grid1 <- expand.grid(Xseq1, Yseq1)
      KDE1 <- kde(X = XX, Grid = Grid1, h = h)
      minmappingpath<-matrix(0,nrow(minadjpath),3)
      minmappingpath[,1:2]<-minadjpath
      for(m in 1:nr){
        l<-which(abs(minadjpath[m,1]-Grid1[,1])+abs(minadjpath[m,2]-Grid1[,2])==min(abs(minadjpath[m,1]-Grid1[,1])+abs(minadjpath[m,2]-Grid1[,2])))
        minmappingpath[m,3]<-KDE1[l]
      }
      
      points(trans3d(mappingpoints[i,1],mappingpoints[i,2],minmappingpath[1,3],pmat = res),pch = 19, cex = 1,col="cyan")
      points(trans3d(correpoint[1],correpoint[2],minmappingpath[nrow(minmappingpath),3],pmat = res),pch = 19, cex = 1,col="cyan")
      lines(trans3d(minmappingpath[,1], minmappingpath[,2], minmappingpath[,3], pmat = res), col = "cyan", lwd=3)
      }
    dev.off()
    densitypath<-c(densitypath,minmappingpath=list(minadjpaths))
    }
  
  #######################################################################                
  # Return value                
  # return the density path, which contains the estimated density of the sample points(densityKDE),
  #               the two-dimensional coordinates of the density peaks(KDEdensitypeaks), 
  #               the geodesic distance between the density peaks(dis), 
  #               and the path of the minimum spanning tree on the three-dimensional density surface(minspantreepath).
  #               and the the mapping paths for the points(minadjpaths).
  densitypath<-c(densitypath,densityKDE=list(densityKDE))
  densitypath<-c(densitypath,KDEdensitypeaks=list(KDEdensitypeaks))
  densitypath<-c(densitypath,peaksdistance=list(dis))
  densitypath<-c(densitypath,minspantreepath=list(minspantreepath))

  return(densitypath)
}



###################################################################################### 
# Examples of calling function
# the Example without mapping points
datadensitypath<-DensityPath(XX,k,h,SI_fig_dir)
# the Example with mapping points
datadensitypath<-DensityPath(XX,k,h,SI_fig_dir,X)

