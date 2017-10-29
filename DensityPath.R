

rm(list=ls())
library("TDA")
library("gdistance")
library("raster")
library("vegan")
library("shape")
library("destiny")
library('R.matlab')


path<-('～/densitypath')
pathname<-file.path(path,'testsamples6.mat')
treedata<-readMat(pathname)
X1<-as.matrix(treedata$testsamples6)
XX<-X1[,1:2]

k<-50
h<-0.3

DensityPath<-function(XX,k,h){
  densitypath<-list()
  
  if(ncol(XX)!=2){
    XX<-DiffusionMap(XX,k=2)
  }

  TreeKDE <- clusterTree(XX, k = k, h = h, density = "kde",
                         printProgress = FALSE)
  densityKDE<-TreeKDE$density
 

  idKDE<-setdiff(TreeKDE$id,TreeKDE$parent)
  numKDEleaves<-length(idKDE)
  l<-1
  KDEdensitypeaks<-matrix(1,numKDEleaves,2)
  for (i in idKDE){
    clusterkde<-TreeKDE$DataPoints[[i]]
    KDEdensitypeaks[l,]<-XX[clusterkde[which.max(densityKDE[clusterkde])],]
    l<-l+1
  }
  
  
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
  
  ###################################################
  KDE <- kde(X = XX, Grid = Grid, h = h)
  

  r <- raster(nrows=numy, ncols=numx, xmn=xmin, xmx=xmax, ymn=ymin, ymx=ymax,crs="+proj=utm +units=m")
  r[] <- KDE
  T <- transition(r, function(x) mean(x), 8)
  # 1/mean: reciprocal to get permeability
  T <- geoCorrection(T)
  C <-KDEdensitypeaks
  D <-KDEdensitypeaks
  dis<- costDistance(T, C, D)
  spanningtree <- spantree(dis)
  
  pdf("output.pdf")
  par(mfrow = c(2,2))
  par(mai = c(0.7,0.7,0.3,0.4))
  #2维投影
  plot(XX, pch = 19, cex = 0.6, main = "2D mapping of single cell points",xlab="",ylab="")
  
  #3维密度曲面
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
  
  #离散的density clusters
  plot(1:5,1:5,xlim=c(xmin,xmax),ylim=c(ymin,ymax),type = "n",main = "Density clusters",xlab="",ylab="")
  c<-1
  for (i in idKDE){
    points(matrix(XX[TreeKDE[["DataPoints"]][[i]], ], ncol = 2), col = c,pch = 19, cex = 0.6)
    c<-c+1
  }
  
  #density path
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
  
  densitypath<-c(densitypath,densityKDE=list(densityKDE))
  densitypath<-c(densitypath,KDEdensitypeaks=list(KDEdensitypeaks))
  densitypath<-c(densitypath,peaksdistance=list(dis))
  densitypath<-c(densitypath,minspantreepath=list(minspantreepath))

  return(densitypath)
}


datadesitypath<-DensityPath(XX,k,h)







