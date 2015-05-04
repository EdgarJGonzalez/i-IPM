s1 <- read.delim("evals_PCA1_PCA2_scaled_params.xls")
library("rgl")
s1m <- as.matrix(s1)
surface3d(1:1000,1:1000,s1m)
image(s1m,useRaster=TRUE,col=gray((0:50)/50))
image(s1m[seq(1,1000,by=2),seq(1,1000,by=2)],
      useRaster=TRUE,col=gray((0:50)/50))
contour(s1m,add=TRUE,levels=c(300,310,320,330),col=2,drawlabels=FALSE)
topvals <- head(sort(s1m),200)
toplocs <- which(s1m %in% topvals)
xvec <- (1:1000)/1000
lines(xvec,apply(s1m,1,which.min)/1000,col=4)
          
## figure out x,y coordinates corresponding to array index
xvals <- toplocs %% 1000
yvals <- toplocs %/% 1000
points(xvals/1000,yvals/1000,col=4)               

x <- read.delim("pars_s5_C1.xls")


surface3d(1:1000,1:1000,s1m)

