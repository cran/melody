
plot.sg <- function(x, blanks=TRUE, u.ranges=TRUE, u.yzxs=TRUE, label.cex=.7, xlab="Time", ylab="Frequency", ...){

  image(x=0:x$width, y=0:x$height, z=t(x$x), xlab=xlab, ylab=ylab, axes=F, ...)
  blank.clmns <- rep(x$blank.rls$values, times=x$blank.rls$lengths)

  blank.clmn.idx = (1:ncol(x$x))[blank.clmns]  #for plotting!!

  if(blanks)
    points((1:x$width)[blank.clmn.idx]-.5, rep(1,length(blank.clmn.idx)), pch="|")
  
  if('n' %in% names(x)){
   if(x$n > 0){
    if(u.ranges){  
     points(x$u.xavg -.5, x$u.ymax, pch="-")
     points(x$u.xavg -.5, x$u.ymin, pch="-")
     points(x$u.xavg -.5, x$u.ymid, pch="|") #unitNmbrs
    }
  
    if('u.xavg' %in% names(x)){
     if(u.yzxs)
      for(i in 1:length(x$units))
       points( (x$u.xavg[i]-x$u.width[i]/2 +.5):(x$u.xavg[i]+x$u.width[i]/2-.5), x$u.ymean.l[[i]], pch=".") 

      #text(x$u.xavg - .5, 4 , col=rainbow(x$n, start=0,end=.8),cex=label.cex)
      mtext((1:x$n), side=1,  at=x$u.xavg - .5 , col=rainbow(x$n, start=0,end=.8),cex=label.cex)
     } 
    }
  }
}



graph <- function(sg, label.cex=.7, ...){
        #bmp(file=paste(plotPath,'/',vocFileName,"-matchGraph.bmp",sep=""))#height)
        ## at some point fold rainbow colors and labels into one "colabels" named vector
        plot.network(sg$graph,vertex.cex=sqrt(sg$u.intensity)/sqrt(mean(sg$u.intensity)),
                         vertex.col=rainbow(sg$n,start=0,end=.8),
                         label.col=rainbow(sg$n,start=0,end=.8),
                         displaylabels=T, boxed.labels=F, label.cex=label.cex, ...
         )
}


