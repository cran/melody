    
.nodes2edges <- function(n)
 (((n^2)-n)/2)



## calculate the graph distances between parameters (that start with u, not followed by an x and don't end with a plural 's')
.calcEdgeDistances <- function(sg, u.params = names(sg)[grep('^u\\.[^x][^\\.]+[^s]$', names(sg))]){
#c('u.width','u.ymax','u.ymin','u.ymid','u.ylen','u.ymean','u.intensity','u.tone','u.interval')){
    
  d.params <- c(gsub('u\\.','d.', u.params), 'd.shape')  
  for(d.param in d.params)
    sg[[d.param]] <- matrix(NA, sg$n, sg$n);
  
  bg <- sg$bg.threshold

  for(i in 1:sg$n){    
    for(j in i:sg$n){
      
      if(sg$u.width[i] > sg$u.width[j]){
        #nrrwMtrx <-   sg$units[[j]];  wideMtrx   <- sg$units[[i]]
        nrrwMtrxTF <-   sg$units[[j]]<bg;  wideMtrxTF   <- sg$units[[i]]<bg
        nrrwWdth <-   sg$u.width[j];  wideWdth   <- sg$u.width[i];
        nrrwintensity <- sg$u.intensity[j]; wideintensity <- sg$u.intensity[i];
      }else{
        #nrrwMtrx   <- sg$units[[i]]; wideMtrx   <- sg$units[[j]]
        nrrwMtrxTF   <- sg$units[[i]]<bg; wideMtrxTF   <- sg$units[[j]]<bg
        nrrwWdth   <- sg$u.width[i]; wideWdth   <- sg$u.width[j];
        nrrwintensity <- sg$u.intensity[i];wideintensity <- sg$u.intensity[j];
      }  
      diff <- list()
      for(u.param in u.params){
         param <- gsub('^u\\.','', u.param)
         diff[[param]] <- abs(sg[[u.param]][j] - sg[[u.param]][i])
      }
      longLen <-  with(sg, max(u.ylen[j],u.ylen[i]))

      overlapSumVect <- vector("numeric",length=diff[['width']]+1); #diffVect
      for(k in 1:(diff[['width']]+1)){
        #diffVect[k] <- sum(abs(nrrwMtrx-wideMtrx[,k:(nrrwWdth+k-1)]))
        overlapSumVect[k] <- sum((nrrwMtrxTF & wideMtrxTF[,k:(nrrwWdth+k-1)]))
      }
    
      sg$hdenom.ymax <- quantile(sg$u.ymax, seq(0,1,.1))[9]*2 #mean(sg$u.ymax, na.rm=TRUE)

      sg$d.intensity[i,j] <- diff[['intensity']]/wideintensity
      sg$d.shape[i,j] <- 1-(max(overlapSumVect)/sum(nrrwMtrxTF))
      sg$d.width[i,j] <- diff[['width']] / wideWdth 
      sg$d.ymax[i,j] <- diff[['ymax']]   / sg$hdenom.ymax #sg$height
      sg$d.ymean[i,j] <- diff[['ymean']] / sg$hdenom.ymax #sg$height #?longLen
      sg$d.ymin[i,j] <- diff[['ymin']]   / sg$hdenom.ymax #sg$height
      sg$d.ylen[i,j] <- diff[['ylen']]   / longLen
      sg$d.ymid[i,j] <- diff[['ymid']]   / sg$hdenom.ymax #sg$height #?longLen
      sg$d.pctbg[i,j] <- diff[['pctbg']] 
      sg$d.ymRpuH[i,j] <- diff[['ymRpuH']]
      sg$d.lmpv[i,j] <- diff[['lmpv']]      
 
    } 
  }

  for(dparam in d.params){
     diag(sg[[dparam]]) <- NA
     if(all(is.na(sg[[dparam]])))
      warning(paste("All of the distances for",dparam,"were NA. Undefined distance metric for a new 'u' param?")) 
  }
  return(sg)

}


.edgequery <- function(sg, n1, n2, dls=.dist.limits, fail=TRUE){
  #diagnostic for edge failure (connection or disconnection)
  for(dl in names(dls)){
    if(paste('u',dl,sep='.') %in% names(sg)){
      n1v <- sg[[paste('u',dl,sep='.')]][n1]
      n2v <- sg[[paste('u',dl,sep='.')]][n2]
      details <- paste('(',round(n1v,2),'-',round(n2v,2),'=)')
    }else{
      details <- ''
    }
    ev <- sg[[paste('d',dl,sep='.')]][n1,n2]
    if(!fail | (ev > dls[dl]))
      print(paste(dl,':',details,round(ev,3),'>',dls[dl]))
  }
}


.dist.limits <- nv(
c(    .5,         .6,    .46,    .07,   .15,    .4,    .2,    .2,    .55,     .3),
c('shape','intensity','width','ymean','ymid','ylen','ymin','ymax','ymRpuH','pctbg')
)

.binaryEdge <- function(sg, method=c('limits','avg'), 
                        dist.limits, 
                        dist.weights=rep(1,length(.dist.limits)), dist.avg=.15){  

  method <- match.arg(method)

  d.params <- names(sg)[grep('^d\\.',names(sg))]
  names(dist.limits) <- paste('d',names(dist.limits), sep='.')

  dp.in.dln <- d.params %in% names(dist.limits) 
  if(! all(dp.in.dln))
    stop(paste('every distance comparison must have a limit defined: ', 
                paste(d.params[!dp.in.dln],collapse=',')))

  if(method =='limits'){
    binary <- list()
    for(d.param in d.params)
       binary[[d.param]] <- sg[[d.param]] < dist.limits[d.param]

    sg$edge.mtrx <- Reduce('&',binary[d.params])

  }else if (method=='avg'){
    
    #edge.mtrx.sum <- Reduce("+",sg[d.params])
    #edge.mtrx.mean <- edge.mtrx.sum / length(d.params) 
    require(abind)      
    edge.cube <- abind(sg[d.params], along=3)
    edge.mtrx.mean <- apply(edge.cube, 1:2, weighted.mean, dist.weights)

    sg$edge.mtrx <- edge.mtrx.mean < dist.avg

  }
  return(sg)

}


.plotEdgeParams <- function(VSG, rnss, dl =.dist.limits){
    ## this takes the list of lists of vocalization spectrograms [VSG] and helps you to
    ## tune the clustering paramters in 'dist.limits()'
    VSG.ss <- VSG[names(VSG) %in% rnss]
    #dparams <-  names(VSG.ss[[1]][[1]])[grep('^d\\.',names(VSG.ss[[1]][[1]]))]

    adp <- list()
    par(mfrow=c(2,1 +length(dl)/2),mar=c(2,1,1,0))

    for(d in names(dl)){
    adp[[d]] <- unlist(lapply(VSG.ss, function(x)
                          unlist(lapply(x, function(y){
                                           as.vector(y[[paste('d',d,sep='.')]])
                                         } ))))
    print(d)
    hist(adp[[d]], main=d)
    abline(v=dl[d],col='red')
    }

}
        
.calcNetStats <- function(sg){ 
    sg$e.ct <- sum(sg$edge.mtrx, na.rm=TRUE)
    sg$c.density <- sg$e.ct / sg$e

    sg$edge.mtrx[is.na(sg$edge.mtrx)] <- FALSE
    
    # leverage the network package
    require(network);require(sna)
    sg$graph <- network(symmetrize(sg$edge.mtrx), directed=F)
    unit.clust.membership <- component.dist(sg$graph)$membership
    unit.clust.degrees <- degree(sg$graph, cmode="indegree")
    cluster.sizes <- component.dist(sg$graph)$csize;     
    cluster.maxedgect <- .nodes2edges(cluster.sizes)
    sg$c.ct <- components(sg$graph) 
    #sg$c.ct = weighted.mean(sg$u.pctbg, sg$u.intensity, na.rm=TRUE)

    cluster.energies <- sapply(split(sg$u.intensity, unit.clust.membership), mean);

     ## are isolates 1 repeat or zero or NA?
    cluster.degree.avgs <- sapply(split(unit.clust.degrees, unit.clust.membership), mean);
    cluster.degree.pcts <- cluster.degree.avgs/(cluster.sizes-1)
    cluster.degree.pcts[is.na(cluster.degree.pcts)] <- 0

    sg$c.degavg <- mean(cluster.degree.avgs, na.rm=TRUE)
    sg$c.degpct <- mean(cluster.degree.pcts, na.rm=TRUE)
    sg$c.nreppct <- mean(unit.clust.degrees > 0, na.rm=TRUE)

    sg$c.degavg.ew <- weighted.mean(cluster.degree.avgs, cluster.energies, na.rm=TRUE)
    sg$c.degpct.ew <- weighted.mean(cluster.degree.pcts, cluster.energies, na.rm=TRUE)
    sg$repetition <- sg$c.nreppct.ew <- weighted.mean(unit.clust.degrees > 0, sg$u.intensity, na.rm=TRUE)
    

    #weighted cluster count
    sg$c.ct.ew <- sum(cluster.energies/max(cluster.energies))

    return(sg)
}

cluster <- function(sg, method=c('limits','avg'), dist.limits = nv(
c(    .5,         .6,    .46,    .07,   .15,    .4,    .2,    .2,    .55,     .3),
c('shape','intensity','width','ymean','ymid','ylen','ymin','ymax','ymRpuH','pctbg')
), 
dist.avg=.15, rep.score=c('c.nreppct','c.degpct','cedge.dens'), syl.score=c('c.ct'), intensity.weighted=TRUE){

  method <- match.arg(method)
  rep.score <- match.arg(rep.score)
  syl.score <- match.arg(syl.score)

  if(sg$n > 1){
 
    sg <- .calcEdgeDistances(sg) 
    sg <- .binaryEdge(sg, method=method, dist.limits=dist.limits, dist.avg=dist.avg)
    sg <- .calcNetStats(sg)    

  }else{
      sg$c.ct <- 1          # nClusters
      sg$c.degavg <- 0         # aveClustDegree
      sg$cdist.mtrx <- NA     # distMtrx
      sg$cdist.avg <- Inf     # aveDiff
      sg$cedge.dens <- 0
      sg$c.degpct <- 0
      sg$c.nreppct <- 0       # prcntCmpntsRptd
      sg$c.degpct.ew <- 0
      sg$c.nreppct.ew <- 0   # prcntintensityRptd
      sg$repetition <- 0
      sg$syllable_ct <- 1
      sg$c.ct <- 1
      sg$c.ct.ew <- 1
                              # aveClustSize
                              # aveReps.group
                              # aveDegree.group
                              # aveSaturation
                              # wAveSaturation
    if(sg$n == 0){
      sg$c.ct <- 0
    }
  }   

  sg$repetition  <- sg[[ paste(rep.score, c('','.ew')[intensity.weighted+1],sep='') ]]
  sg$syllable_ct <- sg[[ paste(syl.score, c('','.ew')[intensity.weighted+1],sep='') ]]

  return(sg)
}     
      


#.rti <- function(sg){
#(widthDiff > (2*nrrwWdth)) | (intensityDiff > (2.5*nrrwintensity)) | (yLenDiff > (.3*longLen)) | (yMidDiff > (.3* sg$height))){
#}

