tone <- function(sg, method=c('pctbg','iqsub.pctbg','yDpuH'), intensity.weighted=TRUE, w2h.min=0.1){  #'iqsub.pctbg',

  method <- match.arg(method)#, several.ok=TRUE)

  if(sg$n == 0){
    sg$tone <- sg$pctbg <- 0
    return(sg)
  }

  ## if the height is ~10 times bigger than the width.. penalize it as atonal 
  if(w2h.min < 1){   
     if('u.ylen' %in% names(sg)){
        w2h.p <- sg$u.w2hs <-  sg$u.width / sg$u.ylen
        w2h.p <- w2h.p * 1/w2h.min
        w2h.p[sg$u.w2hs > w2h.min] <- 1
     }else{
        warning('w2h penalty set but no unitYstats available. ignoring')
        w2h.p <- 1
     }
  }
  sg$u.pctbg <- sapply(1:sg$n, function(i) 
                     sum(sg$units[[i]] > sg$bg.threshold) 
                   / prod(dim(sg$units[[i]]))
                       )
  
  sg$pctbg <- mean(sg$u.pctbg * w2h.p, na.rm=TRUE)

  if('u.ymean.l' %in% names(sg)){
   sg$u.ymdeltas.l <- sapply(sg$u.ymean.l, function(x) abs(diff(x)) ) #absdiff?
   sg$u.ymdltavgs <- sapply(sg$u.ymdeltas.l, mean, na.rm=TRUE)
  ##y delta per Range: measures randomness
   sg$yDpuH <- sg$u.ymdltavgs / sg$u.ylen 
  }

  #average the number of (inner-quartile x sub y.max) pixels above the bg.threshold & divide by y.max 
  u.iqsub.pctbg <- sapply(1:sg$n,function(i){
                            xrange <- do.call(":",as.list(quantile(1:sg$u.width[i])[c(2,4)]));
                            sum(sg$units[[i]][, xrange] > sg$bg.threshold) /(sg$height*length(xrange))
                            }
                         )
  sg$iqsub.pctbg <- mean(u.iqsub.pctbg * w2h.p, na.rm=TRUE)

  if('u.intensity' %in% names(sg)){
    sg$pctbg.ew = weighted.mean(sg$u.pctbg * w2h.p, sg$u.intensity, na.rm=TRUE)
    sg$iqsub.pctbg.ew = weighted.mean(sg$iqsub.pctbg * w2h.p, sg$u.intensity, na.rm=TRUE)

  }else{
    if(intensity.weighted){
      warning('intensity weighting is not possible until after unitstats is run to obtain unit energies.\
                 setting "tone" to the unweighted mean')
      intensity.weighted <- FALSE
    }
  }

  sg$tone <- sg[[paste(method, c('','.ew')[intensity.weighted+1],sep='')]]
             #do.call('prod', lapply(method, function(mt)  <snip> ))

  return(sg)
}

interval <- function(sg, method=c('ymRpuH'), intensity.weighted=TRUE){

  reqd.unit.stats <- c('u.intensity','u.ymean.l','u.ylen')
  if(!all(reqd.unit.stats %in% names(sg)))
      warning('calculating interval is not possible until after unitstats is run')
  
  method <- match.arg(method)

  if(sg$n == 0){
    sg$interval <- 0
    return(sg)
  }

  ## check for inervalic (sloping) content 
  #VyRangePerHeights
  sg$u.ymean.rngs <- sapply(sg$u.ymean.l, function(x) abs(diff(range(x, na.rm=TRUE)))) #yMeanRanges

  sg$u.ymRpuH <- sg$u.ymean.rngs / sg$u.ylen ; 
  sg$u.ymRpuH[is.na(sg$u.ymRpuH)] <- 0    #measures inerval
    
  sg$ymRpuH = mean(sg$u.ymRpuH, na.rm=TRUE)
  sg$ymRpuH.ew = weighted.mean(sg$u.ymRpuH, sg$u.intensity, na.rm=TRUE)

  sg$interval <- sg[[paste(method,c('','.ew')[intensity.weighted+1],sep='')]]

  return(sg)
}

#intensity.weighted.means <- function(sg)

.filter <- function(sg){
  # not used
  tone <- sg$TpLpH  < .5                 
  vert <- sg$ymRpuH  > .35                
  grad <- sg$yDpT  < 1.75                   
  grad <- sg$yDpRpW  < 1                   

  sg$keep <-  tone & vert & grad
  return(sg)
}


