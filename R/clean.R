#remove background noise

#affineTransform(){
#}
#collapseHarmonics(){
#}

threshold <- function(sg, type=c('bw','div'), pct.max=0.89){
  sg$gray.dev <- sqrt(var(as.numeric(sg$x)))
  sg$gray.mean <- mean(sg$x)
  sg$gray.max <- max(sg$x)
  sg$bg.threshold <- sg$gray.max * pct.max # 227/255  #BlkWhtThreshold
  return(sg)
}


partition <- function(sg, ms=1, ds=1.62, lambda=5){

  sg$blank.rls <- .findBlanks(sg$x, ms=ms, ds=ds, bg.th=sg$bg.threshold, lambda=lambda)
  
  sg$units <- .splitOnBlanks(sg$x, blank.rls=sg$blank.rls)

  sg$n <- length(sg$units)
  sg$e <- .nodes2edges(sg$n)  
  return(sg)
}

.findBlanks <- function(x, ms=1, ds=1.62, bg.th=227, lambda=5, bgrp=.2, bgcp=.5){   

  if(ds == 0 & ms == 0){
    db <- mb <- 0 # this isn't totally necessary but faster
  }else{
    row.dev.v <- apply(x,1,function(r) sqrt(var(r)))
    db <- row.dev.v * ds

    ## create an exponentially decaying (low) background removal curve
    efnc <- dexp((1:nrow(x))/(nrow(x)), rate=lambda)

    ## create a multiplier to assess how noisy the bottom of the spec is
    ## dbp stands for dark bottom percent
    bgrc <- round(bgrp*nrow(x))
    dbp <- sum(x[1:bgrc, ] < bg.th) / (ncol(x) * bgrc)
    dbp <- dbp/.40 #anything greater than 40 percent dark bottom
    # should get full effect of exponential mean threshold blank finding
    dbp[dbp > 1] <- 1

    row.avg.v <- apply(x,1,mean)
    row.avg.thresh.diff.v <- bg.th -row.avg.v
    mb <- row.avg.thresh.diff.v * (efnc/max(efnc)) * ms * dbp 
  }
  
  blanks <- apply(x,2,function(c) 
    !any(c < (bg.th - (db + mb ))))

  ## run length encode the blanks
  return(rle(blanks))

}

.splitOnBlanks <- function(X, blank.rls){
  blank.rls.cs <- cumsum(blank.rls$lengths)
  blank.rls.cs[blank.rls$values] <- NA
  ## now do the actual partitioning, dropping the blank (NA) clmns 
  ## (t matrix rotation enables horizontally acting 'split()')
  units <- split(t(X), f=factor(rep(blank.rls.cs, times=blank.rls$length)))

  ## need a check in here or above in '.findBlanks' to check for empty-ish units (nothing below 228

  #reformat each unit into a matrix and calculate the height of each component 
  lapply(units,function(x) matrix(x,nrow=nrow(X),byrow=T))  # byrow flips back upright

}


unitYstats <- function(sg, harmonics=FALSE){

  if(harmonics){  #if we just want to use a single harmonic (the fundamental) for calculating y stats
    u.harmonics <- lapply(sg$units, function(u) 
         .splitOnBlanks(t(u), blank.rls=.findBlanks(t(u), ms=0, ds=0, bg.th=sg$bg.threshold)))
    ## find fundamental based on max intensity
    u.harm.intens.max <- lapply(u.harmonics, function(uhs) max(sapply(uhs, function(h) sum(sg$gray.max-h)))) 
    u.harm.intens.med <- lapply(u.harmonics, function(uhs) median(sapply(uhs, function(h) sum(sg$gray.max-h)))) 
    u.harm.widths <- lapply(u.harmonics, function(uhs) lapply(uhs, function(h) 
                                                           sum(!is.na((1:nrow(h))[apply(h,1,min) < sg$bg.threshold])))) 
    u.harm.intensities <- lapply(u.harmonics, function(uhs) sapply(uhs, function(h) sum(sg$gray.max-h)))  
    u.harm.intense.pct <- lapply(1:length(u.harmonics), function(i)  u.harm.intensities[[i]] / u.harm.intens.max[[i]])
    ## filter out extremely weak narrow or low harmonics
    u.harmonics <- lapply(1:length(u.harmonics), function(i)  u.harmonics[[i]][(u.harm.intense.pct[[i]] > .01) 
                                                                             & (u.harm.widths[[i]]      > (.85 * sg$u.width[i]))
                                                                             & as.numeric(names(u.harmonics[[i]])) > 5
                                                                              ])
    ## make and index of the units that appear to be truely harmonic
    u.harm.idx <- (1:length(u.harmonics))[sapply(u.harmonics, length) > 0]
 
    ## set up default units and ymax lists (null non-harmonic hypothesis)
    u.fund.ymax <- sapply(sg$units, nrow) 
    units <- sg$units                     
    ## replace the harmonic units with just their first fundamental (reject the null)
    u.fund.ymax[u.harm.idx] <- sapply(u.harm.idx, function(i) as.numeric(names(u.harmonics[[i]])[1]))  
    units[u.harm.idx] <- lapply(u.harm.idx, function(i) t(u.harmonics[[i]][[1]]))  
    
    u.ybtms <- unlist(u.fund.ymax) - sapply(units, nrow)
  }

  if(!harmonics){
    units <- sg$units
    u.ybtms <- rep(0, length(units))
  }

  sg$u.ys <- lapply(units,function(x) (1:nrow(x))[apply(x,1,min) < sg$bg.threshold])    
  sg$u.ys <- lapply(1:length(units), function(i) sg$u.ys[[i]] + u.ybtms[i]) # add bottoms back on if necessary
  sg$u.ys[sapply(sg$u.ys, length)==0] <- NA
  sg$u.ymin <- sapply(sg$u.ys, min) #+ u.ybtms
  sg$u.ymax <- sapply(sg$u.ys, max) #+ u.ybtms     

  sg$u.ymid <- (sg$u.ymin + sg$u.ymax)/2  #+ u.ybtms
  sg$u.ylen <- sg$u.ymax - sg$u.ymin 

  sg$u.ymean.l <- lapply(units,function(x) 
                       apply(x,2,function(z) {
                        dark <- z<sg$bg.threshold;
                        if(any(dark))
                         weighted.mean((1:nrow(x))[dark],(sg$gray.max - z[dark]))
                        else
                          NA
                        }
                      )) 
  sg$u.ymean.l <- lapply(1:length(units), function(i) sg$u.ymean.l[[i]] + u.ybtms[i])  
  sg$u.ymean <- sapply(sg$u.ymean.l, mean, na.rm=TRUE) #+ u.ybtms

  return(sg)

}

unitstats <- function(sg, harmonics=FALSE){


  sg$u.width <- sapply(sg$units, ncol)  # oldWidths  # for plotting!
  sg$u.xavg <- as.numeric(names(sg$units)) - sg$u.width/2 
 
  sg <- clean(sg)

  sg <- unitYstats(sg, harmonics=harmonics)

  sg$u.intensity <- sapply(sg$units,function(x) sum(-(x- sg$gray.max)))
  sg$intensity <- sum(sg$u.intensity)

  return(sg)
}



clean <- function(sg, wlim=c(3,sg$width), ylim=c(0,1), ylen=3, imin=0.01){ #, ilim.sd=3){ #, wpmmin=.05){

#  widx <- sg$u.width >= wpmmin * max(sg$u.width)
  wide <- sg$u.width >= wlim[1]                   
  naro <- sg$u.width <= wlim[2]              
 
  dark <- sapply(sg$units, function(x) any(x < sg$bg.threshold))

  unit.stats.avail <- ('u.ylen' %in% names(sg))

  if(unit.stats.avail){
    tall <- sg$u.ylen >= ylen               # tallEnough                    
    high <- sg$u.ymin >= ylim[1] *sg$height # highEnough         
    low  <- sg$u.ymax <= ylim[2] *sg$height # lowEnough #0.999       
    loud <- sg$u.intensity > imin * sg$intensity  # strongEnough
  }else{
    tall <- high <- low <- loud <- rep(TRUE, sg$n)
  }

  keep <- k <- wide & naro & tall & high & low & loud & dark

  # attempt to exclude large outliers that likely didn't get split up well by partition step
  #if(unit.stats.avail & (sum(keep) >= 5)){
  #  ## because 'quiet' is relativistic, we do that last after the initial filter
  #  quiet <- sg$u.intensity[k] < (median(sg$u.intensity[k]) + (ilim.sd * sd(sg$u.intensity[k])))
  #  keep[k] <- keep[k] & quiet
  #}
  #sg$keep <- keep

  ## update any slots that start with a "u" 
  for(u in grep('^u',names(sg), value=TRUE))
    sg[[u]] <- sg[[u]][keep]
  
  if('u.intensity' %in% names(sg))
    sg$intensity <- sum(sg$u.intensity)


  sg$n <- length(sg$units)
  sg$e <- .nodes2edges(sg$n)  

  if(length(sg$units) ==0)
    warning('over-cleaned: no units left!')    

  return(sg)
}



