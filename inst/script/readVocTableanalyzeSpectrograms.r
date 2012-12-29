#path <- "H:\\My Docs\\research\\music\\music evolution\\spectroAnalysis\\code"

verbose <- F
requireInterval <- T  #require a filter for units which are not sloped enough
useDiffThresholding <- T

baseDIR <-  "~/primateVocs/"

source("utilFunctions.r")
library(network) #for plotting the blob comparison diffmatrix clustering 
library(sna)     #for the "symetrize function
library(zoo)     #for finding local minima to split the clustered blobs
library(scatterplot3d)  #for 3D plotting... but no text labeling

vocs <- read.table("priVocTable.txt",sep="\t",quote="",as.is=T,header=T)
vocs <- vocs[!(is.na(vocs$fig)|is.na(vocs$number)),] #remove the voc types that have no spectros

figRagArray <- sapply(vocs$figure,function(x) strsplit(as.character(x),",")) #split up the comma delimited list of figures for each voc type
spectrosPerVocVect <- unlist(lapply(figRagArray,length)) # count the number of spectrograms per vocalization type

totalNvocs <- dim(vocs)[1]
   
gifPath <- paste(baseDIR,"vocs",SLASH,sep="")
rSzdPath <- paste(baseDIR,"rSzd",SLASH,sep="")
pgmPath <- paste(baseDIR,"pgm",SLASH,sep="")
plotPath <- paste(baseDIR,"plots",SLASH,sep="")

###########################################################
###########################################################

#########          START MAIN LOOP         ################

###########################################################
###########################################################
#main program
VrepMatches <- VnClusters <- VnRepGroups <- VprcntCmpntsRptd <- VprcntEnergyRptd <- VpctTotVocEnergryCnsdrd <- VaveClustSize <- VaveClustDegree <- VaveReps.group <- VaveDegree.group  <- VaveSaturation <- VwAveSaturation  <- VaveIntervalVar <- VwAveIntervalVar <- VaveMaxPrecision <- spectrosPerVoc <- VunitDiffThreshold <- VaveDiff <- VpctTotVocEnergryRptd <- VyRangePerHeights <- as.vector(rep(NA,totalNvocs))



#loop through each vocalization type
#nmbrVect <-c(671,665,637,636,634,628,621,513,512,474,464,460,278,273,250,251,244,239,238,237,235,230,202,85,84,43,42,39,38,35)
for(vocNmbr in (1:totalNvocs)){  #nmbrVect ){   #
  #vocNmbr <- 238
  if(verbose)cat(paste("\nVoc",vocNmbr,":"))
   auth <- year <- vocs$ref[vocNmbr]
   auth <- gsub("^\([a-zA-Z]*\),[^~]*","\\1",auth)
   year <- gsub(".*(\([0-9]{4}\)).*","\\1",year)
   genus <- vocs$genus[vocNmbr]
   species <- vocs$species[vocNmbr]
   vocNumber <- vocs$number[vocNmbr]
   vocName <- vocs$name[vocNmbr]
   
   figures <- figRagArray[[vocNmbr]]  #vocs$figure[vocNmbr]
   duration <- vocs$duration[vocNmbr]
   frequency <- vocs$frequency[vocNmbr]
   newHeight <- round(log(1+frequency)*100) # round(frequency*20)  #this is just used for imageMagick's "convert -resize function"
   #if(!is.na(vocs$height) | vocs$height > newHeight) vocs$height <- newHeight  #not yet available
   #newWidth <- #round(duration*20) 
   #if(!is.na(vocs$width) | vocs$width > newWidth) vocs$width <- newWidth  #not yet available
   
   
   vocFileName <- paste(genus,species,vocNumber,vocName,paste(auth,year,sep=""),sep="-")
   #read in the spectrogram(s) for each vocalization type (joining spectros if there is more than one)
   PGMmtrx <- NULL #matrix()
  spectrosPerVoc[vocNmbr] <- 0
   for(figure in figures){
     if(spectrosPerVoc[vocNmbr] == 0){  #this is for a red seperator line inbetween spectro instances
       red <- 255
     }else{
       red <- 0
     }
     spectroFileName <- paste(vocFileName,figure,sep="-"); if(verbose)cat(paste(spectroFileName, ":"))
     #system(paste("ls ",pgmPath,spectroFileName,".gif",sep="")) #this is just a test line to make sure the names are kosher
     #system(paste("convert ",gifPath,spectroFileName,".gif"," -resize x",newHeight," ",rSzdPath,spectroFileName,".gif",sep=""))
     #system(paste("convert -compress none ",rSzdPath,spectroFileName,".gif"," ",pgmPath,spectroFileName,".pgm",sep=""))
     PGMmtrx <- cbind(PGMmtrx,red,255,readPGM(paste(pgmPath,spectroFileName,".pgm",sep="")),255)
     spectrosPerVoc[vocNmbr] <- spectrosPerVoc[vocNmbr] + 1
   }
 #} #this bracket used when converting all of the gifs (voc directory) into pgms with log height scaling
      ##################################################
      ##### direct PGM split window frame analysis #####
      ##################################################
    #  par(mfrow=c(2,2))
      meanGray <- mean(PGMmtrx)
      devGray <- sqrt(var(as.numeric(PGMmtrx)))
      maxGray <- max(PGMmtrx)
      BlkWhtThreshold <- 0.89*maxGray  #227/255
      devThreshold <- .4*maxGray   #100/255
      #now reformat the matrix so that it doesn't break the splitting algorythm below
      #PGMmtrx <- cbind(PGMmtrx,255)#do this by adding a light col @ begin & end of matrix (begining added above)#not needed anymore

      height <- dim(PGMmtrx)[1]
      width <- dim(PGMmtrx)[2]
      #print(voc); #print(table(apply(PGMmtrx,2,sum)))#REMOVE LATER
      rowDevianceVect <- apply(PGMmtrx,1,function(x) sqrt(var(x)))
      #devMax <- .5*maxGray  #100 #max(rowDevianceVect)
      #print(paste("highest deviance =",max(rowDevianceVect)," width:",width))
      rowAverageVect <- apply(PGMmtrx,1,mean)
      rowAveThresDiffVect <- BlkWhtThreshold-rowAverageVect
      #rowAveThresCrctdDiffVect <- (rowAveThresDiffVect-min(rowAveThresDiffVect))
    #  blankClmnVect <- apply(PGMmtrx,2,function(x) sum(x<BlkWhtThreshold)<2)#(meanGray-((.1*maxGray)/2)))<2)#!any(n<200)) #use histogram here too? upper bound? can't cross < 30?
      #frctnMaxRowDeviance <- -((devMax-rowDevianceVect)/devMax-1)
   
     ms <-  3 #  2.75
     vs <- .6 # 0.65
      blankClmnVect <- apply(PGMmtrx,2,function(x) !any(x<(BlkWhtThreshold-(rowAveThresDiffVect/ms + rowDevianceVect/vs))))
      #blankClmnVect <- apply(PGMmtrx,2,function(x) !any(x<(BlkWhtThreshold-(rowAveThresDiffVect)(frctnMaxRowDeviance))))#devThreshold-fraction
      #blankClmnVect <- apply(PGMmtrx,2,function(x) !any(x<(BlkWhtThreshold-(devThreshold*frctnMaxRowDeviance))))#devThreshold-fraction...
      blankClmnVect[length(blankClmnVect)] <- T #make sure that the last col is always blank (should i do this for first too?)
     # hist(PGMmtrx,breaks=30)
      #the column is blank if its not the case that any of the pixels are stronger  than x (200??)
      
    nBlankCols <- sum(as.numeric(blankClmnVect))      

           

      blankColsPos <- (1:width)[blankClmnVect]
      if(!(max(blankColsPos)==width)) blankColsPos <- c(blankColsPos,width) #in case the last column isn't blank (need to do for 1st too?)
    
      
      repTimes <- c(blankColsPos,0) - c(0,blankColsPos) - 1  #subtract the positions from each other to get the distances between positions
      repTimes <- c(repTimes[-c(1,length(repTimes))],0) #remove first and last and add a 0 to the end (last col is always blank)
      unitGroupings <- rep(blankColsPos,times=repTimes)#seq(1:nBlankCols)
      unitNmbrs <- oldUnitPositions <- unique(unitGroupings)
      PGMmtrxNoBlank <- PGMmtrx[,!(blankClmnVect)]
       #next is a sanity check. What caused it has been fixed by adding two white cols @ begin & end of matrix (above)
      if(!dim(PGMmtrxNoBlank)[2]==sum(repTimes))print("non-blank Matrix width doesn't match regrouping vector length")

#add recursive cutting based on vertical midpoint differences

      PGMmtrxNoBlankGrpList <- list()
      initialSplitNmbr <- length(unitNmbrs)
   if(verbose)cat(paste(" spliting into",initialSplitNmbr,"units -"))
      if(initialSplitNmbr>1){
       PGMmtrxNoBlankGrpList <- split(t(PGMmtrxNoBlank),unitGroupings)# I t'd because split is designed for splitting DFs by rows
      }else{
       PGMmtrxNoBlankGrpList[[1]] <- PGMmtrxNoBlank  #in the case that there is only one isolated repeat
       VrepMatches[vocNmbr] <- repMatches <- 1
      }
   if(verbose)cat(" calculating unit params -")

      whyExclude <- ""
  
      widths <- oldWidths <- sapply(PGMmtrxNoBlankGrpList,length)/height
      #now make sure these measurements are okay
      wideEnough <- widths >= 3                    ; whyExclude <- paste(whyExclude,c("n\n","")[wideEnough+1],sep="")
      nrrwEnough <- widths <= 1*width             ; whyExclude <- paste(whyExclude,c("w\n","")[nrrwEnough+1],sep="")   #(mean(widths)*5) 

      #reformat each unit into a matrix and calculate the height of each component 
      PGMmtrxNoBlankGrpMtrxList <- lapply(PGMmtrxNoBlankGrpList,function(x) matrix(x,nrow=height,byrow=T))

      yMins <- oldYmins <- sapply(PGMmtrxNoBlankGrpMtrxList,function(x) min((1:height)[apply(x,1,min)<BlkWhtThreshold]))      
      yMaxs <- oldYmaxs <- sapply(PGMmtrxNoBlankGrpMtrxList,function(x) max((1:height)[apply(x,1,min)<BlkWhtThreshold]))
  VaveMaxPrecision[vocNmbr] <- exp(-sqrt(var(yMaxs)))  #used to translate the 0 to inf var values to the 1 to 0 precision values
      yMids <- (yMins+yMaxs)/2
      yLens <- yMaxs-yMins
      #now make sure these measurements are okay
       tallEnough <- yLens >= 3                     ; whyExclude <- paste(whyExclude,c("s\n","")[tallEnough+1],sep="")
      highEnough <- yMins >= (0*height)            ; whyExclude <- paste(whyExclude,c("l\n","")[highEnough+1],sep="")  #not used
      lowEnough <- yMaxs <= (.999*height)             ; whyExclude <- paste(whyExclude,c("h\n","")[lowEnough+1],sep="") #not used

      totEnergy <- sum(-(PGMmtrxNoBlank-maxGray))
      energies <- sapply(PGMmtrxNoBlankGrpList,function(x) sum(-(x-maxGray)))
      strongEnough <- energies > (.01*totEnergy)  ; whyExclude <-  paste(whyExclude,c("w\n","")[strongEnough+1],sep="") #not used

    if(requireInterval){
      if(verbose)cat(" interval detection -") 

     # wghtdYcolVectLists <- sapply(PGMmtrxNoBlankGrpMtrxList,function(x) apply(x,2,function(z) rep((1:height)[z<BlkWhtThreshold],-(z[z<BlkWhtThreshold]-maxGray))))                   
     # defunctyMeansVectors <- sapply(wghtdYcolVectLists,function(x) sapply(x,mean))

      yMeansVectors <- sapply(PGMmtrxNoBlankGrpMtrxList,function(x) apply(x,2,function(z) weighted.mean((1:height)[z<BlkWhtThreshold],-(z[z<BlkWhtThreshold]-maxGray))))
      #yMedianVectors <- sapply(wghtdYcolVectLists,function(x) sapply(x,function(z) median))
      splitHarmonicGroups <- function(z){
         diffVector <- c((c(z,maxGray)- c(0,z))[-c(1,length(z)+1)],2)
         breakPoints <- diffVector>1
         #breakPoints <- breakPoints | c(breakPoints[-length(breakPoints)])
         nHarmonics <- sum(breakPoints)
         harmonics <- list()
         if(nHarmonics>1){
           breakPointFactor <- cut(z,breaks=z[breakPoints],include.lowest=T)
           harmonics <- split(z,breakPointFactor)
         }else{
           harmonics[[1]] <- z
         }
         harmonics
       }
    #  wghtdYcolHarmonicLists <- sapply(wghtdYcolVectLists, function(x) sapply(x, function(y) splitHarmonicGroups(y)))
       wghtdYcolHarmonicLists <- sapply(PGMmtrxNoBlankGrpMtrxList,function(x) apply(x,2,function(z) splitHarmonicGroups((1:height)[z<BlkWhtThreshold])))

    yThicknesses <- sapply(wghtdYcolHarmonicLists,function(u) mean(sapply(u,function(cl) mean(sapply(cl,function(h) max(h)-min(h)+1)))))
    aveNharmonics <- sapply(wghtdYcolHarmonicLists,function(u) mean(sapply(u,function(cl) length(cl))))

      
      yMeanDeltasVectors <- sapply(yMeansVectors,function(x) abs(c(x,0)-c(0,x))[-c(1,length(x))])
      
      yMeanRanges <- sapply(yMeansVectors,function(x) range(x)[2]-range(x)[1])
      aveYmeanDeltas <- sapply(yMeanDeltasVectors,mean)

      #whyExclude <- paste(whyExclude,paste("t:",round(yThicknesses*frequency/height,1),"\n",sep=""),"") #not really an exclusion factor
      yRpH <- yMeanRanges/yLens ; yRpH[is.na(yRpH)] <- 0                 #measures inerval
      #yRpTh <- yMeanRanges/yThicknesses ; yRpTh[is.na(yRpTh)] <- 0     
      yDpT <- aveYmeanDeltas/yThicknesses#yRpTh;                        #yDelta per Range   #measures randomness
     # TpLpH <- yThicknesses/(yLens/aveNharmonics)                               #this should be a measure of how tonal something is
      #yDpRpW <- aveYmeanDeltas/yRpTh/widths;                        #yDelta per Range

      VyRangePerHeights[vocNmbr] <- mean(yRpH)  #need to RAISE yRpH to above .3.. especially when harmonic flatenning comes around
    # harmonicEnough <- TpLpH  < .5                 ; whyExclude <- paste(whyExclude,paste("h:",round(TpLpH,1),"\n",sep=""),"")#[vertEnough+1])
      vertEnough <- yRpH  > .25                    ; whyExclude <- paste(whyExclude,paste("v:",round(yRpH,2),"\n",sep=""),"")#[vertEnough+1])
      gradEnough <- yDpT  < 1.75                    ; whyExclude <- paste(whyExclude,paste("g:",round(yDpT,1),"\n",sep=""),"")#[gradEnough+1]) #10 for avYmD/yR
#      gradEnough <- yDpRpW  < 1                    ; whyExclude <- paste(whyExclude,paste("\nv:\n",round(yDpRpW,2),"\n",sep="")) #10 for avYmD/yRpTh
      #the above is largely dependent upon wether or not there are harmonics present...  devide by number of harmonics?
     }else{
       #slopedEnough <- T
       vertEnough <- T
       gradEnough <- T   
     }
      
    if(verbose)cat(" outlier removal -")
      #THROW OUT OUTLIERS BY SIZE OR WIDTH
      unitIsAppropriateToUse <- (wideEnough & nrrwEnough & tallEnough & highEnough & lowEnough & strongEnough & vertEnough &gradEnough) #slopedEnough   
      PGMmtrxNoBlankGrpList <- PGMmtrxNoBlankGrpList[unitIsAppropriateToUse]
      unitNmbrs <- unitNmbrs[unitIsAppropriateToUse]
      PGMmtrxNoBlankGrpMtrxList <- PGMmtrxNoBlankGrpMtrxList[unitIsAppropriateToUse]
      widths <- widths[unitIsAppropriateToUse]
      yMins <- yMins[unitIsAppropriateToUse]
      yMaxs <- yMaxs[unitIsAppropriateToUse]
      yMids <- yMids[unitIsAppropriateToUse]
      yLens <- yLens[unitIsAppropriateToUse]
      energies <- energies[unitIsAppropriateToUse]
   
      nComponents <- length(PGMmtrxNoBlankGrpList)
      #xlim=c(0,duration),ylim=c(0,frequency)

     VpctTotVocEnergryCnsdrd[vocNmbr] <- sum(energies)/totEnergy
 

  
      bmp(file=paste(plotPath,vocFileName,"-chopPlot.bmp",sep=""))#,width=width+200,height=height+200)
      image(x=0:width,y=0:height,z=t(PGMmtrx),xlab="time, s",ylab="frequency, kHz",axes=F)#,asp=height/width)
        title(paste("#",vocNmbr,":",genus,species,"\n call",vocNumber,paste("\"",vocName,"\"",sep="")),
              sub =paste(auth,paste("(",year,")",sep=""),paste(figures,collapse=",")),
              cex.main = 1.2,   font.main= 2, col.main= "black",
              cex.sub = 0.75, font.sub = 4, col.sub = "black")
      axis(1,at = seq(0,width,by = width/duration),labels = seq(0,floor(duration)))   
      axis(2,at = seq(0,height,by= height/frequency),labels = seq(0,floor(frequency)))
   
   
      points((1:width)[blankClmnVect]-.5,rep(1,nBlankCols),pch="|")
     
      points(oldUnitPositions+(.5*oldWidths-.5),oldYmins,pch="-")
      points(oldUnitPositions+(.5*oldWidths-.5),oldYmaxs,pch="-")
      points(unitNmbrs+(.5*widths-.5),yMids,pch="|")

      oldUnits <- length(oldUnitPositions)
      if(oldUnits > 1){
        YmeansXcoord <- NULL
        YmeansYcoord <- NULL
        for(i in 1:oldUnits){
          YmeansXcoord <- c(YmeansXcoord,oldUnitPositions[i]:(oldUnitPositions[i]+oldWidths[i]-1))
          YmeansYcoord <- c(YmeansYcoord,yMeansVectors[[i]])
        }
        points(YmeansXcoord,YmeansYcoord,pch=".")
      }
      stagger <- (((1:length(oldUnitPositions))%%10)+.5)*.10*height  #cascade the heights down a tenth at a time 10 to 1
      text(oldUnitPositions+(.5*oldWidths-.5),stagger,labels=whyExclude,cex=.4)
    
      if(nComponents>0){
        text(unitNmbrs+(.5*widths-.5),4,col=rainbow(length(unitNmbrs),start=0,end=.8),cex=.6)
      }
      dev.off()
      

   
      if(nComponents>1){

        if(verbose)cat(" matching units (by parameters") 
      diffMtrx <- matrix(NA,nComponents,nComponents);
      #PGMmtrxNoBlankGrpMtrxList <- PGMmtrxNoBlankGrpMtrxList[widths]
      
      for(i in 1:nComponents){
        
        for(j in i:nComponents){
          nrrwMtrx <- PGMmtrxNoBlankGrpMtrxList[[i]] #start off assuming that i is more narrow
          wideMtrx <- PGMmtrxNoBlankGrpMtrxList[[j]]
          if(widths[i] > widths[j]){
            temp <- nrrwMtrx;nrrwMtrx <- wideMtrx; wideMtrx <- temp; #swap if necessary
            nrrwWdth <- widths[j];wideWdth <- widths[i];
            nrrwEnergy <- energies[j];wideEnergy <- energies[i];
          }else{
            nrrwWdth <- widths[i];wideWdth <- widths[j];
            nrrwEnergy <- energies[i];wideEnergy <- energies[j];
          }
          widthDiff <- wideWdth-nrrwWdth
          energyDiff <- abs(energies[j]-energies[i])#wideEnergy-nrrwEnergy)
          yMaxDiff <- abs(yMaxs[j]-yMaxs[i])#wideyMax-nrrwyMax)
          yMinDiff <- abs(yMins[j]-yMins[i])#wideyMin-nrrwyMin)
          yLenDiff <- abs(yLens[j]-yLens[i])
          yMidDiff <- abs(yMids[j]-yMids[i])
          longLen <- max(yLens[i],yLens[j])
          diffVect <- vector("numeric",length=widthDiff+1);
#          if((widthDiff > (2.5*nrrwWdth)) | (energyDiff > (2*nrrwEnergy)) | (yLenDiff > (.25*longLen)) | (yMidDiff > (.2*height))){#(yMaxDiff > (.27*longLen)) | (yMinDiff > (.27*longLen)) ){
           if((widthDiff > (2*nrrwWdth)) | (energyDiff > (2.5*nrrwEnergy)) | (yLenDiff > (.3*longLen)) | (yMidDiff > (.3*height))){#(yMaxDiff > (.27*longLen)) | (yMinDiff > (.27*longLen)) ){

            diffMtrx[i,j] <- NA
          }else{
            for(k in 1:(widthDiff+1)){  
              diffVect[k] <- sum(abs(nrrwMtrx-wideMtrx[,k:(nrrwWdth+k-1)]))
              #image(t(wideMtrx[,k:(nrrwWdth+k-1)]))
              #image(t(nrrwMtrx))
              #diffSquare <- nrrwMtrx-wideMtrx[,k:(nrrwWdth+k-1)]
              #image(t(diffSquare))
            }
            diffMtrx[i,j] <- min(diffVect)/nrrwEnergy #(height*nrrwWidth)
          }
        } 
      }
      diag(diffMtrx) <- NA

      
      VaveDiff[vocNmbr] <- mean(diffMtrx,na.rm=T)
        if(verbose)cat(" & by subtraction) -")   
      
      defaultGrySclThreshold <- 1#/maxGrayValue

    if(useDiffThresholding){
      maxGrySclThreshold <- 1.5#/maxGrayValue      what about 1.7??->     Cebuella-pygmaea-11-Alerting.Call.B-Pola1975-chopPlot.bmp  
      minGrySclThreshold <- .5#/maxGrayValue
      xMin <- defaultGrySclThreshold  #important re-initialization (so that we don't reuse the xMin from the previous loop)
      #par(mfrow=c(2,2))
      if(sum(!is.na(diffMtrx))>2){
      #hist(diffMtrx,breaks=20)
        diffDensity <- density(diffMtrx,na.rm=T)
        f <- function(x) x[2] > max(x[-2]) 
        minMaxVect <- rapply(zoo(rank(diffDensity$y, ties = "first")), 3, function(x) f(x) - f(-x)) # uses zoo library
        diffDensXY <- cbind(diffDensity$x,diffDensity$y)
        minMtrxXY <- diffDensXY[minMaxVect==-1,]  #only one of the  mins if there are multiple matches
        if(length(minMtrxXY)>2){
          xMin <- minMtrxXY[minMtrxXY[,2]==min(minMtrxXY[,2]),1][1]  #of all the local minima.. find the lowest one. #if several =0 take first
        }else{
          xMin <- minMtrxXY[1]
        }
        
        maxMtrxXY <- diffDensXY[minMaxVect== 1,]

        bmp(file=paste(plotPath,vocFileName,"-matchDist.bmp",sep=""))#height)
        plot(diffDensity,main=paste("Distribution of Distance \n Between Repeated Component Candidates in,\n",genus,species,"call",vocNumber,"\"",vocName,"\""),xlab="Difference, gray values / pixel")
        
        if(!is.na(xMin)){
            abline(v=xMin); xMinIsNumber <- T
            text(c(((xMin-min(diffDensXY[,1]))/2),((xMin+max(diffDensXY[,1]))/2)),.01,c("matches\n(repeats)","mismatches\n(non-repeats)"))
            dev.off()
            VunitDiffThreshold[vocNmbr] <- xMin    #for use later to find an appropriate average threshold (instead of the current 13)
        }else{
          xMin <- defaultGrySclThreshold  #default empirically determined #grayScaleValues (out of maxGray) per pixel to count as a shape match
          xMinIsNumber <- F
          dev.off()
          VunitDiffThreshold[vocNmbr] <- NA    #for use later to find an appropriate average threshold (instead of the current 13)
        }
      if(xMinIsNumber & ((xMin> maxGrySclThreshold)|(xMin < minGrySclThreshold))){
        xMin <- defaultGrySclThreshold
        VunitDiffThreshold[vocNmbr] <- NA
      }
      }
      #print(paste("xMin =", xMin))
    }else{
      xMin <- 1
    }#end if useDiffThresholding
      
  GSpPXthresh <- xMin #use a histogram???? a distribution?  minimum between 5 and 25 
      VrepMatches[vocNmbr] <- repMatches <- sum(as.numeric(diffMtrx<GSpPXthresh),na.rm=T)
  #  if(is.na(repMatches))
   #   repMatches <- 0      #in cases where there is only one unit in the graph
      


      if(repMatches==0 ){
VrepMatches[vocNmbr] <- VnClusters[vocNmbr] <- VnRepGroups[vocNmbr] <- VprcntCmpntsRptd[vocNmbr] <- VprcntEnergyRptd[vocNmbr] <- VaveClustSize[vocNmbr] <- VaveClustDegree[vocNmbr] <- VaveReps.group[vocNmbr] <- VaveDegree.group [vocNmbr] <- VaveSaturation[vocNmbr] <- VwAveSaturation[vocNmbr] <- 0
VaveIntervalVar[vocNmbr] <- VwAveIntervalVar[vocNmbr] <- 0  #arbitrarily low.. labeled incorrectly now.. should be precision vector

        }else{
 if(verbose)cat(" calculating stats -")
          
      net <- network(symmetrize(diffMtrx<GSpPXthresh),directed=F)

     
      clustMembershipVect <- component.dist(net)$membership
      degreeVect <- degree(net,cmode="indegree")
      
      clustSizeVect <- component.dist(net)$csize;     repVect <- clustSizeVect   #are isolates 1 repeat or zero or NA?
      aveClustDegreeVect <- sapply(split(degreeVect,clustMembershipVect),mean); aveRepGroupDegreeVect <- aveClustDegreeVect
      aveClustEnergyVect <- sapply(split(energies,clustMembershipVect),sum);  #should this be mean? change crashes tho.. just a weight?

      VnClusters[vocNmbr] <- nClusters <- components(net) #includes isolates   #repGroups don't include isolates
      isRepGroup <- clustSizeVect>1
      VnRepGroups[vocNmbr] <- nRepGroups <- sum(isRepGroup) #only counts clusters with 2 or more repeats
      VprcntCmpntsRptd[vocNmbr] <- prcntCmpntsRptd <- nRepGroups/nClusters
      VprcntEnergyRptd[vocNmbr] <- prcntEnergyRptd <- sum(sapply(split(energies,clustMembershipVect),sum)[isRepGroup])/sum(energies)
      
      repVect[!isRepGroup] <- NA
      aveRepGroupDegreeVect[!isRepGroup] <- NA
      VaveReps.group[vocNmbr] <- aveReps.group <- mean(repVect,na.rm=T)
      VaveDegree.group[vocNmbr] <- aveDegree.group <- mean(aveRepGroupDegreeVect,na.rm=T)
      VaveClustSize[vocNmbr] <- aveClustSize <- mean(clustSizeVect)
      VaveClustDegree[vocNmbr] <- aveClustDegree <- mean(aveClustDegreeVect)
      
      saturationVect <- aveRepGroupDegreeVect/choose(repVect,2)
      VaveSaturation[vocNmbr] <- aveSaturation <- mean(saturationVect,na.rm=T)
      VwAveSaturation[vocNmbr] <- wAveSaturation <- weighted.mean(saturationVect,aveClustEnergyVect,na.rm=T) # degree of repetition is weighted by energy put in to them
      
      nodeNumbers <- sapply(split(1:nComponents,clustMembershipVect),paste,collapse=",")
      
      
      #yMinVarVect <- sapply(split(yMins,clustMembershipVect),var)     # temporary until i figure out a way to reject entirely zero min reps
      #yMaxVarVect <- sapply(split(yMaxs,clustMembershipVect),var)
      yMaxPrecision.repGroup <- sapply(split(yMaxs,clustMembershipVect),function(x) exp(-var(x,na.rm=T)))
      #yIntervalVarMtrx <- rbind(yMinVarVect,yMaxVarVect)               # temporary until i figure out a way to reject entirely zero min reps
      yIntervalVarVect <- yMaxPrecision.repGroup  #apply(yIntervalVarMtrx,2,mean)  # temporary until i figure out a way to reject entirely zero min reps
      clusterStatsDF <- data.frame(Cluster.Size=clustSizeVect,Repeats=repVect,Average.Degree=aveRepGroupDegreeVect,Interval.Variance=yIntervalVarVect,Nodes.In.Cluster=nodeNumbers,saturation=saturationVect)

      VaveIntervalVar[vocNmbr] <- aveIntervalVar <- mean(yIntervalVarVect,na.rm=T)
      VwAveIntervalVar[vocNmbr] <- wAveIntervalVar <- weighted.mean(yIntervalVarVect,saturationVect,na.rm=T) #precision is weighted by confidence in actual repetition   #not used
       if(verbose)cat(" graphing -")
       bmp(file=paste(plotPath,vocFileName,"-matchGraph.bmp",sep=""))#height)
      plot.network(net,vertex.cex=energies/mean(energies),vertex.col=rainbow(nComponents,start=0,end=.8),displaylabels=T,boxed.labels=F,label.cex=.65,label.col=rainbow(nComponents,start=0,end=.8),main=paste("Graph of Repeat Groupings in\n",genus,species,"call",vocNumber,"\"",vocName,"\"","\n (paired if difference <",round(GSpPXthresh,1),"gray values/pixel)"),sub=paste(nRepGroups,"Repetition Groups (",round(VprcntEnergyRptd[vocNmbr]*100),"% of graph",round(VprcntEnergyRptd[vocNmbr]*VpctTotVocEnergryCnsdrd[vocNmbr]*100),"% of voc)with\n ",round(aveReps.group,1),"reps and",round(aveDegree.group,1),"degree (average) per group \n (Energy Weighted Mean Saturation of",round(wAveSaturation*100,1),"%) &\na Saturation Weighted Mean Interval Variance of",round(wAveIntervalVar,1),"\n (this graph represents",round(VpctTotVocEnergryCnsdrd[vocNmbr]*100),"% of the total vocalization)"))
       dev.off()
         }
    }else{
      if(verbose)cat("only one or fewer component found... no matches")
VrepMatches[vocNmbr] <- VnClusters[vocNmbr] <- VnRepGroups[vocNmbr] <- VprcntCmpntsRptd[vocNmbr] <- VprcntEnergyRptd[vocNmbr] <- VaveClustSize[vocNmbr] <- VaveClustDegree[vocNmbr] <- VaveReps.group[vocNmbr] <- VaveDegree.group [vocNmbr] <- VaveSaturation[vocNmbr] <- VwAveSaturation[vocNmbr] <- 0
VaveIntervalVar[vocNmbr] <- VwAveIntervalVar[vocNmbr] <- 0  #arbitrarily low.. labeled incorrectly now.. should be precision vector
    }
     if(verbose)cat("\n")
}#end loop through vocalizations


pctUnitsUsed <- VpctTotVocEnergryCnsdrd/mean(VpctTotVocEnergryCnsdrd,na.rm=T)
pctUsedUnitsGrouped <- VprcntEnergyRptd/mean(VprcntEnergyRptd,na.rm=T)
groupageWeight <- (VpctTotVocEnergryCnsdrd * VprcntEnergyRptd)   #UnitsUsed*pctUsedUnitsGrouped


wideVocs <- data.frame(vocs,VrepMatches, VnClusters, VnRepGroups, VprcntCmpntsRptd, VprcntEnergyRptd,VaveClustSize, VaveClustDegree, VaveReps.group, VaveDegree.group , VaveSaturation, VwAveSaturation , VaveIntervalVar, VwAveIntervalVar,spectrosPerVoc, VunitDiffThreshold, VaveDiff,VpctTotVocEnergryCnsdrd,VyRangePerHeights, VaveMaxPrecision,groupageWeight)

write.table(wideVocs,"priVocTable-gramParsed.txt",quote=F,sep="\t")


repMtrx <- data.frame(vocs$prRep,vocs$prMelody,VrepMatches, VnClusters, VnRepGroups, VprcntCmpntsRptd,VprcntEnergyRptd, VaveClustSize, VaveClustDegree, VaveReps.group, VaveDegree.group , VaveSaturation, VwAveSaturation , VaveIntervalVar, VwAveIntervalVar,VaveDiff,VyRangePerHeights)



#plot(repMtrx)
par(mfrow=c(1,2))
#clustDegMod <- lm(wideVocs$prMelody~wideVocs$VaveClustDegree + wideVocs$isTrill + wideVocs$duration + wideVocs$spectrosPerVoc + groupageWeight);summary(clustDegMod)
#plot(wideVocs$prMelody,log(wideVocs$VaveClustDegree),col=2,cex=groupageWeight))  #my personal subjective plot (best measure of how closely flat repeated
#text(wideVocs$prMelody,wideVocs$VaveClustDegree,wideVocs$genus,cex=groupageWeight))
#abline(clustDegMod) 

#clustSizeMod <- lm(wideVocs$prMelody~wideVocs$VaveClustSize + wideVocs$isTrill + wideVocs$duration + wideVocs$spectrosPerVoc + groupageWeight);summary(clustSizeMod)
#plot(wideVocs$prMelody,log(wideVocs$VaveClustSize),col=2,cex=groupageWeight))              #best with interval detection included (pre harmonic collapse or trill busting)
#text(wideVocs$prMelody,wideVocs$VaveClustSize,wideVocs$genus,cex=groupageWeight))
#abline(clustSizeMod)

#intrvlVarMod <- lm(wideVocs$prMelody ~ wideVocs$VwAveIntervalVar + wideVocs$isTrill + wideVocs$duration + wideVocs$spectrosPerVoc + groupageWeight);summary(intrvlVarMod)
#plot(wideVocs$prMelody,log(wideVocs$VwAveIntervalVar),col=2,cex=groupageWeight))
#text(wideVocs$prMelody,wideVocs$VwAveIntervalVar,wideVocs$genus,cex=groupageWeight))
#abline(intrvlVarMod)


intrvlVarMod <- lm(wideVocs$prMelody ~ wideVocs$VaveMaxPrecision + wideVocs$isTrill + wideVocs$duration + wideVocs$spectrosPerVoc + groupageWeight);summary(intrvlVarMod)
plot(wideVocs$prMelody,log(wideVocs$VaveMaxPrecision),col=2,cex=groupageWeight)
#text(wideVocs$prMelody,log(wideVocs$VaveMaxPrecision),wideVocs$genus,cex=groupageWeight))
abline(intrvlVarMod)



#diffMod <- lm(wideVocs$prMelody~wideVocs$VaveDiff + wideVocs$isTrill + wideVocs$duration + wideVocs$spectrosPerVoc  + groupageWeight);summary(diffMod)
#plot(wideVocs$prMelody,wideVocs$VaveDiff,col=2,cex=groupageWeight))
#text(wideVocs$prMelody,wideVocs$VaveDiff,wideVocs$genus,cex=groupageWeight))
#abline(diffMod)


diffMod <- lm(wideVocs$prMelody~wideVocs$VnRepGroups + wideVocs$isTrill + wideVocs$duration + wideVocs$spectrosPerVoc + groupageWeight);summary(diffMod)
plot(wideVocs$prMelody,wideVocs$VnRepGroups*wideVocs$groupageWeight,col=0)
text(wideVocs$prMelody,wideVocs$VnRepGroups*wideVocs$groupageWeight,c("/","T")[wideVocs$isTrill+1])
#text(wideVocs$prMelody,wideVocs$VnRepGroups,wideVocs$genus,cex=groupageWeight)
abline(diffMod)

summary(lm(wideVocs$prMelody~  VnRepGroups + groupageWeight +VaveMaxPrecision + wideVocs$duration))# + wideVocs$isTrill))


library(stats)
repMeldCor <- cor(repMtrx,use = "pairwise.complete.obs")[1:2,];repMeldCor
#apply(repMeldCor,2,sum)

