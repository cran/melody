readPGM <- function(file){
 
  if (is.character(file)) {
    file <- file(file, "rt")
    on.exit(close(file))#, new = TRUE)
  }
  if (!inherits(file, "connection")) 
    stop("argument 'file' must be a character string or connection")
  if (!isOpen(file)) {
    open(file, "rt")
    on.exit(close(file))#, new = TRUE)
  }
  
  id <- readLines(file,1) #first 2 lines
  temp <- readLines(file,1) #first 2 lines
  if(grepl("^#", temp))  # second line was comment
    hw.line <- readLines(file,1)
  else
    hw.line <- temp     
  #need to make this optional ... if there is a # sign in front
  lineChars <- strsplit(hw.line, " ")[[1]]
  lineVect <- lineChars[lineChars != ""]
  w <- as.numeric(lineVect[1])
  h <- as.numeric(lineVect[2])
  maxVal <- readLines(file,1) #max gray value
  vect <- scan(file,what="numeric",nmax=w*h) 
  matrx <- matrix(as.numeric(vect),nrow=h,ncol=w,byrow=T)
  #matrxI <- -(matrx-255)
  matrx <-matrx[h:1,]
  matrx  #directoryPath <- baseDir
}

