
\name{sg-class}
\alias{sg-class}
\title{Spectrogram Class}
\description{
An R class for storing a spectrogram matrix and associated statistics on units and clustering.
}

\section{Slots/List Components}{
Objects of this class contain (at least) the following list components:
  
  \code{sg}: numeric matrix containing the grayscale spectrographic values.
  
  \code{height}: height of the spectrogram matrix (sg).
  \code{width}: width of the spectrogram matrix (sg).
  
}

\section{Methods}{
This class inherits directly from class \code{list} so any operation appropriate for lists will work on objects of this class.
\code{sg} objects also have a \code{show} method.
}

\keyword{classes}
