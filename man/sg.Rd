
\name{sg}
\alias{sg}
\title{Create a Spectrogram Object}
\description{
Create an S3 Spectrogram object.   This function is a simple instantiation function that performs checks and populates a properly sized 2-D matrix of grayscale values. 
}
\usage{
sg(x, ...)
}
\arguments{
\item{x}{a matrix representing a grayscale spectrogram}
\item{...}{additional parameters passed on to list}
}
\value{
a spectrogram object
}
\examples{

pgm <- readPGM(system.file('extdata','calls','Tarsius-spectrum-Duet.female-Nietsch1998-2b2.pgm', package='melody'))
spectrogram.object <- sg(pgm)

str(spectrogram.object)
plot(spectrogram.object)

}

