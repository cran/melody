
\name{threshold}
\alias{threshold}
\title{Calculate Background Grayscale Level}
\description{
This function calculates the threshold level between foreground and background and assigns it to a slot in the spectrogram object.
}
\usage{
threshold(sg, type=c('bw','div'), pct.max=0.89)
}
\arguments{
\item{sg}{a spectrogram object.}
\item{type}{the type of thresholding calculation to perform}
\item{pct.max}{percentage of the maximum gray value to use for thresholding.}
}
\value{
sg object
}
\examples{

s <- sg(readPGM(system.file('extdata','calls','Tarsius-spectrum-Duet.female-Nietsch1998-2b2.pgm', package='melody')))

names(s)
s <- threshold(s, pct.max = 0.89)
names(s)

print(s$bg.threshold)


}


