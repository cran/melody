
\name{partition}
\alias{partition}
\title{Partition a Spectrogram}
\description{
Split a spectrogram matrix at relatively quiet breakpoints along the temporal axis.
}
\usage{
partition(sg, ms=1, ds=1.62, lambda=5)
}
\arguments{
\item{sg}{a spectrogram object.}
\item{ms}{mean scalar}
\item{ds}{deviance scalar}
\item{lambda}{exponential parameter determining spectrogram bottom (low-frequency) subtraction factor}
}
\value{
a partitioned spectrogram object
}
\examples{

load(system.file('extdata','calls','Tarsius-spectrum-Duet.f.RData', package='melody'))
ls()

par(mfrow=c(1,2))
plot(s, main='before')
s <- partition(s)
plot(s, main='after')

}

