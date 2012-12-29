\name{clean}
\alias{clean}
\title{Clean a Spectrogram}
\description{
This function provides automated cleaning of outlier (too brief, too quiet, etc.) spectrographic units.
}
\usage{
clean(sg, wlim=c(3,sg$width), ylim=c(0,1), ylen=3, imin=0.01)
}
\arguments{
\item{sg}{a spectrogram object.}
\item{wlim}{min and max limit on width}
\item{ylim}{min and max limit on y}
\item{ylen}{min y length (unit height)}
\item{imin}{min intensity of unit}
}
\value{
a cleaned spectrogram object
}
\examples{

load(system.file('extdata','calls','Tarsius-spectrum-Duet.f-split.clean.RData', package='melody'))

length(s$units)
s <- clean(s)
length(s$units)
}
