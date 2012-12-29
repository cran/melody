
\name{graph}
\alias{graph}
\title{Plot Spectrogram as a Graph}
\description{
Plot a mathematical graph showing the clustering patterns of the units and grouping different syllable types.
}
\usage{
graph(sg, label.cex = 0.7, ...)
}
\arguments{
\item{sg}{a spectrogram object.}
\item{label.cex}{character expansion labels for graph nodes}
\item{...}{additional parameters passed to underlying network:::plot}
}
\value{
a plot of the mathematical graph representation of units of the spectrogram
}
\examples{

load(system.file('extdata','calls','Tarsius-spectrum-Duet.f-split.clean.clst.RData', package='melody'))

par(mfrow=c(1,2))
plot(s)
graph(s)

}
