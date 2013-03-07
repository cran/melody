
\name{plot.sg}
\alias{plot.sg}
\alias{spectrogram}
\title{Plot a Spectrogram Object}
\description{
This function plots a two dimensional spectrogram object as a heatmap image with default labels of time (x) versus frequency (y).  Depending on which unit stats are available in the object, they will be plotted as annotative markings on the units.  Additionally, a colored sequence of numbers is plotted below the units for later reference in a 'graph' (clustering) plot. 
}
\usage{
\method{plot}{sg}(x, blanks=TRUE, u.ranges=TRUE, u.yzxs=TRUE, label.cex = 0.7, xlab='Time', ylab='Frequency', ...)
}
\arguments{
\item{x}{a spectrogram object}
\item{blanks}{boolean specifiying if the blank columns should be indicated with a black "rug" along the bottom of the plot}
\item{u.ranges}{boolean specifying if the unit min and max should be indicated with pairs of "-"s }
\item{u.yzxs}{boolean specifying if the per column y averages should be indicated as a string of "."s}
\item{label.cex}{character expansion factor for unit labels}
\item{xlab}{x axis label}
\item{ylab}{y axis label}
\item{...}{additional parameters passed on to "plot"}
}
\value{
a spectrogram plot
}
\examples{

load(system.file('extdata','calls','Tarsius-spectrum-Duet.f.RData', package='melody'))

par(mfrow=c(1,2))
plot(s, main='before')
s <- partition(s)
s <- unitstats(s)
plot(s, main='after')


}
