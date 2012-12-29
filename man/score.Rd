
\name{score}
\alias{tone}
\alias{interval}
\title{Score Tonality and Interval of a Spectrogram}
\description{
'pctbg': percent of pixels in units that are considered 'background'
'ymRpuH' y mean rise per unit height
}
\usage{
tone(sg, method=c('pctbg', "iqsub.pctbg", 'yDpuH'), intensity.weighted=TRUE, w2h.min=0.1)
interval(sg, method=c('ymRpuH'), intensity.weighted=TRUE)
}
\arguments{
\item{sg}{a spectrogram object.}
\item{method}{the method to use for scoring: pctbg (default) or yDpuH.}
\item{intensity.weighted}{boolean specifying if individual unit scores should be weighted by their intensities in determining the overall tonality for the spectrogram}
\item{w2h.min}{minimum threshold for the ratio of width to height of a unit before it gets penalized as atonal}
}
\value{
a spectrogram with tone and/or interval scores
}
\examples{

load(system.file('extdata','calls','Tarsius-spectrum-Duet.f-split.clean.RData', package='melody'))

s <- tone(s, method='pctbg', intensity.weighted=TRUE)
print(s$tone)

s <- interval(s, method='ymRpuH', intensity.weighted=TRUE)
print(s$interval)


}

