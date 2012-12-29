
\name{unitstats}
\alias{unitstats}
\alias{unitYstats}
\title{Calculate Spectrogram Unit Stats}
\description{
Calculate per-unit statistics such as height, width, mean y-value, mid x-value, mid y-value, interval, tone, etc.  
}
\usage{
unitYstats(sg, harmonics=FALSE)
unitstats(sg, harmonics=FALSE)
}
\arguments{
\item{sg}{a spectrogram object.}
\item{harmonics}{boolean for cleaning out upper harmonics and to estimate stats from just the fundamental frequency }
}
\value{
spectrogram object with stats
}
\examples{

load(system.file('extdata','calls','Tarsius-spectrum-Duet.f-split.RData', package='melody'))
s <- unitstats(s)
s <- tone(s)
s <- unitYstats(s, harmonics=c(FALSE,TRUE)[(s$tone > .75) +1])

unit.param.names <- grep('^u', names(s), value=TRUE)
print(sapply(unit.param.names, function(x) s[[x]]))


}
