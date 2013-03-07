\name{cluster}
\alias{cluster}
\title{Cluster Spectrographic Units}
\description{
This funtion matches similar units into clusters: simultaneously calculating the number of syllables (a cluster count) and repetition (the average degree per cluster).  Both measures have the option of being weighted by the amplitude of the nodes involved in the calculation of each. 
}
\usage{
cluster(sg, method=c('limits','avg'), dist.limits = nv(
c(   0.5,        0.6,   0.46, 0.07,  0.15,  0.4,   0.2,   0.2,   0.55,    0.3),
c('shape','intensity','width','ymean','ymid','ylen','ymin','ymax','ymRpuH','pctbg')
), 
dist.avg=.15, rep.score=c('c.nreppct','c.degpct','cedge.dens'), syl.score=c('c.ct'), intensity.weighted=TRUE)
}
\arguments{
\item{sg}{a spectrogram object.}
\item{method}{which edge determining method to use}
\item{dist.limits}{distance limits used to determine associations (or "edges") between spectrographic units (or "nodes")}
\item{dist.avg}{for method == 'avg' the single numeric threshold value determining edges between unit nodes}
\item{rep.score}{scoring method used to determine unit repetition. default is 'cnreppct' }
\item{syl.score}{scoring method used to determine unit variability. default is 'c.ct' or 'cluster count'}
\item{intensity.weighted}{boolean specifying if scores should be weighted by unit intensities}
}
\value{
a clustered spectrogram object
}
\examples{

load(system.file('extdata','calls','Tarsius-spectrum-Duet.f-split.clean.RData', package='melody'))

## match the units of the spectrgram to each other, clustering like shapes
s <- cluster(s, method='limits')  

print(s$repetition)  # a score between 0 and 1
print(s$syllable_ct)  # a count of different syllables (unit types)

}
