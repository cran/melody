
\name{readPGM}
\alias{readPGM}
\title{Read a PGM-3 Formatted Image File.}
\description{
Reads a PGM version 3 formatted text image file into R as a two-dimensional matrix.
}
\usage{
readPGM(file)
}
\arguments{
\item{file}{file path to read}
}
\value{
an image matrix of grayscale values
}
\examples{

tarsier.song <- system.file('extdata','calls','Tarsius-spectrum-Duet.female-Nietsch1998-2b2.pgm', package='melody')
spectrogram.matrix <- readPGM(tarsier.song)

dim(spectrogram.matrix)
image(t(spectrogram.matrix))

}


