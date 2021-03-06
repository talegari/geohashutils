% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{gh_cover}
\alias{gh_cover}
\title{Get a set of geohashes covering a set of object}
\usage{
gh_cover(sfs, precision = NULL)
}
\arguments{
\item{sfs}{(list of sf objects)}

\item{precision}{(integer) precision on the covering geohashes}
}
\value{
(list) of covering geohashes
}
\description{
Get a set of geohashes covering a set of object
}
\examples{
fname = system.file("shape/nc.shp", package="sf")
temp = sf::st_read(fname)
temp
gu = geohashTools::gh_to_sf(gh_cover(temp, 4)[[1]])
mapview::mapview(list(gu, temp))
# comparse with 'geohashTools'
mapview::mapview(geohashTools::gh_covering(temp, 4))
}
