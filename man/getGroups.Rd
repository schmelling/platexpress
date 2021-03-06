% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/platereader.R
\name{getGroups}
\alias{getGroups}
\title{group wells by experiment annotations (in plate map file)}
\usage{
getGroups(plate, by = "medium", verb = TRUE)
}
\arguments{
\item{by}{a list of column IDs of the plate layout}
}
\value{
Returns a list of well IDs for the identified grouping. This list
can be used in viewGroups(data,groups) to summarize data for these groups.
}
\description{
group wells by experiment annotations (in plate map file)
}
\details{
Calculates the distinct groups from the plate layout by the selected
experimental parameters.
}
\seealso{
\code{\link{readPlateMap}}, \code{\link{viewGroups}}
}

