% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/main_functions.r
\name{fit.cd19.sub}
\alias{fit.cd19.sub}
\title{Fit model for CD19 B cell subtypes for whole blood}
\usage{
fit.cd19.sub(data, fit.main, cpg.list = cpg.list.cd19)
}
\arguments{
\item{data}{A data frame}

\item{fit.main}{A numeric matrix}

\item{cpg.list}{A character}
}
\value{
A numeric matrix
}
\description{
Fit model for CD19 B cell subtypes for whole blood
}
\seealso{
Other Model fitting: \code{\link{MethylDeconBloodSubtypes}},
  \code{\link{fit.cd19.sorted}},
  \code{\link{fit.cd4.sorted}},
  \code{\link{fit.cd8.sub.sorted}},
  \code{\link{fit.cd8.sub}}, \code{\link{fit.main.nosum1}},
  \code{\link{fit.stage2b}},
  \code{\link{get.cpg.list.cd19.sub}}
}

