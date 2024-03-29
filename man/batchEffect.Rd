% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/batchEffect.R
\name{batchEffect}
\alias{batchEffect}
\title{Analyse possible batch effect}
\usage{
batchEffect(project, automataID)
}
\value{
A CSV file with proportions of TSS; a CSV filr with proportions of sample type: a
CSV file with proportions of TSS and sample type
}
\description{
Compares proportion of different TSS and sample type for each LPD group to see if there
is a significant diffence between them.
}
