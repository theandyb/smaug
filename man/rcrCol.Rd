% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bed_funs.r
\name{rcrCol}
\alias{rcrCol}
\title{Get recombination rate at each site}
\usage{
rcrCol(sites, rcrfile)
}
\arguments{
\item{sites}{A dataframe with 10 columns (CHR, POS, Sequence, Dummy Variables for the 6 mutation types, and read depth)}

\item{rcrfile}{Path to file containing recombination rates}
}
\value{
Vector containing recombination rates in order matching input dataset
}
\description{
Get recombination rate at each site
}
