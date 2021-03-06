% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/helper_functions.R
\name{ldobf}
\alias{ldobf}
\title{Lan-DeMets O'Brien-Fleming alpha-spending function}
\usage{
ldobf(t, alpha_one_sided = 0.025)
}
\arguments{
\item{t}{Information fraction. Vector of numbers between 0 and 1.}

\item{alpha_one_sided}{One-sided alpha level.}
}
\value{
A vector of cumulative alpha spend.
}
\description{
Find the cumulative alpha-spend at various information times according to Lan-DeMets OBF function.
}
