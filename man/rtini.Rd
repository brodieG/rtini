% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rtini.R
\docType{package}
\name{rtini}
\alias{rtini}
\alias{rtini-package}
\title{Right Triangle Irregular Network Surface Approximation}
\description{
A vectorized implementation of the Right-Triangulated Irregular Networks
(RTIN) surface approximation algorithm.  This implementation generates meshes
that are triangulated subsets of maps, subject to error tolerances.  The
vertices in the reduced mesh are vertices in the original, and they are
arranged to form right angle triangles.
}
\details{
The workhorse functions are \code{\link[=rtini_error]{rtini_error()}} and \code{\link[=rtini_mesh]{rtini_mesh()}}.  The former
computes the surface approximation error at all levels of approximation, and
the latter extract mesh triangles subject to an error tolerance.

This is package is experimental, lightly tested, and not intended for
production use.
}
