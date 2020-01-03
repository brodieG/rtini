#' Right Triangle Irregular Network Surface Approximation
#'
#' A vectorized implementation of the Right-Triangulated Irregular Networks
#' (RTIN) surface approximation algorithm.  This implementation generates meshes
#' that are triangulated subsets of maps, subject to error tolerances.  The
#' vertices in the reduced mesh are vertices in the original, and they are
#' arranged to form right angle triangles.
#'
#' The workhorse functions are [rtini_error()] and [rtini_mesh()].  The former
#' computes the surface approximation error at all levels of approximation, and
#' the latter extract mesh triangles subject to an error tolerance.
#'
#' This is package is experimental, lightly tested, and not intended for
#' production use.
#'
#' @import vetr
#' @name rtini
#' @docType package

NULL

