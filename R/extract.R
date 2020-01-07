## Copyright (C) 2020  Brodie Gaslam
##
## This file is part of "rtini - Right Triangle Irregular Network Surface
## Approximation"
##
## This program is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 2 of the License, or
## (at your option) any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## Go to <https://www.r-project.org/Licenses/GPL-2> for a copy of the license.

#' Extract A Triangle Mesh With Approximation Error Below Tolerance
#'
#' Given a surface approximation error matrix as computed with [rtini_error()],
#' extract a triangle mesh that approximates the original heightmap with error
#' less than or equal to a specified tolerance.
#'
#' @param error a numeric matrix as produced by [rtini_error()]
#' @param tol scalar positive numeric, the error tolerance
#' @return a list of numeric matrices.  Each column in a matrix corresponds to
#'   the "linearized" vertex coordinates of a triangle.  See examples for how to
#'   convert to x-y coordinates.  Each matrix corresponds to a different
#'   approximation layer.
#' @export
#' @seealso [rtini_error()]
#' @examples
#' err <- rtini_error(volcano)
#' mesh <- rtini_extract(err, tol=10)
#'
#' ## Convert linearized indices to x/y in [0,1]
#' tris <- rbind(do.call(cbind, mesh), NA)
#' plot.new()
#' polygon(
#'   x=((tris - 1) %/% nrow(err)) / (ncol(err) - 1),
#'   y=((tris - 1) %% nrow(err)) / (nrow(err) - 1),
#'   col='grey90'
#' )
#'
#' ## One benefit of the linearized indices is that you can
#' ## index directly into the height and error matrices.  Here
#' ## We use that to show the errors at every point that was not
#' ## plotted versus volcano height.
#' plot(volcano[-unlist(mesh)], err[-unlist(mesh)])

rtini_extract <- function(error, tol) {
  vetr(matrix(numeric(), 0, 0) && dim(.) > 2 && dim(.) %% 2 > 0, NUM.1.POS)
  nr <- nrow(error)
  nc <- ncol(error)
  layers <- floor(min(log2(c(nr, nc) - 1L)))
  res <- vector('list', 2L * layers + 1L)
  id.dat <- id.pass <- id.fail <- seed_ids(0L, 0L, 0L, FALSE)
  trow <- tcol <- 0L

  for(i in seq_len(layers)) {
    # seed initial and additional triangles
    m <- as.integer(2^(layers - i))
    trow.p <- trow * 2L
    tcol.p <- tcol * 2L
    trow <- ((nr - 1L) %/% (m * 2L))
    tcol <- ((nc - 1L) %/% (m * 2L))
    row.ex <- trow.p < trow
    col.ex <- tcol.p < tcol
    seed.h <- seed.v <- replicate(4L, integer(), simplify=FALSE)
    if(i == 1L) {
      if(trow < tcol) seed.h <- seed_ids(0L, m, tcol, TRUE)
      else seed.v <- seed_ids(0L, m, trow, FALSE)
    } else {
      if(row.ex) seed.h <- seed_ids(trow.p * m * 2L, m, tcol.p + col.ex, TRUE)
      if(col.ex) seed.v <- seed_ids(tcol.p * m * 2L, m, trow.p, FALSE)
    }
    id.dat[] <- Map(c, id.dat, seed.h, seed.v)

    # diag first, then axis
    for(j in 1:2) {
      ids <- id.dat[['x','tar']] * nr + id.dat[['y','tar']] + 1L
      pass <- error[ids] <= tol
      wpass <- which(pass)
      wfail <- which(!pass)

      id.pass[] <- lapply(id.dat, '[', wpass)
      id.fail[] <- lapply(id.dat, '[', wfail)

      res[[i * 2L - (j == 1L)]] <- extract_tris(id.pass, nr)

      id.dat <- next_children(id.fail)
  } }
  # Any remaining failures must be split to lowest level

  res[[2L * layers + 1L]] <- extract_tris(id.dat, nr)
  res
}


# helper functions for extraction

extract_tris <- function(id, nr) {
  dx <- id[['x','tar']] - id[['x','par']]
  dy <- id[['y','tar']] - id[['y','par']]

  res.x <- rbind(id[['x','par']], id[['x','tar']] + dy, id[['x','tar']] - dy)
  res.y <- rbind(id[['y','par']], id[['y','tar']] - dx, id[['y','tar']] + dx)

  res.x * nr + res.y + 1L
}
next_children <- function(id) {
  dx <- (id[['x','tar']] - id[['x','par']])/2L
  dy <- (id[['y','tar']] - id[['y','par']])/2L
  dx_p_dy <- dx + dy
  dx_m_dy <- dx - dy

  id[,'par'] <- lapply(id[,'tar'], rep, 2L)
  id[,'tar'] <- list(
      c(id[['x','tar']] - dx_p_dy, id[['x','tar']] - dx_m_dy),
      c(id[['y','tar']] + dx_m_dy, id[['y','tar']] - dx_p_dy)
  )
  id
}
seed_ids <- function(offset, m, length, hrz) {
  mids <- seq(m, length.out=length, by=m * 2L)
  a <- rep(mids, 2L)
  b <- rep(m + offset, 2L * length)
  c <- c(mids - m, mids + m)
  o <- rep(c(m, -m), length.out=length)
  d <- b + c(o, -o)

  matrix(
    if(hrz) list(a, b, c, d) else list(b, a, d, c),
    2L, 2L, dimnames=list(c('x', 'y'), c('tar', 'par'))
  )
}
