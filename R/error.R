#' Compute RTIN Surface Approximation Error on a Height Map
#'
#' A vectorized adaptation of \href{https://observablehq.com/@mourner/martin-real-time-rtin-terrain-mesh}{Vladimir Agafonkin's MARTINI implementation}
#' of the \href{https://www.cs.ubc.ca/~will/papers/rtin.pdf}{Right Triangulated Irregular Networks} surface approximation
#' algorithm.  Accepts a broader set of input matrix dimensions than MARTINI.
#'
#' @param map numeric matrix with at least three rows and three columns, and
#'   with an odd number of each.  It is expected that the input is finite and
#'   does not contain NAs.
#' @return numeric matrix of the same dimension as `map` with the
#'   maximum of the rtin approximation error at each point and the equivalent
#'   approximation at the child points.
#' @seealso [rtin_extract()]
#' @export
#' @examples
#' err <- rtini_error(volcano)
#' rtin_extract(err, tol=10)

rtini_error <- function(map) {
  vetr(matrix(numeric(), 0, 0) && dim(.) > 2 && dim(.) %% 2 > 0)

  # - Helper Funs --------------------------------------------------------------
  .get_child_err <- function(ids.mid, which, type) {
    if(identical(type, 'axis')) {          # square sides
      col.off <- c(-1L, 1L, 1L, -1L)
      row.off <- c(-1L, -1L, 1L, 1L)
    } else if (identical(type, 'diag')) {  # diagonals
      col.off <- c(0L, 2L, 0L, -2L)
      row.off <- c(-2L, 0L, 2L, 0L)
    } else stop("bad input")
    offset <- (row.off * mult) %/% 4L + nr * (col.off * mult) %/% 4L
    lapply(
      seq_along(offset)[which], function(i) errors[ids.mid + offset[i]]
  ) }
  .get_par_err <- function(ids, offsets)
    abs(map[ids] - (map[ids + offsets[1L]] + map[ids + offsets[2L]]) / 2)

  nr <- nrow(map)
  nc <- ncol(map)
  layers <- floor(min(log2(c(nr, nc) - 1L)))
  errors <- array(0, dim=dim(map))

  # Need to add errors at the smallest seam between differnt size adjoining
  # tiles to make sure the larger ones are broken up.  What's here currently is
  # too aggressive

  if((r.extra <- (nr - 1L) %% 2L^layers)) {
    while(r.extra - 2^(floor(log2(r.extra))))
      r.extra <- r.extra - 2^(floor(log2(r.extra)))
    errors[nr - r.extra, seq(r.extra + 1L, nc - 1L, by=r.extra)] <- Inf
  }
  if((c.extra <- (nc - 1L) %% 2L^layers)) {
    while(c.extra - 2^(floor(log2(c.extra))))
      c.extra <- c.extra - 2^(floor(log2(c.extra)))
    errors[seq(c.extra + 1L, nr - 1L, by=c.extra), nc - c.extra] <- Inf
  }

  for(i in seq_len(layers)) {
    mult <- as.integer(2^i)
    mhalf <- mult %/% 2L
    tile.r <- ((nr - 1L) %/% mult)
    tile.c <- ((nc - 1L) %/% mult)

    # - Axis (vertical/horizontal) ---------------------------------------------

    # v in var names designates vertical, h horizontal
    col.vn <- tile.c + 1L
    col.hn <- tile.c
    col.v <- seq(1L, length.out=tile.c + 1L, by=mult)
    col.h <- seq(2L, length.out=tile.c, by=mult)
    ids.v <- seq(mhalf + 1L, tile.r * mult, mult)
    ids.h <- seq(nr * mhalf + 1L, nr * mhalf + 1L + tile.r * mult, by=mult)
    v.len <- length(ids.v)
    h.len <- length(ids.h)

    # Non-square cases where there is no outside right/bot
    v.ex <- tile.c * mult + mhalf + 1 > nc
    h.ex <- tile.r * mult + mhalf + 1 > nr

    # Compute IDs, distinguish b/w inside vs. outside
    ids <- list(
      v.out.l=ids.v,
      v.in=ids.v + rep_each((col.v[-c(1L,col.vn * v.ex)] - 1L) * nr, v.len),
      v.out.r=if(v.ex) ids.v + (col.v[col.vn] - 1L) * nr,
      h.out.t=ids.h[1L] + (seq_len(col.hn) - 1L) * nr * mult,
      h.in=ids.h[-c(1L, h.len * h.ex)] +
        rep_each((seq_len(col.hn)-1L) * nr * mult, length(ids.h) - 1L - h.ex),
      h.out.b=if(h.ex) ids.h[h.len] + (seq_len(col.hn) - 1L) * nr * mult
    )
    # Errors
    err.child <- if(i > 1L) {
      which <- list(2:3, 1:4, c(1L,4L), 3:4, 1:4, 1:2)
      Map(.get_child_err, ids, which, 'axis')
    }
    err.par <- c(
      Map(.get_par_err, ids[1:3], list(c(-mhalf, mhalf))),
      Map(.get_par_err, ids[4:6], list(c(-mhalf * nr, mhalf * nr)))
    )
    for(j in seq_along(ids)) {
      errors[ids[[j]]] <- do.call(
        pmax, c(list(na.rm=TRUE, errors[ids[[j]]]), err.par[j], err.child[[j]])
    ) }
    # - Diagonals --------------------------------------------------------------

    ids.raw <- seq(mhalf + nr * mhalf + 1L, length.out=tile.r, by=mult) +
      matrix((seq_len(tile.c) - 1L) * mult * nr, tile.r, tile.c, byrow=TRUE)
    which.sw <- xor((col(ids.raw) %% 2L), (row(ids.raw) %% 2L))
    ids <- list(nw=ids.raw[!which.sw], sw=ids.raw[which.sw])

    err.child <- Map(.get_child_err, ids, list(1:4), 'diag')
    err.par <- Map(
      .get_par_err, ids,
      list((mhalf * nr + mhalf) * c(1L, -1L), (mhalf * nr - mhalf) * c(1L, -1L))
    )
    for(j in seq_along(ids)) {
      errors[ids[[j]]] <- do.call(
        pmax, c(list(na.rm=TRUE, errors[ids[[j]]]), err.par[j], err.child[[j]])
    ) }
  }
  errors
}


rep_each <- function(x, each) rep(x, rep(each, length(x)))
