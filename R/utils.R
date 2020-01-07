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

#' Convert Linearized Coordinates Into x/y/z
#'
#' [rtini_extract()] produces "linearized" coordinates that are useful for
#' directly indexing back into the elevation map, but can be clumsy for
#' plotting.  "Linearized" coordinates are offsets that can be used to access a
#' matrix or array as if it were a dimensionless vector.  This function
#' converts those coordinates to x/y/z values.
#'
#' @export
#' @param mesh list of numeric matrices as produced by [rtini_extract()].
#' @param map numeric matrix elevation map used to produce the mesh.
#' @param normalize FALSE unimplemented
#' @return list with three components, one for each dimension.  Each component
#'   is a three row matrix in which each column represents one triangle.
#' @examples
#' err <- rtini_error(volcano)
#' mesh <- rtini_extract(err, 10)
#' coord <- id_to_coord(mesh, volcano)

id_to_coord <- function(mesh, map, normalize=FALSE) {
  if(!isFALSE(normalize)) stop('`normalize` not supported yet')

  tris <- do.call(cbind, mesh)
  list(
    x=((tris - 1) %/% nrow(map)),
    y=((tris - 1) %% nrow(map)),
    z=map[unlist(mesh)]
  )
}

