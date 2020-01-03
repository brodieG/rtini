<!-- README.md is generated from README.Rmd. Please edit that file

library(rmarkdown)
render('README.Rmd', output_format=md_document())

-->
rtini - Right Triangulated Irregular Network Surface Approximation
==================================================================

Overview
--------

A vectorized adaption of [Vladimir Agafonkin’s
MARTINI](https://observablehq.com/@mourner/martin-real-time-rtin-terrain-mesh)
implementation of the Evans etal. [Right-Triangulated Irregular Networks
(RTIN)](https://www.cs.ubc.ca/~will/papers/rtin.pdf) surface
approximation algorithm. This implementation supports a broader range of
dimensions and should also be faster in many use cases[1].

The workhorse functions are `rtini_error()` and `rtini_mesh()`. The
former computes the surface approximation error at all levels of
approximation, and the latter extract mesh triangles subject to an error
tolerance.

This is package is experimental, lightly tested, and not intended for
production use. It is unlikely this package will ever be posted on CRAN,
but if you are interested in the functionality and know of a GPL
licensed home for it do let me know.

    library(rtini)

    err <- rtini_error(volcano)

    library(rgl)
    par3d(windowRect=c(100,100,800,400))
    mfrow3d(nr=1, nc=3, sharedMouse=TRUE)
    for(tol in c(2, 10, 30)){
      next3d()
      ids <- rtini_extract(err, tol)
      xyz <- id_to_coord(ids, volcano)
      triangles3d(xyz, color='grey50')
    }

![](extra/figures/README-rgl-shot.png)

This implementation works for any elevation map with odd number of rows
and columns, but works best when the dimensions are (2^k + 1, n \* 2^k +
1) as otherwise many tiles must be split just to get them to conform
with each other, as can be seen here with `volcano` along the bottom and
right edges.

Installation
------------

This package is only available on github. Use your favorite package
installer or:

    f.dl <- tempfile()
    f.uz <- tempfile()
    github.url <- 'https://github.com/brodieG/rtini/archive/master.zip'
    download.file(github.url, f.dl)
    unzip(f.dl, exdir=f.uz)
    install.packages(file.path(f.uz, 'rtini-master'), repos=NULL, type='source')
    unlink(c(f.dl, f.uz))

Acknowledgments
---------------

-   Vladimir Agafonkin for the implementation of and [wonderful
    post](https://observablehq.com/@mourner/martin-real-time-rtin-terrain-mesh)
    on the topic.
-   William Evans, David Kirkpatrick, and Gregg Townsend for the
    [original paper](https://www.cs.ubc.ca/~will/papers/rtin.pdf).
-   Daniel Adler, Duncan Murdoch, etal. for `rgl` with which we made the
    image in this README.
-   R Core for developing and maintaining such a wonderful language.

[1] Error calculation should be generally faster, and mesh extraction
should be comparable except when the reduced mesh has very large numbers
of triangles, in which case the JavaScript implementation is faster.
