#' @name is_integerish
#' @title is a number convertible to integer
#' @description is a number convertible to integer
#' @param x input
#' @details from hadley's assertthat
#' @return (flag)
#'
is_integerish = function(x){
  res = is.integer(x) || (is.numeric(x) && all(x == trunc(x)) &&
      !is.na(x))
  return(res)
  }

#' @name gh_is_valid
#' @title Is a string is a valid geohash
#' @description Vectorized check if a string is a valid geohash
#' @param ghs (character vector) Input strings
#' @return (Logical vector) with same length as input
#' @examples
#' to_test = c("tdxw1", "notgeohash")
#' gh_is_valid(to_test)
#' @export
#'
gh_is_valid = function(ghs){

  stopifnot(is.character(ghs))
  valid = strsplit("0123456789bcdefghjkmnpqrstuvwxyz", "")[[1]]

  res = vapply(ghs
               , function(x) all(strsplit(x, "")[[1]] %in% valid)
               , logical(1L)
               , USE.NAMES = FALSE
               )

  return(res)
}

#' @name gh_expand
#' @title Get one level deeper geohashes inside a geohash
#' @description Vectorized listing of all geohashes of one level lower
#'   resolution
#' @param ghs (character vector) Input geohashes
#' @return (List) of character vectors containing one level deeper geohashes
#' @examples
#' to_test = c("tdxw1", "tdr1u")
#' gh_expand(to_test)
#' @export
#'
gh_expand = function(ghs){

  stopifnot(gh_is_valid(ghs))
  valid = strsplit("0123456789bcdefghjkmnpqrstuvwxyz", "")[[1]]
  res   = lapply(ghs, function(x) paste0(x, valid))
  names(res) = ghs

  return(res)

}

#' @name gh_has_precision
#' @title Check precision of geohashes
#' @description Vectorized checking of precision of geohashes
#' @param ghs (character vector) Input geohashes
#' @param precision (integer vector) Either a single precision value or a vector
#'   of length equal to geohashes
#' @return TRUE or FALSE When precision has length 1. Else, a logical vector of
#'   length equal to number of geohashes
#' @examples
#' gh_has_precision(c("tdxw1", "tdr1"), 5)
#' #gh_has_precision(c("tdxw1", "tdr1"), precision = c(5, 4))
#' @export
#'
gh_has_precision = function(ghs, precision){

  stopifnot(gh_is_valid(ghs))
  stopifnot(is_integerish(precision))
  len_precision = length(precision)
  stopifnot(len_precision == 1L || len_precision == length(ghs))

  res = (nchar(ghs) == precision)
  return(res)
}

#' @name gh_pairwise_distances
#' @title Compute pairwise distance between a pair of geohash vectors
#' @description If one of the vector is a singleton and the other has length >
#'   1, then the singleton vector is extended to that length of the longer one
#' @param gh_1 (character vector) first geohash vector
#' @param gh_2 (character vector) second geohash vector
#' @param ... arguments to be passed to geodist::geodist()
#' @return (numeric vector) Vector of pairwise distances with length equal to
#'   the longest input vector
#' @examples
#' gh_pairwise_distances(c("tdrw1", "tdrw2"), c("tdrw3", "tdrw4"))
#' @export
#'
gh_pairwise_distances = function(gh_1, gh_2, ...){

  stopifnot(gh_is_valid(gh_1))
  stopifnot(gh_is_valid(gh_2))

  len_1 = length(gh_1)
  len_2 = length(gh_2)
  # extend the singleton to the size of the larger one
  if (len_1 != len_2) {

      stopifnot(min(len_1, len_2) == 1)
      if (which.min(c(len_1, len_2)) == 1){
          gh_1 = rep(gh_1, max(len_1, len_2))
      } else {
          gh_2 = rep(gh_2, max(len_1, len_2))
      }
  }

  mat_1 = do.call(cbind, geohashTools::gh_decode(gh_1))
  mat_2 = do.call(cbind, geohashTools::gh_decode(gh_2))

  res = geodist::geodist(mat_1, mat_2, paired = TRUE, sequential = FALSE, ...)
  return(res)
}

#' @name gh_distances
#' @title Compute distances between a pair of geohash vectors
#' @description If one of the vector is a singleton and the other has length >
#'   1, then the singleton vector is extended to that length of the longer one
#' @param gh_1 (character vector) first geohash vector
#' @param gh_2 (character vector) second geohash vector
#' @param ... arguments to be passed to geodist::geodist()
#' @return (numeric vector) matrix of cartesian product of distances with
#'   dimension of length(gh_1) X length(gh_2)
#' @examples
#' gh_distances(c("tdrw1", "tdrw2"), c("tdrw3", "tdrw4"))
#' @export
#'
gh_distances = function(gh_1, gh_2, ...){

  stopifnot(gh_is_valid(gh_1))
  stopifnot(gh_is_valid(gh_2))

  mat_1 = do.call(cbind, geohashTools::gh_decode(gh_1))
  mat_2 = do.call(cbind, geohashTools::gh_decode(gh_2))

  res = geodist::geodist(mat_1
                         , mat_2
                         , paired = FALSE
                         , sequential = FALSE
                         ,  ...
                         )
  dimnames(res) = list(gh_1, gh_2)
  return(res)
}

#' @name gh_dist
#' @title Compute dist object for a set of geohashes
#' @description set of geohashes is expected to be unique
#' @param ghs (character vector) geohash vector
#' @param ... arguments to be passed to geodist::geodist()
#' @return dist object
#' @examples
#' gh_dist(c("tdrw1", "tdrw2", "tdrw3"))
#' @export
#'
gh_dist = function(ghs, ...){

  stopifnot(length(ghs) == length(unique(ghs)))

  geo_mat = do.call(cbind, geohashTools::gh_decode(ghs))
  res     = geodist::geodist(geo_mat, ...)
  dimnames(res) = list(ghs, ghs)
  res = stats::as.dist(res)

  return(res)
}

#' @name gh_cover_
#' @title Get a set of geohashes covering a sf object
#' @description Get a set of geohashes covering a sf object
#' @param x (sf object)
#' @param precision (integer) precision on the covering geohashes
#' @return (character vector) of covering geohashes
#'
gh_cover_ = function(x, precision = NULL){

  stopifnot(inherits(x, "sf"))

  # get all geohashes in the bounding box
  bb = sf::st_bbox(x)
  delta = 2 * geohashTools::gh_delta(precision)
  grid = expand.grid(
    lat = seq(bb['ymin'], bb['ymax'] + delta[1L], by = delta[1L]),
    lng = seq(bb['xmin'], bb['xmax'] + delta[2L], by = delta[2L])
    )
  ghs = geohashTools::gh_encode(grid[["lat"]], grid[["lng"]], precision)
  ghs = unique(ghs)

  # filter the geohashes which intersect with the sf object
  ghs_sf = geohashTools::gh_to_sf(ghs)
  x = sf::st_transform(x, crs = 4326)
  suppressMessages({
    flags  = (rowSums(sf::st_intersects(ghs_sf, x, sparse = FALSE)) > 0)
  })
  ghs_in = ghs[flags]

  return(ghs_in)
}

#' @name gh_cover
#' @title Get a set of geohashes covering a set of object
#' @description Get a set of geohashes covering a set of object
#' @param sfs (list of sf objects)
#' @param precision (integer) precision on the covering geohashes
#' @return (list) of covering geohashes
#' @examples
#' fname = system.file("shape/nc.shp", package="sf")
#' temp = sf::st_read(fname)
#' temp
#' gu = geohashTools::gh_to_sf(gh_cover(temp, 4)[[1]])
#' mapview::mapview(list(gu, temp))
#' @export
#'
gh_cover = function(sfs, precision = NULL){

  if (inherits(sfs, "sf")) {
    sfs = list(sfs)
  }
  stopifnot(is.list(sfs))
  stopifnot(all(sapply(sfs, function(x) inherits(x, "sf"))))
  stopifnot(is_integerish(precision) && length(precision) == 1L)

  res = lapply(sfs, gh_cover_, precision)
  return(res)
}

#' @name gh_pad_
#' @title Pad a geohash with neighboring one-level smaller geohashes
#' @description Pad a geohash with neighboring one-level smaller geohash
#' @param gh (character vector of length 1) geohash
#' @param self (flag) Whether the input geohash should be included in the output
#' @return (character vector) of padding geohashes
#'
gh_pad_ = function(gh, self = TRUE){

  neighbors          = unlist(geohashTools::gh_neighbours(gh, self = FALSE))
  neighbors_small    = unlist(gh_expand(neighbors), use.names = FALSE)
  neighbors_small_sf = lapply(neighbors_small, geohashTools::gh_to_sf)
  gh_sf              = geohashTools::gh_to_sf(gh)

  # filter the ones that touch the geohash
  suppressMessages({
    touch_flag = vapply(neighbors_small_sf
                        , function(x){
                            sf::st_touches(gh_sf
                                           , x
                                           , sparse = FALSE
                                           )[, 1][1]
                        }
                        , logical(1)
                        , USE.NAMES = FALSE
                        )
  })
  res = neighbors_small[touch_flag]
  if (self) { res = c(gh, res) }
  return(res)
}

#' @name gh_pad
#' @title Pad geohashes with neighboring one-level smaller geohashes
#' @description Pad geohashes with neighboring one-level smaller geohash
#' @param ghs (character vector) geohashes
#' @param self (flag) Whether the input geohash should be included in the output
#' @return (list of character vectors) of padding geohashes
#' @examples
#' gh_pad(c("tdr1w", "tdr1q"))
#' @export
#'
gh_pad = function(ghs, self = TRUE){

  stopifnot(all(gh_is_valid(ghs)))
  stopifnot(is.logical(self) && length(self) == 1L)
  res = lapply(ghs, gh_pad_)
  names(res) = ghs

  return(res)
}

#' @name gh_envelope_bbox
#' @title Single bounding box for a set of geohashes
#' @description as sf bbox object
#' @param ghs (character vector) of geohashes
#' @param crs one of (i) character: a string accepted by GDAL, (ii) integer, a
#'   valid EPSG value (numeric), or (iii) an object of class crs
#' @return object of class 'bbox'
#' @examples
#' gh_envelope_bbox(c("tdr1w", "tdr1q"))
#' @export
#'
gh_envelope_bbox = function(ghs, crs = 4236){

  stopifnot(gh_is_valid(ghs))

  ne = simplify2array(geohashTools::gh_decode(ghs, coord_loc = "ne"))
  sw = simplify2array(geohashTools::gh_decode(ghs, coord_loc = "sw"))

  x_max = max(ne[, 2])
  x_min = min(sw[, 2])
  y_max = max(ne[, 1])
  y_min = min(sw[, 1])

  bbox = sf::st_bbox(c(xmin = x_min, xmax = x_max, ymax = y_max, ymin = y_min)
                     , crs = sf::st_crs(crs)
                     )
  return(bbox)
}

#' @name gh_bbox
#' @title Get bounding boxes
#' @description Vectorized function to obtain sf bbox class object for a set of
#'   geohashes
#' @param ghs (character vector) of geohashes
#' @param crs one of (i) character: a string accepted by GDAL, (ii) integer, a
#'   valid EPSG value (numeric), or (iii) an object of class crs
#' @return list of 'bbox'es
#' @examples
#' gh_bbox(c("tdr1w", "tdr1q"))
#' @export
#'
gh_bbox = function(ghs, crs = 4236){

  stopifnot(gh_is_valid(ghs))
  stopifnot(is_integerish(crs) & length(crs) == 1)

  bounds =
    t(cbind(simplify2array(geohashTools::gh_decode(ghs, coord_loc = "ne")),
            simplify2array(geohashTools::gh_decode(ghs, coord_loc = "sw"))
            )
      )
  rownames(bounds) = c("ymax", "xmax", "ymin", "xmin")
  crs_object = sf::st_crs(crs)

  res = lapply(
    1:ncol(bounds)
    , function(bs){
        do.call(sf::st_bbox, list(obj = bounds[, bs], crs = crs_object))
      }
    )
  names(res) = ghs

  return(res)
}

#' @name gh_is_neighbor_
#' @title check neighborhood relationship btween a geohash pair
#' @description check neighborhood relationship btween a geohash pair
#' @param gh_1 first geohash
#' @param gh_2 second geohash
#' @return (flag)
#'
gh_is_neighbor_ = function(gh_1, gh_2){

  gh_1_sf = geohashTools::gh_to_sf(gh_1)
  gh_2_sf = geohashTools::gh_to_sf(gh_2)

  touches  = sf::st_touches(gh_1_sf, gh_2_sf, sparse = FALSE)[1,1]
  contains = sf::st_contains(gh_1_sf, gh_2_sf, sparse = FALSE)[1,1]

  res = touches && !contains
  return(res)
}

#' @name gh_is_neighbor
#' @title check neighborhood relationship between geohash vectors
#' @description vectorized over geohash vectors
#' @param gh_1 first geohash
#' @param gh_2 second geohash
#' @return (logical vector) with length equal to longer of the two geohash
#'   vectors
#' @examples
#' gh_is_neighbor(c("tdrw1", "tdrw5"), c("tdrw2","tdrw3"))
#' @export
#'
gh_is_neighbor = function(gh_1, gh_2){

  stopifnot(all(gh_is_valid(gh_1)))
  stopifnot(all(gh_is_valid(gh_2)))

  suppressMessages({
    res = Vectorize(gh_is_neighbor_, USE.NAMES = FALSE)(gh_1, gh_2)
  })

  return(res)
}

#' @name gh_is_interior_
#' @title check whether gh_2 is in the proper interior of gh_1
#' @description check whether gh_2 is in the proper interior of gh_1
#' @param gh_1 first geohash
#' @param gh_2 second geohash
#' @return (flag)
#'
gh_is_interior_ = function(gh_1, gh_2){

  gh_1_sf = geohashTools::gh_to_sf(gh_1)
  gh_2_sf = geohashTools::gh_to_sf(gh_2)

  res = sf::st_contains_properly(gh_1_sf, gh_2_sf, sparse = FALSE)[1,1]
  return(res)
}

#' @name gh_is_interior
#' @title check whether gh_2 is in the proper interior of gh_1
#' @description vectorized check whether gh_2 is in the proper interior of gh_1
#' @param gh_1 first geohash vector
#' @param gh_2 second geohash vector
#' @return (logical vector) with length equal to longer of the two geohash
#'   vectors
#' @examples
#' gh_is_interior("tdr1w", c("tdr1wd", "tdr1w1"))
#' mapview::mapview(geohashTools::gh_to_sf(c("tdr1w", "tdr1wd", "tdr1w1")))
#' @export
#'
gh_is_interior = function(gh_1, gh_2){

  stopifnot(all(gh_is_valid(gh_1)))
  stopifnot(all(gh_is_valid(gh_2)))

  suppressMessages({
    res = Vectorize(gh_is_interior_, USE.NAMES = FALSE)(gh_1, gh_2)
  })
  return(res)
}

#' @name gh_hull
#' @title convex hull of geohashes
#' @description Obtain sf POINT geometry of convex hull covering the nodes of a set of geohashes
#' @param ghs (character vector) of geohashes
#' @return Object of class "sfc_POINT"
#' @examples
#' temp = c("tdr1w", "tdr3qt", "tdr1qw")
#' hull_sfc = gh_hull(temp)
#' mapview::mapview(list(hull_sfc, geohashTools::gh_to_sf(temp)))
#' @export
#'
gh_hull = function(ghs){

  stopifnot(all(gh_is_valid(ghs)))

  geom_layer = geohashTools::gh_to_sf(ghs)$geometry
  points = sf::st_union(sf::st_cast(geom_layer, to = "MULTIPOINT"))
  res = sf::st_convex_hull(points)

  return(res)
}
