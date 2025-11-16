#' Compute distance from points to one habitat feature
#'
#' @param pts_sf sf POINT.
#' @param habitat_sf sf object (POINT / LINESTRING / POLYGON, etc.).
#'
#' @return Numeric vector of minimum distances for each point.
#' @export
compute_distance_to_feature <- function(pts_sf, habitat_sf) {
  pts_sf     <- sf::st_as_sf(pts_sf)
  habitat_sf <- sf::st_as_sf(habitat_sf)

  # Ensure both are in the same CRS
  habitat_sf <- sf::st_transform(habitat_sf, sf::st_crs(pts_sf))

  # st_distance returns a matrix (units), rows = pts, cols = habitat features
  dmat <- sf::st_distance(pts_sf, habitat_sf)

  # Take the minimum distance to any habitat feature for each point
  dmin <- apply(dmat, 1, min)

  as.numeric(dmin)
}

#' Compute distances from occurrence and background points to multiple habitats
#'
#' @description
#' Given occurrence points, background points, and a named list of habitat
#' features, compute the minimum distance from each point to each habitat.
#'
#' @param occ_sf sf POINT, occurrence locations.
#' @param bg_sf sf POINT, background locations.
#' @param habitats Named list of sf objects, e.g.
#'   \code{list(roads = roads_sf, rivers = rivers_sf)}.
#'
#' @return list with two data.frames:
#'   \itemize{
#'     \item \code{occ}: columns are habitats, rows are occurrence points;
#'     \item \code{bg}:  columns are habitats, rows are background points.
#'   }
#' @export
compute_all_distances_generic <- function(occ_sf, bg_sf, habitats) {
  occ_sf <- sf::st_as_sf(occ_sf)
  bg_sf  <- sf::st_as_sf(bg_sf)

  if (is.null(names(habitats)) || any(names(habitats) == "")) {
    stop("Argument 'habitats' must be a named list, e.g. list(roads = roads_sf, rivers = rivers_sf).")
  }

  occ_list <- list()
  bg_list  <- list()

  for (hname in names(habitats)) {
    hab <- habitats[[hname]]
    hab_sf <- sf::st_as_sf(hab)
    hab_sf <- sf::st_transform(hab_sf, sf::st_crs(occ_sf))

    occ_list[[hname]] <- compute_distance_to_feature(occ_sf, hab_sf)
    bg_list[[hname]]  <- compute_distance_to_feature(bg_sf,  hab_sf)
  }

  list(
    occ = as.data.frame(occ_list),
    bg  = as.data.frame(bg_list)
  )
}
