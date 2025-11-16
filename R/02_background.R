#' Generate background points for invasive species analysis
#'
#' @description
#' Generate background points representing available environment.
#'
#' Logic:
#' \itemize{
#'   \item If \code{study_area} is provided:
#'     \itemize{
#'       \item method = "random": sample uniformly inside study_area;
#'       \item method = "restricted": sample inside a buffer around occurrences
#'             (optionally intersected with study_area).
#'     }
#'   \item If \code{study_area} is NULL:
#'     \itemize{
#'       \item method = "random": sample inside an expanded bounding box of occurrences;
#'       \item method = "restricted": sample inside a buffer around occurrences.
#'     }
#' }
#'
#' @param occ_sf sf POINT. Invasive species occurrence locations.
#' @param study_area Optional polygon (sf or terra SpatVector). If NULL, a
#'   bounding box around occurrences is used when method = "random".
#' @param n_bg Number of background points.
#' @param method One of \code{c("random","restricted","target_group")}.
#' @param buffer_dist Buffer distance around occurrences (same units as CRS),
#'   required when \code{method = "restricted"}.
#' @param bbox_expand Proportion to expand the occurrence bounding box in each
#'   direction when \code{study_area = NULL} and \code{method = "random"}.
#' @param target_points Optional sf POINT for \code{method = "target_group"}.
#'
#' @return sf POINT with background locations.
#' @export
generate_background_points <- function(occ_sf,
                                       study_area = NULL,
                                       n_bg = 10000,
                                       method = c("random", "restricted", "target_group"),
                                       buffer_dist = NULL,
                                       bbox_expand = 0.1,
                                       target_points = NULL) {
  method <- match.arg(method)

  # 1. Ensure CRS consistency ----------------------------------------------
  occ_sf <- sf::st_as_sf(occ_sf)

  if (!is.null(study_area)) {
    study_area <- sf::st_as_sf(study_area)
    study_area <- sf::st_transform(study_area, sf::st_crs(occ_sf))
  }

  # Convert occurrences to terra for geometry operations and sampling
  occ_v <- terra::vect(occ_sf)

  # 2. Determine sampling area (area_v) ------------------------------------
  area_v <- NULL

  if (!is.null(study_area)) {
    area_v <- terra::vect(study_area)
  } else if (method == "random") {
    # Use an expanded bounding box around occurrences
    ext <- terra::ext(occ_v)
    dx  <- ext$xmax - ext$xmin
    dy  <- ext$ymax - ext$ymin

    ext_expanded <- terra::ext(
      ext$xmin - dx * bbox_expand,
      ext$xmax + dx * bbox_expand,
      ext$ymin - dy * bbox_expand,
      ext$ymax + dy * bbox_expand
    )

    area_v <- terra::as.polygons(ext_expanded)
  }

  # 3. Generate background points by method --------------------------------
  if (method == "random") {
    if (is.null(area_v)) {
      stop("method = 'random' but no sampling area could be determined.")
    }

    # IMPORTANT: no 'as.points' argument here (terra no longer supports it)
    bg_v <- terra::spatSample(
      x      = area_v,
      size   = n_bg,
      method = "random"
    )

  } else if (method == "restricted") {
    if (is.null(buffer_dist)) {
      stop("method = 'restricted' requires 'buffer_dist'.")
    }

    buf_v <- terra::buffer(occ_v, width = buffer_dist)

    if (!is.null(area_v)) {
      buf_v <- terra::intersect(buf_v, area_v)
    }

    bg_v <- terra::spatSample(
      x      = buf_v,
      size   = n_bg,
      method = "random"
    )

  } else if (method == "target_group") {
    if (is.null(target_points)) {
      stop("method = 'target_group' requires 'target_points'.")
    }

    tg_sf <- sf::st_as_sf(target_points)
    tg_v  <- terra::vect(tg_sf)

    if (nrow(tg_v) > n_bg) {
      idx  <- sample(seq_len(nrow(tg_v)), n_bg)
      bg_v <- tg_v[idx, ]
    } else {
      bg_v <- tg_v
    }
  }

  # 4. Convert back to sf ---------------------------------------------------
  bg_sf <- sf::st_as_sf(bg_v)
  sf::st_set_crs(bg_sf, sf::st_crs(occ_sf))
}
