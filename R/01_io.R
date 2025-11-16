#' Read occurrence points of an invasive species
#'
#' @description
#' Read longitude and latitude coordinates from a CSV file and convert them
#' into an sf POINT object. This function is typically used to import
#' invasive species occurrence records for spatial analysis.
#'
#' @param file_path Character string. Path to the CSV file containing
#'   occurrence records (e.g., "data/occ.csv").
#' @param lon_col Character string. Name of the longitude column.
#'   Default: "lon".
#' @param lat_col Character string. Name of the latitude column.
#'   Default: "lat".
#' @param crs_string Character string specifying the coordinate reference
#'   system (CRS) of the coordinates. Default: "EPSG:4326".
#'
#' @return An sf POINT object representing the occurrence locations.
#' @export
read_occurrence_points <- function(file_path,
                                   lon_col = "lon",
                                   lat_col = "lat",
                                   crs_string = "EPSG:4326") {
  df <- readr::read_csv(file_path, show_col_types = FALSE)

  if (!all(c(lon_col, lat_col) %in% names(df))) {
    stop("Longitude or latitude column not found. Please check lon_col / lat_col.")
  }

  sf::st_as_sf(
    df,
    coords = c(lon_col, lat_col),
    crs = crs_string
  )
}


#' Read a habitat–distance frequency table (e.g., a.csv)
#'
#' @description
#' Read a formatted CSV file containing distance–frequency information for
#' sampling points and occurrence points with respect to multiple habitat
#' types. The function renames verbose column names in the input file to the
#' concise internal naming convention used by the WAN package.
#'
#' The input file must contain:
#'   - a distance column
#'   - paired frequency columns for sampling and occurrence for each habitat
#'
#' These will be renamed into:
#'   distance, sampling_<habitat>, occurrence_<habitat>
#'
#' @param file_path Character string. Path to the CSV file (e.g., "a.csv").
#'
#' @return A data.frame containing standardized columns:
#'   distance,
#'   sampling_roads, occurrence_roads,
#'   sampling_sand, occurrence_sand,
#'   sampling_rivers, occurrence_rivers,
#'   sampling_stations, occurrence_stations.
#'
#' @export
read_habitat_data <- function(file_path) {
  data <- readr::read_csv(file_path, show_col_types = FALSE)
  message("Loaded data dimensions: ", paste(dim(data), collapse = " x "))

  data <- data %>%
    dplyr::rename(
      distance =
        "Distance from frequently distributed habitats (km)",

      sampling_roads =
        "Frequency of sampling points with the increase of the distance from roads  (%)",
      occurrence_roads =
        "Frequency of occurrence points with the increase of the distance from roads (%)",

      sampling_sand =
        "Frequency of sampling points with the increase of the distance from construction sand distribution centers (%)",
      occurrence_sand =
        "Frequency of occurrence points with the increase of the distance from construction sand distribution centers (%)",

      sampling_rivers =
        "Frequency of sampling points with the increase of the distance from rivers  (%)",
      occurrence_rivers =
        "Frequency of occurrence points with the increase of the distance from rivers(%)",

      sampling_stations =
        "Frequency of sampling points with the increase of the distance from petrol service stations (%)",
      occurrence_stations =
        "Frequency of occurrence points with the increase of the distance from petrol service stations (%)"
    )

  data
}
