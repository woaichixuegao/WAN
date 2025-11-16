#' One-click invasive species habitat preference analysis (from coordinates)
#'
#' @description
#' Starting from occurrence point coordinates (longitude/latitude), this
#' function:
#' \enumerate{
#'   \item reads occurrence data from a CSV file,
#'   \item generates background points,
#'   \item computes distances from points to multiple habitat layers,
#'   \item builds distance–frequency tables.
#' }
#'
#' The resulting objects can be used for further analysis and plotting,
#' such as decay curves, trend fitting (LM/GAM), and statistical tests.
#'
#' If \code{breaks = NULL}, distance class breakpoints are generated
#' automatically based on the distances of both occurrence and background
#' points (controlled by \code{n_bins} and \code{method}).
#'
#' @param occ_file Character string. Path to the CSV file containing
#'   occurrence coordinates.
#' @param lon_col,lat_col Character strings. Column names for longitude
#'   and latitude in \code{occ_file}. Defaults are \code{"lon"} and \code{"lat"}.
#' @param habitats A named list of habitat layers, e.g.
#'   \code{list(roads = roads_sf, rivers = rivers_sf)}.
#'   Each element should be an \code{sf} or \code{SpatVector} object.
#' @param study_area Optional polygon (sf or terra SpatVector). If \code{NULL},
#'   the background sampling area is derived from the occurrence bounding box
#'   (optionally expanded via \code{bbox_expand}).
#' @param breaks Optional numeric vector of distance breakpoints. If \code{NULL},
#'   breaks are generated automatically.
#' @param n_bins Integer. Desired number of bins when \code{breaks = NULL}.
#' @param method Character string specifying the method used to generate distance
#'   breakpoints when \code{breaks = NULL}:
#'   \code{"pretty"}, \code{"equal"}, or \code{"quantile"}.
#' @param n_bg Integer. Number of background points to generate.
#' @param bg_method Character string specifying the background sampling method:
#'   \code{"random"}, \code{"restricted"}, \code{"target_group"}, or \code{"bias"}.
#' @param buffer_dist Numeric. Buffer distance (in CRS units) used when
#'   \code{bg_method = "restricted"}.
#' @param bbox_expand Numeric. Expansion factor applied to the occurrence
#'   bounding box when \code{study_area = NULL} and \code{bg_method = "random"}.
#' @param target_points Optional sf POINT object used when
#'   \code{bg_method = "target_group"}.
#'
#' @return A list with components:
#' \itemize{
#'   \item \code{freq_long}: long-format distance–frequency table
#'         (\code{habitat, distance, sampling, occurrence}).
#'   \item \code{freq_wide}: wide-format table with columns
#'         \code{distance, sampling_<habitat>, occurrence_<habitat>}.
#'   \item \code{occ_sf}: sf POINT object of occurrence locations.
#'   \item \code{bg_sf}: sf POINT object of background locations.
#' }
#' @export
run_invasion_habitat_pipeline_generic <- function(occ_file,
                                                  lon_col = "lon",
                                                  lat_col = "lat",
                                                  habitats,
                                                  study_area = NULL,
                                                  breaks = NULL,
                                                  n_bins = 6,
                                                  method = c("pretty", "equal", "quantile"),
                                                  n_bg = 10000,
                                                  bg_method = c("random", "restricted", "target_group", "bias"),
                                                  buffer_dist = 50000,
                                                  bbox_expand = 0.1,
                                                  target_points = NULL) {
  method    <- match.arg(method)
  bg_method <- match.arg(bg_method)

  # 1. Read occurrence points ----------------------------------------------
  occ_sf <- read_occurrence_points(
    file_path  = occ_file,
    lon_col    = lon_col,
    lat_col    = lat_col,
    crs_string = "EPSG:4326"
  )

  # 2. Generate background points ------------------------------------------
  bg_sf <- generate_background_points(
    occ_sf        = occ_sf,
    study_area    = study_area,
    n_bg          = n_bg,
    method        = bg_method,
    buffer_dist   = buffer_dist,
    bbox_expand   = bbox_expand,
    target_points = target_points
  )

  # 3. Compute distances to each habitat -----------------------------------
  dist_list <- compute_all_distances_generic(
    occ_sf   = occ_sf,
    bg_sf    = bg_sf,
    habitats = habitats
  )

  # 4. Build frequency tables (auto or manual breaks) ----------------------
  freq_long <- build_frequency_table_generic_long(
    dist_occ = dist_list$occ,
    dist_bg  = dist_list$bg,
    breaks   = breaks,
    n_bins   = n_bins,
    method   = method
  )

  freq_wide <- build_frequency_table_generic_wide(freq_long)

  list(
    freq_long = freq_long,
    freq_wide = freq_wide,
    occ_sf    = occ_sf,
    bg_sf     = bg_sf
  )
}


#' Run habitat preference analysis from a frequency table
#'
#' @description
#' Given a pre-assembled wide-format frequency table (containing
#' \code{distance}, \code{sampling_*}, and \code{occurrence_*} columns),
#' this function:
#' \enumerate{
#'   \item reshapes it into long format,
#'   \item fits exponential decay curves for multiple habitats,
#'   \item produces a multi-panel decay plot,
#'   \item computes statistical summaries of habitat preference (decay-based),
#'   \item optionally fits LM/GAM trends and produces trend plots and summaries.
#' }
#'
#' @param freq_wide A wide-format data.frame with columns:
#'   \code{distance}, \code{sampling_<habitat>}, \code{occurrence_<habitat>}.
#' @param trend_methods Optional character vector specifying which trend
#'   models to fit for frequency–distance relationships:
#'   any of \code{"lm"} and/or \code{"gam"}.
#'   If \code{NULL} (default), no trend models are fitted.
#' @param lm_degree Integer. Polynomial degree for LM trends (default 1 = linear),
#'   used when \code{"lm"} is included in \code{trend_methods}.
#' @param gam_k Integer. Basis dimension for GAM smooth trends (default 5),
#'   used when \code{"gam"} is included in \code{trend_methods}.
#' @param make_trend_plots Logical. If \code{TRUE} and \code{trend_methods}
#'   is not \code{NULL}, multi-panel LM/GAM trend plots are produced.
#'
#' @return A list with components:
#' \itemize{
#'   \item \code{decay_plot}: patchwork object with multi-habitat exponential
#'         decay curves and confidence intervals.
#'   \item \code{decay_summary}: data.frame with decay-based statistics
#'         per habitat (e.g. decay rates, p-values, ecological interpretation).
#'   \item \code{trend_plot}: patchwork object with LM/GAM trends (or \code{NULL}
#'         if no trend plots were requested).
#'   \item \code{trend_summary}: data.frame with LM/GAM trend statistics
#'         per habitat (or \code{NULL} if no trend models were fitted).
#' }
#' @export
run_from_frequency_table <- function(freq_wide,
                                     trend_methods   = NULL,
                                     lm_degree       = 1,
                                     gam_k           = 5,
                                     make_trend_plots = TRUE) {
  # 1. Long-format for plotting --------------------------------------------
  freq_long <- freq_wide %>%
    tidyr::pivot_longer(
      cols      = -distance,
      names_to  = c(".value", "habitat"),
      names_sep = "_"
    )

  # 2. Exponential decay plot + summary ------------------------------------
  decay_plot    <- plot_multi_habitat_decay(freq_long)
  decay_summary <- generate_statistical_summary_multi(freq_wide)

  # 3. Optional LM/GAM trend models ----------------------------------------
  trend_plot    <- NULL
  trend_summary <- NULL

  if (!is.null(trend_methods)) {
    trend_methods <- intersect(trend_methods, c("lm", "gam"))

    if (length(trend_methods) > 0) {
      # Trend summary table (uses LM/GAM fitting and model statistics)
      trend_summary <- generate_trend_summary_multi(
        freq_wide = freq_wide,
        methods   = trend_methods,
        degree    = lm_degree,
        k         = gam_k
      )

      # Trend plots (LM/GAM curves + 95% CI ribbons)
      if (isTRUE(make_trend_plots)) {
        trend_plot <- plot_multi_habitat_trends(
          freq_long = freq_long,
          methods   = trend_methods,
          degree    = lm_degree,
          k         = gam_k,
          ncol      = 2
        )
      }
    }
  }

  list(
    decay_plot    = decay_plot,
    decay_summary = decay_summary,
    trend_plot    = trend_plot,
    trend_summary = trend_summary
  )
}


#' Run habitat preference analysis from a frequency table file
#'
#' @description
#' Wrapper function to run the habitat preference analysis directly from
#' a CSV frequency table file (e.g. an “a.csv”-style input).
#'
#' Internally, it:
#' \enumerate{
#'   \item reads the frequency table using \code{read_habitat_data()},
#'   \item extracts the relevant columns into a wide-format table,
#'   \item runs \code{run_from_frequency_table()},
#'   \item saves the multi-panel decay plot as an image,
#'   \item optionally saves the LM/GAM trend plot as an image,
#'   \item saves the decay-based statistical summary as a CSV file.
#' }
#'
#' @param file_path Character string. Path to the frequency table CSV file.
#' @param save_decay_plot_path Character string. File path for saving the
#'   exponential decay plot (default: \code{"four_panel_plot_simple.png"}).
#' @param save_trend_plot_path Optional character string. File path for saving
#'   the LM/GAM trend plot (e.g. \code{"trend_plot_lm_gam.png"}). If \code{NULL}
#'   (default), the trend plot is not saved.
#' @param save_summary_path Character string. File path for saving the
#'   decay-based statistical summary CSV (default:
#'   \code{"habitat_preference_summary.csv"}).
#' @param trend_methods Optional character vector specifying which trend
#'   models to fit and (optionally) plot; passed to
#'   \code{run_from_frequency_table()}.
#' @param lm_degree Integer. Polynomial degree for LM trends (default 1).
#' @param gam_k Integer. Basis dimension for GAM smooth trends (default 5).
#'
#' @return A list with components:
#' \itemize{
#'   \item \code{decay_plot}: multi-habitat exponential decay plot.
#'   \item \code{decay_summary}: decay-based statistical summary.
#'   \item \code{trend_plot}: LM/GAM trend plot (or \code{NULL}).
#'   \item \code{trend_summary}: LM/GAM trend summary (or \code{NULL}).
#' }
#' @export
run_from_frequency_file <- function(file_path,
                                    save_decay_plot_path = "four_panel_plot_simple.png",
                                    save_trend_plot_path = NULL,
                                    save_summary_path    = "habitat_preference_summary.csv",
                                    trend_methods        = NULL,
                                    lm_degree            = 1,
                                    gam_k                = 5) {
  # 1. Read frequency table -------------------------------------------------
  data <- read_habitat_data(file_path)

  # 2. Extract wide-format table with standard columns ----------------------
  freq_wide <- data %>%
    dplyr::select(
      distance,
      sampling_roads,    occurrence_roads,
      sampling_sand,     occurrence_sand,
      sampling_rivers,   occurrence_rivers,
      sampling_stations, occurrence_stations
    )

  # 3. Run full analysis (decay + optional trends) --------------------------
  res <- run_from_frequency_table(
    freq_wide        = freq_wide,
    trend_methods    = trend_methods,
    lm_degree        = lm_degree,
    gam_k            = gam_k,
    make_trend_plots = !is.null(save_trend_plot_path)
  )

  # 4. Save decay plot ------------------------------------------------------
  ggplot2::ggsave(
    filename = save_decay_plot_path,
    plot     = res$decay_plot,
    width    = 10,
    height   = 8,
    dpi      = 300,
    bg       = "white"
  )

  # 5. Optionally save trend plot -------------------------------------------
  if (!is.null(save_trend_plot_path) && !is.null(res$trend_plot)) {
    ggplot2::ggsave(
      filename = save_trend_plot_path,
      plot     = res$trend_plot,
      width    = 10,
      height   = 8,
      dpi      = 300,
      bg       = "white"
    )
  }

  # 6. Save decay-based statistical summary ---------------------------------
  utils::write.csv(
    res$decay_summary,
    save_summary_path,
    row.names    = FALSE,
    fileEncoding = "UTF-8"
  )

  message("Decay plot saved to: ", save_decay_plot_path)
  if (!is.null(save_trend_plot_path) && !is.null(res$trend_plot)) {
    message("Trend plot saved to: ", save_trend_plot_path)
  }
  message("Decay-based statistical summary saved to: ", save_summary_path)

  res
}
