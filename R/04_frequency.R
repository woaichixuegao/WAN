#' Automatically generate distance class breakpoints
#'
#' @description
#' Automatically generate a set of meaningful breakpoints for grouping
#' distance values of occurrence points and background points.
#' The function merges all distance values (occurrence + background),
#' determines the global minimum and maximum, and computes breakpoints
#' using one of the following methods:
#'
#' \itemize{
#'   \item \code{pretty}: Use \code{pretty()} to compute “nice” human-readable breaks.
#'   \item \code{equal}: Use \code{seq()} to create equal-width bins.
#'   \item \code{quantile}: Use quantiles to ensure roughly equal sample size per bin.
#' }
#'
#' The function forces the lower bound to be ≤ 0 so that distance zero is always included.
#'
#' @param dist_occ A data.frame or numeric vector of distances from occurrence points.
#' @param dist_bg A data.frame or numeric vector of distances from background points (optional).
#' @param n_bins Integer. Desired number of bins (default = 6, producing ~7 breakpoints).
#' @param method Character string specifying the breakpoint method:
#'   \code{"pretty"}, \code{"equal"}, or \code{"quantile"}.
#'
#' @return A numeric vector of breakpoints sorted from small to large, suitable
#'   for use in \code{cut()}.
#' @export
auto_distance_breaks <- function(dist_occ,
                                 dist_bg = NULL,
                                 n_bins = 6,
                                 method = c("pretty", "equal", "quantile")) {
  method <- match.arg(method)

  get_vals <- function(x) {
    if (is.null(x)) return(numeric(0))
    if (is.data.frame(x)) {
      as.numeric(unlist(x, use.names = FALSE))
    } else {
      as.numeric(x)
    }
  }

  vals <- c(get_vals(dist_occ), get_vals(dist_bg))
  vals <- vals[is.finite(vals) & !is.na(vals)]

  if (length(vals) == 0) {
    stop("auto_distance_breaks(): No valid distance values found (all NA or empty).")
  }

  rng <- range(vals)
  rng[1] <- min(0, rng[1])   # ensure lower bound covers 0

  # Handle degenerate cases where range is zero
  if (rng[1] == rng[2]) {
    if (rng[2] == 0) {
      rng[2] <- 1
    } else {
      rng[1] <- 0
    }
  }

  # Choose method
  if (method == "pretty") {
    brks <- pretty(rng, n = n_bins)
  } else if (method == "equal") {
    brks <- seq(rng[1], rng[2], length.out = n_bins + 1)
  } else if (method == "quantile") {
    probs <- seq(0, 1, length.out = n_bins + 1)
    brks <- stats::quantile(vals, probs = probs, na.rm = TRUE, names = FALSE)
    brks[1] <- min(0, brks[1])
  }

  brks <- sort(unique(brks))
  if (length(brks) < 2) {
    brks <- seq(rng[1], rng[2], length.out = 2)
  }

  brks
}


#' Compute frequency of distances within defined bins
#'
#' @param dist_vec Numeric vector of distances.
#' @param breaks Numeric vector of class breakpoints.
#'
#' @return A data.frame with columns:
#'   \code{distance} (bin midpoint),
#'   \code{frequency} (percentage within each bin).
#' @export
compute_distance_frequency <- function(dist_vec, breaks) {
  cut_dist <- cut(dist_vec,
                  breaks = breaks,
                  include.lowest = TRUE,
                  right = FALSE)

  tab <- table(cut_dist)
  freq_pct <- 100 * as.numeric(tab) / sum(tab)

  mids <- (head(breaks, -1) + tail(breaks, -1)) / 2

  data.frame(
    distance = mids,
    frequency = freq_pct
  )
}


#' Build long-format frequency table for any number of habitats
#'
#' @description
#' Bin distances of occurrence points and background points into classes and
#' compute frequency (%) for each habitat and each distance class.
#'
#' If \code{breaks = NULL}, the breakpoints are generated automatically
#' using \code{auto_distance_breaks()}.
#'
#' @param dist_occ A data.frame where each column is a habitat and each value
#'   is the distance from an occurrence point to that habitat.
#' @param dist_bg A data.frame where each column is a habitat and each value
#'   is the distance from a background point to that habitat.
#' @param breaks Numeric vector of class breakpoints. If NULL, breaks are
#'   generated automatically.
#' @param n_bins Number of bins to generate when \code{breaks = NULL}.
#' @param method Method for automatic break generation:
#'   \code{"pretty"}, \code{"equal"}, or \code{"quantile"}.
#'
#' @return A long-format data.frame with columns:
#'   \code{habitat}, \code{distance}, \code{sampling}, \code{occurrence}.
#' @export
build_frequency_table_generic_long <- function(dist_occ,
                                               dist_bg,
                                               breaks = NULL,
                                               n_bins = 6,
                                               method = c("pretty", "equal", "quantile")) {
  method <- match.arg(method)

  # If no breaks supplied, generate automatically
  if (is.null(breaks)) {
    breaks <- auto_distance_breaks(
      dist_occ = dist_occ,
      dist_bg  = dist_bg,
      n_bins   = n_bins,
      method   = method
    )
  }

  habitats <- intersect(names(dist_occ), names(dist_bg))
  mids <- (head(breaks, -1) + tail(breaks, -1)) / 2

  all_list <- list()

  for (hname in habitats) {
    fo <- compute_distance_frequency(dist_occ[[hname]], breaks)
    fb <- compute_distance_frequency(dist_bg[[hname]],  breaks)

    all_list[[hname]] <- data.frame(
      habitat    = hname,
      distance   = mids,
      sampling   = fb$frequency,
      occurrence = fo$frequency
    )
  }

  do.call(rbind, all_list)
}


#' Convert long-format frequency table to wide format
#'
#' @param freq_long A long-format table with columns:
#'   \code{habitat}, \code{distance}, \code{sampling}, \code{occurrence}.
#'
#' @return A wide-format data.frame with columns:
#'   \code{distance},
#'   \code{sampling_<habitat>},
#'   \code{occurrence_<habitat>}.
#' @export
build_frequency_table_generic_wide <- function(freq_long) {
  freq_long %>%
    tidyr::pivot_wider(
      id_cols = distance,
      names_from = habitat,
      values_from = c(sampling, occurrence),
      names_glue = "{.value}_{habitat}"
    )
}


#' Build frequency tables directly from uploaded distance CSV files
#'
#' @description
#' This function is used when the user already has precomputed distance values
#' (e.g., from GIS software). The occurrence-distance file and background-distance
#' file must contain the same set of habitat distance columns.
#'
#' The function:
#'   1. Reads the two files,
#'   2. Extracts the specified habitat distance columns,
#'   3. Computes distance–frequency distributions,
#'   4. Returns both long and wide format tables.
#'
#' If \code{breaks = NULL}, breakpoints are generated automatically.
#'
#' @param occ_file Character string. Path to the occurrence distance CSV file.
#' @param bg_file Character string. Path to the background distance CSV file.
#' @param habitat_cols Character vector listing the column names representing
#'   habitat distance variables.
#' @param breaks Optional numeric vector of class breakpoints. If NULL, breaks
#'   are generated automatically.
#' @param n_bins Number of bins when breaks are generated automatically.
#' @param method Breakpoint method: \code{"pretty"}, \code{"equal"}, \code{"quantile"}.
#'
#' @return A list containing:
#'   \itemize{
#'     \item \code{freq_long}: long-format frequency table.
#'     \item \code{freq_wide}: wide-format frequency table.
#'     \item \code{breaks}: the breakpoints used.
#'   }
#' @export
build_frequency_from_uploaded_files <- function(occ_file,
                                                bg_file,
                                                habitat_cols,
                                                breaks = NULL,
                                                n_bins = 6,
                                                method = c("pretty", "equal", "quantile")) {
  method <- match.arg(method)

  occ_raw <- readr::read_csv(occ_file, show_col_types = FALSE)
  bg_raw  <- readr::read_csv(bg_file,  show_col_types = FALSE)

  # Check for missing columns
  missing_occ <- setdiff(habitat_cols, names(occ_raw))
  missing_bg  <- setdiff(habitat_cols, names(bg_raw))

  if (length(missing_occ) > 0) {
    stop("Missing habitat distance columns in occurrence file: ",
         paste(missing_occ, collapse = ", "))
  }
  if (length(missing_bg) > 0) {
    stop("Missing habitat distance columns in background file: ",
         paste(missing_bg, collapse = ", "))
  }

  # Extract distance columns
  dist_occ <- occ_raw[, habitat_cols, drop = FALSE]
  dist_bg  <- bg_raw[,  habitat_cols, drop = FALSE]

  # Build frequency tables
  freq_long <- build_frequency_table_generic_long(
    dist_occ  = dist_occ,
    dist_bg   = dist_bg,
    breaks    = breaks,
    n_bins    = n_bins,
    method    = method
  )

  freq_wide <- build_frequency_table_generic_wide(freq_long)

  # Return breaks if they were generated internally
  if (is.null(breaks)) {
    breaks <- auto_distance_breaks(
      dist_occ = dist_occ,
      dist_bg  = dist_bg,
      n_bins   = n_bins,
      method   = method
    )
  }

  list(
    freq_long = freq_long,
    freq_wide = freq_wide,
    breaks    = breaks
  )
}
