# 05_trend_models.R
# Trend fitting and summary for frequency–distance relationships
# (Least-squares LM and GAM, with 95% confidence intervals)

#' Fit a least-squares trend (linear or polynomial) for frequency vs distance
#'
#' @description
#' Fit a simple least-squares model using \code{lm()}:
#' \itemize{
#'   \item degree = 1: linear model, frequency ~ distance
#'   \item degree > 1: polynomial model, frequency ~ poly(distance, degree)
#' }
#'
#' The function returns fitted values and approximate 95% confidence
#' intervals for frequency (in %). Predictions and intervals are truncated
#' to the range 0–100.
#'
#' @param distance Numeric vector of distances (e.g. bin midpoints).
#' @param frequency Numeric vector of frequencies (%).
#' @param degree Integer, polynomial degree (1 = linear).
#' @param x_pred Optional numeric vector of distances
#'   where predictions should be made. If NULL, an equally spaced grid
#'   is created between the minimum and maximum observed distance.
#'
#' @return A list with components:
#'   \itemize{
#'     \item \code{model}: fitted \code{lm} object;
#'     \item \code{data}: data.frame with non-NA input data;
#'     \item \code{pred}: data.frame with columns:
#'       \code{distance}, \code{fit}, \code{lower}, \code{upper}
#'       (approx. 95\% confidence intervals, truncated to the range 0–100).
#'   }
#'   Returns \code{NULL} if there are not enough points to fit the model.
#' @export
fit_lm_trend <- function(distance,
                         frequency,
                         degree = 1,
                         x_pred = NULL) {
  df <- data.frame(
    distance  = distance,
    frequency = frequency
  )
  df <- df[is.finite(df$distance) & is.finite(df$frequency), , drop = FALSE]

  if (nrow(df) < degree + 2) {
    return(NULL)
  }

  if (degree == 1) {
    form <- stats::as.formula("frequency ~ distance")
  } else {
    form <- stats::as.formula(
      paste0("frequency ~ poly(distance, ", degree, ", raw = TRUE)")
    )
  }

  fit <- stats::lm(form, data = df)

  # Prediction grid
  if (is.null(x_pred)) {
    xr <- range(df$distance, na.rm = TRUE)
    x_pred <- seq(xr[1], xr[2], length.out = 200)
  }

  newdf <- data.frame(distance = x_pred)
  pred  <- stats::predict(fit, newdata = newdf, se.fit = TRUE)

  crit  <- 1.96  # ~95% CI
  fit_y <- as.numeric(pred$fit)
  se_y  <- as.numeric(pred$se.fit)

  lower_raw <- fit_y - crit * se_y
  upper_raw <- fit_y + crit * se_y

  # Truncate to 0–100 because frequency is a percentage
  fit_y  <- pmin(pmax(fit_y,  0), 100)
  lower  <- pmin(pmax(lower_raw, 0), 100)
  upper  <- pmin(pmax(upper_raw, 0), 100)

  pred_df <- data.frame(
    distance = x_pred,
    fit      = fit_y,
    lower    = lower,
    upper    = upper
  )

  list(
    model = fit,
    data  = df,
    pred  = pred_df
  )
}

#' Fit a GAM smooth trend for frequency vs distance
#'
#' @description
#' Fit a generalized additive model using \code{mgcv::gam}:
#' \deqn{frequency = s(distance)}
#' with a Gaussian family. This captures smooth, non-linear trends
#' in frequency vs distance.
#'
#' The function returns fitted values and approximate 95% confidence
#' intervals for frequency (in %). Predictions and intervals are truncated
#' to the range 0–100.
#'
#' @param distance Numeric vector of distances.
#' @param frequency Numeric vector of frequencies (%).
#' @param k Basis dimension for the smooth (passed to \code{s(distance, k = k)}).
#' @param x_pred Optional numeric vector of distances for predictions.
#'
#' @return A list with components:
#'   \itemize{
#'     \item \code{model}: fitted \code{gam} object;
#'     \item \code{data}: data.frame with non-NA input;
#'     \item \code{pred}: data.frame with columns:
#'       \code{distance}, \code{fit}, \code{lower}, \code{upper}
#'       (approx. 95\% confidence intervals, truncated to the range 0–100).
#'   }
#'   Returns \code{NULL} if there are not enough points to fit the smooth.
#' @export
fit_gam_trend <- function(distance,
                          frequency,
                          k = 5,
                          x_pred = NULL) {
  df <- data.frame(
    distance  = distance,
    frequency = frequency
  )
  df <- df[is.finite(df$distance) & is.finite(df$frequency), , drop = FALSE]

  if (nrow(df) < (k + 2)) {
    # not enough data to justify a smooth with k basis functions
    return(NULL)
  }

  # IMPORTANT: use s(), not mgcv::s(), inside the formula
  fit <- mgcv::gam(
    frequency ~ s(distance, k = k),
    data   = df,
    family = stats::gaussian()
  )

  if (is.null(x_pred)) {
    xr <- range(df$distance, na.rm = TRUE)
    x_pred <- seq(xr[1], xr[2], length.out = 200)
  }

  newdf <- data.frame(distance = x_pred)
  pred  <- stats::predict(fit, newdata = newdf, se.fit = TRUE)

  crit  <- 1.96  # ~95% CI
  fit_y <- as.numeric(pred$fit)
  se_y  <- as.numeric(pred$se.fit)

  lower_raw <- fit_y - crit * se_y
  upper_raw <- fit_y + crit * se_y

  # Truncate to 0–100 because frequency is a percentage
  fit_y  <- pmin(pmax(fit_y,  0), 100)
  lower  <- pmin(pmax(lower_raw, 0), 100)
  upper  <- pmin(pmax(upper_raw, 0), 100)

  pred_df <- data.frame(
    distance = x_pred,
    fit      = fit_y,
    lower    = lower,
    upper    = upper
  )

  list(
    model = fit,
    data  = df,
    pred  = pred_df
  )
}

#' Fit LM and/or GAM trends for one habitat (sampling vs occurrence)
#'
#' @description
#' Convenience helper to fit trends for both sampling (background)
#' and occurrence frequencies for a single habitat. Each fitted model
#' carries predicted values and 95% confidence intervals via
#' \code{fit_lm_trend()} and \code{fit_gam_trend()}.
#'
#' @param df_habitat data.frame with columns:
#'   \code{distance}, \code{sampling}, \code{occurrence}.
#' @param methods Character vector, any of \code{c("lm", "gam")}.
#' @param degree Polynomial degree for LM (default 1).
#' @param k Basis dimension for GAM (default 5).
#'
#' @return A list with elements:
#'   \itemize{
#'     \item \code{sampling_lm}, \code{occurrence_lm} (if "lm" in methods);
#'     \item \code{sampling_gam}, \code{occurrence_gam} (if "gam" in methods).
#'   }
#'   Each element may be \code{NULL} if the corresponding model could not be fitted.
#' @export
fit_trends_for_habitat <- function(df_habitat,
                                   methods = c("lm", "gam"),
                                   degree = 1,
                                   k = 5) {
  methods <- intersect(methods, c("lm", "gam"))
  if (length(methods) == 0) {
    stop("methods must include at least one of 'lm' or 'gam'.")
  }

  distance   <- df_habitat$distance
  sampling   <- df_habitat$sampling
  occurrence <- df_habitat$occurrence

  res <- list(
    sampling_lm    = NULL,
    occurrence_lm  = NULL,
    sampling_gam   = NULL,
    occurrence_gam = NULL
  )

  if ("lm" %in% methods) {
    res$sampling_lm   <- fit_lm_trend(distance, sampling,   degree = degree)
    res$occurrence_lm <- fit_lm_trend(distance, occurrence, degree = degree)
  }

  if ("gam" %in% methods) {
    res$sampling_gam   <- fit_gam_trend(distance, sampling,   k = k)
    res$occurrence_gam <- fit_gam_trend(distance, occurrence, k = k)
  }

  res
}

#' Summarize LM and GAM trends for a single habitat
#'
#' @description
#' Internal helper: given a data.frame with columns
#' \code{distance}, \code{sampling}, \code{occurrence}, fit
#' LM and/or GAM trends (via \code{fit_trends_for_habitat()}) and
#' extract model statistics (R², AIC, explained deviance, etc.).
#' Confidence intervals are handled inside the fitted objects'
#' \code{pred} components and are not repeated here.
#'
#' @param df_habitat data.frame with columns:
#'   \itemize{
#'     \item \code{distance}
#'     \item \code{sampling}
#'     \item \code{occurrence}
#'   }
#' @param habitat_name Character, habitat name for labeling.
#' @param methods Character vector, subset of \code{c("lm", "gam")}.
#' @param degree Polynomial degree for LM (default 1).
#' @param k Basis dimension for GAM smooth (default 5).
#'
#' @return A one-row data.frame with summary statistics for this habitat.
#' @keywords internal
summarize_trends_for_habitat <- function(df_habitat,
                                         habitat_name,
                                         methods = c("lm", "gam"),
                                         degree = 1,
                                         k = 5) {
  methods <- intersect(methods, c("lm", "gam"))
  if (length(methods) == 0) {
    stop("methods must include at least one of 'lm' or 'gam'.")
  }

  # Fit trends
  trend_res <- fit_trends_for_habitat(
    df_habitat = df_habitat,
    methods    = methods,
    degree     = degree,
    k          = k
  )

  # Initialize result row with NA values
  out <- data.frame(
    habitat = habitat_name,

    # LM (sampling)
    lm_sampling_degree       = NA_real_,
    lm_sampling_R2           = NA_real_,
    lm_sampling_adjR2        = NA_real_,
    lm_sampling_AIC          = NA_real_,

    # LM (occurrence)
    lm_occurrence_degree     = NA_real_,
    lm_occurrence_R2         = NA_real_,
    lm_occurrence_adjR2      = NA_real_,
    lm_occurrence_AIC        = NA_real_,

    # GAM (sampling)
    gam_sampling_k           = NA_real_,
    gam_sampling_expl_dev    = NA_real_,
    gam_sampling_AIC         = NA_real_,

    # GAM (occurrence)
    gam_occurrence_k         = NA_real_,
    gam_occurrence_expl_dev  = NA_real_,
    gam_occurrence_AIC       = NA_real_,

    stringsAsFactors = FALSE
  )

  ## --- LM summaries -------------------------------------------------------
  if ("lm" %in% methods) {
    # sampling
    if (!is.null(trend_res$sampling_lm)) {
      fit_lm_s <- trend_res$sampling_lm$model
      sm_s     <- summary(fit_lm_s)

      out$lm_sampling_degree <- degree
      out$lm_sampling_R2     <- unname(sm_s$r.squared)
      out$lm_sampling_adjR2  <- unname(sm_s$adj.r.squared)
      out$lm_sampling_AIC    <- stats::AIC(fit_lm_s)
    }

    # occurrence
    if (!is.null(trend_res$occurrence_lm)) {
      fit_lm_o <- trend_res$occurrence_lm$model
      sm_o     <- summary(fit_lm_o)

      out$lm_occurrence_degree <- degree
      out$lm_occurrence_R2     <- unname(sm_o$r.squared)
      out$lm_occurrence_adjR2  <- unname(sm_o$adj.r.squared)
      out$lm_occurrence_AIC    <- stats::AIC(fit_lm_o)
    }
  }

  ## --- GAM summaries ------------------------------------------------------
  if ("gam" %in% methods) {
    # sampling
    if (!is.null(trend_res$sampling_gam)) {
      fit_gam_s <- trend_res$sampling_gam$model
      sm_gs     <- mgcv::summary.gam(fit_gam_s)

      # mgcv::summary.gam reports dev.expl as a fraction (0–1)
      expl_dev_s <- sm_gs$dev.expl * 100

      out$gam_sampling_k         <- k
      out$gam_sampling_expl_dev  <- expl_dev_s
      out$gam_sampling_AIC       <- stats::AIC(fit_gam_s)
    }

    # occurrence
    if (!is.null(trend_res$occurrence_gam)) {
      fit_gam_o <- trend_res$occurrence_gam$model
      sm_go     <- mgcv::summary.gam(fit_gam_o)
      expl_dev_o <- sm_go$dev.expl * 100

      out$gam_occurrence_k        <- k
      out$gam_occurrence_expl_dev <- expl_dev_o
      out$gam_occurrence_AIC      <- stats::AIC(fit_gam_o)
    }
  }

  out
}

#' Generate LM/GAM trend summary for multiple habitats (generic)
#'
#' @description
#' For each habitat, fit least-squares (LM) and/or GAM smooth trends
#' of frequency vs distance for:
#' \itemize{
#'   \item sampling frequencies (background points);
#'   \item occurrence frequencies (occurrence points).
#' }
#'
#' For LM, the summary includes:
#' \itemize{
#'   \item polynomial degree;
#'   \item R² and adjusted R²;
#'   \item AIC.
#' }
#'
#' For GAM, the summary includes:
#' \itemize{
#'   \item basis dimension \code{k};
#'   \item explained deviance (%);
#'   \item AIC.
#' }
#'
#' The function returns a wide-format table with one row per habitat.
#' You can choose LM only, GAM only, or both via the \code{methods} argument.
#'
#' @param freq_wide A wide-format frequency table, typically produced by
#'   \code{build_frequency_table_generic_wide()}, with columns:
#'   \itemize{
#'     \item \code{distance}
#'     \item \code{sampling_<habitat>}
#'     \item \code{occurrence_<habitat>}
#'   }
#' @param methods Character vector, subset of \code{c("lm", "gam")}.
#'   Defaults to both: \code{c("lm", "gam")}.
#' @param degree Polynomial degree for LM (default 1).
#' @param k Basis dimension for GAM smooth (default 5).
#'
#' @return A data.frame with one row per habitat and columns summarizing
#'   LM and/or GAM fit quality for sampling and occurrence trends.
#' @export
generate_trend_summary_multi <- function(freq_wide,
                                         methods = c("lm", "gam"),
                                         degree = 1,
                                         k = 5) {
  methods <- intersect(methods, c("lm", "gam"))
  if (length(methods) == 0) {
    stop("methods must include at least one of 'lm' or 'gam'.")
  }

  all_cols      <- names(freq_wide)
  sampling_cols <- grep("^sampling_", all_cols, value = TRUE)
  habitats      <- sub("^sampling_", "", sampling_cols)

  out_list <- list()

  for (hname in habitats) {
    s_col <- paste0("sampling_", hname)
    o_col <- paste0("occurrence_", hname)

    if (!all(c(s_col, o_col, "distance") %in% all_cols)) {
      next
    }

    df_h <- data.frame(
      distance   = freq_wide$distance,
      sampling   = freq_wide[[s_col]],
      occurrence = freq_wide[[o_col]]
    )

    row_h <- summarize_trends_for_habitat(
      df_habitat   = df_h,
      habitat_name = hname,
      methods      = methods,
      degree       = degree,
      k            = k
    )

    out_list[[hname]] <- row_h
  }

  if (length(out_list) == 0) {
    warning("No habitats found with both sampling_ and occurrence_ columns.")
    return(NULL)
  }

  do.call(rbind, out_list)
}

#' LM-only trend summary (least-squares option, separated)
#'
#' @description
#' A convenience wrapper that runs \code{generate_trend_summary_multi()}
#' using only LM (least-squares) trends. This keeps LM analysis
#' completely separate from GAM, while still using 95% confidence
#' intervals internally for plotting via \code{fit_lm_trend()}.
#'
#' @param freq_wide A wide-format frequency table, typically produced by
#'   \code{build_frequency_table_generic_wide()}, with columns:
#'   \code{distance}, \code{sampling_<habitat>}, \code{occurrence_<habitat>}.
#' @param degree Polynomial degree for LM (default 1 = linear).
#'
#' @return A data.frame with one row per habitat summarizing LM trends
#'   (R², adjusted R², AIC, and degree).
#' @export
generate_lm_trend_summary_multi <- function(freq_wide,
                                            degree = 1) {
  generate_trend_summary_multi(
    freq_wide = freq_wide,
    methods   = "lm",
    degree    = degree,
    k         = 5  # ignored for LM
  )
}

#' GAM-only trend summary (smooth option, separated)
#'
#' @description
#' A convenience wrapper that runs \code{generate_trend_summary_multi()}
#' using only GAM trends. This keeps GAM analysis completely
#' separate from LM, while still using 95% confidence intervals
#' internally for plotting via \code{fit_gam_trend()}.
#'
#' @param freq_wide A wide-format frequency table, typically produced by
#'   \code{build_frequency_table_generic_wide()}, with columns:
#'   \code{distance}, \code{sampling_<habitat>}, \code{occurrence_<habitat>}.
#' @param k Basis dimension for GAM smooth (default 5).
#'
#' @return A data.frame with one row per habitat summarizing GAM trends
#'   (explained deviance %, AIC, and k).
#' @export
generate_gam_trend_summary_multi <- function(freq_wide,
                                             k = 5) {
  generate_trend_summary_multi(
    freq_wide = freq_wide,
    methods   = "gam",
    degree    = 1, # ignored for GAM
    k         = k
  )
}
