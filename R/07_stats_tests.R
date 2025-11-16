#' Bootstrap analysis of decay rate (including zeros)
#'
#' @description
#' Perform a bootstrap analysis of the exponential decay rate parameter \code{b}
#' using the fitted exponential decay model. Zero distances and zero
#' frequencies are allowed in the input.
#'
#' @param distance Numeric vector of distances.
#' @param frequency Numeric vector of frequencies.
#' @param n_iterations Integer. Number of bootstrap resamples (default = 1000).
#' @param ci_level Confidence level for the bootstrap confidence interval
#'   (default = 0.95).
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{b_bootstrap}: numeric vector with bootstrap estimates of \code{b}.
#'   \item \code{params}: estimated parameters from the original fit.
#'   \item \code{r_squared}: \eqn{R^2} of the original fit.
#'   \item \code{conf_interval}: bootstrap confidence interval for \code{b}.
#'   \item \code{data}: cleaned data used for fitting.
#' }
#'   Returns \code{NULL} if the initial fit fails or no valid bootstrap
#'   iteration can be obtained.
#' @export
bootstrap_decay_analysis_with_zeros <- function(distance,
                                                frequency,
                                                n_iterations = 1000,
                                                ci_level = 0.95) {
  fit_result <- fit_exponential_decay(distance, frequency)
  if (!fit_result$success) return(NULL)

  n <- nrow(fit_result$data)
  b_bootstrap <- numeric(n_iterations)
  valid_iterations <- 0

  for (i in seq_len(n_iterations)) {
    idx <- sample(seq_len(n), n, replace = TRUE)
    data_bs <- fit_result$data[idx, ]

    tryCatch({
      fit_bs <- minpack.lm::nlsLM(
        frequency ~ exponential_decay(distance, N0, b),
        data = data_bs,
        start = list(N0 = max(data_bs$frequency), b = 0.01),
        control = minpack.lm::nls.lm.control(maxiter = 500)
      )
      valid_iterations <- valid_iterations + 1
      b_bootstrap[valid_iterations] <- stats::coef(fit_bs)["b"]
    }, error = function(e) {})
  }

  if (valid_iterations == 0) return(NULL)

  b_bootstrap <- b_bootstrap[seq_len(valid_iterations)]
  alpha <- 1 - ci_level
  conf_interval <- stats::quantile(
    b_bootstrap,
    c(alpha / 2, 1 - alpha / 2),
    na.rm = TRUE
  )

  list(
    b_bootstrap   = b_bootstrap,
    params        = fit_result$params,
    r_squared     = fit_result$r_squared,
    conf_interval = conf_interval,
    data          = fit_result$data
  )
}

#' Permutation test comparing decay rates between two groups
#'
#' @description
#' Perform a permutation test to compare the exponential decay rate \code{b}
#' between two groups (e.g. sampling vs. occurrence). The null hypothesis
#' is that the two groups share the same decay rate.
#'
#' @param distance_sampling Numeric vector of distances for the sampling group.
#' @param freq_sampling Numeric vector of frequencies for the sampling group.
#' @param distance_occurrence Numeric vector of distances for the occurrence group.
#' @param freq_occurrence Numeric vector of frequencies for the occurrence group.
#' @param n_perm Integer. Number of permutations (default = 2000).
#'
#' @return A list with elements:
#' \itemize{
#'   \item \code{p_value}: permutation p-value for the difference in \code{b}.
#'   \item \code{delta_b_obs}: observed difference \code{b_occurrence - b_sampling}.
#'   \item \code{delta_perm}: vector of permuted differences in \code{b}.
#'   \item \code{valid_iter}: number of permutations with successful model fits.
#'   \item \code{message}: text message describing the result.
#' }
#'   If data are insufficient or fitting fails, \code{p_value} and/or
#'   \code{delta_b_obs} may be \code{NA}.
#' @export
permutation_test_decay_rate <- function(
    distance_sampling, freq_sampling,
    distance_occurrence, freq_occurrence,
    n_perm = 2000
) {
  df_sampling <- data.frame(
    distance  = distance_sampling,
    frequency = freq_sampling
  ) %>%
    dplyr::filter(!is.na(distance), !is.na(frequency))

  df_occurrence <- data.frame(
    distance  = distance_occurrence,
    frequency = freq_occurrence
  ) %>%
    dplyr::filter(!is.na(distance), !is.na(frequency))

  if (nrow(df_sampling) < 3 || nrow(df_occurrence) < 3) {
    return(list(
      p_value     = NA,
      delta_b_obs = NA,
      message     = "Not enough data points to perform permutation test."
    ))
  }

  fit_samp <- fit_exponential_decay(df_sampling$distance, df_sampling$frequency)
  fit_occ  <- fit_exponential_decay(df_occurrence$distance, df_occurrence$frequency)

  if (!fit_samp$success || !fit_occ$success) {
    return(list(
      p_value     = NA,
      delta_b_obs = NA,
      message     = "Model fitting failed on original data, permutation test not performed."
    ))
  }

  b_samp <- fit_samp$params["b"]
  b_occ  <- fit_occ$params["b"]
  delta_b_obs <- as.numeric(b_occ - b_samp)

  df_all <- rbind(
    data.frame(distance = df_sampling$distance,
               frequency = df_sampling$frequency,
               group = "sampling"),
    data.frame(distance = df_occurrence$distance,
               frequency = df_occurrence$frequency,
               group = "occurrence")
  )

  n_samp <- nrow(df_sampling)
  n_occ  <- nrow(df_occurrence)

  delta_perm <- numeric(n_perm)
  valid_iter <- 0

  for (i in seq_len(n_perm)) {
    perm_idx <- sample(seq_len(nrow(df_all)))
    group_perm <- rep(NA_character_, nrow(df_all))
    group_perm[perm_idx[seq_len(n_samp)]] <- "sampling"
    group_perm[perm_idx[(n_samp + 1):(n_samp + n_occ)]] <- "occurrence"

    df_all$group_perm <- group_perm

    df_s <- df_all[df_all$group_perm == "sampling", ]
    df_o <- df_all[df_all$group_perm == "occurrence", ]

    perm_fit_s <- fit_exponential_decay(df_s$distance, df_s$frequency)
    perm_fit_o <- fit_exponential_decay(df_o$distance, df_o$frequency)

    if (perm_fit_s$success && perm_fit_o$success) {
      valid_iter <- valid_iter + 1
      delta_perm[valid_iter] <- as.numeric(
        perm_fit_o$params["b"] - perm_fit_s$params["b"]
      )
    }
  }

  if (valid_iter == 0) {
    return(list(
      p_value     = NA,
      delta_b_obs = delta_b_obs,
      message     = "All permutations failed to produce valid model fits."
    ))
  }

  delta_perm <- delta_perm[seq_len(valid_iter)]
  p_value <- mean(delta_perm >= delta_b_obs)

  list(
    p_value      = p_value,
    delta_b_obs  = delta_b_obs,
    delta_perm   = delta_perm,
    valid_iter   = valid_iter,
    message      = "Permutation test completed."
  )
}

#' Kolmogorov–Smirnov test comparing two distance distributions
#'
#' @description
#' Perform a Kolmogorov–Smirnov (K-S) test to compare the empirical
#' distributions of distances between two groups (e.g. sampling vs. occurrence).
#'
#' @param distance_sampling Numeric vector of distances for the sampling group.
#' @param distance_occurrence Numeric vector of distances for the occurrence group.
#'
#' @return A list with elements:
#' \itemize{
#'   \item \code{ks_p_value}: p-value of the K-S test.
#'   \item \code{ks_D}: K-S test statistic.
#'   \item \code{message}: text message describing the result.
#' }
#'   Returns \code{NA} values if there are fewer than 5 non-NA distances in either group.
#' @export
ks_test_for_habitat <- function(distance_sampling,
                                distance_occurrence) {
  sampling_dist   <- distance_sampling[!is.na(distance_sampling)]
  occurrence_dist <- distance_occurrence[!is.na(distance_occurrence)]

  if (length(sampling_dist) < 5 || length(occurrence_dist) < 5) {
    return(list(
      ks_p_value = NA,
      ks_D       = NA,
      message    = "Not enough data points to perform K-S test."
    ))
  }

  ks_res <- stats::ks.test(sampling_dist, occurrence_dist)

  list(
    ks_p_value = ks_res$p.value,
    ks_D       = ks_res$statistic,
    message    = "K-S test completed."
  )
}

#' Compare decay patterns for a single habitat (Bootstrap + permutation + K-S)
#'
#' @description
#' For a single habitat, compare the decay rates and distance distributions
#' between sampling and occurrence frequencies using:
#'   - bootstrap analysis of the exponential decay rate \code{b},
#'   - a permutation test of \code{b},
#'   - and a K-S test on distance distributions.
#'
#' @param distance Numeric vector of distances (aligned with \code{freq_sampling}
#'   and \code{freq_occurrence}).
#' @param freq_sampling Numeric vector of sampling frequencies (e.g., by distance bin).
#' @param freq_occurrence Numeric vector of occurrence frequencies.
#' @param n_bootstrap Integer. Number of bootstrap iterations (default = 1000).
#' @param n_perm Integer. Number of permutations for the permutation test
#'   (default = 2000).
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{sampling_analysis}: bootstrap result for the sampling group.
#'   \item \code{occurrence_analysis}: bootstrap result for the occurrence group.
#'   \item \code{p_value_bootstrap}: p-value based on the bootstrap comparison.
#'   \item \code{p_value_permutation}: permutation test p-value for \code{b}.
#'   \item \code{delta_b_obs}: observed difference in \code{b}.
#'   \item \code{ks_p_value}: K-S test p-value.
#'   \item \code{ks_D}: K-S test statistic.
#' }
#' @export
compare_decay_one_habitat <- function(distance,
                                      freq_sampling,
                                      freq_occurrence,
                                      n_bootstrap = 1000,
                                      n_perm = 2000) {
  sampling_analysis <- bootstrap_decay_analysis_with_zeros(
    distance, freq_sampling, n_iterations = n_bootstrap
  )
  occurrence_analysis <- bootstrap_decay_analysis_with_zeros(
    distance, freq_occurrence, n_iterations = n_bootstrap
  )

  p_value_boot <- NULL
  if (!is.null(sampling_analysis) && !is.null(occurrence_analysis)) {
    p_value_boot <- mean(
      occurrence_analysis$b_bootstrap > sampling_analysis$b_bootstrap,
      na.rm = TRUE
    )
  }

  perm_res <- permutation_test_decay_rate(
    distance_sampling   = distance,
    freq_sampling       = freq_sampling,
    distance_occurrence = distance,
    freq_occurrence     = freq_occurrence,
    n_perm = n_perm
  )

  ks_res <- ks_test_for_habitat(
    distance_sampling   = distance,
    distance_occurrence = distance
  )

  list(
    sampling_analysis    = sampling_analysis,
    occurrence_analysis  = occurrence_analysis,
    p_value_bootstrap    = p_value_boot,
    p_value_permutation  = perm_res$p_value,
    delta_b_obs          = perm_res$delta_b_obs,
    ks_p_value           = ks_res$ks_p_value,
    ks_D                 = ks_res$ks_D
  )
}

#' Multi-habitat statistical summary of decay rates and distance distributions
#'
#' @description
#' For each habitat, summarize the decay rate parameters and statistical
#' tests comparing sampling vs. occurrence patterns. This function expects a
#' wide-format frequency table produced by
#' \code{build_frequency_table_generic_wide()}, containing columns:
#' \code{distance}, \code{sampling_*}, \code{occurrence_*}.
#'
#' @param freq_wide A wide-format data.frame with columns:
#'   \code{distance}, \code{sampling_<habitat>}, \code{occurrence_<habitat>}.
#'
#' @return A data.frame with one row per habitat and columns:
#' \itemize{
#'   \item \code{habitat}
#'   \item \code{b_sampling}: decay rate for sampling frequencies.
#'   \item \code{b_occurrence}: decay rate for occurrence frequencies.
#'   \item \code{p_bootstrap}: bootstrap-based p-value.
#'   \item \code{p_permutation}: permutation test p-value.
#'   \item \code{p_ks}: K-S test p-value for distance distributions.
#'   \item \code{interpretation}: textual ecological interpretation.
#' }
#' @export
generate_statistical_summary_multi <- function(freq_wide) {
  all_cols      <- names(freq_wide)
  sampling_cols <- grep("^sampling_", all_cols, value = TRUE)
  habitats      <- sub("^sampling_", "", sampling_cols)

  summary_df <- data.frame(
    habitat       = character(),
    b_sampling    = numeric(),
    b_occurrence  = numeric(),
    p_bootstrap   = numeric(),
    p_permutation = numeric(),
    p_ks          = numeric(),
    interpretation = character(),
    stringsAsFactors = FALSE
  )

  for (hname in habitats) {
    s_col <- paste0("sampling_",   hname)
    o_col <- paste0("occurrence_", hname)

    if (!all(c(s_col, o_col) %in% all_cols)) next

    dist  <- freq_wide$distance
    fsamp <- freq_wide[[s_col]]
    focc  <- freq_wide[[o_col]]

    res_h <- compare_decay_one_habitat(
      distance        = dist,
      freq_sampling   = fsamp,
      freq_occurrence = focc
    )

    sampling_b <- if (!is.null(res_h$sampling_analysis))
      round(res_h$sampling_analysis$params["b"], 4) else NA
    occurrence_b <- if (!is.null(res_h$occurrence_analysis))
      round(res_h$occurrence_analysis$params["b"], 4) else NA

    p_boot <- if (!is.null(res_h$p_value_bootstrap))
      round(res_h$p_value_bootstrap, 4) else NA
    p_perm <- if (!is.null(res_h$p_value_permutation))
      round(res_h$p_value_permutation, 4) else NA
    p_ks   <- if (!is.null(res_h$ks_p_value))
      round(res_h$ks_p_value, 4) else NA

    interpretation <- if (!is.na(p_perm) && p_perm < 0.05) {
      "Permutation test indicates a significantly higher decay rate for occurrence frequencies, suggesting a preference for this habitat."
    } else if (!is.na(p_boot) && p_boot < 0.05) {
      "Bootstrap analysis indicates a significant difference in decay rates between sampling and occurrence."
    } else if (!is.na(p_ks) && p_ks < 0.05) {
      "K-S test indicates a significant difference in distance distributions between sampling and occurrence."
    } else {
      "No strong statistical evidence of habitat preference; the effect may be weak or sample size may be insufficient."
    }

    summary_df <- rbind(summary_df, data.frame(
      habitat        = hname,
      b_sampling     = sampling_b,
      b_occurrence   = occurrence_b,
      p_bootstrap    = p_boot,
      p_permutation  = p_perm,
      p_ks           = p_ks,
      interpretation = interpretation,
      stringsAsFactors = FALSE
    ))
  }

  summary_df
}
