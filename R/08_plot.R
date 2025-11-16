############################################################
# Decay plots for multiple habitats
############################################################

#' Plot exponential decay for a single habitat
#'
#' @description
#' Plot sampling and occurrence frequencies vs distance for a single habitat,
#' together with fitted exponential decay curves and confidence intervals.
#'
#' @param df_habitat A data.frame with columns:
#'   \code{distance}, \code{sampling}, \code{occurrence}.
#' @param habitat_name Character string used as plot title.
#'
#' @return A ggplot object.
#' @export
plot_single_habitat_decay <- function(df_habitat,
                                      habitat_name) {

  # Fit exponential decay models
  fit_samp <- fit_exponential_decay(df_habitat$distance, df_habitat$sampling)
  fit_occ  <- fit_exponential_decay(df_habitat$distance, df_habitat$occurrence)

  # Prediction grid
  x_range <- seq(
    min(df_habitat$distance, na.rm = TRUE),
    max(df_habitat$distance, na.rm = TRUE),
    length.out = 200
  )

  ci_samp <- if (isTRUE(fit_samp$success))
    calculate_confidence_intervals(fit_samp, x_range) else NULL
  ci_occ  <- if (isTRUE(fit_occ$success))
    calculate_confidence_intervals(fit_occ,  x_range) else NULL

  p <- ggplot2::ggplot(df_habitat) +
    ggplot2::geom_point(
      ggplot2::aes(x = distance, y = sampling, color = "Sampling"),
      size = 1.8, alpha = 0.8
    ) +
    ggplot2::geom_point(
      ggplot2::aes(x = distance, y = occurrence, color = "Occurrence"),
      size = 1.8, alpha = 0.8
    )

  if (!is.null(ci_samp)) {
    p <- p +
      ggplot2::geom_line(
        data = ci_samp,
        ggplot2::aes(x = distance, y = fit, color = "Sampling"),
        linewidth = 0.9
      ) +
      ggplot2::geom_ribbon(
        data = ci_samp,
        ggplot2::aes(x = distance, ymin = lower, ymax = upper),
        inherit.aes = FALSE,
        fill = "#1f77b4", alpha = 0.18
      )
  }

  if (!is.null(ci_occ)) {
    p <- p +
      ggplot2::geom_line(
        data = ci_occ,
        ggplot2::aes(x = distance, y = fit, color = "Occurrence"),
        linewidth = 0.9
      ) +
      ggplot2::geom_ribbon(
        data = ci_occ,
        ggplot2::aes(x = distance, ymin = lower, ymax = upper),
        inherit.aes = FALSE,
        fill = "#ff7f0e", alpha = 0.18
      )
  }

  p +
    ggplot2::scale_color_manual(
      values = c("Sampling" = "#1f77b4", "Occurrence" = "#ff7f0e")
    ) +
    ggplot2::labs(
      x = "Distance (km)",
      y = "Frequency (%)",
      color = "",
      title = habitat_name
    ) +
    ggplot2::theme_bw(base_size = 10) +
    ggplot2::theme(
      legend.position  = "bottom",
      panel.grid.major = ggplot2::element_line(color = "grey90", linewidth = 0.2),
      panel.grid.minor = ggplot2::element_blank(),
      plot.title       = ggplot2::element_text(face = "bold", size = 11)
    )
}


#' Plot exponential decay for multiple habitats (multi-panel)
#'
#' @description
#' Generate a multi-panel figure showing exponential decay curves and
#' confidence intervals for multiple habitats.
#'
#' @param freq_long A data.frame with columns:
#'   \code{habitat}, \code{distance}, \code{sampling}, \code{occurrence}.
#' @param ncol Integer. Number of columns in the multi-panel layout (default 2).
#'
#' @return A patchwork object with multiple panels.
#' @export
plot_multi_habitat_decay <- function(freq_long, ncol = 2) {
  habitats <- unique(freq_long$habitat)

  plot_list <- lapply(seq_along(habitats), function(i) {
    hname <- habitats[i]
    df_h  <- freq_long[freq_long$habitat == hname, , drop = FALSE]

    letter <- letters[i]
    title  <- paste0("(", letter, ") ", hname)

    plot_single_habitat_decay(df_habitat = df_h, habitat_name = title)
  })

  patchwork::wrap_plots(plot_list, ncol = ncol) +
    patchwork::plot_annotation(
      theme = ggplot2::theme(
        plot.margin = ggplot2::margin(8, 8, 8, 8)
      )
    )
}

############################################################
# LM / GAM trend plots for single and multiple habitats
############################################################

#' Plot LM/GAM trends for a single habitat
#'
#' @description
#' Plot least-squares (LM) and/or GAM smooth trends of frequency vs distance
#' for both sampling (background) and occurrence frequencies.
#'
#' Trend curves and confidence ribbons use the predictions produced by
#' \code{fit_trends_for_habitat()}, which internally calls
#' \code{fit_lm_trend()} and/or \code{fit_gam_trend()} and returns
#' approximate 95\% confidence intervals.
#'
#' @param df_habitat data.frame with columns:
#'   \code{distance}, \code{sampling}, \code{occurrence}.
#' @param habitat_name Character, used as plot title.
#' @param methods Character vector, any of "lm", "gam".
#'   \itemize{
#'     \item "lm": least-squares linear/polynomial trend;
#'     \item "gam": smooth non-linear trend (GAM).
#'   }
#' @param degree Polynomial degree for LM (default 1 = linear).
#' @param k Basis dimension for GAM smooth (default 5).
#'
#' @return A ggplot object.
#' @export
plot_single_habitat_trend <- function(df_habitat,
                                      habitat_name,
                                      methods = c("lm", "gam"),
                                      degree = 1,
                                      k = 5) {
  methods <- intersect(methods, c("lm", "gam"))
  if (length(methods) == 0) {
    stop("methods must include at least one of 'lm' or 'gam'.")
  }

  plot_data <- df_habitat

  # Fit trends (uses 05_trend_models.R implementations)
  trend_res <- fit_trends_for_habitat(
    df_habitat = plot_data,
    methods    = methods,
    degree     = degree,
    k          = k
  )

  # Base scatter plot (points)
  p <- ggplot2::ggplot(plot_data) +
    ggplot2::geom_point(
      ggplot2::aes(x = distance, y = sampling),
      color = "#1f77b4", size = 1.5, alpha = 0.7
    ) +
    ggplot2::geom_point(
      ggplot2::aes(x = distance, y = occurrence),
      color = "#ff7f0e", size = 1.5, alpha = 0.7
    )

  # Add LM trends (if fitted)
  if ("lm" %in% methods) {
    if (!is.null(trend_res$sampling_lm)) {
      df_lm_s <- trend_res$sampling_lm$pred
      p <- p +
        ggplot2::geom_line(
          data = df_lm_s,
          ggplot2::aes(x = distance, y = fit),
          inherit.aes = FALSE,
          color = "#1f77b4",
          linewidth = 0.8,
          linetype  = "dashed"
        ) +
        ggplot2::geom_ribbon(
          data = df_lm_s,
          ggplot2::aes(x = distance, ymin = lower, ymax = upper),
          inherit.aes = FALSE,
          fill = "#1f77b4", alpha = 0.15
        )
    }
    if (!is.null(trend_res$occurrence_lm)) {
      df_lm_o <- trend_res$occurrence_lm$pred
      p <- p +
        ggplot2::geom_line(
          data = df_lm_o,
          ggplot2::aes(x = distance, y = fit),
          inherit.aes = FALSE,
          color = "#ff7f0e",
          linewidth = 0.8,
          linetype  = "dashed"
        ) +
        ggplot2::geom_ribbon(
          data = df_lm_o,
          ggplot2::aes(x = distance, ymin = lower, ymax = upper),
          inherit.aes = FALSE,
          fill = "#ff7f0e", alpha = 0.15
        )
    }
  }

  # Add GAM trends (if fitted)
  if ("gam" %in% methods) {
    if (!is.null(trend_res$sampling_gam)) {
      df_gam_s <- trend_res$sampling_gam$pred
      p <- p +
        ggplot2::geom_line(
          data = df_gam_s,
          ggplot2::aes(x = distance, y = fit),
          inherit.aes = FALSE,
          color = "#1f77b4",
          linewidth = 1.2
        ) +
        ggplot2::geom_ribbon(
          data = df_gam_s,
          ggplot2::aes(x = distance, ymin = lower, ymax = upper),
          inherit.aes = FALSE,
          fill = "#1f77b4", alpha = 0.15
        )
    }
    if (!is.null(trend_res$occurrence_gam)) {
      df_gam_o <- trend_res$occurrence_gam$pred
      p <- p +
        ggplot2::geom_line(
          data = df_gam_o,
          ggplot2::aes(x = distance, y = fit),
          inherit.aes = FALSE,
          color = "#ff7f0e",
          linewidth = 1.2
        ) +
        ggplot2::geom_ribbon(
          data = df_gam_o,
          ggplot2::aes(x = distance, ymin = lower, ymax = upper),
          inherit.aes = FALSE,
          fill = "#ff7f0e", alpha = 0.15
        )
    }
  }

  p +
    ggplot2::labs(
      x = "Distance (km)",
      y = "Frequency (%)",
      title = habitat_name,
      subtitle = paste("Trends:", paste(methods, collapse = " + "))
    ) +
    ggplot2::ylim(0, 50) +
    ggplot2::theme_bw(base_size = 10) +
    ggplot2::theme(
      panel.grid.major = ggplot2::element_line(color = "grey90", linewidth = 0.2),
      panel.grid.minor = ggplot2::element_blank(),
      plot.title       = ggplot2::element_text(hjust = 0, face = "bold", size = 12),
      axis.title       = ggplot2::element_text(size = 9),
      axis.text        = ggplot2::element_text(size = 8),
      legend.position  = "none"
    )
}

#' Plot LM/GAM trends for multiple habitats (multi-panel)
#'
#' @description
#' Generate a multi-panel figure showing LM and/or GAM trend curves
#' (with 95\% confidence ribbons) for multiple habitats.
#'
#' @param freq_long data.frame with columns:
#'   \code{habitat}, \code{distance}, \code{sampling}, \code{occurrence}.
#' @param methods Character vector, any of "lm", "gam".
#' @param degree Polynomial degree for LM (default 1).
#' @param k Basis dimension for GAM (default 5).
#' @param ncol Number of columns in the multi-panel layout.
#'
#' @return A patchwork object with multiple habitat panels.
#' @export
plot_multi_habitat_trends <- function(freq_long,
                                      methods = c("lm", "gam"),
                                      degree = 1,
                                      k = 5,
                                      ncol = 2) {
  methods <- intersect(methods, c("lm", "gam"))
  if (length(methods) == 0) {
    stop("methods must include at least one of 'lm' or 'gam'.")
  }

  habitats <- unique(freq_long$habitat)
  plots <- vector("list", length(habitats))

  for (i in seq_along(habitats)) {
    hname <- habitats[i]
    df_h  <- freq_long[freq_long$habitat == hname, , drop = FALSE]

    letter <- letters[i]
    title  <- paste0("(", letter, ") ", hname)

    plots[[i]] <- plot_single_habitat_trend(
      df_habitat   = df_h,
      habitat_name = title,
      methods      = methods,
      degree       = degree,
      k            = k
    )
  }

  patchwork::wrap_plots(plots, ncol = ncol) +
    patchwork::plot_annotation(
      theme = ggplot2::theme(
        plot.margin = ggplot2::margin(10, 10, 10, 10)
      )
    )
}
