#' Exponential decay model
#'
#' @param x Numeric vector of distances.
#' @param N0 Intercept parameter (frequency at distance 0).
#' @param b Decay rate parameter.
#'
#' @return Numeric vector of predicted frequency values.
#' @export
exponential_decay <- function(x, N0, b) {
  N0 * exp(-b * x)
}

#' Fit an exponential decay model
#'
#' @description
#' Fit an exponential decay model describing frequency as a function of
#' distance using non-linear least squares (\code{nlsLM}).
#'
#' The model has the form:
#' \deqn{frequency = N0 * exp(-b * distance)}
#'
#' @param distance Numeric vector of distances.
#' @param frequency Numeric vector of frequencies (percentages, \code{0–100}).
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{params}: estimated parameters \code{N0} and \code{b}.
#'   \item \code{r_squared}: coefficient of determination (\eqn{R^2}).
#'   \item \code{fit}: the fitted \code{nlsLM} model object.
#'   \item \code{data}: cleaned data used in the fit.
#'   \item \code{success}: logical, whether the fit succeeded.
#'   \item \code{message}: error message (if the fit failed).
#' }
#' @export
fit_exponential_decay <- function(distance, frequency) {
  valid_data <- data.frame(
    distance  = distance,
    frequency = frequency
  ) %>%
    dplyr::filter(!is.na(frequency), !is.na(distance))

  if (nrow(valid_data) < 3) {
    return(list(success = FALSE, message = "Not enough data points to fit the model."))
  }

  res <- tryCatch({
    fit <- minpack.lm::nlsLM(
      frequency ~ exponential_decay(distance, N0, b),
      data = valid_data,
      start = list(N0 = max(valid_data$frequency), b = 0.01),
      control = minpack.lm::nls.lm.control(maxiter = 500)
    )

    y_pred   <- stats::predict(fit)
    ss_res   <- sum((valid_data$frequency - y_pred)^2)
    ss_tot   <- sum((valid_data$frequency - mean(valid_data$frequency))^2)
    r_squared <- 1 - (ss_res / ss_tot)

    list(
      params   = stats::coef(fit),
      r_squared = r_squared,
      fit      = fit,
      data     = valid_data,
      success  = TRUE
    )
  }, error = function(e) {
    list(success = FALSE, message = paste("Model fitting failed:", e$message))
  })

  res
}

#' Calculate confidence intervals for the exponential decay fit
#'
#' @param fit_result A result list returned by \code{fit_exponential_decay()}.
#' @param x_range Numeric vector of distances at which predictions should be made.
#'
#' @return A data.frame with columns:
#'   \code{distance}, \code{fit}, \code{upper}, \code{lower};
#'   or \code{NULL} if the fit was not successful.
#' @export
calculate_confidence_intervals <- function(fit_result, x_range) {
  if (!fit_result$success) return(NULL)

  x_pred <- x_range
  y_pred <- exponential_decay(
    x_pred,
    fit_result$params["N0"],
    fit_result$params["b"]
  )

  residuals <- fit_result$data$frequency - stats::predict(fit_result$fit)
  se        <- stats::sd(residuals)

  upper <- y_pred + 1.96 * se
  lower <- pmax(y_pred - 1.96 * se, 0)

  data.frame(
    distance = x_pred,
    fit      = y_pred,
    upper    = upper,
    lower    = lower
  )
}

#' Fit a linear model: freq = a + c * distance
#'
#' @param distance Numeric vector of distances.
#' @param frequency Numeric vector of frequencies.
#'
#' @return A list with elements \code{model}, \code{fit}, \code{AIC}, and \code{data},
#'   or \code{NULL} if there are fewer than 3 valid points.
#' @export
fit_linear_model <- function(distance, frequency) {
  df <- data.frame(distance = distance, frequency = frequency) %>%
    dplyr::filter(!is.na(distance), !is.na(frequency))

  if (nrow(df) < 3) return(NULL)

  fit <- stats::lm(frequency ~ distance, data = df)

  list(
    model = "linear",
    fit   = fit,
    AIC   = stats::AIC(fit),
    data  = df
  )
}

#' Fit a power-law model: freq = a * distance^(-b)
#'
#' @param distance Numeric vector of distances.
#' @param frequency Numeric vector of frequencies.
#'
#' @return A list with elements \code{model}, \code{fit}, \code{AIC},
#'   \code{params}, and \code{data}, or \code{NULL} if fewer than 3 valid
#'   positive points are available.
#' @export
fit_power_model <- function(distance, frequency) {
  df <- data.frame(distance = distance, frequency = frequency) %>%
    dplyr::filter(!is.na(distance), !is.na(frequency),
                  distance > 0, frequency > 0)

  if (nrow(df) < 3) return(NULL)

  fit <- stats::lm(log(frequency) ~ log(distance), data = df)
  a_hat <- exp(stats::coef(fit)[1])
  b_hat <- -stats::coef(fit)[2]

  list(
    model  = "power",
    fit    = fit,
    AIC    = stats::AIC(fit),
    params = c(a = a_hat, b = b_hat),
    data   = df
  )
}

#' Fit a Gaussian model: freq = a * exp(-(distance / b)^2)
#'
#' @param distance Numeric vector of distances.
#' @param frequency Numeric vector of frequencies.
#'
#' @return A list with elements \code{model}, \code{fit}, \code{AIC},
#'   \code{params}, and \code{data}, or \code{NULL} if fitting fails
#'   or there are fewer than 3 valid points.
#' @export
fit_gaussian_model <- function(distance, frequency) {
  df <- data.frame(distance = distance, frequency = frequency) %>%
    dplyr::filter(!is.na(distance), !is.na(frequency))

  if (nrow(df) < 3) return(NULL)

  res <- tryCatch({
    fit <- minpack.lm::nlsLM(
      frequency ~ a * exp(-(distance / b)^2),
      data = df,
      start = list(
        a = max(df$frequency, na.rm = TRUE),
        b = diff(range(df$distance)) / 3
      ),
      control = minpack.lm::nls.lm.control(maxiter = 500)
    )
    list(
      model  = "gaussian",
      fit    = fit,
      AIC    = stats::AIC(fit),
      params = stats::coef(fit),
      data   = df
    )
  }, error = function(e) NULL)

  res
}

#' Compare multiple decay models using AIC
#'
#' @description
#' Fit several candidate models (exponential, linear, power-law, Gaussian)
#' to the same distance–frequency data and compare them using the Akaike
#' Information Criterion (AIC).
#'
#' @param distance Numeric vector of distances.
#' @param frequency Numeric vector of frequencies.
#'
#' @return A list with elements:
#' \itemize{
#'   \item \code{details}: a named list of model fit results.
#'   \item \code{AIC_table}: a data.frame with columns \code{model} and \code{AIC},
#'         sorted by AIC (ascending).
#'   \item \code{best_model}: character string with the name of the model
#'         that has the lowest AIC.
#' }
#'   Returns \code{NULL} if no model could be fitted.
#' @export
compare_models_AIC <- function(distance, frequency) {
  res <- list()

  # Exponential decay model
  exp_fit <- fit_exponential_decay(distance, frequency)
  if (exp_fit$success) {
    res[["exponential"]] <- list(
      model = "exponential",
      fit   = exp_fit$fit,
      AIC   = stats::AIC(exp_fit$fit),
      extra = exp_fit
    )
  }

  # Linear model
  lin <- fit_linear_model(distance, frequency)
  if (!is.null(lin)) res[["linear"]] <- lin

  # Power-law model
  pow <- fit_power_model(distance, frequency)
  if (!is.null(pow)) res[["power"]] <- pow

  # Gaussian model
  gau <- fit_gaussian_model(distance, frequency)
  if (!is.null(gau)) res[["gaussian"]] <- gau

  if (length(res) == 0) return(NULL)

  # Construct AIC comparison table
  aic_df <- data.frame(
    model = sapply(res, function(x) x$model),
    AIC   = sapply(res, function(x) x$AIC)
  ) %>%
    dplyr::arrange(.data$AIC)

  list(
    details    = res,
    AIC_table  = aic_df,
    best_model = aic_df$model[1]
  )
}
