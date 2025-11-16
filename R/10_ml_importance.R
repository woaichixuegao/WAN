#' Prepare machine-learning input data (point-based)
#'
#' @description
#' Combine occurrence and background distances to multiple habitats (and
#' optionally environmental variables) into a single \code{data.frame}
#' suitable for machine-learning models such as Random Forest.
#'
#' The output contains:
#' \itemize{
#'   \item \code{presence = 1} for occurrence points,
#'   \item \code{presence = 0} for background points.
#' }
#'
#' @param dist_occ A data.frame of distances from occurrence points to
#'   multiple habitats (one column per habitat).
#' @param dist_bg A data.frame of distances from background points to
#'   multiple habitats.
#' @param env_occ Optional data.frame of environmental variables for
#'   occurrence points (must have the same number of rows as \code{dist_occ}).
#'   Examples: mean annual temperature (MAT), mean annual precipitation (MAP).
#' @param env_bg Optional data.frame of environmental variables for
#'   background points (same number of rows as \code{dist_bg}).
#'
#' @return A combined data.frame with a \code{presence} column and all
#'   predictor variables.
#' @export
prepare_ml_data <- function(dist_occ,
                            dist_bg,
                            env_occ = NULL,
                            env_bg  = NULL) {

  occ_df <- dist_occ
  bg_df  <- dist_bg

  if (!is.null(env_occ)) {
    occ_df <- cbind(occ_df, env_occ)
  }
  if (!is.null(env_bg)) {
    bg_df <- cbind(bg_df, env_bg)
  }

  occ_df$presence <- 1L
  bg_df$presence  <- 0L

  all_df <- rbind(occ_df, bg_df)
  all_df$presence <- factor(all_df$presence, levels = c(0, 1))

  all_df
}


#' Random Forest: Ranking importance of habitat and environmental predictors
#'
#' @description
#' Fit a Random Forest model (\code{randomForest} package) using:
#'
#' \deqn{presence ~ predictor\_1 + predictor\_2 + \cdots}
#'
#' The function returns the model object and a neat importance table showing
#' MeanDecreaseGini for each predictor (including habitat distances and
#' optional climate variables such as MAT and MAP).
#'
#' @param ml_data A data.frame created by \code{prepare_ml_data()}, containing
#'   a factor \code{presence} column and multiple predictor columns.
#' @param n_trees Integer. Number of trees to grow (default = 500).
#'
#' @return A list containing:
#'   \itemize{
#'     \item \code{model}: the fitted \code{randomForest} object,
#'     \item \code{importance}: a data.frame ranked by MeanDecreaseGini.
#'   }
#'
#' @export
fit_random_forest_importance <- function(ml_data,
                                         n_trees = 500) {

  if (!"presence" %in% names(ml_data)) {
    stop("ml_data must contain a 'presence' column. Use prepare_ml_data() first.")
  }

  pred_cols <- setdiff(names(ml_data), "presence")

  formula_rf <- stats::as.formula(
    paste("presence ~", paste(pred_cols, collapse = " + "))
  )

  rf_model <- randomForest::randomForest(
    formula_rf,
    data       = ml_data,
    ntree      = n_trees,
    importance = TRUE
  )

  imp <- randomForest::importance(rf_model)

  imp_df <- data.frame(
    variable          = rownames(imp),
    MeanDecreaseGini  = imp[, "MeanDecreaseGini"],
    row.names         = NULL
  ) %>%
    dplyr::arrange(dplyr::desc(.data$MeanDecreaseGini))

  list(
    model      = rf_model,
    importance = imp_df
  )
}
