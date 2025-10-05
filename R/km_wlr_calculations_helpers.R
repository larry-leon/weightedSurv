
#' Format p-value for display
#'
#' Formats a p-value for display, showing "<eps" for small values.
#'
#' @param pval Numeric p-value.
#' @param eps Threshold for small p-values.
#' @param digits Number of digits to display.
#' @return Formatted p-value as character.
#' @export

format_pval <- function(pval, eps = 0.001, digits = 3) {
  if (is.na(pval)) return(NA)
  if (pval < eps) return(paste0("<", eps))
  format(round(pval, digits), nsmall = digits)
}

#' Validate required columns in a data frame
#'
#' Checks that all required columns are present in a data frame.
#'
#' @param df Data frame to check.
#' @param required_cols Character vector of required column names.
#' @return NULL if all columns present, otherwise error.
#' @export

validate_input <- function(df, required_cols) {
  missing <- setdiff(required_cols, names(df))
  if (length(missing) > 0) stop(paste("Missing required columns:", paste(missing, collapse = ", ")))
  invisible(NULL)
}


#' Weighted counting process
#'
#' Computes the weighted count of events up to a specified time.
#'
#' @param x Time point.
#' @param y Vector of event/censoring times.
#' @param w Weights (default 1).
#' @return Weighted count of events up to time x.
#' @export

count_weighted <- function(x, y, w = rep(1, length(y))) {
  sum(w * (y <= x))
}

#' Weighted risk set
#'
#' Computes the weighted number at risk at a specified time.
#'
#' @param x Time point.
#' @param y Vector of event/censoring times.
#' @param w Weights (default 1).
#' @return Weighted number at risk at time x.
#' @export

risk_weighted <- function(x, y, w = rep(1, length(y))) {
  sum(w * (y >= x))
}

#' Kaplan-Meier quantile calculation
#'
#' Calculates the quantile time for a Kaplan-Meier curve.
#'
#' @param time_points Vector of time points.
#' @param survival_probs Vector of survival probabilities.
#' @param qprob Quantile probability (default 0.5).
#' @param type Calculation type (midpoint or min).
#' @return Estimated quantile time.
#' @importFrom stats quantile
#' @export

kmq_calculations <- function(time_points, survival_probs, qprob = 0.5, type = "midpoint") {
  tq2 <- suppressWarnings(min(time_points[which(survival_probs <= qprob)]))
  loc.tq2 <- which(time_points == tq2)
  # jump point prior to quant
  qjt1 <- suppressWarnings(min(survival_probs[survival_probs > qprob]))
  tq1 <- suppressWarnings(min(time_points[which(survival_probs == qjt1)]))
  tq.hat <- tq2
  mid_flag <- !is.na(qjt1) && (round(qjt1, 12) == qprob)
  if (type == "midpoint" && mid_flag) {
    tq.hat <- tq1 + (tq2 - tq1) / 2
  }
  if (is.infinite(tq.hat) || is.na(tq.hat)) tq.hat <- NA
  return(tq.hat)
}

#' Kaplan-Meier quantile and confidence interval
#'
#' Calculates the quantile and confidence interval for a Kaplan-Meier curve.
#'
#' @param time_points Vector of time points.
#' @param survival_probs Vector of survival probabilities.
#' @param se_probs Standard errors of survival probabilities.
#' @param qprob Quantile probability (default 0.5).
#' @param type Calculation type (midpoint or min).
#' @param conf_level Confidence level (default 0.95).
#' @return List with quantile and confidence interval.
#' @importFrom stats quantile
#' @export

km_quantile <- function(time_points, survival_probs, se_probs = NULL, qprob = 0.5, type = c("midpoint","min"), conf_level = 0.95) {
  type <- match.arg(type)
  z <- qnorm(1 - (1 - conf_level) / 2)
  qhat <- kmq_calculations(time_points = time_points, survival_probs = survival_probs, qprob = qprob, type = type)
  qhat[is.infinite(qhat)] <- NA
  qhat_lower <- qhat_upper <- NA
  if (!is.null(se_probs)) {
    # log transform for CI
    lower_probs <- exp(log(survival_probs) - z * se_probs / survival_probs)
    upper_probs <- exp(log(survival_probs) + z * se_probs / survival_probs)
    qhat_lower <- kmq_calculations(time_points = time_points, survival_probs = lower_probs, qprob = qprob, type = type)
    qhat_lower[is.infinite(qhat_lower)] <- NA
    qhat_upper <- kmq_calculations(time_points = time_points, survival_probs = upper_probs, qprob = qprob, type = type)
    qhat_upper[is.infinite(qhat_upper)] <- NA
  }
  return(list(qhat = qhat, lower = qhat_lower, upper = qhat_upper))
}

#' Table of KM quantiles for two groups
#'
#' Returns a data frame of quantiles and confidence intervals for two groups.
#'
#' @param time_points Vector of time points.
#' @param surv0 Survival probabilities for group 0.
#' @param se0 Standard errors for group 0.
#' @param surv1 Survival probabilities for group 1.
#' @param se1 Standard errors for group 1.
#' @param arms Group labels.
#' @param qprob Quantile probability.
#' @param type Calculation type.
#' @param conf_level Confidence level.
#' @return Data frame of quantiles and CIs for each group.
#' @export

km_quantile_table <- function(time_points, surv0, se0, surv1, se1, arms = c("treat", "control"), qprob = 0.5, type = c("midpoint","min"), conf_level = 0.95) {
  type <- match.arg(type)
  kmq0 <- km_quantile(time_points = time_points, survival_probs = surv0, se_probs = se0, qprob = qprob, type = type, conf_level = conf_level)
  df0 <- data.frame(group = arms[2], quantile = kmq0$qhat, lower = kmq0$lower, upper = kmq0$upper)
  kmq1 <- km_quantile(time_points = time_points, survival_probs = surv1, se_probs = se1, qprob = qprob, type = type, conf_level = conf_level)
  df1 <- data.frame(group = arms[1], quantile = kmq1$qhat, lower = kmq1$lower, upper = kmq1$upper)
  quantiles_df <- rbind(df1, df0)
  return(quantiles_df)
}


#' Weighted Log-Rank and Difference Estimate at a Specified Time
#'
#' Computes the weighted log-rank statistic, its variance, the difference in survival at a specified time (`tzero`),
#' the variance of the difference, their covariance, and correlation, using flexible time-dependent weights.
#' The weighting scheme is selected via the \code{scheme} argument and is calculated using \code{wt.rg.S}.
#'
#' @param dfcounting List output from \code{df_counting} containing risk sets, event counts, and survival estimates.
#' @param rho Numeric weighting parameter (used for \"fh\" and \"custom_code\" schemes).
#' @param gamma Numeric weighting parameter (used for \"fh\" and \"custom_code\" schemes).
#' @param tzero Time point at which to evaluate the difference in survival (default: 24).
#' @param scheme Character string specifying weighting scheme. One of:
#'   \"fh\" (Fleming-Harrington), \"schemper\", \"XO\", \"MB\", \"custom_time\", or \"custom_code\".
#' @param Scensor Optional numeric vector of censoring probabilities (for \"schemper\" scheme).
#' @param Ybar Optional numeric vector of risk set sizes (for \"XO\" scheme).
#' @param tpoints Optional numeric vector of time points (for \"custom_time\" and \"MB\" schemes).
#' @param t.tau Optional time cutoff for \"custom_time\" weights.
#' @param w0.tau, w1.tau Weights before/after t.tau (for \"custom_time\" scheme).
#' @param mb_tstar Optional time for MB weights.
#'
#' @return A list with elements:
#'   \item{lr}{Weighted log-rank test statistic.}
#'   \item{sig2_lr}{Variance of the log-rank statistic.}
#'   \item{dhat}{Difference in survival at \code{tzero}.}
#'   \item{cov_wlr_dhat}{Covariance between log-rank and difference at \code{tzero}.}
#'   \item{sig2_dhat}{Variance of the difference at \code{tzero}.}
#'   \item{cor_wlr_dhat}{Correlation between log-rank and difference at \code{tzero}.}
#'
#' @details
#' The weighting scheme is selected via the \code{scheme} argument and calculated using \code{wt.rg.S}.
#' Supports standard Fleming-Harrington, Schemper, XO (Xu & O'Quigley), MB (Maggir-Burman), custom time-based, and custom code weights.
#'
#' @export

wlr_dhat_estimates <- function(dfcounting,
                               scheme = "fh", scheme_params = list(rho = 0, gamma = 0), tzero = NULL
) {
  # Check for required columns
  required_cols <- c("at_points", "nbar0", "nbar1", "ybar0", "ybar1", "survP")
  if(scheme == "schemper") requires_cols <- c(required_cols,"survG")
  missing_cols <- setdiff(required_cols, names(dfcounting))
  if (length(missing_cols) > 0) {
    stop("df_weights is missing required columns: ", paste(missing_cols, collapse = ", "))
  }

  at_points <- dfcounting$at_points
  nbar0 <- dfcounting$nbar0
  nbar1 <- dfcounting$nbar1
  ybar0 <- dfcounting$ybar0
  ybar1 <- dfcounting$ybar1
  S.pool <- dfcounting$survP
  G.pool <- dfcounting$survG

  # weights
  dfwlr <- data.frame(at_points, S.pool, G.pool, Ybar = ybar0 + ybar1)

  wt_rg <- get_validated_weights(dfwlr, scheme, scheme_params)

  S1 <- dfcounting$surv1
  S0 <- dfcounting$surv0

  dN0 <- diff(c(0, nbar0))
  dN1 <- diff(c(0, nbar1))
  dN <- dN0 + dN1
  ybar <- ybar0 + ybar1

  if(!is.null(tzero)){
    loc_tzero <- which.max(at_points > tzero)
    if (at_points[loc_tzero] <= tzero & at_points[loc_tzero + 1] > tzero) {
      dhat_tzero <- S1[loc_tzero] - S0[loc_tzero]
    } else {
      dhat_tzero <- S1[loc_tzero - 1] - S0[loc_tzero - 1]
    }
    Sp_tzero <- S.pool[loc_tzero]
  }

  K <- ifelse(ybar > 0, wt_rg * (ybar0 * ybar1) / ybar, 0.0)
  drisk0 <- sum(ifelse(ybar0 > 0, (K / ybar0) * dN0, 0.0))
  drisk1 <- sum(ifelse(ybar1 > 0, (K / ybar1) * dN1, 0.0))
  lr <- drisk0 - drisk1
  h0 <- ifelse(ybar0 == 0, 0, (K^2 / ybar0))
  h1 <- ifelse(ybar1 == 0, 0, (K^2 / ybar1))
  dJ <- ifelse(ybar == 1, 0, (dN - 1) / (ybar - 1))
  dL <- ifelse(ybar == 0, 0, dN / ybar)
  sig2_lr <- sum((h0 + h1) * (1 - dJ) * dL)

  if(!is.null(tzero)){
    w_tzero <- wt_rg * ifelse(at_points <= tzero, 1, 0)
    w_integral_t0 <- sum(w_tzero * (1 - dJ) * dL)
    cov_wlr_dhat <- Sp_tzero * w_integral_t0
    h2 <- ifelse(ybar0 * ybar1 > 0, (ybar / (ybar0 * ybar1)), 0)
    h2 <- h2 * ifelse(at_points <= tzero, 1, 0)
    sig2_dhat <- (Sp_tzero^2) * sum(h2 * (1 - dJ) * dL)
    # Correlation between log-rank statistic and dhat at time tzero
    cor_wlr_dhat <- cov_wlr_dhat / (sqrt(sig2_lr) * sqrt(sig2_dhat))

    ans <- list(
      score = lr, sig2.score = sig2_lr, z.score = lr / sqrt(sig2_lr), dhat = dhat_tzero,
      cov_wlr_dhat = cov_wlr_dhat, sig2_dhat = sig2_dhat, cor_wlr_dhat = cor_wlr_dhat
    )
  } else{
    ans <- list(
      score = lr, sig2.score  = sig2_lr, z.score = lr / sqrt(sig2_lr)
    )
  }
  return(ans)
}




#' Kaplan-Meier Survival Estimates and Variance
#'
#' Computes Kaplan-Meier survival estimates and their variances given risk and event counts.
#'
#' @param ybar Vector of risk set sizes at each time point.
#' @param nbar Vector of event counts at each time point.
#' @param sig2w_multiplier Optional vector for variance calculation. If NULL, calculated internally.
#' @return List with survival estimates and variances.
#' @export

KM_estimates <- function(ybar, nbar, sig2w_multiplier = NULL){
  dN <- diff(c(0, nbar))
  dN_risk <- ifelse(ybar > 0, dN / ybar, 0.0)
  S_KM <- cumprod(1 - dN_risk)
  if(is.null(sig2w_multiplier)){
    sig2w_multiplier  <- ifelse(ybar > 0 & ybar > dN, dN / (ybar * (ybar - dN)), 0.0)
  }
  var_KM <- (S_KM^2) * cumsum(sig2w_multiplier)
  list(S_KM = S_KM, sig2_KM = var_KM)
}

#' Event and Risk Matrices for Survival Analysis
#'
#' Constructs matrices indicating event and risk status for each subject at specified time points.
#'
#' @param U Vector of observed times (e.g., time-to-event).
#' @param at.points Vector of time points at which to evaluate events and risk.
#' @return A list with event and risk matrices.
#' @export

get_event_risk_matrices <- function(U, at.points) {
  event_mat <- outer(U, at.points, FUN = "<=")
  risk_mat  <- outer(U, at.points, FUN = ">=")
  list(event_mat = event_mat, risk_mat = risk_mat)
}



#' Resampling Survival Curves for Confidence Bands
#'
#' Performs resampling to generate survival curves for a group, used for constructing confidence bands.
#'
#' @param U Vector of observed times.
#' @param W Vector of weights.
#' @param D Vector of event indicators (0/1).
#' @param at.points Vector of time points for evaluation.
#' @param draws.band Number of resampling draws.
#' @param surv Vector of survival estimates.
#' @param G_draws Matrix of random draws for resampling.
#' @return Matrix of resampled survival curves.
#' @export

resampling_survival <- function(U, W, D, at.points, draws.band, surv, G_draws) {
  mats <- get_event_risk_matrices(U, at.points)
  risk_w <- colSums(mats$risk_mat * W)
  counting_star_all <- t(mats$event_mat * W) %*% (D * G_draws)
  dN_star_all <- apply(counting_star_all, 2, function(x) diff(c(0, x)))
  drisk_star <- sweep(dN_star_all, 1, risk_w, "/")
  drisk_star[is.infinite(drisk_star) | is.nan(drisk_star)] <- 0
  surv_star <- (-1) * surv * apply(drisk_star, 2, cumsum)
  return(surv_star)
}

#' KM difference between groups
#'
#' Calculates the difference in KM curves between two groups, with confidence intervals and optional resampling bands.
#'
#' @param df Data frame with survival data.
#' @param tte.name Name of time-to-event variable.
#' @param event.name Name of event indicator variable.
#' @param treat.name Name of treatment group variable.
#' @param weight.name Name of weights variable.
#' @param at_points Time points for calculation.
#' @param alpha Significance level.
#' @param seedstart Random seed.
#' @param draws Number of resampling draws.
#' @param risk.points Risk points for calculation.
#' @param draws.band Number of draws for simultaneous bands.
#' @param tau.seq Step size for tau sequence.
#' @param qtau Quantile for tau range.
#' @param show_resamples Logical; whether to plot resamples.
#' @param modify_tau Logical; restrict time range for bands.
#' @return Data frame with KM differences and confidence intervals.
#' @importFrom survival Surv
#' @importFrom stats quantile
#' @importFrom graphics matplot
#' @export

KM_diff <- function(df, tte.name = "tte", event.name = "event", treat.name = "treat", weight.name=NULL, at_points = sort(df[[tte.name]]), alpha = 0.05, seedstart = 8316951, draws = 0,
                    risk.points, draws.band = 0, tau.seq = 0.25, qtau = 0.025, show_resamples = TRUE, modify_tau = FALSE) {

  required_cols <- c(tte.name, event.name, treat.name)
  missing_cols <- setdiff(required_cols, names(df))
  if (length(missing_cols) > 0) {
    stop(paste("Missing required columns in df:", paste(missing_cols, collapse = ", ")))
  }

  # Check treatment coding
  if (!all(df[[treat.name]] %in% c(0, 1))) {
    stop("Treatment must be numerical indicator: 0=control, 1=experimental")
  }

  # Check event coding
  if (!all(df[[event.name]] %in% c(0, 1))) {
    stop("Event must be binary (0/1).")
  }

  # Check for NA values
  if (any(is.na(df[[tte.name]]))) warning("NA values found in time-to-event column.")
  if (any(is.na(df[[event.name]]))) warning("NA values found in event column.")
  if (any(is.na(df[[treat.name]]))) warning("NA values found in treatment column.")

  tfixed <- aeqSurv(Surv(df[[tte.name]],df[[event.name]]))
  time<- tfixed[,"time"]
  delta <- tfixed[,"status"]
  z <- df[[treat.name]]
  wgt <- if (!is.null(weight.name)) df[[weight.name]] else rep(1, length(time))

  if (!all(z %in% c(0, 1))) stop("Treatment must be numerical indicator: 0=control, 1=experimental")

  if (is.unsorted(time)) {
    ord <- order(time)
    time <- time[ord]
    delta <- delta[ord]
    z <- z[ord]
    wgt <- wgt[ord]
  }

  # For simultaneous bands restrict time range if modify_tau (otherwise align with RMST 'max(tau)')
  if(draws.band > 0 && modify_tau){
    taus <- quantile(time[delta ==1], c(qtau, 1-qtau))
    at_points<-seq(taus[1],taus[2], by =tau.seq)
    riskp <- risk.points[which(risk.points <= taus[2])]
    at_points <- sort(unique(c(at_points, riskp)))
  }

  # Treatment arm
  group_data <- extract_group_data(time, delta, wgt, z, group = 1)
  risk_event <- calculate_risk_event_counts(group_data$U, group_data$D, group_data$W, at_points = at_points, draws = draws, seedstart = seedstart)
  temp <- KM_estimates(ybar = risk_event$ybar, nbar = risk_event$nbar, sig2w_multiplier = risk_event$sig2w_multiplier)
  risk_event1 <- risk_event
  group_data1 <- group_data
  surv1 <- temp$S_KM
  sig2_surv1 <- temp$sig2_KM

  # Treatment arm
  group_data <- extract_group_data(time, delta, wgt, z, group = 0)
  risk_event <- calculate_risk_event_counts(group_data$U, group_data$D, group_data$W, at_points = at_points, draws = draws, seedstart = seedstart)
  temp <- KM_estimates(ybar = risk_event$ybar, nbar = risk_event$nbar, sig2w_multiplier = risk_event$sig2w_multiplier)
  risk_event0 <- risk_event
  group_data0 <- group_data
  surv0 <- temp$S_KM
  sig2_surv0 <- temp$sig2_KM

  dhat <- surv1 - surv0
  sig2_dhat <- sig2_surv0 + sig2_surv1

  surv0_star <- surv1_star <- dhat_star <- NULL
  c_alpha_band <- sb_lower <- sb_upper <- NULL
  Zdhat_star <- NULL
  if(draws.band > 0){
    # draws > 0
    set.seed(seedstart)

    # Control group
    n0 <- length(group_data0$U)
    G0.draws <- matrix(rnorm(draws.band * n0), ncol = draws.band)
    surv0_star <- resampling_survival(group_data0$U, group_data0$W, group_data0$D, at_points, draws.band, surv0, G0.draws)
    # Treatment group
    n1 <- length(group_data1$U)
    G1.draws <- matrix(rnorm(draws.band * n1), ncol = draws.band)
    surv1_star <- resampling_survival(group_data1$U, group_data1$W, group_data1$D, at_points, draws.band, surv1, G1.draws)

    dhat_star <- (surv1_star - surv0_star)

    Zdhat_star <- dhat_star / sqrt(sig2_dhat)

    # simultaneous band
    sups <- apply(abs(Zdhat_star), 2, max, na.rm = TRUE)
    c_alpha_band <- quantile(sups,c(0.95))
    # Show first 20
    if(show_resamples){
      matplot(at_points, Zdhat_star[,c(1:20)], type="s", lty=2, lwd = 1, xlab="time", ylab = "Normalized survival differences (1st 20)",
              main = sprintf("c_alpha (simul. band): %.2f", c_alpha_band)
      )
    }
    # simulataneous band
    sb_lower <- dhat - c_alpha_band * sqrt(sig2_dhat)
    sb_upper <- dhat + c_alpha_band * sqrt(sig2_dhat)
  }
  # Standard point-wise CIs
  c_alpha <- qnorm(1 - alpha / 2)
  lower <- dhat - c_alpha * sqrt(sig2_dhat)
  upper <- dhat + c_alpha * sqrt(sig2_dhat)
  list(
    at_points = at_points, surv0 = surv0, sig2_surv0 = sig2_surv0,
    surv1 = surv1, sig2_surv1 = sig2_surv1, dhat = dhat, sig2_dhat = sig2_dhat,
    lower = lower, upper = upper,
    dhat_star = dhat_star,
    Zdhat_star = Zdhat_star,
    surv0_star = surv0_star, surv1_star = surv1_star,
    c_alpha_band = c_alpha_band, sb_lower = sb_lower, sb_upper = sb_upper
  )
}


#' Weighted log-rank cumulative statistics
#'
#' Calculates cumulative weighted log-rank statistics for survival data.
#'
#' @param df_weights Data frame with weights and risk/event counts.
#' @param scheme Weighting scheme.
#' @param scheme_params List of scheme parameters.
#' @param return_cumulative Logical; whether to return cumulative statistics.
#' @return List with score and variance.
#' @export

wlr_cumulative <- function(df_weights, scheme, scheme_params = list(rho = 0, gamma = 0), return_cumulative = FALSE) {
# Check for required columns
  required_cols <- c("at_points", "nbar0", "nbar1", "ybar0", "ybar1", "survP")
  if(scheme == "schemper") requires_cols <- c(required_cols,"survG")
  missing_cols <- setdiff(required_cols, names(df_weights))
  if (length(missing_cols) > 0) {
    stop("df_weights is missing required columns: ", paste(missing_cols, collapse = ", "))
  }

  at_points <- df_weights$at_points
  nbar0 <- df_weights$nbar0
  nbar1 <- df_weights$nbar1
  ybar0 <- df_weights$ybar0
  ybar1 <- df_weights$ybar1
  S.pool <- df_weights$survP
  G.pool <- df_weights$survG

  # weights
  dfwlr <- data.frame(at_points, S.pool, G.pool, Ybar = ybar0 + ybar1)
  wt_rg <- get_validated_weights(dfwlr, scheme, scheme_params)

  dN0 <- diff(c(0, nbar0))
  dN1 <- diff(c(0, nbar1))
  dN <- dN0 + dN1
  ybar <- ybar0 + ybar1
  K <- ifelse(ybar > 0, wt_rg * (ybar0 * ybar1) / ybar, 0.0)
  h0 <- ifelse(ybar0 == 0, 0, (K^2 / ybar0))
  h1 <- ifelse(ybar1 == 0, 0, (K^2 / ybar1))
  dJ <- ifelse(ybar == 1, 0, (dN - 1) / (ybar - 1))
  dL <- ifelse(ybar == 0, 0, dN / ybar)

  if(return_cumulative){
    drisk0 <- cumsum(ifelse(ybar0 > 0, (K / ybar0) * dN0, 0.0))
    drisk1 <- cumsum(ifelse(ybar1 > 0, (K / ybar1) * dN1, 0.0))
    sig2_lr <- cumsum((h0 + h1) * (1 - dJ) * dL)
  } else {
    drisk0 <- sum(ifelse(ybar0 > 0, (K / ybar0) * dN0, 0.0))
    drisk1 <- sum(ifelse(ybar1 > 0, (K / ybar1) * dN1, 0.0))
    sig2_lr <- sum((h0 + h1) * (1 - dJ) * dL)
  }
  lr <- drisk0 - drisk1
  z <- lr / sqrt(sig2_lr)
  list(z.score = z, time = at_points, score = lr, sig2.score = sig2_lr)
}

