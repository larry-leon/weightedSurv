# ---- Utility Functions ----
#' Safe execution wrapper
#' @export

safe_run <- function(expr) {
  tryCatch(expr, error = function(e) {
    message("Error: ", e$message)
    NULL
  })
}


#' Get df_counting
#' @export

get_dfcounting <- function(df, tte.name, event.name, treat.name, arms, by.risk=12, cox.digits=3, lr.digits=3,
                           qprob=0.50, strata.name=NULL, weight.name=NULL, check.KM = TRUE, scheme = "fh", scheme_params = list(rho = 0, gamma = 0), draws = 0, seedstart = 8316951, check.seKM = FALSE) {
  safe_run({
    dfcount <- df_counting(
      df=df, tte.name=tte.name, event.name=event.name, treat.name=treat.name,
      arms=arms, by.risk=by.risk, cox.digits=cox.digits, lr.digits=lr.digits,
      qprob=qprob, strata.name=strata.name, weight.name=weight.name, check.KM = check.KM, check.seKM = check.seKM,
      draws = draws, seedstart=seedstart, scheme = scheme, scheme_params = scheme_params
    )
    return(dfcount)
  })
}

#' Checking results
#' @export

check_results <- function(dfcount){
  zlr_sq  <- with(dfcount,lr^2/sig2_lr)
  zCox_sq <-  with(dfcount,z.score^2)
  cat(sprintf("zlr_sq=%.6f, logrank=%.6f, zCox_sq=%.6f\n", zlr_sq, dfcount$logrank_results$chisq, zCox_sq))
}


#' Plot Kaplan-Meier curves
#' @export

plot_km <- function(df, tte.name, event.name, treat.name, weights=NULL, ...) {
  safe_run({
    surv_obj <- Surv(df[[tte.name]], df[[event.name]])
    formula <- as.formula(paste("surv_obj ~", treat.name))
    if (!is.null(weights)) {
      km_fit <- survfit(formula, data=df, weights=df[[weights]])
    } else {
      km_fit <- survfit(formula, data=df)
    }
    plot(km_fit, mark.time=TRUE, ...)
    invisible(km_fit)
  })
}

#' Plot weighted KM using custom function
#' @export

plot_weighted_km <- function(dfcount, ...) {
  safe_run({
    KM_plot_2sample_weighted_counting(dfcount=dfcount, ...
    )
  })
}


# Helper: Extract Group Data
#' Extract time, event, and weight data for a group
#'
#' @param time Numeric vector of times
#' @param delta Numeric vector of event indicators (1=event, 0=censored)
#' @param wgt Numeric vector of weights
#' @param z Numeric vector of group indicators
#' @param group Value of group to extract (default 1)
#' @return List with U (times), D (events), W (weights)
#' @export

extract_group_data <- function(time, delta, wgt, z, group = 1) {
  list(
    U = time[z == group],
    D = delta[z == group],
    W = wgt[z == group]
  )
}

# Helper: Calculate Risk and Event Counts
#' Calculate risk set and event counts at time points
#'
#' @param U Numeric vector of times for group
#' @param D Numeric vector of event indicators for group
#' @param W Numeric vector of weights for group
#' @param at_points Numeric vector of time points
#' @return List with ybar (risk set counts), nbar (event counts)
#' @export

calculate_risk_event_counts <- function(U, D, W, at_points, draws = 0, seedstart = 816951) {
  ybar <- colSums(outer(U, at_points, FUN = ">=") * W)
  nbar <- colSums(outer(U[D == 1], at_points, FUN = "<=") * W[D == 1])
  dN <- diff(c(0, nbar))
  # For un-weighted (all weights equal) return standard variance term
  if(length(unique(W)) == 1){
    # Greenwood
    sig2w_multiplier <- ifelse(ybar > 0 & ybar > dN, dN / (ybar * (ybar-dN)), 0.0)
    # Alternative with dN / (ybar^2)
  }
  if(length(unique(W)) > 1){
    n <- length(U)
    event_mat <- outer(U, at_points, FUN = "<=")
    risk_mat  <- outer(U, at_points, FUN = ">=")
    risk_w <- colSums(risk_mat *  W)
    result <- switch(
      as.integer(draws > 0) + 1,
      {
        # draws == 0
        counting <- colSums(event_mat * (D * W))
        dN_w <- diff(c(0, counting))
        dJ <- ifelse(risk_w == 1, 0, (dN_w - 1) / (risk_w - 1))
        dL <- ifelse(risk_w == 0, 0, dN_w / risk_w)
        h2 <- ifelse(risk_w > 0, (1 / (risk_w)), 0)
        sig2w_multiplier <- (h2 * (dN_w - dL))^2
        sig2w_multiplier
      },
      {
        # draws > 0
        set.seed(seedstart)
        G.draws <- matrix(rnorm(draws * n), ncol = draws)
        counting_star_all <- t(event_mat * W) %*% (D * G.draws)
        dN_star_all <- apply(counting_star_all, 2, function(x) diff(c(0, x)))
        drisk_star <- sweep(dN_star_all, 1, risk_w, "/")
        drisk_star[is.infinite(drisk_star) | is.nan(drisk_star)] <- 0
        sig2w_multiplier <- apply(drisk_star, 1, var)
        sig2w_multiplier
      }
    )
  }
  list(ybar = ybar, nbar = nbar, sig2w_multiplier = sig2w_multiplier)
}

# Helper : Get Censoring and Event Times
#' Get censoring and event times and their indices
#'
#' @param time Numeric vector of times
#' @param delta Numeric vector of event indicators
#' @param z Numeric vector of group indicators
#' @param group Value of group to extract
#' @param censoring_allmarks Logical; if FALSE, remove events from censored
#' @param at_points Numeric vector of time points
#' @return List with cens (censored times), ev (event times), idx_cens, idx_ev, idx_ev_full
#' @export

get_censoring_and_events <- function(time, delta, z, group, censoring_allmarks, at_points) {
  cens <- time[z == group & delta == 0]
  ev <- sort(unique(time[z == group & delta == 1]))
  if (!censoring_allmarks) cens <- setdiff(cens, ev)
  idx_cens <- match(cens, at_points)
  idx_ev <- match(ev, at_points)
  ev <- c(ev, max(time[z == group]))
  idx_ev_full <- match(ev, at_points)
  list(
    cens = cens,
    ev = ev,
    idx_cens = idx_cens,
    idx_ev = idx_ev,
    idx_ev_full = idx_ev_full
  )
}

# Helper : Get Risk Points
#' Get risk set counts at specified risk points
#'
#' @param ybar Numeric vector of risk set counts
#' @param risk_points Numeric vector of risk points
#' @param at_points Numeric vector of time points
#' @return Numeric vector of risk set counts at risk points
#' @export

get_riskpoints <- function(ybar, risk_points, at_points) {
  ybar[match(risk_points, at_points)]
}
