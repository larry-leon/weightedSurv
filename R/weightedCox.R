#' Confidence interval for Cox model estimate
#'
#' Calculates the confidence interval for a Cox model log hazard ratio estimate.
#'
#' @param bhat Estimated log hazard ratio.
#' @param sig_bhat Standard error of the estimate.
#' @param alpha Significance level (default: 0.05).
#' @param verbose Logical; if TRUE, prints interval.
#' @return Data frame with log hazard ratio, standard error, hazard ratio, lower and upper confidence limits.
#' @importFrom stats qnorm
#' @export

ci_cox  <- function(bhat, sig_bhat, alpha = 0.05, verbose = FALSE) {
  z <- qnorm(1 - alpha / 2)
  bhat_lower <- bhat - z * sig_bhat
  bhat_upper <- bhat + z * sig_bhat
  hr <- exp(bhat)
  lower <- exp(bhat_lower)
  upper <- exp(bhat_upper)
  if (verbose) {
    cat(sprintf("Hazard Ratio (HR): %.3f\n", hr))
    cat(sprintf("95%% CI: [%.3f, %.3f]\n", lower, upper))
  }
  result <- data.frame(
    beta = bhat,
    sig_bhat = sig_bhat,
    hr = hr,
    lower = lower,
    upper = upper
  )
  return(result)
}



#' Weighted event count for Cox model
#'
#' Calculates the weighted event count for the Cox model with rho-gamma weights.
#'
#' @param x Numeric vector of time points.
#' @param error Numeric vector of error terms.
#' @param delta Numeric vector of event indicators.
#' @param weight Numeric vector of weights (default: 1).
#' @return Numeric value of weighted event count.
#' @export

N_rhogamma <- function(x, error, delta, weight = 1) {
  sum(weight * delta * (error <= x))
}

#' Root-finding for Cox score function
#'
#' Finds the root of the Cox model score equation for weighted log-rank statistics.
#'
#' @param time Numeric vector of event times.
#' @param delta Numeric vector of event indicators.
#' @param z Numeric vector of group indicators.
#' @param w_hat Numeric vector of weights.
#' @param wt_rg Numeric vector of rho-gamma weights.
#' @return List with root and additional information, or NA if not found.
#' @importFrom stats uniroot
#' @export

find_cox_root <- function(time, delta, z, w_hat, wt_rg) {
  tryCatch(
    uniroot(f = cox_score_rhogamma, interval = c(-15, 15), extendInt = "yes", tol = 1e-10,
            time = time, delta = delta, z = z, w_hat = w_hat, wt_rg = wt_rg),
    error = function(e) NA
  )
}




#' Score calculation for weighted Cox model
#'
#' Calculates the score and variance for the weighted Cox model.
#'
#' @param ybar1 Numeric vector of event counts for group 1.
#' @param ybar0 Numeric vector of event counts for group 0.
#' @param dN1 Numeric vector of event increments for group 1.
#' @param dN0 Numeric vector of event increments for group 0.
#' @param wt_rg Numeric vector of rho-gamma weights.
#' @return List with score, variance, information, and weights.
#' @export

score_calculation <- function(ybar1, ybar0, dN1, dN0, wt_rg){
dN <- dN1 + dN0
num <- wt_rg * ybar1 * ybar0
den <- ybar0 + ybar1
K <- ifelse(den > 0, num / den, 0.0)

drisk1 <- ifelse(ybar1 > 0, dN1 / ybar1, 0.0)
drisk0 <- ifelse(ybar0 > 0, dN0 / ybar0, 0.0)
score <- sum(K * (drisk0 - drisk1))

h1 <- ifelse(ybar1 > 0, (K^2 / ybar1), 0.0)
h2 <- ifelse(ybar0 > 0, (K^2 / ybar0), 0.0)
temp <- c(den - 1)
ybar_mod <- ifelse(temp < 1, 1, temp)
dH1 <- ifelse(ybar_mod > 0, (dN-1) / ybar_mod, 0.0)
dH2 <- ifelse(den > 0, dN / den, 0.0)
sig2s <- (h1+h2)*(1-dH1)*dH2
sig2U <- sum(sig2s)
i_bhat <- sum(ifelse(den > 0, (num / (den^2)) * (dN0 + dN1), 0.0))
sig2_beta_asy <- sig2U/i_bhat^2
return(list(score = score, sig2U = sig2U, sig2_beta_asy = sig2_beta_asy, i_bhat = i_bhat, K_wt_rg = K))
}


#' Cox score with rho-gamma weights
#'
#' Calculates the Cox score statistic for weighted log-rank tests.
#'
#' @param beta Log hazard ratio parameter.
#' @param time Numeric vector of event times.
#' @param delta Numeric vector of event indicators.
#' @param z Numeric vector of group indicators.
#' @param w_hat Numeric vector of weights.
#' @param wt_rg Numeric vector of rho-gamma weights.
#' @param score_only Logical; if TRUE, returns only the score.
#' @return Numeric value of the score or a list with additional results.
#' @export

cox_score_rhogamma <- function(beta, time, delta, z, w_hat = rep(1,length(time)), wt_rg = rep(1,length(time)), score_only = TRUE) {
  at_points <- time
  tt0 <- time[z == 0]
  dd0 <- delta[z == 0]
  w0_hat <- w_hat[z == 0]
  ybar0 <- colSums(outer(tt0, at_points, FUN = ">=") * w0_hat)
  event_mat0 <- outer(tt0[dd0 == 1], at_points, FUN = "<=") * w0_hat[dd0 == 1]
  counting0 <- colSums(event_mat0)
  dN0 <- diff(c(0, counting0))
  tt1 <- time[z == 1]
  dd1 <- delta[z == 1]
  w1_hat <- w_hat[z == 1]
  ybar1 <- colSums(outer(tt1, at_points, FUN = ">=") * w1_hat * exp(beta))
  event_mat1 <- outer(tt1[dd1 == 1], at_points, FUN = "<=") * w1_hat[dd1 == 1]
  counting1 <- colSums(event_mat1)
  dN1 <- diff(c(0, counting1))

  score_stats <- score_calculation(ybar1 = ybar1, ybar0 = ybar0, dN1 = dN1, dN0 = dN0, wt_rg = wt_rg)

  score <- score_stats$score

  if(score_only){
    return(score)
  }
  else{
  sig2U <- score_stats$sig2U
  sig2_beta_asy <- score_stats$sig2_beta_asy
  K_wt_rg <- score_stats$K_wt_rg
  i_bhat <- score_stats$i_bhat
  return(list(score = score, sig2_score = sig2U, sig2_beta_asy = sig2_beta_asy, K_wt_rg = K_wt_rg, i_bhat = i_bhat))
  }
}

#' Weighted Cox Model with Rho-Gamma Weights
#'
#' Fits a weighted Cox proportional hazards model using flexible time-dependent
#' weights (e.g., Fleming-Harrington, Magirr-Burman). Supports resampling-based
#' inference for variance estimation and bias correction.
#'
#' @param dfcount List; output from \code{\link{df_counting}} containing counting
#'   process data including time, delta, z, w_hat, survP_all, survG_all, etc.
#' @param scheme Character; weighting scheme. See \code{\link{wt.rg.S}} for options.
#'   Default: "fh".
#' @param scheme_params List; parameters for the selected scheme. Default:
#'   \code{list(rho = 0, gamma = 0.5)}. Required parameters depend on scheme:
#'   \itemize{
#'     \item "fh": \code{rho}, \code{gamma}
#'     \item "MB": \code{mb_tstar}
#'     \item "custom_time": \code{t.tau}, \code{w0.tau}, \code{w1.tau}
#'   }
#' @param draws Integer; number of resampling draws for variance estimation and
#'   bias correction. If 0, only asymptotic inference is performed. Default: 0.
#' @param alpha Numeric; significance level for confidence intervals. Default: 0.05.
#' @param verbose Logical; whether to print detailed output. Default: FALSE.
#' @param lr.digits Integer; number of decimal places for formatted output. Default: 4.
#'
#' @return A list containing:
#' \describe{
#'   \item{fit}{List with fitted model components:
#'     \itemize{
#'       \item \code{bhat}: Estimated log hazard ratio
#'       \item \code{sig_bhat_asy}: Asymptotic standard error
#'       \item \code{u.zero}: Score statistic at beta=0 (log-rank)
#'       \item \code{z.score}: Standardized score statistic
#'       \item \code{sig2_score}: Variance of score statistic
#'       \item \code{wt_rg}: Vector of time-dependent weights
#'       \item \code{bhat_debiased}: Bias-corrected estimate (if draws > 0)
#'       \item \code{sig_bhat_star}: Resampling-based standard error (if draws > 0)
#'     }
#'   }
#'   \item{hr_ci_asy}{Data frame with asymptotic HR and CI}
#'   \item{hr_ci_star}{Data frame with resampling-based HR and CI (if draws > 0)}
#'   \item{cox_text_asy}{Formatted string with HR and asymptotic CI}
#'   \item{cox_text_star}{Formatted string with HR and resampling CI (if draws > 0)}
#'   \item{z.score_debiased}{Bias-corrected score statistic (if draws > 0)}
#'   \item{zlogrank_text}{Formatted log-rank test result}
#' }
#'
#' @details
#' This function solves the weighted Cox partial likelihood score equation:
#' \deqn{U(\beta) = \sum_i w_i K_i \left(\frac{dN_0}{Y_0} - \frac{dN_1}{Y_1}\right) = 0}
#'
#' where \eqn{K_i} are time-dependent weights and \eqn{Y_j, dN_j} are risk sets
#' and event counts for group j.
#'
#' When \code{draws > 0}, the function performs resampling to:
#' \itemize{
#'   \item Estimate finite-sample variance (more accurate than asymptotic)
#'   \item Compute bias correction for \eqn{\hat{\beta}}
#'   \item Provide improved confidence intervals for small samples
#' }
#'
#' The score test at \eqn{\beta=0} corresponds to the weighted log-rank test.
#'
#' @note The treatment variable in \code{dfcount} must be coded as 0=control, 1=experimental.
#'
#' @examples
#' \dontrun{
#' # First get counting process data
#' library(survival)
#' data(veteran)
#' veteran$treat <- as.numeric(veteran$trt) - 1
#'
#' dfcount <- df_counting(
#'   df = veteran,
#'   tte.name = "time",
#'   event.name = "status",
#'   treat.name = "treat"
#' )
#'
#' # Fit weighted Cox model with FH(0,0.5) weights
#' fit <- cox_rhogamma(
#'   dfcount = dfcount,
#'   scheme = "fh",
#'   scheme_params = list(rho = 0, gamma = 0.5),
#'   draws = 1000,
#'   verbose = TRUE
#' )
#'
#' print(fit$cox_text_star)  # Resampling-based CI
#' print(fit$zlogrank_text)  # Weighted log-rank test
#'
#' # Compare asymptotic and resampling CIs
#' print(fit$hr_ci_asy)
#' print(fit$hr_ci_star)
#' }
#'
#' @seealso
#' \code{\link{df_counting}} for preprocessing
#' \code{\link{wt.rg.S}} for weighting schemes
#' \code{\link{cox_score_rhogamma}} for score function
#'
#' @references
#' Magirr, D. and Burman, C. F. (2019). Modestly weighted logrank tests.
#' Statistics in Medicine, 38(20), 3782-3790.
#'
#' @family survival_analysis
#' @family weighted_tests
#' @importFrom stats confint var pnorm
#' @export

cox_rhogamma <- function(dfcount, scheme = "fh", scheme_params = list(rho = 0, gamma = 0.5), draws = 0, alpha = 0.05, verbose = FALSE, lr.digits = 4) {
# Extract weights for Cox estimation to align with weights functions

  df_weights_cox <- data.frame(at_points = dfcount$time, S.pool = dfcount$survP_all, G.pool = dfcount$survG_all, Ybar = dfcount$ybar_all)

  wt_rg <- get_validated_weights(df_weights_cox, scheme, scheme_params)

  ans <- list()
  # Extract variables
  time   <- dfcount$time
  delta  <- dfcount$delta
  z      <- dfcount$z
  w_hat  <- dfcount$w_hat
  S.pool <- dfcount$survP_all

  # Input checks
  stopifnot(is.numeric(time), is.numeric(delta), is.numeric(z), is.numeric(w_hat), is.numeric(S.pool))
  n <- length(time)
  n0 <- sum(z == 0)
  n1 <- sum(z == 1)
  if (n0 + n1 != n) stop("z must be a (0/1) treatment group indicator")


  # Find root of score function
  get_Cox <- find_cox_root(time, delta, z, w_hat, wt_rg)

  if (is.na(get_Cox$root)) {
    warning("Root finding failed.")
    return(list(bhat = NA, u.beta = NA, u.zero = NA, status = "fail"))
  }

  bhat_rhogamma <- get_Cox$root

  # Score test: U(beta=0)
  temp <- cox_score_rhogamma(beta = 0, time = time, delta = delta, w_hat = w_hat, z = z, wt_rg = wt_rg, score_only = FALSE)
  u.zero <- temp$score
  z.score <- u.zero / sqrt(temp$sig2_score)
  i_zero <- temp$i_bhat
  K_zero <- temp$K_wt_rg
  sig2_score <- temp$sig2_score
  rm("temp")

  ans$z.score <- z.score

  temp <- cox_score_rhogamma(beta = bhat_rhogamma, time = time, delta = delta, w_hat = w_hat, z = z, wt_rg = wt_rg, score_only = FALSE)
  u.beta <- temp$score
  sig_bhat_asy <- sqrt(temp$sig2_beta_asy)
  i_bhat <- temp$i_bhat
  K_wt_rg <- temp$K_wt_rg
  rm("temp")

  pval <- 1 - pnorm(z.score)
  ans$zlogrank_text <- paste0("logrank (1-sided) p = ", format_pval(pval, eps = 0.001, digits = lr.digits))
  if(verbose) cat("z-statistic: ", ans$zlogrank_text, "\n")

  # CI based on asymptotic SE
  hr_ci_asy <- ci_cox(bhat = bhat_rhogamma, sig_bhat = sig_bhat_asy, alpha = alpha, verbose = verbose)
  ans$hr_ci_asy <- hr_ci_asy

  fit_rhogamma <- list(
    bhat = bhat_rhogamma,
    sig_bhat_asy = sig_bhat_asy,
    u.beta = u.beta,
    u.zero = u.zero,
    z.score = z.score,
    sig2_score = sig2_score,
    status = "ok",
    wt_rg = wt_rg,
    time = time,
    delta = delta,
    z = z,
    w_hat = w_hat
  )

  ans$fit <- fit_rhogamma

  if(draws > 0){

  get_resamples <- cox_rhogamma_resample(fit_rhogamma = fit_rhogamma, i_bhat = i_bhat, K_wt_rg = K_wt_rg,
                                         i_zero = i_zero, K_zero = K_zero, draws = draws
                                         )
  ans$fit_resamples <- get_resamples
  sig_bhat_star <- get_resamples$sig_bhat_star
  ans$fit$sig_bhat_star <- sig_bhat_star
  # De-biased hr
  bhat_debiased <- fit_rhogamma$bhat - mean(get_resamples$bhat_center_star, na.rm = TRUE)
  ans$fit$bhat_debiased <- bhat_debiased
  ans$fit$wald_debiased1 <- bhat_debiased / sig_bhat_asy
  ans$fit$wald_debiased2 <- bhat_debiased / sig_bhat_star
  ans$fit$wald <- fit_rhogamma$bhat / sig_bhat_asy
  # CI based on de-biased and resampled SEs
  hr_ci_star <- ci_cox(bhat = bhat_debiased, sig_bhat = sig_bhat_star, alpha = alpha, verbose = verbose)
  ans$hr_ci_star <- hr_ci_star
   # De-biased score
  sig2_score_star <- var(get_resamples$score_star, na.rm = TRUE)
  ans$sig2_score_star <- sig2_score_star
  u.zero_debiased <- u.zero - mean(get_resamples$score_star_null, na.rm = TRUE)
  ans$z.score_debiased <- u.zero_debiased / sqrt(sig2_score)
  #ans$z.score_debiased <- u.zero / sqrt(sig2_score_star)

  cox_text_star <- paste0("HR = ", round(hr_ci_asy$hr, lr.digits),
                          " (", round(hr_ci_star$lower, lr.digits), ", ", round(hr_ci_star$upper, lr.digits), ")")

  ans$cox_text_star <- cox_text_star


  }

cox_text_asy <- paste0("HR = ", round(hr_ci_asy$hr, lr.digits),
                     " (", round(hr_ci_asy$lower, lr.digits), ", ", round(hr_ci_asy$upper, lr.digits), ")")

ans$cox_text_asy <- cox_text_asy


ans
}



#' Resampling for Weighted Cox Model (rho, gamma)
#'
#' Performs resampling to estimate uncertainty for the weighted Cox model (rho, gamma).
#'
#' @param fit_rhogamma List with fitted Cox model results.
#' @param i_bhat Information at estimated beta.
#' @param K_wt_rg Weights at estimated beta.
#' @param i_zero Information at beta=0.
#' @param K_zero Weights at beta=0.
#' @param G1.draws Optional: pre-generated random draws for groups.
#' @param G0.draws Optional: pre-generated random draws for groups.
#' @param draws Number of resampling iterations (default: 100).
#' @param seedstart Random seed for reproducibility (default: 8316951).
#' @return List with resampling results (score, beta, standard error, etc.).
#' @importFrom stats var rnorm
#' @export

cox_rhogamma_resample <- function(fit_rhogamma, i_bhat, K_wt_rg, i_zero, K_zero, G1.draws = NULL, G0.draws = NULL,
                                  draws = 100, seedstart=8316951
) {

  bhat <- fit_rhogamma$bhat
  sig_bhat_asy <- fit_rhogamma$sig_bhat_asy
  time <- fit_rhogamma$time
  delta <- fit_rhogamma$delta
  z <- fit_rhogamma$z
  w_hat <- fit_rhogamma$w_hat
  wt_rg <- fit_rhogamma$wt_rg
  at_points <- time

  stopifnot(is.numeric(bhat), is.numeric(time), is.numeric(delta), is.numeric(z), is.numeric(w_hat), is.numeric(wt_rg))

  n <- length(time)
  n0 <- sum(z == 0)
  n1 <- sum(z == 1)
  if (n0 + n1 != n) stop("z must be a (0/1) treatment group indicator")

  if (is.null(G0.draws) && is.null(G1.draws) && draws > 0) {
    set.seed(seedstart)
    G0.draws <- matrix(rnorm(draws * n0), ncol = draws)
    G1.draws <- matrix(rnorm(draws * n1), ncol = draws)
  }
  idx0 <- which(z == 0)
  idx1 <- which(z == 1)
  y0 <- time[idx0]; d0 <- delta[idx0]
  y1 <- time[idx1]; d1 <- delta[idx1]
  w0_hat <- w_hat[idx0]; w1_hat <- w_hat[idx1]

  event_mat0 <- outer(y0, at_points, FUN = "<=")
  risk_mat0  <- outer(y0, at_points, FUN = ">=")
  event_mat1 <- outer(y1, at_points, FUN = "<=")
  risk_mat1  <- outer(y1, at_points, FUN = ">=")
  risk_z0 <- colSums(risk_mat0 *  w0_hat)
  risk_z1 <- colSums(risk_mat1 *  w1_hat * exp(bhat))


  counting0 <- colSums(event_mat0 * (d0 *  w0_hat))
  counting1 <- colSums(event_mat1 * (d1 *  w1_hat))
  dN_z0 <- diff(c(0, counting0))
  dN_z1 <- diff(c(0, counting1))

      counting0_star_all <- t(event_mat0 *  w0_hat) %*% (d0 * G0.draws)
      counting1_star_all <- t(event_mat1 *  w1_hat) %*% (d1 * G1.draws)
      dN_z0_star_all <- apply(counting0_star_all, 2, function(x) diff(c(0, x)))
      dN_z1_star_all <- apply(counting1_star_all, 2, function(x) diff(c(0, x)))

      drisk1_star <- sweep(dN_z1_star_all, 1, risk_z1, "/")
      drisk1_star[is.infinite(drisk1_star) | is.nan(drisk1_star)] <- 0
      drisk0_star <- sweep(dN_z0_star_all, 1, risk_z0, "/")
      drisk0_star[is.infinite(drisk0_star) | is.nan(drisk0_star)] <- 0

      score_star <- colSums(K_wt_rg * (drisk0_star - drisk1_star))

      bhat_center_star <- score_star / i_bhat
      Z_bstar <- bhat_center_star / sig_bhat_asy
      var_bhat_star <- var(bhat_center_star, na.rm = TRUE)
      sig_bhat_star <- sqrt(var_bhat_star)

       # Logrank stat (replace risk_z1(bhat) with risk_z1_null = risk_z1(bhat=0))
       risk_z1_null <- colSums(risk_mat1 *  w1_hat)
       drisk1_star_null <- sweep(dN_z1_star_all, 1, risk_z1_null, "/")
       drisk1_star_null[is.infinite(drisk1_star_null) | is.nan(drisk1_star_null)] <- 0
       score_star_null <- colSums(K_zero * (drisk0_star - drisk1_star_null))

      ans <- list(
        score_star = score_star,
        score_star_null = score_star_null,
        bhat_center_star = bhat_center_star,
        Z_bstar = Z_bstar,
        sig_bhat_star = sig_bhat_star
      )

return(ans)
  }

