#' Compute Time-Dependent Weights for Survival Analysis
#'
#' Calculates time-dependent weights for survival analysis according to various schemes 
#' (Fleming-Harrington, Schemper, XO, MB, custom).
#'
#' @param S Numeric vector of survival probabilities at each time point.
#' @param scheme Character; weighting scheme. One of: 'fh', 'schemper', 'XO', 'MB', 
#'   'custom_time', 'fh_exp1', 'fh_exp2'.
#' @param rho Numeric; rho parameter for Fleming-Harrington weights. Required if scheme='fh'.
#' @param gamma Numeric; gamma parameter for Fleming-Harrington weights. Required if scheme='fh'.
#' @param Scensor Numeric vector; censoring KM curve for Schemper weights. Must have same length as S.
#' @param Ybar Numeric vector; risk set sizes for XO weights. Must have same length as S.
#' @param tpoints Numeric vector; time points corresponding to S. Required for MB and custom_time schemes.
#' @param t.tau Numeric; cutoff time for custom_time weights.
#' @param w0.tau Numeric; weight before t.tau for custom_time weights.
#' @param w1.tau Numeric; weight after t.tau for custom_time weights.
#' @param mb_tstar Numeric; cutoff time for Magirr-Burman weights.
#' @param details Logical; if TRUE, returns list with weights and diagnostic info. Default: FALSE.
#'
#' @return 
#' If \code{details=FALSE} (default): Numeric vector of weights with length equal to \code{S}.
#' 
#' If \code{details=TRUE}: List containing:
#' \describe{
#'   \item{weights}{Numeric vector of calculated weights}
#'   \item{S}{Input survival probabilities}
#'   \item{S_left}{Left-continuous survival probabilities}
#'   \item{scheme}{Scheme name}
#'   \item{Additional components specific to the scheme}
#' }
#'
#' @details
#' This function implements several weighting schemes for weighted log-rank tests:
#'
#' \describe{
#'   \item{Fleming-Harrington (fh)}{\eqn{w(t) = S(t-)^{\rho} \times (1-S(t-))^{\gamma}}
#'     \itemize{
#'       \item \eqn{\rho=0, \gamma=0}: Standard log-rank (equal weights)
#'       \item \eqn{\rho=0, \gamma=1}: Emphasizes late differences
#'       \item \eqn{\rho=1, \gamma=0}: Emphasizes early differences
#'       \item \eqn{\rho=0.5, \gamma=0.5}: Balanced weighting
#'     }
#'   }
#'   \item{Schemper}{\eqn{w(t) = S(t-)/G(t-)} where G is the censoring distribution.
#'     Upweights times with heavy censoring.
#'   }
#'   \item{Xu-O'Quigley (XO)}{\eqn{w(t) = S(t-)/Y(t)} where Y is risk set size.
#'     Downweights early times with large risk sets.
#'   }
#'   \item{Magirr-Burman (MB)}{\eqn{w(t) = 1/\max(S(t-), S(t^*))} 
#'     Modest downweighting after cutoff time \eqn{t^*}.
#'   }
#'   \item{Custom time}{Step function with weight \eqn{w_0} before \eqn{t^*} 
#'     and \eqn{w_1} after \eqn{t^*}.
#'   }
#' }
#'
#' @note 
#' All weights are calculated using left-continuous survival probabilities \eqn{S(t-)}
#' to ensure consistency with counting process notation.
#'
#' @examples
#' \dontrun{
#' # Generate example survival curve
#' time <- seq(0, 24, by = 0.5)
#' surv <- exp(-0.05 * time)  # Exponential survival
#' 
#' # Fleming-Harrington (0,1) weights
#' w_fh01 <- wt.rg.S(surv, scheme = "fh", rho = 0, gamma = 1, tpoints = time)
#' plot(time, w_fh01, type = "l", main = "FH(0,1) Weights")
#' 
#' # Magirr-Burman weights
#' w_mb <- wt.rg.S(surv, scheme = "MB", mb_tstar = 12, tpoints = time)
#' plot(time, w_mb, type = "l", main = "MB(12) Weights")
#' 
#' # Compare multiple schemes
#' w_lr <- wt.rg.S(surv, scheme = "fh", rho = 0, gamma = 0)  # Log-rank
#' w_fh <- wt.rg.S(surv, scheme = "fh", rho = 0, gamma = 1)  # Late emphasis
#' 
#' plot(time, w_lr, type = "l", ylim = c(0, 2))
#' lines(time, w_fh, col = "blue")
#' legend("topleft", c("Log-rank", "FH(0,1)"), col = 1:2, lty = 1)
#' }
#'
#' @references
#' Fleming, T. R. and Harrington, D. P. (1991). Counting Processes and Survival Analysis. Wiley.
#' 
#' Magirr, D. and Burman, C. F. (2019). Modestly weighted logrank tests. 
#' Statistics in Medicine, 38(20), 3782-3790.
#' 
#' Schemper, M., Wakounig, S., and Heinze, G. (2009). The estimation of average hazard ratios 
#' by weighted Cox regression. Statistics in Medicine, 28(19), 2473-2489.
#'
#' @family survival_analysis
#' @family weighted_tests
#' @export
wt.rg.S <- function(S, scheme = c('fh', 'schemper', 'XO', 'MB', 'custom_time', 'fh_exp1','fh_exp2'),
                    rho = NULL, gamma = NULL, Scensor = NULL, Ybar = NULL, tpoints = NULL,
                    t.tau = NULL, w0.tau = 0, w1.tau = 1, mb_tstar = NULL, details = FALSE) {
  # Function body here (add your original function code)
}