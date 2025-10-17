#' Weighted and Stratified Survival Analysis
#'
#' Performs comprehensive weighted and/or stratified survival analysis, including
#' Cox proportional hazards model, logrank/Fleming-Harrington tests, and calculation
#' of risk/event sets, Kaplan-Meier curves, quantiles, and variance estimates.
#'
#' @param df Data frame containing survival data.
#' @param tte.name Character; name of the time-to-event column in \code{df}.
#' @param event.name Character; name of the event indicator column in \code{df} (1=event, 0=censored).
#' @param treat.name Character; name of the treatment/group column in \code{df} (must be coded as 0=control, 1=experimental).
#' @param weight.name Character or NULL; name of the weights column in \code{df}. If NULL, equal weights are used.
#' @param strata.name Character or NULL; name of the strata column in \code{df} for stratified analysis.
#' @param arms Character vector of length 2; group labels. Default: \code{c("treat","control")}.
#' @param time.zero Numeric; time value to use as zero. Default: 0.
#' @param tpoints.add Numeric vector; additional time points to include in calculations. Default: \code{c(0)}.
#' @param by.risk Numeric; interval for risk set time points. Default: 6.
#' @param time.zero.label Numeric; label for time zero in output. Default: 0.0.
#' @param risk.add Numeric vector or NULL; additional specific risk points to include.
#' @param get.cox Logical; whether to fit Cox proportional hazards model. Default: TRUE.
#' @param cox.digits Integer; number of decimal places for Cox output formatting. Default: 2.
#' @param lr.digits Integer; number of decimal places for logrank output formatting. Default: 2.
#' @param cox.eps Numeric; threshold for Cox p-value formatting (values below shown as "<eps"). Default: 0.001.
#' @param lr.eps Numeric; threshold for logrank p-value formatting. Default: 0.001.
#' @param verbose Logical; whether to print warnings and diagnostic messages. Default: FALSE.
#' @param qprob Numeric in (0,1); quantile probability for KM quantile table. Default: 0.5 (median).
#' @param scheme Character; weighting scheme for logrank/Fleming-Harrington test. 
#'   Options: "fh", "schemper", "XO", "MB", "custom_time", "fh_exp1", "fh_exp2". Default: "fh".
#' @param scheme_params List; parameters for the selected weighting scheme. Default: \code{list(rho = 0, gamma = 0)}.
#'   \itemize{
#'     \item For "fh": \code{rho} and \code{gamma} (Fleming-Harrington parameters)
#'     \item For "MB": \code{mb_tstar} (cutoff time)
#'     \item For "custom_time": \code{t.tau}, \code{w0.tau}, \code{w1.tau}
#'   }
#' @param conf_level Numeric in (0,1); confidence level for quantile intervals. Default: 0.95.
#' @param check.KM Logical; whether to check KM curve validity against \code{survival::survfit}. Default: TRUE.
#' @param check.seKM Logical; whether to check KM standard error estimates. Default: FALSE.
#' @param draws Integer; number of draws for resampling-based variance estimation. Default: 0 (no resampling).
#' @param seedstart Integer; random seed for reproducible resampling. Default: 8316951.
#' @param stop.onerror Logical; whether to stop execution on errors (TRUE) or issue warnings (FALSE). Default: FALSE.
#' @param censoring_allmarks Logical; if FALSE, removes events from censored time points. Default: TRUE.
#'
#' @return A list with the following components:
#' \describe{
#'   \item{cox_results}{List with Cox model results including hazard ratio, confidence interval, p-value, and formatted text}
#'   \item{logrank_results}{List with log-rank test chi-square statistic, p-value, and formatted text}
#'   \item{z.score}{Standardized weighted log-rank test statistic}
#'   \item{at_points}{Vector of all time points used in calculations}
#'   \item{surv0, surv1}{Kaplan-Meier survival estimates for control and treatment groups}
#'   \item{sig2_surv0, sig2_surv1}{Variance estimates for survival curves}
#'   \item{survP}{Pooled survival estimates}
#'   \item{survG}{Censoring distribution estimates}
#'   \item{quantile_results}{Data frame with median survival and confidence intervals by group}
#'   \item{lr, sig2_lr}{Weighted log-rank statistic and its variance}
#'   \item{riskpoints0, riskpoints1}{Risk set counts at specified time points}
#'   \item{z.score_stratified}{Stratified z-score (if stratified analysis)}
#' }
#'
#' @details
#' This function implements a comprehensive survival analysis framework supporting:
#' \itemize{
#'   \item Weighted observations via \code{weight.name}
#'   \item Stratified analysis via \code{strata.name}
#'   \item Multiple weighting schemes for log-rank tests
#'   \item Resampling-based variance estimation
#'   \item Automatic validation against \code{survival} package results
#' }
#'
#' The function performs time-fixing using \code{survival::aeqSurv} to handle tied event times.
#' For stratified analyses, stratum-specific estimates are computed and combined appropriately.
#'
#' @section Weighting Schemes:
#' \describe{
#'   \item{fh}{Fleming-Harrington: \eqn{w(t) = S(t)^{\rho} \times (1-S(t))^{\gamma}}}
#'   \item{MB}{Magirr-Burman: \eqn{w(t) = 1/\max(S(t), S(t^*))}}
#'   \item{schemper}{Schemper: \eqn{w(t) = S(t)/G(t)} where G is censoring distribution}
#'   \item{XO}{Xu-O'Quigley: \eqn{w(t) = S(t)/Y(t)} where Y is risk set size}
#' }
#'
#' @examples
#' \dontrun{
#' # Basic survival analysis
#' library(survival)
#' data(veteran)
#' veteran$treat <- as.numeric(veteran$trt) - 1
#' 
#' result <- df_counting(
#'   df = veteran,
#'   tte.name = "time",
#'   event.name = "status",
#'   treat.name = "treat",
#'   arms = c("Treatment", "Control")
#' )
#' 
#' # Print results
#' print(result$cox_results$cox_text)
#' print(result$zlogrank_text)
#' 
#' # Fleming-Harrington (0,1) weights (emphasizing late differences)
#' result_fh <- df_counting(
#'   df = veteran,
#'   tte.name = "time",
#'   event.name = "status",
#'   treat.name = "treat",
#'   scheme = "fh",
#'   scheme_params = list(rho = 0, gamma = 1)
#' )
#' 
#' # Stratified analysis
#' result_strat <- df_counting(
#'   df = veteran,
#'   tte.name = "time",
#'   event.name = "status",
#'   treat.name = "treat",
#'   strata.name = "celltype"
#' )
#' }
#'
#' @seealso 
#' \code{\link[survival]{coxph}}, \code{\link[survival]{survdiff}}, \code{\link[survival]{survfit}}
#' \code{\link{cox_rhogamma}} for weighted Cox models
#' \code{\link{KM_diff}} for Kaplan-Meier differences
#'
#' @references
#' Fleming, T. R. and Harrington, D. P. (1991). Counting Processes and Survival Analysis. Wiley.
#' 
#' Magirr, D. and Burman, C. F. (2019). Modestly weighted logrank tests. 
#' Statistics in Medicine, 38(20), 3782-3790.
#'
#' @family survival_analysis
#' @family weighted_tests
#' @importFrom stats as.formula pchisq
#' @importFrom survival aeqSurv Surv coxph survfit survdiff
#' @export
df_counting <- function(df, tte.name = "tte", event.name = "event", treat.name = "treat", 
                        weight.name=NULL, strata.name = NULL, arms=c("treat","control"),
                        time.zero=0, tpoints.add=c(0), by.risk=6, time.zero.label = 0.0, 
                        risk.add=NULL, get.cox=TRUE, cox.digits=2, lr.digits=2,
                        cox.eps = 0.001, lr.eps = 0.001, verbose = FALSE, qprob=0.5,
                        scheme = "fh", scheme_params = list(rho = 0, gamma = 0),
                        conf_level = 0.95, check.KM = TRUE, check.seKM = FALSE, 
                        draws = 0, seedstart = 8316951, stop.onerror=FALSE,
                        censoring_allmarks=TRUE) {
  # Function body here (add your original function code)
}