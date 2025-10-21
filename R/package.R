#' Weighted Survival Analysis Tools
#'
#' Provides functions for weighted and stratified survival analysis, including
#' Cox proportional hazards models, weighted log-rank tests, Kaplan-Meier curves,
#' and RMST calculations with various weighting schemes.
#'
#' @section Main Functions:
#' \itemize{
#'   \item \code{\link{df_counting}}: Main function for weighted survival analysis
#'   \item \code{\link{KM_diff}}: Kaplan-Meier difference calculations
#'   \item \code{\link{cox_rhogamma}}: Weighted Cox model with flexible weights
#'   \item \code{\link{plot_weight_schemes}}: Visualize weighting schemes
#'   \item \code{\link{cumulative_rmst_bands}}: RMST analysis with confidence bands
#' }
#'
#' @section Weighting Schemes:
#' The package supports multiple weighting schemes for log-rank tests:
#' \itemize{
#'   \item Fleming-Harrington (fh): Emphasizes early or late differences
#'   \item Magirr-Burman (MB): Modest downweighting after cutoff time
#'   \item Schemper: Adjusts for censoring distribution
#'   \item Xu-O'Quigley (XO): Adjusts for risk set size
#'   \item Custom time-based weights
#' }
#'
#' @section Typical Workflow:
#' \enumerate{
#'   \item Prepare data with treatment (0/1), time-to-event, and event indicator
#'   \item Run \code{df_counting()} for comprehensive analysis
#'   \item Use \code{plot_weighted_km()} to visualize survival curves
#'   \item Calculate survival differences with \code{KM_diff()}
#'   \item Perform weighted Cox regression with \code{cox_rhogamma()}
#'   \item Compute RMST with \code{cumulative_rmst_bands()}
#' }
#'
#' @section Key Features:
#' \itemize{
#'   \item Flexible weighting schemes for hypothesis testing
#'   \item Resampling-based inference for improved small-sample properties
#'   \item Simultaneous confidence bands for survival curves
#'   \item Stratified analysis support
#'   \item Weighted observations
#'   \item Comprehensive diagnostic checks
#' }
#'
#' @examples
#' library(survival)
#' str(veteran)
#' veteran$treat <- as.numeric(veteran$trt) - 1
#'
#' # Basic analysis
#' result <- df_counting(
#'   df = veteran,
#'   tte.name = "time",
#'   event.name = "status",
#'   treat.name = "treat"
#' )
#'
#' # Plot results
#' plot_weighted_km(result)
#'
#' # Weighted log-rank emphasizing late differences
#' result_fh <- df_counting(
#'   df = veteran,
#'   tte.name = "time",
#'   event.name = "status",
#'   treat.name = "treat",
#'   scheme = "fh",
#'   scheme_params = list(rho = 0, gamma = 1)
#' )
#'
#' @references
#' Fleming, T. R. and Harrington, D. P. (1991). Counting Processes and Survival Analysis. Wiley.
#' 
#' Magirr, D. and Burman, C. F. (2019). Modestly weighted logrank tests. 
#' Statistics in Medicine, 38(20), 3782-3790.
#'
#' @keywords internal
"_PACKAGE"