#' Check and load required R packages
#'
#' Checks if required R packages are installed and loads them.
#'
#' @param pkgs Character vector of package names.
#' @return Logical vector indicating which packages are available.
#' @importFrom utils installed.packages
#' @export

check_required_packages <- function(pkgs) {
  missing <- pkgs[!sapply(pkgs, requireNamespace, quietly = TRUE)]
  if (length(missing) > 0) {
    stop("The following required packages are missing: ", paste(missing, collapse = ", "),
         ". Please install them before running this function.", call. = FALSE)
  }
}

#' Set up parallel processing
#'
#' Sets up parallel processing using the specified approach and number of workers.
#'
#' @param approach Character string specifying parallel approach ("multisession", "multicore", "callr").
#' @param workers Number of parallel workers (default: 4).
#' @return Parallel backend object.
#' @export

setup_parallel <- function(approach = c("multisession", "multicore", "callr"), workers = 4) {
  approach <- match.arg(approach)
  library(future)
  if (approach == "multisession") {
    plan(multisession, workers = workers)
    message("Parallel plan: multisession with ", workers, " workers.")
  } else if (approach == "multicore") {
    plan(multicore, workers = workers)
    message("Parallel plan: multicore with ", workers, " workers.")
  } else if (approach == "callr") {
    if (!requireNamespace("future.callr", quietly = TRUE)) {
      stop("The 'future.callr' package is required for the callr approach.")
    }
    plan(callr, workers = workers)
    message("Parallel plan: callr with ", workers, " workers.")
  } else {
    stop("Unknown parallel approach: ", approach)
  }
}

#' Check if true value is covered by confidence interval
#'
#' Checks if a target value is within a given confidence interval.
#'
#' @param target Numeric value to check.
#' @param hr_ci Numeric vector of length 2 (lower and upper confidence interval).
#' @return Logical indicating if target is covered.
#' @export

is_covered <- function(target,hr_ci){
ifelse(hr_ci$lower <= target & hr_ci$upper >= target, 1, 0)
}

#' Run a single simulation analysis for weighted Cox model
#'
#' Runs a single simulation scenario for weighted Cox model analysis.
#'
#' @param scen List or object describing the scenario.
#' @param enroll_rate Enrollment rate.
#' @param dropout_rate Dropout rate.
#' @param fr Follow-up rate.
#' @param delay Delay time (default: 12).
#' @param n_sample Number of samples (default: 698).
#' @param sim_num Simulation number.
#' @param mart_draws Number of martingale draws (default: 300).
#' @param hr_true True hazard ratio.
#' @param seedstart Random seed (default: 8316951).
#' @return List or data frame with simulation results.
#' @export

sim_fn_analysis <- function(scen, enroll_rate, dropout_rate, fr, delay = 12, n_sample = 698, sim_num, mart_draws = 300, hr_true, seedstart = 8316951){
if(is.na(hr_true) | length(hr_true) !=1) stop("Target hazard-ratio hr_true is missing or of length > 1")
res <- data.table()
res$Scenario <- c(scen)

fail_rate <- data.frame(stratum = rep("All", 2 * nrow(fr)),
                          period = rep(fr$Period, 2),
                          treatment = c(rep("control", nrow(fr)), rep("experimental", nrow(fr))),
                          duration = rep(fr$duration, 2),
                          rate = c(fr$rate, fr$rate * fr$hr)
  )

  # Generate a dataset
  set.seed(seedstart + 1000*sim_num)

  dat <- sim_pw_surv(n = n_sample, enroll_rate = enroll_rate,
                     fail_rate = fail_rate, dropout_rate = dropout_rate)

  analysis_data <- cut_data_by_date(dat, 36)
  dfa <- analysis_data
  dfa$treat <- ifelse(dfa$treatment=="experimental",1,0)

  # reverse sign for superiority z-stat
  maxcombopz <- (analysis_data |> maxcombo(rho = c(0,0), gamma = c(0,0.5), return_corr = TRUE))
  res$logrankz <- (-1) * maxcombopz$z[1]
  res$fh05z <- (-1) *maxcombopz$z[2]
  res$maxcombop <- maxcombopz$p_value
  # Seems for wlr weighting that the signs are reversed relative to fh?
  res$fh01z <- (analysis_data |> wlr(weight = fh(rho = 0, gamma =1)))$z
  res$mb6z <- (analysis_data |> wlr(weight = mb(delay = 6, w_max = Inf)))$z
  res$mb12z <- (analysis_data |> wlr(weight = mb(delay = 12, w_max = Inf)))$z
  res$mb16z <- (analysis_data |> wlr(weight = mb(delay = 16, w_max = Inf)))$z


  # Weighted Cox versions
  # Note: turnoff KM and seKM checking so that plots are not called within parallel sessions
  dfcount <- get_dfcounting(df=dfa, tte.name = "tte", event.name = "event", treat.name = "treat", arms = "treatment", by.risk= 6, check.KM = FALSE, check.seKM = FALSE)
  # MB to compare with simtrial
  res$mb6z_mine <- cox_rhogamma(dfcount, scheme = "MB", scheme_params = list(mb_tstar = 6))$z.score
  #res$mb12z_mine <- cox_rhogamma(dfcount = dfcount, scheme = "MB", scheme_params = list(mb_tstar = 12))$z.score
  res$mb16z_mine <- cox_rhogamma(dfcount, scheme = "MB", scheme_params = list(mb_tstar = 16))$z.score
  # Exponential versions
  #res$fhexp1z <- cox_rhogamma(dfcount = dfcount, scheme = "fh_exp1")$z.score
  #res$fhexp2z <- cox_rhogamma(dfcount = dfcount, scheme = "fh_exp2")$z.score

  temp <- cox_rhogamma(dfcount, scheme = "MB", scheme_params = list(mb_tstar = 12), draws = mart_draws)

  res$mb12z_mine <- temp$z.score
  res$mb12z_debiased <- temp$z.score_debiased
  res$mb12_bhat <- with(temp$fit,bhat)
  res$mb12_bhatdebiased <- with(temp,hr_ci_star$beta)
  res$mb12_wald <- with(temp,hr_ci_asy$upper)
  res$mb12_waldstar <- with(temp,hr_ci_star$upper)
  res$mb12_sigbhat <- with(temp$fit,sig_bhat_asy)
  res$mb12_sigbhatstar <- with(temp$fit,sig_bhat_star)
  res$mb12_cover <- with(temp, is_covered(hr_true, hr_ci_asy))
  res$mb12_coverstar <- with(temp, is_covered(hr_true, hr_ci_star))
  rm("temp")


  # Try de-biased versions for tests and estimates
  # log-rank
  temp <- cox_rhogamma(dfcount, scheme = "fh", scheme_params = list(rho = 0, gamma = 0), draws = mart_draws)
  res$fh00z_mine <- temp$z.score
  res$fh00z_debiased <- temp$z.score_debiased
  res$fh00_bhat <- with(temp$fit,bhat)
  res$fh00_bhatdebiased <- with(temp,hr_ci_star$beta)
  res$fh00_wald <- with(temp,hr_ci_asy$upper)
  res$fh00_waldstar <- with(temp,hr_ci_star$upper)
  res$fh00_sigbhat <- with(temp$fit,sig_bhat_asy)
  res$fh00_sigbhatstar <- with(temp$fit,sig_bhat_star)
  res$fh00_cover <- with(temp, is_covered(hr_true, hr_ci_asy))
  res$fh00_coverstar <- with(temp, is_covered(hr_true, hr_ci_star))
  rm("temp")

  temp <- cox_rhogamma(dfcount, scheme = "fh_exp1", draws = mart_draws)
  res$fhe1z <- temp$z.score
  res$fhe1z_debiased <- temp$z.score_debiased
  res$fhe1_bhat <- with(temp$fit,bhat)
  res$fhe1_bhatdebiased <- with(temp,hr_ci_star$beta)
  res$fhe1_wald <- with(temp,hr_ci_asy$upper)
  res$fhe1_waldstar <- with(temp,hr_ci_star$upper)
  res$fhe1_sigbhat <- with(temp$fit,sig_bhat_asy)
  res$fhe1_sigbhatstar <- with(temp$fit,sig_bhat_star)
  res$fhe1_cover <- with(temp, is_covered(hr_true, hr_ci_asy))
  res$fhe1_coverstar <- with(temp, is_covered(hr_true, hr_ci_star))
  rm("temp")

  temp <- cox_rhogamma(dfcount, scheme = "fh_exp2", draws = mart_draws)
  res$fhe2z <- temp$z.score
  res$fhe2z_debiased <- temp$z.score_debiased
  res$fhe2_bhat <- with(temp$fit,bhat)
  res$fhe2_bhatdebiased <- with(temp,hr_ci_star$beta)
  res$fhe2_wald <- with(temp,hr_ci_asy$upper)
  res$fhe2_waldstar <- with(temp,hr_ci_star$upper)
  res$fhe2_sigbhat <- with(temp$fit,sig_bhat_asy)
  res$fhe2_sigbhatstar <- with(temp$fit,sig_bhat_star)
  res$fhe2_cover <- with(temp, is_covered(hr_true, hr_ci_asy))
  res$fhe2_coverstar <- with(temp, is_covered(hr_true, hr_ci_star))
  rm("temp")


  # FH(0,1)
  temp <- cox_rhogamma(dfcount, scheme = "fh", scheme_params = list(rho = 0, gamma = 1), draws = mart_draws)
  res$fh01z_mine <- temp$z.score
  res$fh01z_debiased <- temp$z.score_debiased
  res$fh01_bhat <- with(temp$fit,bhat)
  res$fh01_bhatdebiased <- with(temp,hr_ci_star$beta)
  res$fh01_wald <- with(temp,hr_ci_asy$upper)
  res$fh01_waldstar <- with(temp,hr_ci_star$upper)
  res$fh01_sigbhat <- with(temp$fit,sig_bhat_asy)
  res$fh01_sigbhatstar <- with(temp$fit,sig_bhat_star)
  res$fh01_cover <- with(temp, is_covered(hr_true, hr_ci_asy))
  res$fh01_coverstar <- with(temp, is_covered(hr_true, hr_ci_star))
  rm("temp")

  # FH(0,0.5)
  temp <- cox_rhogamma(dfcount, scheme = "fh", scheme_params = list(rho = 0, gamma = 0.5), draws = mart_draws)
  res$fh05z_mine <- temp$z.score
  res$fh05z_debiased <- temp$z.score_debiased
  res$fh05_bhat <- with(temp$fit,bhat)
  res$fh05_bhatdebiased <- with(temp,hr_ci_star$beta)
  res$fh05_wald <- with(temp,hr_ci_asy$upper)
  res$fh05_waldstar <- with(temp,hr_ci_star$upper)
  res$fh05_sigbhat <- with(temp$fit,sig_bhat_asy)
  res$fh05_sigbhatstar <- with(temp$fit,sig_bhat_star)
  res$fh05_cover <- with(temp, is_covered(hr_true, hr_ci_asy))
  res$fh05_coverstar <- with(temp, is_covered(hr_true, hr_ci_star))
  rm("temp")


  # t6(0,1)
  temp <- cox_rhogamma(dfcount, scheme = "custom_time", scheme_params = list(t.tau = 6, w0.tau = 0, w1.tau = 1), draws = mart_draws)
  res$t601z <- temp$z.score
  res$t601z_debiased <- temp$z.score_debiased
  res$t601_bhat <- with(temp$fit,bhat)
  res$t601_bhatdebiased <- with(temp,hr_ci_star$beta)
  res$t601_wald <- with(temp,hr_ci_asy$upper)
  res$t601_waldstar <- with(temp,hr_ci_star$upper)
  res$t601_sigbhat <- with(temp$fit,sig_bhat_asy)
  res$t601_sigbhatstar <- with(temp$fit,sig_bhat_star)
  res$t601_cover <- with(temp, is_covered(hr_true, hr_ci_asy))
  res$t601_coverstar <- with(temp, is_covered(hr_true, hr_ci_star))
  rm("temp")
  return(as.data.frame(res))
}


#' Run multiple simulations for weighted Cox model
#'
#' Runs multiple simulations in parallel for weighted Cox model analysis.
#'
#' @param n_sim Number of simulations.
#' @param dof_approach Parallelization approach ("callr", etc.).
#' @param num_workers Number of parallel workers (default: 4).
#' @param n_sample Number of samples (default: 698).
#' @param seedstart Random seed (default: 8316951).
#' @param file_togo File path to save results.
#' @param save_results Logical; if TRUE, saves results to file.
#' @param verbose Logical; if TRUE, prints progress.
#' @param mart_draws Number of martingale draws (default: 100).
#' @return List or data frame with all simulation results.
#' @export

get_sims <- function(n_sim, dof_approach = "callr", num_workers = 4,  n_sample = 698, seedstart = 8316951, file_togo = c("results/sims_example_new.RData"), save_results = FALSE, verbose = TRUE, mart_draws = 100){

  required_pkgs <- c("dplyr", "tibble", "foreach", "future", "tictoc", "simtrial", "doFuture")
  if(dof_approach == "callr") required_pkgs <- c(required_pkgs, "future.callr")

  check_required_packages(required_pkgs)

  if(save_results){
    dir_name <- dirname(file_togo)
    dir_dne <- !dir.exists(dir_name)
    if(dir_dne) stop(paste("File location for saving results does not exist (check where you are running from):", dir_name))
  }

  get_setup <- sims_setup()
  fr <- get_setup$fr
  dropout_rate <- get_setup$dropout_rate
  enroll_rate <- get_setup$er

  # Consider the "true hr" to correspond to PH scenario
  fr_PH = fr |> dplyr::filter(Name == "PH")

  hr_PH <- fr_PH[1,]$hr

  setup_parallel(approach = dof_approach, workers = num_workers)

  fr$stratum <- "All"

  tictoc::tic(log = FALSE)
  # Do simulations by scenario
  results_sims <- foreach(scen = 1:6, .combine = 'rbind')%:%
    foreach(sim = 1:n_sim, .combine = 'rbind', .options.future = list(seed = TRUE)) %dofuture% {
      library(simtrial)
      tryCatch({sim_fn_analysis(scen = scen, enroll_rate = enroll_rate, n_sample = n_sample, dropout_rate = dropout_rate, mart_draws = mart_draws, hr_true = hr_PH, sim_num = sim,
                                fr = fr |> dplyr::filter(Scenario == scen))},
               error = function(e) NA)
    }

  invisible(plan(sequential))

  toc_result <- tictoc::toc(log = FALSE)
  elapsed_seconds <- as.numeric(toc_result$toc - toc_result$tic)

  if(verbose) cat("Nsims, timing minutes =", c(n_sim, round(elapsed_seconds / 60, 2)), "\n")


  res_out <- list(get_setup = get_setup, results_sims = results_sims, tminutes = c(elapsed_seconds / 60), thours = c(elapsed_seconds / (60^2)),
                  number_sims = n_sim, hr_target = hr_PH, seedstart = seedstart)

  if(save_results) save(res_out, file = file_togo)

  return(res_out)
}
