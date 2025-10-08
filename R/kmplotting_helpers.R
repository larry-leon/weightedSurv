
#' Plot confidence interval polygon for KM curve
#'
#' Plots a shaded polygon representing the confidence interval for a Kaplan-Meier survival curve.
#'
#' @param x Numeric vector of time points.
#' @param surv Numeric vector of survival probabilities.
#' @param se Numeric vector of standard errors of survival probabilities.
#' @param conf_level Confidence level for interval (default 0.95).
#' @param col Color for the polygon.
#' @importFrom graphics polygon
#' @export

plot_km_confint_polygon <- function(x, surv, se, conf_level, col) {
  z <- qnorm(1 - (1 - conf_level) / 2)
  lower <- exp(log(surv) - z * se / surv)
  upper <- exp(log(surv) + z * se / surv)
  polygon(
    c(x[order(x)], rev(x[order(x)])),
    c(lower[order(x)], rev(upper[order(x)])),
    col = col, border = FALSE
  )
}

#' Plot KM curves for two groups with optional confidence intervals and censoring marks
#'
#' Plots Kaplan-Meier survival curves for two groups, with options for confidence intervals and censoring marks.
#'
#' @param at_points Numeric vector of time points for plotting.
#' @param S0.KM Numeric vector of survival probabilities for group 0.
#' @param idx0 Indices for censoring in group 0.
#' @param idv0 Indices for events in group 0.
#' @param S1.KM Numeric vector of survival probabilities for group 1.
#' @param idx1 Indices for censoring in group 1.
#' @param idv1 Indices for events in group 1.
#' @param col.0 Color for group 0.
#' @param col.1 Color for group 1.
#' @param ltys Line types for groups.
#' @param lwds Line widths for groups.
#' @param Xlab X-axis label.
#' @param Ylab Y-axis label.
#' @param ylim Y-axis limits.
#' @param xlim X-axis limits.
#' @param show.ticks Logical; show censoring marks (default FALSE).
#' @param cens0 Numeric vector of censoring times for group 0.
#' @param risk.points Numeric vector of risk time points.
#' @param risk.points.label Character vector of labels for risk time points.
#' @param cens1 Numeric vector of censoring times for group 1.
#' @param se0.KM Numeric vector of standard errors for group 0.
#' @param se1.KM Numeric vector of standard errors for group 1.
#' @param conf.int Logical; show confidence intervals (default FALSE).
#' @param conf_level Confidence level (default 0.95).
#' @param censor.cex Numeric; censoring mark size (default 1.0).
#' @param ... Additional arguments to plot.
#' @inheritParams KM_plot_2sample_weighted_counting
#' @importFrom graphics plot lines legend title axis box points
#' @export

plot_km_curves_counting <- function(
    at_points, S0.KM, idx0, idv0, S1.KM, idx1, idv1, col.0, col.1, ltys, lwds,
    Xlab, Ylab, ylim, xlim, show.ticks = FALSE, cens0 = NULL, risk.points, risk.points.label,
    cens1 = NULL, se0.KM = NULL, se1.KM = NULL, conf.int = FALSE, conf_level = 0.95,
    censor.cex = 1.0, time.zero = 0, tpoints.add = c(0), ...
) {
  # Input validation
  if (length(ltys) < 2) stop("ltys must have at least two elements")
  if (length(lwds) < 2) stop("lwds must have at least two elements")
  if (missing(col.0) || missing(col.1)) stop("col.0 and col.1 must be specified")

  plot(
    at_points[idv1], S1.KM[idv1], type = "n", ylim = ylim, xlim = xlim,
    lty = ltys[1], col = col.1, lwd = lwds[1], xlab = Xlab, ylab = Ylab, yaxt = "n",
    xaxt = "n", ...
  )
  axis(2, las = 2)
  axis(1, at = risk.points, labels = risk.points.label)

  if (conf.int && !is.null(se0.KM) && !is.null(se1.KM)) {
    plot_km_confint_polygon(at_points[idv0], S0.KM[idv0], se0.KM[idv0], conf_level, "lightgrey")
    plot_km_confint_polygon(at_points[idv1], S1.KM[idv1], se1.KM[idv1], conf_level, "lightblue")
  }

  lines(at_points[idv1], S1.KM[idv1], type = "s", lty = ltys[1], col = col.1, lwd = lwds[1])
  lines(at_points[idv0], S0.KM[idv0], type = "s", lty = ltys[2], col = col.0, lwd = lwds[2])

  if (show.ticks && !is.null(cens0)) {
    points(cens0, S0.KM[idx0], pch = 3, col = col.0, cex = censor.cex)
  }
  if (show.ticks && !is.null(cens1)) {
    points(cens1, S1.KM[idx1], pch = 3, col = col.1, cex = censor.cex)
  }
}


#' Check KM curve for validity
#'
#' Checks that a Kaplan-Meier curve is valid (values in \[0,1\], non-increasing).
#'
#' @param S.KM Numeric vector of survival probabilities.
#' @param group_name Character; name of the group.
#' @param stop_on_error Logical; whether to stop on error (default TRUE).
#' @return None. Stops or warns if invalid.
#' @export

check_km_curve <- function(S.KM, group_name = "Group", stop_on_error = TRUE) {
  if (any(S.KM < 0 | S.KM > 1, na.rm = TRUE)) {
    msg <- paste0(group_name, " KM curve has values outside [0,1].")
    if (stop_on_error) stop(msg) else warning(msg)
  }
  if (any(diff(S.KM) > 0, na.rm = TRUE)) {
    msg <- paste0(group_name, " KM curve is not non-increasing.")
    if (stop_on_error) stop(msg) else warning(msg)
  }
}

#' Add risk table annotation to KM plot
#'
#' Adds risk set counts for two groups to a Kaplan-Meier plot.
#'
#' @param risk.points Numeric vector of risk time points.
#' @param rpoints0 Numeric vector of risk set counts for group 0.
#' @param rpoints1 Numeric vector of risk set counts for group 1.
#' @param col.0 Color for group 0.
#' @param col.1 Color for group 1.
#' @param risk.cex Numeric; text size for risk table.
#' @param ymin Numeric; minimum y value.
#' @param risk_offset Numeric; offset for risk table.
#' @param risk_delta Numeric; delta for risk table.
#' @param y.risk0 Numeric; y position for group 0 risk table.
#' @param y.risk1 Numeric; y position for group 1 risk table.
#' @importFrom graphics text axis
#' @export

add_risk_table <- function(risk.points, rpoints0, rpoints1, col.0, col.1, risk.cex, ymin, risk_offset, risk_delta, y.risk0, y.risk1) {
  text(risk.points, rep(ifelse(is.null(y.risk0), ymin - risk_offset, y.risk0), length(risk.points)),
       round(rpoints0), col = col.0, cex = risk.cex)
  text(risk.points, rep(ifelse(is.null(y.risk1), ymin - risk_offset + risk_delta, y.risk1), length(risk.points)),
       round(rpoints1), col = col.1, cex = risk.cex)
}

#' Add median annotation to KM plot
#'
#' Adds median survival annotation to a Kaplan-Meier plot.
#'
#' @param medians_df Data frame with quantile results.
#' @param med.digits Integer; digits for median.
#' @param med.cex Numeric; text size for median.
#' @param med.font Integer; font for median.
#' @param xmed.fraction Numeric; x position fraction.
#' @param ymed.offset Numeric; y offset.
#' @importFrom graphics text
#' @export

add_median_annotation <- function(medians_df, med.digits, med.cex, med.font, xmed.fraction, ymed.offset) {
  medians_df$label <- paste0(
    format(medians_df$quantile, digits = med.digits), " [",
    format(medians_df$lower, digits = med.digits), ", ",
    format(medians_df$upper, digits = med.digits), "]"
  )
  annotation_text <- paste(
    paste(medians_df$group, medians_df$label, sep = ": "),
    collapse = "\n"
  )

  usr <- par("usr")
  x_pos <- usr[2] * xmed.fraction
  y_pos <- usr[4] - ymed.offset

  text(x = x_pos, y = y_pos, labels = annotation_text, adj = c(0, 1), cex = med.cex, font = med.font)
}

#' Add legends to KM plot
#'
#' Adds legends for Cox model, log-rank test, and arms to a Kaplan-Meier plot.
#'
#' @param dfcount List with results.
#' @param show.cox Logical; show Cox legend.
#' @param cox.cex Numeric; Cox legend size.
#' @param put.legend.cox Character; Cox legend position.
#' @param show.logrank Logical; show logrank legend.
#' @param logrank.cex Numeric; logrank legend size.
#' @param put.legend.lr Character; logrank legend position.
#' @param show_arm_legend Logical; show arm legend.
#' @param arms Character vector of arm labels.
#' @param col.1 Color for group 1.
#' @param col.0 Color for group 0.
#' @param ltys Line types.
#' @param lwds Line widths.
#' @param arm.cex Numeric; arm legend size.
#' @param put.legend.arms Character; arm legend position.
#' @importFrom graphics legend
#' @export

add_legends <- function(dfcount, show.cox, cox.cex, put.legend.cox, show.logrank, logrank.cex, put.legend.lr, show_arm_legend, arms, col.1, col.0, ltys, lwds, arm.cex, put.legend.arms) {
  if (show.cox && !is.null(dfcount$cox_results)) {
    legend(put.legend.cox, legend = dfcount$cox_results$cox_text, cex = cox.cex, bty = "n")
  }
  if (show.logrank && !is.null(dfcount$z.score)) {
    legend(put.legend.lr, legend = dfcount$zlogrank_text, cex = logrank.cex, bty = "n")
  }
  if (show_arm_legend) {
    legend(put.legend.arms, legend = arms, col = c(col.1, col.0), lty = ltys, lwd = lwds, cex = arm.cex, bty = "n")
  }
}


#' Plot Weighted Kaplan-Meier Curves for Two Samples (Counting Process Format)
#'
#' Plots Kaplan-Meier survival curves for two groups using precomputed risk/event counts and survival estimates.
#' Optionally displays confidence intervals, risk tables, median survival annotations, and statistical test results.
#'
#' @param dfcount List containing precomputed survival data.
#' @param show.cox Logical; show Cox model results.
#' @param cox.cex Numeric; text size for Cox annotation.
#' @param show.logrank Logical; show log-rank test results.
#' @param logrank.cex Numeric; text size for log-rank annotation.
#' @param cox.eps Numeric; small values for Cox calculations.
#' @param lr.eps Numeric; small values for log-rank calculations.
#' @param show_arm_legend Logical; show arm legend.
#' @param arms Character vector of arm labels.
#' @param put.legend.arms Character; legend positions.
#' @param put.legend.cox Character; legend positions.
#' @param put.legend.lr Character; legend positions.
#' @param stop.onerror Logical; stop on KM curve errors.
#' @param check.KM Logical; check KM curve validity.
#' @param lr.digits Integer; digits for test results.
#' @param cox.digits Integer; digits for test results.
#' @param tpoints.add Numeric; additional time points for risk table.
#' @param by.risk Numeric; interval for risk table time points.
#' @param Xlab Character; axis labels.
#' @param Ylab Character; axis labels.
#' @param col.0 Color for control curve.
#' @param col.1 Color for treatment curve.
#' @param show.med Logical; annotate median survival.
#' @param med.digits Median annotation settings.
#' @param med.font Median annotation settings.
#' @param med.cex Median annotation settings.
#' @param ymed.offset Median annotation settings.
#' @param xmed.fraction Median annotation settings.
#' @param conf.int Logical; plot confidence intervals.
#' @param conf_level Numeric; confidence level for intervals.
#' @param choose_ylim Logical; auto-select y-axis limits.
#' @param arm.cex Numeric; text size for arm legend.
#' @param quant Numeric; quantile for annotation.
#' @param qlabel character; label for annotation.
#' @param risk.cex Numeric; text size for risk table.
#' @param ltys Integer; line types for curves.
#' @param lwds Integer; line widths for curves.
#' @param censor.mark.all Logical; mark all censored times.
#' @param censor.cex Numeric; size of censor marks.
#' @param show.ticks Logical; show axis ticks.
#' @param risk.set Logical; display risk table.
#' @param ymin Additional graphical and calculation parameters.
#' @param ymax Additional graphical and calculation parameters.
#' @param ymin.del Additional graphical and calculation parameters.
#' @param ymin2 Additional graphical and calculation parameters.
#' @param risk_offset Additional graphical and calculation parameters.
#' @param risk_delta Additional graphical and calculation parameters.
#' @param y.risk0 Additional graphical and calculation parameters.
#' @param y.risk1 Additional graphical and calculation parameters.
#' @param show.Y.axis Additional graphical and calculation parameters.
#' @param cex_Yaxis Additional graphical and calculation parameters.
#' @param add.segment Additional graphical and calculation parameters.
#' @param risk.add Additional graphical and calculation parameters.
#' @param xmin Additional graphical and calculation parameters.
#' @param xmax Additional graphical and calculation parameters.
#' @param x.truncate Additional graphical and calculation parameters.
#' @param time.zero Additional graphical and calculation parameters.
#' @param prob.points Additional graphical and calculation parameters.
#' @return Invisibly returns NULL. Used for plotting side effects.
#' @importFrom graphics plot lines legend title axis box abline par points
#' @export

KM_plot_2sample_weighted_counting <- function(
    dfcount, show.cox = TRUE, cox.cex = 0.725, show.logrank = FALSE,
    logrank.cex = 0.725, cox.eps = 0.001, lr.eps = 0.001, show_arm_legend = TRUE,
    arms = c("treat", "control"), put.legend.arms = "left", stop.onerror = TRUE,
    check.KM = TRUE, put.legend.cox = "topright", put.legend.lr = "topleft",
    lr.digits = 2, cox.digits = 2, tpoints.add = c(0), by.risk = NULL, Xlab = "time",
    Ylab = "proportion surviving", col.0 = "black", col.1 = "blue", show.med = TRUE,
    med.digits = 2, med.font = 4, conf.int = FALSE, conf_level = 0.95,
    choose_ylim = FALSE, arm.cex = 0.7,
    quant = 0.5, qlabel = "median =", med.cex = 0.725,
    ymed.offset = 0.10, xmed.fraction = 0.80,
    risk.cex = 0.725,
    ltys = c(1, 1), lwds = c(1, 1), censor.mark.all = TRUE, censor.cex = 0.5,
    show.ticks = TRUE, risk.set = TRUE, ymin = 0, ymax = 1, ymin.del = 0.035,
    ymin2 = NULL, risk_offset = 0.125, risk_delta = 0.05, y.risk0 = NULL,
    show.Y.axis = TRUE, cex_Yaxis = 1, y.risk1 = NULL, add.segment = FALSE,
    risk.add = NULL, xmin = 0, xmax = NULL, x.truncate = NULL, time.zero = 0.0,
    prob.points = NULL
) {
  # Validate input
  stopifnot(is.list(dfcount))
  required_fields <- c("at_points", "risk.points", "risk.points.label", "idv0", "surv0", "sig2_surv0", "idx0", "cens0", "riskpoints0",
                       "idv1", "surv1", "sig2_surv1", "idx1", "cens1", "riskpoints1")
  missing_fields <- setdiff(required_fields, names(dfcount))
  if (length(missing_fields) > 0) stop("Missing fields in dfcount: ", paste(missing_fields, collapse = ", "))

  at_points <- dfcount$at_points
  risk.points <- dfcount$risk.points
  risk.points.label <- dfcount$risk.points.label

  idv0 <- dfcount$idv0
  S0.KM <- dfcount$surv0
  se0.KM <- sqrt(dfcount$sig2_surv0)
  idx0 <- dfcount$idx0
  cens0 <- dfcount$cens0
  rpoints0 <- dfcount$riskpoints0

  idv1 <- dfcount$idv1
  S1.KM <- dfcount$surv1
  se1.KM <- sqrt(dfcount$sig2_surv1)
  idx1 <- dfcount$idx1
  cens1 <- dfcount$cens1
  rpoints1 <- dfcount$riskpoints1

  # KM curve checks
  if (check.KM) {
    check_km_curve(S0.KM, "Group 0", stop_on_error = stop.onerror)
    check_km_curve(S1.KM, "Group 1", stop_on_error = stop.onerror)
  }

  plot_km_curves_counting(
    at_points, S0.KM, idx0, idv0, S1.KM, idx1, idv1, col.0, col.1, ltys, lwds, Xlab, Ylab,
    ylim = c(ifelse(is.null(ymin2), ymin - risk_offset, ymin2), ymax),
    risk.points = risk.points, risk.points.label = risk.points.label,
    xlim = c(xmin, ifelse(is.null(xmax), max(c(at_points)), xmax)),
    show.ticks = show.ticks, cens0 = cens0, cens1 = cens1,
    censor.cex = censor.cex, conf.int = conf.int, conf_level = conf_level, se1.KM = se1.KM, se0.KM = se0.KM
  )

  if (risk.set) {
    add_risk_table(risk.points, rpoints0, rpoints1, col.0, col.1, risk.cex, ymin, risk_offset, risk_delta, y.risk0, y.risk1)
  }

  abline(h = ymin - ymin.del, lty = 1, col = 1)
  add_legends(dfcount, show.cox, cox.cex, put.legend.cox, show.logrank, logrank.cex, put.legend.lr, show_arm_legend, arms, col.1, col.0, ltys, lwds, arm.cex, put.legend.arms)

  if (show.med && !is.null(dfcount$quantile_results)) {
    add_median_annotation(dfcount$quantile_results, med.digits, med.cex, med.font, xmed.fraction, ymed.offset)
  }
}


#' Plot Kaplan-Meier Survival Difference Curves with Subgroups and Confidence Bands
#'
#' Plots the difference in Kaplan-Meier survival curves between two groups, optionally including simultaneous confidence bands and subgroup curves. Also displays risk tables for the overall population and specified subgroups.
#'
#' @param df Data frame containing survival data.
#' @param tte.name Name of the time-to-event column.
#' @param event.name Name of the event indicator column (0/1).
#' @param treat.name Name of the treatment group column (0/1).
#' @param weight.name Optional name of the weights column.
#' @param sg_labels Character vector of subgroup definitions (as logical expressions).
#' @param ltype Line type for curves (default: "s").
#' @param lty Line style for curves (default: 1).
#' @param draws Number of draws for resampling (default: 20).
#' @param lwd Line width for curves (default: 2).
#' @param sg_colors Colors for subgroup curves.
#' @param color Color for confidence band polygon (default: "lightgrey").
#' @param ymax.pad Padding for y-axis limits.
#' @param ymin.pad Padding for y-axis limits.
#' @param taus Vector for time truncation (default: c(-Inf, Inf)).
#' @param yseq_length Number of y-axis ticks (default: 5).
#' @param cex_Yaxis Text size for axis.
#' @param risk_cex Text size for risk table.
#' @param by.risk Interval for risk table time points (default: 6).
#' @param risk.add Additional time points for risk table.
#' @param xmax Additional graphical and calculation parameters.
#' @param ymin Additional graphical and calculation parameters.
#' @param ymax Additional graphical and calculation parameters.
#' @param ymin.del Additional graphical and calculation parameters.
#' @param y.risk1 Additional graphical and calculation parameters.
#' @param y.risk2 Additional graphical and calculation parameters.
#' @param ymin2 Additional graphical and calculation parameters.
#' @param risk_offset Additional graphical and calculation parameters.
#' @param risk.pad Additional graphical and calculation parameters.
#' @param risk_delta Additional graphical and calculation parameters.
#' @param tau_add Additional graphical and calculation parameters.
#' @param time.zero.pad Additional graphical and calculation parameters.
#' @param time.zero.label Additional graphical and calculation parameters.
#' @param xlabel Additional graphical and calculation parameters.
#' @param ylabel Additional graphical and calculation parameters.
#' @param Maxtau Additional graphical and calculation parameters.
#' @param seedstart Additional graphical and calculation parameters.
#' @param ylim Additional graphical and calculation parameters.
#' @param draws.band Number of draws for simultaneous confidence bands (default: 20).
#' @param qtau Quantile for time range in simultaneous bands (default: 0.025).
#' @param show_resamples Logical; whether to plot resampled curves (default: FALSE).
#' @param modify_tau Logical; restrict time range for bands.
#' @return (Invisible) list containing KM_diff results, time points, subgroup curves, risk tables, and confidence intervals.
#' @importFrom graphics plot lines legend title axis box polygon matplot
#' @export

plotKM.band_subgroups <- function(
    df, tte.name, event.name, treat.name, weight.name = NULL,
    sg_labels = NULL,
    ltype = "s", lty = 1, draws = 20, lwd = 2,
    sg_colors = NULL, color="lightgrey",
    ymax.pad = 0.01, ymin.pad = -0.01,
    taus = c(-Inf, Inf), yseq_length = 5, cex_Yaxis = 0.8, risk_cex = 0.8,
    by.risk = 6, risk.add = NULL, xmax = NULL, ymin = NULL, ymax = NULL, ymin.del = 0.035,
    y.risk1 = NULL, y.risk2 = NULL, ymin2 = NULL, risk_offset = NULL, risk.pad = 0.01,
    risk_delta = 0.0275, tau_add = NULL, time.zero.pad = 0, time.zero.label = 0.0,
    xlabel = NULL, ylabel = NULL, Maxtau = NULL, seedstart = 8316951,
    ylim = NULL, draws.band = 20, qtau = 0.025, show_resamples = FALSE, modify_tau = FALSE
) {
  # Input checks

  required_cols <- c(tte.name, event.name, treat.name, weight.name)
  missing_cols <- setdiff(required_cols, names(df))
  if (length(missing_cols) > 0) {
    stop(paste("Missing required columns in df:", paste(missing_cols, collapse = ", ")))
  }

  if (length(sg_labels) != length(sg_colors)) stop("SG labels and colors do not match")
  if (is.null(risk_offset)) risk_offset <- (1 + length(sg_labels)) * risk_delta
  if (is.null(ylabel)) ylabel <- expression(delta(t) == hat(S)[1](t) - hat(S)[0](t))
  if (is.null(xlabel)) xlabel <- "Months"

  if (!is.numeric(df[[tte.name]])) stop("Time-to-event column must be numeric.")
  if (!all(df[[event.name]] %in% c(0, 1))) stop("Event column must be binary (0/1).")
  if (!all(df[[treat.name]] %in% c(0, 1))) stop("Treatment column must be binary (0/1).")
  if (!is.null(weight.name) && !all(df[[weight.name]] >= 0)) stop("Weights must be non-negative.")


  # ITT
  Y <- df[[tte.name]]
  E <- df[[event.name]]
  Treat <- df[[treat.name]]
  maxe_0 <- max(Y[E == 1 & Treat == 0])
  maxe_1 <- max(Y[E == 1 & Treat == 1])
  max_tau <- min(c(maxe_0, maxe_1))
  if (!is.null(Maxtau)) max_tau <- Maxtau
  if (!is.null(tau_add)){
  if(tau_add > max_tau){ warning("May not want to estimate differences beyond max_tau (recommend tau_add = NULL or below max_tau)")
  }
    max_tau <- max(c(max_tau, tau_add))
  }

# For risk-set display want to include max_tau_risk
  if(max_tau < 1){
    max_tau_risk <- ceiling(max_tau)}
  else { max_tau_risk <- floor(max_tau)
  }
  # Add by.risk points
  riskpoints <- seq(0, max_tau, by = by.risk)

  # Truncate the observed timepoints at max_tau
  atpoints <- Y[Y <= max_tau]
  at_points <- sort(unique(c(atpoints, max_tau, tau_add, riskpoints)))

  fit <- KM_diff(
    df = df, tte.name = tte.name, event.name = event.name, weight.name = weight.name,
    treat.name = treat.name, at_points = at_points, alpha = 0.05, risk.points = riskpoints,
    modify_tau = modify_tau,
    draws = draws, seedstart = seedstart, draws.band = draws.band, qtau = qtau, show_resamples = show_resamples
  )

  # Add max_tau to fit
  fit$max_tau <- max_tau
  fit$risk.points <- riskpoints
  at_points <- fit$at_points
  dhat <- fit$dhat
  # pointwise CIs
  l0_pw <- fit$lower
  u0_pw <- fit$upper
  # simultaneous bands (draws.band > 0)
  l0_sb <- fit$sb_lower
  u0_sb <- fit$sb_upper

  risk.points <- round(seq(time.zero.label, max(at_points), by = by.risk))
  risk.points <- sort(unique(c(risk.points, risk.add)))
  risk.points <- c(time.zero.label - time.zero.pad, risk.points)

  risk.points <- risk.points[which(risk.points <= max(fit$at_points))]
  # include max_tau_risk in risk-set display
  risk.points <- unique(c(risk.points, max_tau_risk))

  # Total risk for ITT
  risk0 <- colSums(outer(Y, risk.points, FUN = ">="))
  risk0 <- round(risk0)
  # Subgroups (SG) defined by sg_labels
  Dsg_mat <- NULL
  Rsg_mat <- NULL
  S0_mat <- NULL
  S1_mat <- NULL

  if(length(sg_labels) > 0){
  sg_flags <- lapply(sg_labels, function(expr) with(df, eval(parse(text = expr))))
  # Preallocate matrices
  Dsg_mat <- matrix(NA, nrow = length(at_points), ncol = length(sg_labels))
  Rsg_mat <- matrix(NA, nrow = length(risk0), ncol = length(sg_labels))
  S0_mat <- matrix(NA, nrow = length(at_points), ncol = length(sg_labels))
  S1_mat <- matrix(NA, nrow = length(at_points), ncol = length(sg_labels))

  # Here we only want point estimates (SE's not considered for subgroup plot)
  for (i in seq_along(sg_flags)) {
    df_sg <- df[sg_flags[[i]], ]
    Y_sg <- df_sg[[tte.name]]
    E_sg <- df_sg[[event.name]]
    Treat_sg <- df_sg[[treat.name]]
    res <- KM_diff(
      df = df_sg, tte.name = tte.name, event.name = event.name, weight.name = weight.name,
      treat.name = treat.name, at_points = at_points, alpha = 0.05, risk.points = risk.points,
      draws = 0, draws.band = 0
    )

    Dsg_mat[, i] <- res$dhat
    rr <- colSums(outer(Y_sg, risk.points, FUN = ">="))
    Rsg_mat[, i] <- round(rr)
    S0_mat[, i] <- res$surv0
    S1_mat[, i] <- res$surv1
  }
  }

  x <- at_points
  mean.value <- dhat
  lower <- l0_pw
  upper <- u0_pw

  dhats_all <- dhat
  if(length(sg_labels)>0){
  dhats_all <- cbind(dhat,Dsg_mat)
  }

  if (is.null(ymax)) {
    ymax <- max(dhats_all, na.rm = TRUE) + ymax.pad
    ymax <- round(max(c(u0_pw, u0_sb, ymax), na.rm = TRUE), 2)
  }
  if (is.null(ymin)) {
    ymin <- min(dhats_all, na.rm = TRUE) - ymin.pad
    ymin <- round(min(c(l0_pw, l0_sb,ymin), na.rm = TRUE), 2)
  }
  if (is.null(ymin2)) ymin2 <- ymin - risk_offset

  yrisks <- seq(ymin2, ymin - (ymin.del + risk.pad), length.out = 1 + length(sg_labels))
  yrisks <- sort(yrisks, decreasing = TRUE)

  if (is.null(ylim)) ylim <- c(ymin2, ymax)

  plot(
    x[order(x)], mean.value[order(x)], type = "n", axes = FALSE, xlab = xlabel, lty = lty,
    ylab = ylabel, ylim = ylim, cex.lab = cex_Yaxis
  )
  if(draws.band ==0){
  polygon(
    c(x[order(x)], rev(x[order(x)])),
    c(l0_pw[order(x)], rev(u0_pw[order(x)])),
    col = color, border = FALSE
  )
  }

  if(draws.band > 0){
    polygon(
      c(x[order(x)], rev(x[order(x)])),
      c(l0_sb[order(x)], rev(u0_sb[order(x)])),
      col = color, border = FALSE
    )
  lines(x[order(x)], l0_pw[order(x)], lty=2, type="s")
  lines(x[order(x)], u0_pw[order(x)], lty=2, type="s")
  }

  lines(x[order(x)], mean.value[order(x)], lty = lty, lwd = lwd, type = ltype)
  abline(h = ymin - ymin.del, lty = 1, col = "black")
  abline(h = time.zero.label, lty = 2, col = "black", lwd = 0.5)

  if(length(sg_labels)>0){

    matplot(
      x, Dsg_mat, type = ltype, lty = lty, lwd = lwd,
      col = sg_colors, add = TRUE
    )


  }

  d_minmax <- round((ymax-ymin)/10,2)
  by_ypoints <- min(c(d_minmax, 0.2))

  if ((ymax - ymin2) <= 0.5) {
  ypoints <- sort(c(time.zero.label, seq(ymin,ymax, by = by_ypoints)))
  } else {
  ypoints <- seq(ymin, ymax, by = by_ypoints)
  }
  ypoints <- sort(c(time.zero.label, ypoints))
  ypoints <- format(round(ypoints,3), scientific = FALSE)

  axis(2, at = ypoints, cex.axis = cex_Yaxis, las = 1)
  risk.points.label <- as.character(c(time.zero.label, risk.points[-1]))

  axis(1, at = risk.points, labels = risk.points.label, cex.axis = cex_Yaxis)
  box()

  # ITT risk
  text(risk.points, yrisks[1], risk0, col = "black", cex = risk_cex)

  if(length(sg_labels)>0){
   for (sg in seq_along(sg_labels)) {
    text(risk.points, yrisks[sg + 1], Rsg_mat[, sg], col = sg_colors[sg], cex = risk_cex)
  }
    }
  invisible(list( fit_itt = fit,
    xpoints = at_points, Dhat_subgroups = Dsg_mat, s0_subgroups = S0_mat, s1_subgroups = S1_mat,
    rpoints = risk.points, Risk_subgroups = Rsg_mat, mean = mean.value, lower = lower, upper = upper
  ))
}
