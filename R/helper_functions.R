#' Lan-DeMets O'Brien-Fleming alpha-spending function
#'
#' Find the cumulative alpha-spend at various information times according to Lan-DeMets OBF function.
#'
#' @param t Information fraction. Vector of numbers between 0 and 1.
#' @param alpha_one_sided One-sided alpha level.
#' @return A vector of cumulative alpha spend.
#' @export
#'
ldobf <- function(t, alpha_one_sided = 0.025) 1 - 1 * pnorm(qnorm(1 - alpha_one_sided) / sqrt(t))



#' Critical value at interim analysis
#'
#' Find the critical value at the interim analysis of a two-stage design
#'
#' @param var_u_int The observed variance of the weighted log-rank statistic at the interim.
#' @param design The design object, created using \code{two_stage_design}
#' @param alpha_spend_f The alpha-spending function. Default is the Lan-DeMets O'Brien-Fleming function.
#' @param alpha_one_sided One-sided alpha level.
#' @return The critical value for the interim z-statistic.
#' @export
#'
crit_1_of_2 <- function(var_u_int,
                        design,
                        alpha_spend_f = ldobf,
                        alpha_one_sided = 0.025){

  qnorm(alpha_spend_f(var_u_int / design$var_u[2], alpha_one_sided))

}
#' Critical value at final analysis
#'
#' Find the critical value at the final analysis of a two-stage design
#'
#' @param var_u_int The observed variance of the weighted log-rank statistic at the interim.
#' @param var_u_final The observed variance of the weighted log-rank statistic at the final analysis.
#' @param design The design object, created using \code{two_stage_design}
#' @param alpha_spend_f The alpha-spending function. Default is the Lan-DeMets O'Brien-Fleming function.
#' @param alpha_one_sided One-sided alpha level.
#' @return The critical value for the final z-statistic.
#' @export
#'
crit_2_of_2 <- function(var_u_int,
                        var_u_final,
                        design,
                        alpha_spend_f = ldobf,
                        alpha_one_sided = 0.025){

  crit_1 <- qnorm(alpha_spend_f(var_u_int / design$var_u[2], alpha_one_sided))


  -uniroot(find_crit_2, c(0,100),
           crit_1 = -crit_1,
           info_frac = var_u_int / var_u_final,
           alpha_one_sided = alpha_one_sided)$root

}


#' Critical value at first interim analysis
#'
#' Find the critical value at the first interim analysis of a three-stage design
#'
#' @param var_u_int_1 The observed variance of the weighted log-rank statistic at the first interim.
#' @param design The design object, created using \code{three_stage_design}
#' @param alpha_spend_f The alpha-spending function. Default is the Lan-DeMets O'Brien-Fleming function.
#' @param alpha_one_sided One-sided alpha level.
#' @return The critical value for the first interim z-statistic.
#' @export
#'
crit_1_of_3 <- function(var_u_int_1,
                        design,
                        alpha_spend_f = ldobf,
                        alpha_one_sided = 0.025){

  qnorm(alpha_spend_f(var_u_int_1 / design$var_u[3], alpha_one_sided))

}


#' Critical value at second interim analysis
#'
#' Find the critical value at the second interim analysis of a three-stage design
#'
#' @param var_u_int_2 The observed variance of the weighted log-rank statistic at the second interim.
#' @param var_u_int_1 The observed variance of the weighted log-rank statistic at the first interim.
#' @param crit_1 The critical value for the z-statistic at the first interim, found via \code{crit_1_of_3}.
#' @param design The design object, created using \code{three_stage_design}
#' @param alpha_spend_f The alpha-spending function. Default is the Lan-DeMets O'Brien-Fleming function.
#' @param alpha_one_sided One-sided alpha level.
#' @return The critical value for the second interim z-statistic.
#' @export
#'
crit_2_of_3 <- function(var_u_int_2,
                        var_u_int_1,
                        crit_1,
                        design,
                        alpha_spend_f = ldobf,
                        alpha_one_sided = 0.025){

  alpha_spend_2 <- alpha_spend_f(var_u_int_2 / design$var_u[3], alpha_one_sided)

  -uniroot(find_crit_2, c(0,100),
           crit_1 = -crit_1,
           info_frac = var_u_int_1 / var_u_int_2,
           alpha_one_sided = alpha_spend_2)$root

}


#' Critical value at final analysis
#'
#' Find the critical value at the final analysis of a three-stage design
#'
#' @param var_u_final The observed variance of the weighted log-rank statistic at the final analysis.
#' @param var_u_int_2 The observed variance of the weighted log-rank statistic at the second interim.
#' @param var_u_int_1 The observed variance of the weighted log-rank statistic at the first interim.
#' @param crit_1 The critical value for the z-statistic at the first interim, found via \code{crit_1_of_3}.
#' @param crit_2 The critical value for the z-statistic at the second interim, found via \code{crit_2_of_3}.
#' @param design The design object, created using \code{three_stage_design}
#' @param alpha_spend_f The alpha-spending function. Default is the Lan-DeMets O'Brien-Fleming function.
#' @param alpha_one_sided One-sided alpha level.
#' @return The critical value for the final z-statistic.
#' @export
#'
crit_3_of_3 <- function(var_u_final,
                        var_u_int_2,
                        var_u_int_1,
                        crit_1,
                        crit_2,
                        design,
                        alpha_spend_f = ldobf,
                        alpha_one_sided = 0.025){




  -uniroot(find_crit_3, c(0,100),
           crit_1 = -crit_1,
           crit_2 = -crit_2,
           info_frac_1 = var_u_int_1 / var_u_final,
           info_frac_2 = var_u_int_2 / var_u_final,
           alpha_one_sided = alpha_one_sided)$root

}
##############################################
#############################################
find_crit_3 <- function(x, crit_1, crit_2, info_frac_1, info_frac_2, alpha_one_sided = 0.025){

  1 - mvtnorm::pmvnorm(lower = c(-Inf, -Inf, -Inf),
                       upper = c(crit_1, crit_2, x),
                       sigma = matrix(c(1, sqrt(info_frac_1 / info_frac_2), sqrt(info_frac_1),
                                        sqrt(info_frac_1 / info_frac_2), 1, sqrt(info_frac_2),
                                        sqrt(info_frac_1), sqrt(info_frac_2), 1), nrow = 3))[1] - alpha_one_sided

}

############################################
find_crit_2 <- function(x, crit_1, info_frac, alpha_one_sided = 0.025){
  1 - mvtnorm::pmvnorm(lower = c(-Inf, -Inf),
                       upper = c(crit_1, x),
                       sigma = matrix(c(1, sqrt(info_frac),
                                        sqrt(info_frac), 1), nrow = 2))[1] - alpha_one_sided
}


############################################
#' Stage-wise p-value at 2nd analysis
#'
#' Find the stage-wise one-sided value, when a trial has crossed efficacy boundary at the second analysis.
#'
#' @param crit_1 Critical value at first analysis. Expecting a negative number.
#' @param z_2 The observed z-statistic at the second analysis..
#' @param v_1 The observed variance of the u-statistic at the first analysis.
#' @param v_2 The observed variance of the u-statistic at the second analysis.
#' @export
#'

stagewise_p_2 <- function(crit_1,
                          z_2,
                          v_1,v_2){

  info_frac <- v_1 / v_2

  1 - mvtnorm::pmvnorm(lower = c(-Inf, -Inf),
                       upper = c(-crit_1, -z_2),
                       sigma = matrix(c(1, sqrt(info_frac),
                                        sqrt(info_frac), 1), nrow = 2))[1]

}

############################################
#' Stage-wise p-value at 3rd analysis
#'
#' Find the stage-wise one-sided value, when a trial has crossed efficacy boundary at the third analysis.
#'
#' @param crit_1 Critical value at first analysis. Expecting a negative number.
#' @param crit_2 Critical value at second analysis. Expecting a negative number.
#' @param z_3 The observed z-statistic at the third analysis..
#' @param v_1 The observed variance of the u-statistic at the first analysis.
#' @param v_2 The observed variance of the u-statistic at the second analysis.
#' @param v_3 The observed variance of the u-statistic at the third analysis.
#' @export
#'
stagewise_p_3 <- function(crit_1,
                          crit_2,
                          z_3,
                          v_1,v_2,v_3){

  info_frac_1 <- v_1 / v_3
  info_frac_2 <- v_2 / v_3

  1 - mvtnorm::pmvnorm(lower = c(-Inf, -Inf, -Inf),
                       upper = c(-crit_1, -crit_2, -z_3),
                       sigma = matrix(c(1, sqrt(info_frac_1 / info_frac_2), sqrt(info_frac_1),
                                        sqrt(info_frac_1 / info_frac_2), 1, sqrt(info_frac_2),
                                        sqrt(info_frac_1), sqrt(info_frac_2), 1), nrow = 3))[1]

}



####### simulate from a piece-wise exponential distribution
t_piecewise_exp <- function(n = 10,
                            change_points = c(6, 12),
                            lambdas = c(log(2) / 9, log(2) / 9, log(2) / 9)){

  t_lim <- matrix(rep(c(diff(c(0, change_points)), Inf), each = n), nrow = n)
  t_sep <- do.call(cbind, purrr::map(lambdas, rexp, n = n))
  which_cells <- t(apply(t_sep < t_lim, 1, function(x){
    rep(c(T,F), c(min(which(x)), length(x) - min(which(x))))
  } ))
  rowSums(pmin(t_sep, t_lim) * which_cells)
}
##########################################################


#' Simulate uncensored data
#'
#' Simulate data as if there is no administrative censoring, i.e., all events observed.
#'
#' @param var_u_int The observed variance of the weighted log-rank statistic at the interim.
#' @param model The piecewise hazard model.
#'   A list containing the \code{change_points} and \code{lambdas}.
#' @param recruitment List of recruitment information.
#'   Containing \enumerate{
#'                 \item Sample size on control, \code{n_0}
#'                 \item Sample size on experimental, \code{n_1}
#'                 \item Recruitment period, \code{r_period}
#'                 \item Recruitment parameter for power model, \code{k}
#'               }
#' @return A data frame containing simulated data.
#' @export
#'

sim_t_uncensored <- function(model,
                             recruitment){

  rec_0 <- recruitment$r_period * runif(recruitment$n_0) ^ (1 / recruitment$k)
  rec_1 <- recruitment$r_period * runif(recruitment$n_1) ^ (1 / recruitment$k)

  time_0 <- t_piecewise_exp(recruitment$n_0, model$change_points, model$lambdas_0)
  time_1 <- t_piecewise_exp(recruitment$n_1, model$change_points, model$lambdas_1)

  data.frame(time = c(time_0, time_1),
             rec = c(rec_0, rec_1),
             group = rep(c("control", "experimental"), c(recruitment$n_0, recruitment$n_1)))

}

#' Apply data cut-off
#'
#' Apply a data cut-off to an uncensored data set.
#'
#' @param df Data frame with uncensored data. Generated from \code{sim_t_uncensored}
#' @param dco Data cut-off time. In time-units since start of the trial.
#' @param events May be specified instead of \code{dco}. The number of events that triggers data cut-off.
#' @return A data frame containing simulated data. Snapshot taken at data cut-off.
#' @export
#'


apply_dco <- function(df,
                      dco = NULL,
                      events = NULL){

  if (is.null(dco) && is.null(events)) stop("Must specify either dco or events")

  df$cal_time <- df$time + df$rec

  if (is.null(dco)){
    dco <- sort(df$cal_time)[events]
  }

  df_dco <- df[df$rec < dco, ]
  df_dco$event <-  df_dco$cal_time <= dco
  df_dco$time <- pmin(df_dco$time, dco - df_dco$rec)
  df_dco$dco <- dco

  df_dco

}

