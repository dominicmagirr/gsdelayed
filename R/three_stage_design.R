#' Three-stage design
#'
#' Three-stage design for a modestly-weighted log-rank test
#'
#' @param t_star Parameter of the modestly-weighted log-rank test. Setting t_star=0 corresponds to a standard log-rank test.
#' @param model The piecewise hazard model.
#'   A list containing the \code{change_points} and \code{lambdas}.
#' @param recruitment List of recruitment information.
#'   Containing \enumerate{
#'                 \item Sample size on control, \code{n_0}
#'                 \item Sample size on experimental, \code{n_1}
#'                 \item Recruitment period, \code{r_period}
#'                 \item Recruitment parameter for power model, \code{k}
#'               }
#' @param dco_int_1 The time of the first interim analysis in time-units since start of the trial.
#' @param dco_int_2 The time of the second interim analysis in time-units since start of the trial.
#' @param dco_final The time of the final analysis in time-units since start of the trial.
#' @param events_int_1 May be specified instead of \code{dco_int_1}. The number of events that triggers the first interim analysis.
#' @param events_int_2 May be specified instead of \code{dco_int_2}. The number of events that triggers the second interim analysis.
#' @param events_final May be specified instead of \code{dco_final}. The number of events that triggers the final analysis.
#' @param alpha_one_sided One-sided alpha level.
#' @param alpha_spending_f The alpha-spending function. Default is the Lan-DeMets O'Brien-Fleming function.
#' @param length_t Number of cutpoints to use when approximating the distribution of the MWLRT-statistic. Default is 18. Can be increased for greater accuracy.
#' @param{rho} rho parameter in a Fleming-Harrington test. Default is NULL. Only used if F-H test used instead of MWLRT.
#' @param{gamma} gamma parameter in a Fleming-Harrington test. Default is NULL. Only used if F-H test used instead of MWLRT.
#' @return A list describing the design.
#' @export


three_stage_design <- function(t_star = NULL,
                               model,
                               recruitment,
                               dco_int_1,
                               dco_int_2,
                               dco_final,
                               events_int_1 = NULL,
                               events_int_2 = NULL,
                               events_final = NULL,
                               alpha_one_sided = 0.025,
                               alpha_spend_f = ldobf,
                               length_t = 18,
                               rho = NULL,
                               gamma = NULL){

  if (all(is.null(c(t_star, rho, gamma)))) stop("Either t_star or rho, gamma must be specified")
  if (is.null(t_star) && is.null(rho)) stop("rho and gamma must be specified")
  if (is.null(t_star) && is.null(gamma)) stop("rho and gamma must be specified")

  if (is.null(dco_int_1) && is.null(events_int_1)) stop("Either dco_int_1 or events_int_1 must be specified.")
  if (is.null(dco_int_2) && is.null(events_int_2)) stop("Either dco_int_2 or events_int_2 must be specified.")
  if (is.null(dco_final) && is.null(events_final)) stop("Either dco_final or events_final must be specified.")
  if (is.null(dco_int_1)){

    dco_int_1 <- expected_dco_two_arm(total_events = events_int_1,
                                      recruitment = recruitment,
                                      model = model)
  }
  if (is.null(dco_int_2)){

    dco_int_2 <- expected_dco_two_arm(total_events = events_int_2,
                                      recruitment = recruitment,
                                      model = model)
  }
  if (is.null(dco_final)){

    dco_final <- expected_dco_two_arm(total_events = events_final,
                                      recruitment = recruitment,
                                      model = model)
  }
  #####################################

  final_analysis <- ncp_power(t_star = t_star,
                              rho = rho,
                              gamma = gamma,
                              model = model,
                              recruitment = recruitment,
                              dco = dco_final,
                              length_t = length_t)


  ######################################

  interim_analysis_2 <- ncp_power(t_star = t_star,
                                  rho = rho,
                                  gamma = gamma,
                                  model = model,
                                  recruitment = recruitment,
                                  dco = dco_int_2,
                                  length_t = length_t)

  ######################################

  interim_analysis_1 <- ncp_power(t_star = t_star,
                                  rho = rho,
                                  gamma = gamma,
                                  model = model,
                                  recruitment = recruitment,
                                  dco = dco_int_1,
                                  length_t = length_t)


  ######################################
  info_frac_2 <- interim_analysis_2$var_u / final_analysis$var_u
  info_frac_1 <- interim_analysis_1$var_u / final_analysis$var_u
  ######################################


  cumulative_spend <- alpha_spend_f(c(info_frac_1, info_frac_2, 1),
                                    alpha_one_sided = alpha_one_sided)



  #######################################
  crit_1 <- qnorm(1 - cumulative_spend[1])
  crit_2 <- uniroot(find_crit_2,
                    c(0, 10),
                    crit_1 = crit_1,
                    info_frac = interim_analysis_1$var_u / interim_analysis_2$var_u,
                    alpha_one_sided = cumulative_spend[2])$root


  crit_3 <- uniroot(find_crit_3, c(0, 10),
                    crit_1 = crit_1,
                    crit_2 = crit_2,
                    info_frac_1 = info_frac_1,
                    info_frac_2 = info_frac_2,
                    alpha_one_sided = alpha_one_sided)$root
  ## check power
  ## possibly loop through again...

  overall_power <- 1 - mvtnorm::pmvnorm(lower = c(-Inf, -Inf, -Inf),
                                        upper = c(crit_1, crit_2, crit_3),
                                        mean = -c(interim_analysis_1$ncp,
                                                  interim_analysis_2$ncp,
                                                  final_analysis$ncp),
                                        sigma = matrix(c(1, sqrt(info_frac_1 / info_frac_2), sqrt(info_frac_1),
                                                         sqrt(info_frac_1 / info_frac_2), 1, sqrt(info_frac_2),
                                                         sqrt(info_frac_1), sqrt(info_frac_2), 1), nrow = 3))[1]


  ### Expected duration of the trial
  p_first_int_stop_alt <- pnorm(-crit_1, mean = interim_analysis_1$ncp)
  p_first_int_stop_null_approx <- pnorm(-crit_1, mean = 0)


  p_first_or_second_int_stop_alt <- 1 - mvtnorm::pmvnorm(lower = c(-Inf, -Inf),
                                                         upper = c(crit_1, crit_2),
                                                         mean = -c(interim_analysis_1$ncp, interim_analysis_2$ncp),
                                                         sigma = matrix(c(1, sqrt(info_frac_1 / info_frac_2),
                                                                          sqrt(info_frac_1 / info_frac_2), 1), nrow = 2))[1]

  p_first_or_second_int_stop_null_approx <- 1 - mvtnorm::pmvnorm(lower = c(-Inf, -Inf),
                                                                 upper = c(crit_1, crit_2),
                                                                 mean =  c(0, 0),
                                                                 sigma = matrix(c(1, sqrt(info_frac_1 / info_frac_2),
                                                                                  sqrt(info_frac_1 / info_frac_2), 1), nrow = 2))[1]


  p_second_int_stop_alt = p_first_or_second_int_stop_alt - p_first_int_stop_alt
  p_second_int_stop_null_approx = p_first_or_second_int_stop_null_approx - p_first_int_stop_null_approx


  expected_t_alt <- dco_int_1 * p_first_int_stop_alt +
                    dco_int_2 * p_second_int_stop_alt +
                    dco_final * (1 - p_first_or_second_int_stop_alt)


  expected_t_null_approx <- dco_int_1 * p_first_int_stop_null_approx +
                            dco_int_2 * p_second_int_stop_null_approx +
                            dco_final * (1 - p_first_or_second_int_stop_null_approx)

  #######################################################
  list(critical_values = -c(crit_1, crit_2, crit_3),
       p_first_int_stop = c(under_alternative = p_first_int_stop_alt,
                            under_null_approx = p_first_int_stop_null_approx),
       p_second_int_stop = c(under_alternative = p_second_int_stop_alt,
                             under_null_approx = p_second_int_stop_null_approx),
       expected_t = c(under_alternative = expected_t_alt,
                      under_null_approx = expected_t_null_approx),
       overall_power = overall_power,
       n_events = c(interim_1 = interim_analysis_1$total_events,
                           interim_2 = interim_analysis_2$total_events,
                           final = final_analysis$total_events),
       var_u = c(interim_1 = interim_analysis_1$var_u,
                 interim_2 = interim_analysis_2$var_u,
                 final = final_analysis$var_u),
       ncp_z = c(interim_1 = interim_analysis_1$ncp,
                 interim_2 = interim_analysis_2$ncp,
                 final = final_analysis$ncp),
       t_star = t_star,
       rho = rho,
       gamma = gamma,
       model = model,
       recruitment = recruitment,
       dco = c(dco_int_1, dco_int_2, dco_final),
       alpha_spend_f = alpha_spend_f,
       alpha_one_sided = alpha_one_sided)

}

