#' Single-stage design
#'
#' Single-stage design for a modestly-weighted log-rank test
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
#' @param dco_final The time of the final analysis in time-units since start of the trial.
#' @param events_final May be specified instead of \code{dco_final}. The number of events that triggers the final analysis.
#' @param alpha_one_sided One-sided alpha level.
#' @param length_t Number of cutpoints to use when approximating the distribution of the MWLRT-statistic. Default is 18. Can be increased for greater accuracy.
#' @return A list describing the design.
#' @export

single_stage_design <- function(t_star,
                                model,
                                recruitment,
                                dco_final,
                                events_final = NULL,
                                alpha_one_sided = 0.025,
                                length_t = 18){

  if (is.null(dco_final) && is.null(events_final)) stop("Either dco_final or events_final must be specified.")

  if (is.null(dco_final)){

    dco_final <- expected_dco_two_arm(total_events = events_final,
                                      recruitment = recruitment,
                                      model = model)
  }


  final_analysis <- ncp_power(t_star = t_star,
                              model = model,
                              recruitment = recruitment,
                              dco = dco_final,
                              length_t = length_t)



  crit_1 <- qnorm(alpha_one_sided)

  overall_power <- pnorm(crit_1, mean = final_analysis$ncp)

  list(critical_values = crit_1,
       expected_t = dco_final,
       overall_power = overall_power,
       expected_events = final_analysis$total_events,
       var_u = final_analysis$var_u,
       t_star = t_star,
       model = model,
       recruitment = recruitment)

}




