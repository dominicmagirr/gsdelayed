############### model extend #######################

model_extend <- function(model,
                         max_t,
                         length_t){

  change_points_plus <- unique(sort(c(seq(0,
                                          max_t,
                                          length.out = length_t),
                                      model$change_points)))

  which_lambda <- purrr::map_dbl(change_points_plus,
                                 function(x) sum(x >= model$change_points)) + 1


  list(change_points = change_points_plus[-1],
       lambdas_0 = model$lambdas_0[which_lambda],
       lambdas_1 = model$lambdas_1[which_lambda])


}


########## \bar{S} ############

s_bar <- function(t, recruitment, model){

  a_0 <- recruitment$n_0 / (recruitment$n_0 + recruitment$n_1)
  a_1 <- 1 - a_0

  s_0 <- surv_pieces_simple(t, model$change_points, model$lambdas_0)
  s_1 <- surv_pieces_simple(t, model$change_points, model$lambdas_1)

  a_0 * s_0 + a_1 * s_1

}

########## w ##################

w <- function(t, t_star, recruitment, model){

  s_bar(pmin(t, t_star), recruitment, model) ^ (-1)

}


########## w ##################

w_fh <- function(t, rho, gamma, recruitment, model){

  s_bar(t, recruitment, model) ^ rho * (1 - s_bar(t, recruitment, model)) ^ gamma

}

######### ncp ############

#' Find the approximate non-centrality parameter and power of a modestWLRT.
#'
#' \code{ncp_power} returns the approximate non-centrality parameter and power of a modestWLRT for a range
#' of possible values of t*.
#' @param{t_star} A vector. A range of possible values for t*.
#' @param{rho} rho parameter in a Fleming-Harrington test. Default is NULL. Only used if F-H test used instead of MWLRT.
#' @param{gamma} gamma parameter in a Fleming-Harrington test. Default is NULL. Only used if F-H test used instead of MWLRT.
#' @param{model} A piecewise constant hazard model.
#'   A list containing the \code{change_points}; the rates \code{lambdas_0} on the control arm;
#'   and the rates \code{lambdas_1} on the treatment arm.
#' @param{recruitment} List of recruitment information.
#'   Containing \enumerate{
#'                 \item Sample size on control, \code{n_0}
#'                 \item Sample size on treatment, \code{n_1}
#'                 \item Recruitment period, \code{r_period}
#'                 \item Recruitment parameter for power model, \code{k}
#'               }
#' @param{dco} Time of data cut-off.
#' @param{length_t} Number of cutpoints to use when approximating the distribution of the MWLRT-statistic. Default is 18. Can be increased for greater accuracy.
#' @param{alpha_one_sided} One-sided alpha level.
#' @return A list containing
#' \enumerate{
#'                 \item non-centrality parameter corresponding to each t* \code{ncp}
#'                 \item power corresponding to each t* \code{power}
#'           }
#' @export

ncp_power <- function(t_star = NULL,
                      rho = NULL,
                      gamma = NULL,
                      model,
                      recruitment,
                      dco,
                      length_t,
                      alpha_one_sided = 0.025){



  if (!is.null(t_star)){
    return(ncp_power_t_star(t_star = t_star,
                            model = model,
                            recruitment = recruitment,
                            dco = dco,
                            length_t = length_t,
                            alpha_one_sided = alpha_one_sided))
  }
  else{
    return(ncp_power_fh(rho = rho,
                        gamma = gamma,
                        model = model,
                        recruitment = recruitment,
                        dco = dco,
                        length_t = length_t,
                        alpha_one_sided = alpha_one_sided))
  }


}



#########################################################################


ncp_power_t_star <- function(t_star,
                             model,
                             recruitment,
                             dco,
                             length_t,
                             alpha_one_sided = 0.025){

  R <- recruitment$n_1 / recruitment$n_0

  model_e <- model_extend(model, max_t = dco, length_t = length_t)

  events_info <- expected_events_two_arm(dco = dco,
                                         recruitment = recruitment,
                                         model = model_e,
                                         total_only = FALSE)

  prop_events <- events_info$prop_events
  total_events <- events_info$total_events

  mid_t <- c(0, model_e$change_points[-length(model_e$change_points)]) + diff(c(0, model_e$change_points)) / 2

  ncp <- numeric(length(t_star))
  power <- numeric(length(t_star))
  var_u <- numeric(length(t_star))
  e_u <- numeric(length(t_star))

  for (i in seq_along(ncp)){

    weights_t_star <- purrr::map_dbl(mid_t,
                                     w,
                                     t_star = t_star[i],
                                     recruitment = recruitment,
                                     model = model_e)


    weights_t_star <- c(weights_t_star, 0)

    log_hr <- log(model_e$lambdas_1 / model_e$lambdas_0)


    num <- sum(weights_t_star * log_hr * prop_events)
    denom <- sqrt(sum(weights_t_star ^ 2 * prop_events))

    ncp[i] <- num / denom * sqrt(total_events * R / (R + 1) ^ 2)
    power[i] <- pnorm(qnorm(alpha_one_sided),
                      mean = num / denom * sqrt(total_events * R / (R + 1) ^ 2))

    e_u[i] <- sum(weights_t_star * log_hr * prop_events * total_events * R / (R + 1) ^ 2)
    var_u[i] <- sum(weights_t_star ^ 2 * prop_events * total_events * R / (R + 1) ^ 2)

  }
  list(ncp = ncp,
       power = power,
       e_u = e_u,
       var_u = var_u,
       total_events = total_events)
}


#########################################################################


ncp_power_fh <- function(rho,
                         gamma,
                         model,
                         recruitment,
                         dco,
                         length_t,
                         alpha_one_sided = 0.025){

  R <- recruitment$n_1 / recruitment$n_0

  model_e <- model_extend(model, max_t = dco, length_t = length_t)

  events_info <- expected_events_two_arm(dco = dco,
                                         recruitment = recruitment,
                                         model = model_e,
                                         total_only = FALSE)

  prop_events <- events_info$prop_events
  total_events <- events_info$total_events

  mid_t <- c(0, model_e$change_points[-length(model_e$change_points)]) + diff(c(0, model_e$change_points)) / 2

  #### need to adjust the w function to allow F-H

  weights_fh <- purrr::map_dbl(mid_t,
                               w_fh,
                               rho = rho,
                               gamma = gamma,
                               recruitment = recruitment,
                               model = model_e)

  weights_fh <- c(weights_fh, 0)

  log_hr <- log(model_e$lambdas_1 / model_e$lambdas_0)


  num <- sum(weights_fh * log_hr * prop_events)
  denom <- sqrt(sum(weights_fh ^ 2 * prop_events))

  ncp <- num / denom * sqrt(total_events * R / (R + 1) ^ 2)
  power <- pnorm(qnorm(alpha_one_sided),
                    mean = num / denom * sqrt(total_events * R / (R + 1) ^ 2))

  e_u <- sum(weights_fh * log_hr * prop_events * total_events * R / (R + 1) ^ 2)
  var_u <- sum(weights_fh ^ 2 * prop_events * total_events * R / (R + 1) ^ 2)


  list(ncp = ncp,
       power = power,
       e_u = e_u,
       var_u = var_u,
       total_events = total_events)
}
