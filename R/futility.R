#' Futility two-stage design
#'
#' Create a two-stage design by adding a futility boundary to an existing two-stage design
#'
#' @param design The design object, created using \code{two_stage_design}
#' @param hr_threshold The threshold for futility on the (average) HR scale
#' @return The design object, with futility boundary and updated power/expected duration.
#' @export



futility_two_stage <- function(design,
                               hr_threshold){


  design$futility_z_value <- log(hr_threshold) * sqrt(design$var_u[1])
  design$p_early_stop_futility <- c(under_alternative = 1 - pnorm(design$futility_z_value, mean = design$ncp_z[1]),
                                    under_null_approx = 1 - pnorm(design$futility_z_value))


  info_frac <- design$var_u[1] / design$var_u[2]
  p2 <- mvtnorm::pmvnorm(lower = c(design$critical_values[1], -Inf),
                         upper = c(design$futility_z_value, design$critical_values[2]),
                         mean = design$ncp_z,
                         sigma = matrix(c(1, sqrt(info_frac),
                                          sqrt(info_frac), 1), nrow = 2))[1]

  design$overall_power <- design$p_early_stop[1] + p2

  p_either_stop <- design$p_early_stop + design$p_early_stop_futility
  design$expected_t <- design$dco[1] * p_either_stop + design$dco[2] * (1 - p_either_stop)

  design

}


#' Futility three-stage design
#'
#' Create a three-stage design by adding a futility boundary to an existing three-stage design
#'
#' @param design The design object, created using \code{three_stage_design}
#' @param hr_threshold_1 The threshold for futility at the first interim, on the (average) HR scale
#' @param hr_threshold_d The threshold for futility at the second interim, on the (average) HR scale
#' @return The design object, with futility boundary and updated power/expected duration.
#' @export



futility_three_stage <- function(design,
                                 hr_threshold_1 = 1.2,
                                 hr_threshold_2 = 1.1){


  design$futility_z_value <- log(c(hr_threshold_1, hr_threshold_2)) * sqrt(design$var_u[1:2])


  ########### 1st interim futility ##################
  design$p_first_int_futility <- c(under_alternative = 1 - pnorm(design$futility_z_value[1], mean = design$ncp_z[1]),
                                   under_null_approx = 1 - pnorm(design$futility_z_value[1]))

  ########### 2nd interim efficacy ##################
  info_frac_1 <- design$var_u[1] / design$var_u[3]
  p2 <- mvtnorm::pmvnorm(lower = c(design$critical_values[1], -Inf),
                         upper = c(design$futility_z_value[1], design$critical_values[2]),
                         mean = design$ncp_z[1:2],
                         sigma = matrix(c(1, sqrt(info_frac_1),
                                          sqrt(info_frac_1), 1), nrow = 2))[1]

  p2_0 <- mvtnorm::pmvnorm(lower = c(design$critical_values[1], -Inf),
                           upper = c(design$futility_z_value[1], design$critical_values[2]),
                           mean = c(0,0),
                           sigma = matrix(c(1, sqrt(info_frac_1),
                                            sqrt(info_frac_1), 1), nrow = 2))[1]

  design$p_second_int_stop <- c(under_alternative = p2,
                                under_null_approx = p2_0)


  ########### 2nd interim futility ##################

  p2_f <- mvtnorm::pmvnorm(lower = c(design$critical_values[1], design$futility_z_value[2]),
                         upper = c(design$futility_z_value[1], Inf),
                         mean = design$ncp_z[1:2],
                         sigma = matrix(c(1, sqrt(info_frac_1),
                                          sqrt(info_frac_1), 1), nrow = 2))[1]

  p2_0_f <- mvtnorm::pmvnorm(lower = c(design$critical_values[1], design$futility_z_value[2]),
                           upper = c(design$futility_z_value[1], Inf),
                           mean = c(0,0),
                           sigma = matrix(c(1, sqrt(info_frac_1),
                                            sqrt(info_frac_1), 1), nrow = 2))[1]

  design$p_second_int_futility <- c(under_alternative = p2_f,
                                    under_null_approx = p2_0_f)

  ########### Final power ############################
  info_frac_2 <- design$var_u[2] / design$var_u[3]
  p3 <- mvtnorm::pmvnorm(lower = c(design$critical_values[1:2], -Inf),
                         upper = c(design$futility_z_value, design$critical_values[3]),
                         mean = design$ncp_z,
                         sigma = matrix(c(1, sqrt(info_frac_1 / info_frac_2), sqrt(info_frac_1),
                                          sqrt(info_frac_1 / info_frac_2), 1, sqrt(info_frac_2),
                                          sqrt(info_frac_1), sqrt(info_frac_2), 1), nrow = 3))[1]

  design$overall_power <- design$p_first_int_stop[1] + p2 + p3

  ########## Expected t #############################


  p_either_stop_1 <- design$p_first_int_stop + design$p_first_int_futility
  p_either_stop_2 <- design$p_second_int_stop + design$p_second_int_futility
  design$expected_t <- design$dco[1] * p_either_stop_1 +
                       design$dco[2] * p_either_stop_2 +
                       design$dco[3] * (1 - p_either_stop_1 - p_either_stop_2)

  design

}


# #######################################################################
# #######################################################################
# t_seq <- seq(0,36, length.out = 100)
# s_c <- purrr::map_dbl(t_seq, surv_pieces_simple, change_points = c(6), lambdas = log(2) / c(9, 9))
# s_e <- purrr::map_dbl(t_seq, surv_pieces_simple, change_points = c(6), lambdas = log(2) / c(9, 16))
#
# plot(t_seq, s_c, ylim = c(0,1), type = "l",
#      xlab = "time (months)", ylab = "survival",
#      main = "Alternative hypothesis in an immuno-oncology RCT")
# points(t_seq, s_e, lty = 2, type = "l")
# legend("topright", c("experimental", "control"), lty = 2:1)
# #######################################################################
# model <- list(change_points = c(6),
#               lambdas_0 = c(log(2) / 9, log(2) / 9),
#               lambdas_1 = c(log(2) / 9, log(2) / 16))
#
# recruitment <- list(n_0 = 225, n_1 = 225, r_period = 12, k = 1)
# #######################################################################
# d2 <- two_stage_design(t_star = 0,
#                        model = model,
#                        recruitment = recruitment,
#                        dco_int = 12,
#                        dco_final = 30,
#                        alpha_spend_f = ldobf,
#                        alpha_one_sided = 0.025)
#
# #####################################################################
# d2f <- futility_two_stage(d2, 1.2)
#
#
# d2
# d2f
#
# #####################################################################
# ### next stage is to write futility three stage...
#
# d3 <- three_stage_design(t_star = 0,
#                          model = model,
#                          recruitment = recruitment,
#                          dco_int_1 = 18,
#                          dco_int_2 = 24,
#                          dco_final = 30,
#                          alpha_spend_f = ldobf,
#                          alpha_one_sided = 0.025)
#
# d3f <- futility_three_stage(d3, 1, 0.95)
# d3f
# d3
