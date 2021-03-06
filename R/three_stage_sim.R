##############################################
three_stage_sim_1 <- function(dummy = 1,
                              design){

  model <- design$model
  recruitment <- design$recruitment
  alpha_spend_f <- design$alpha_spend_f
  alpha_one_sided <- design$alpha_one_sided


  df_uncensored <- sim_t_uncensored(model, recruitment)
  df_interim_1 <- apply_dco(df_uncensored, events = ceiling(design$n_events[1]))
  df_interim_2 <- apply_dco(df_uncensored, events = ceiling(design$n_events[2]))
  df_final <- apply_dco(df_uncensored, events = ceiling(design$n_events[3]))

  if (!is.null(design$t_star)){
    t_star <- design$t_star

    wlrt_interim_1 <- wlrt(df_interim_1,
                           trt_colname = "group",
                           time_colname = "time",
                           event_colname = "event",
                           wlr = "mw",
                           t_star = t_star)


    wlrt_interim_2 <- wlrt(df_interim_2,
                           trt_colname = "group",
                           time_colname = "time",
                           event_colname = "event",
                           wlr = "mw",
                           t_star = t_star)

    wlrt_final <- wlrt(df_final,
                       trt_colname = "group",
                       time_colname = "time",
                       event_colname = "event",
                       wlr = "mw",
                       t_star = t_star)
  }
  else {
    rho <- design$rho
    gamma <- design$gamma

    wlrt_interim_1 <- wlrt(df_interim_1,
                           trt_colname = "group",
                           time_colname = "time",
                           event_colname = "event",
                           wlr = "fh",
                           rho = rho,
                           gamma = gamma)


    wlrt_interim_2 <- wlrt(df_interim_2,
                           trt_colname = "group",
                           time_colname = "time",
                           event_colname = "event",
                           wlr = "fh",
                           rho = rho,
                           gamma = gamma)

    wlrt_final <- wlrt(df_final,
                       trt_colname = "group",
                       time_colname = "time",
                       event_colname = "event",
                       wlr = "fh",
                       rho = rho,
                       gamma = gamma)

  }

  c_1 <- crit_1_of_3(var_u_int_1 = wlrt_interim_1$v_u,
                     design = design,
                     alpha_spend_f = alpha_spend_f,
                     alpha_one_sided = alpha_one_sided)

  c_2 <- crit_2_of_3(var_u_int_2 = wlrt_interim_2$v_u,
                     var_u_int_1 = wlrt_interim_1$v_u,
                     design,
                     alpha_spend_f = alpha_spend_f,
                     alpha_one_sided = alpha_one_sided)

  c_3 <- crit_3_of_3(var_u_final = wlrt_final$v_u,
                     var_u_int_2 = wlrt_interim_2$v_u,
                     var_u_int_1 = wlrt_interim_1$v_u,
                     design,
                     alpha_spend_f = alpha_spend_f,
                     alpha_one_sided = alpha_one_sided)


  data.frame(c_1 = c_1,
             c_2 = c_2,
             c_3 = c_3,
             z_1 = wlrt_interim_1$z,
             z_2 = wlrt_interim_2$z,
             z_3 = wlrt_final$z,
             t_1 = df_interim_1$dco[1],
             t_2 = df_interim_2$dco[1],
             t_3 = df_final$dco[1])

}

#' Three-stage simulation
#'
#' Simulate a three-stage trial, analyzed with a modestly-weighted log-rank test
#'
#' @param n_sims The number of simulations.
#' @param design The design object, created using \code{three_stage_design}
#' @return A data-frame containing: critical values, z-statistics, timing of analyses.
#' @export
#'
three_stage_sim <- function(n_sims = 1, design) purrr::map_df(1:n_sims, three_stage_sim_1, design = design)
