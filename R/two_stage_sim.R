two_stage_sim_1 <- function(dummy = 1,
                            design,
                            model,
                            recruitment){


  if (is.null(model))  {model <- design$model}
  if (is.null(recruitment)) {recruitment <- design$recruitment}


  alpha_spend_f <- design$alpha_spend_f
  alpha_one_sided <- design$alpha_one_sided


  df_uncensored <- sim_t_uncensored(model, recruitment)
  df_interim <- apply_dco(df_uncensored, events = ceiling(design$n_events[1]))
  df_final <- apply_dco(df_uncensored, events = ceiling(design$n_events[2]))

  if (!is.null(design$t_star)){
    t_star <- design$t_star

    wlrt_interim <- wlrt(df_interim,
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

    wlrt_interim <- wlrt(df_interim,
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


  ##########################
  ## allow for under-running
  ##########################


  if (wlrt_interim_1$v_u / design$var_u[2] > 0.99){

    c_1 <- qnorm(alpha_one_sided)

    c_2 <- -Inf

  }
  else {


    c_1 <- crit_1_of_2(var_u_int = wlrt_interim$v_u,
                       design = design,
                       alpha_spend_f = alpha_spend_f,
                       alpha_one_sided = alpha_one_sided)

    c_2 <- crit_2_of_2(var_u_final = wlrt_final$v_u,
                       var_u_int = wlrt_interim$v_u,
                       design = design,
                       alpha_spend_f = alpha_spend_f,
                       alpha_one_sided = alpha_one_sided)

  }

  data.frame(c_1 = c_1,
    c_2 = c_2,
    z_1 = wlrt_interim$z,
    z_2 = wlrt_final$z,
    t_1 = df_interim$dco[1],
    t_2 = df_final$dco[1])

}


#' Two-stage simulation
#'
#' Simulate a two-stage trial, analyzed with a modestly-weighted log-rank test
#'
#' @param n_sims The number of simulations.
#' @param design The design object, created using \code{two_stage_design}
#' @param model Model assumptions. Default is NULL, in which case model assumptions in the design object are used.
#' @param recruitment Recruitment assumptions. Default is NULL, in which case recruitment assumptions in the design object are used.
#' @return A data-frame containing: critical values, z-statistics, timing of analyses.
#' @export

two_stage_sim <- function(n_sims = 1, design, model = NULL, recruitment = NULL){
  purrr::map_df(1:n_sims,
                two_stage_sim_1,
                design = design,
                model = model,
                recruitment = recruitment)
}

