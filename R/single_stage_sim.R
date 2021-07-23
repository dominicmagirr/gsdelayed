single_stage_sim_1 <- function(dummy = 1,
                               design,
                               model,
                               recruitment){


  if (is.null(model))  {model <- design$model}
  if (is.null(recruitment)) {recruitment <- design$recruitment}


  df_uncensored <- sim_t_uncensored(model, recruitment)
  df_final <- apply_dco(df_uncensored, events = ceiling(design$n_events[1]))


  if (!is.null(design$t_star)){
    t_star <- design$t_star


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

    wlrt_final <- wlrt(df_final,
                       trt_colname = "group",
                       time_colname = "time",
                       event_colname = "event",
                       wlr = "fh",
                       rho = rho,
                       gamma = gamma)
  }

  data.frame(c_1 = design$critical_values,
             z_1 = wlrt_final$z,
             t_1 = df_final$dco[1])

}

#' Single-stage simulation
#'
#' Simulate a single-stage trial, analyzed with a modestly-weighted log-rank test
#'
#' @param n_sims The number of simulations.
#' @param design The design object, created using \code{single_stage_design}
#' @param model Model assumptions. Default is NULL, in which case model assumptions in the design object are used.
#' @param recruitment Recruitment assumptions. Default is NULL, in which case recruitment assumptions in the design object are used.
#' @return A data-frame containing: critical value, z-statistic, timing of analysis.
#' @export

single_stage_sim <- function(n_sims = 1, design, model = NULL, recruitment = NULL){
  purrr::map_df(1:n_sims,
                single_stage_sim_1,
                design = design,
                model = model,
                recruitment = recruitment)
}
