#' Time-at-risk intervals.
#'
#' \code{get_t_pieces} returns the time at risk between intervals \code{ts} up to time \code{t}.
#' @param{t} Total time.
#' @param{ts} The breakpoints for intervals, starting at 0, ending at Inf.
#' @return A vector of length one less than \code{ts}, containing time at risk between consecutive breakpoints in \code{ts}.

get_t_pieces = function(t, ts){

  if (length(t) > 1) stop("t must be length 1")

  pmax(0, pmin(t - ts[-length(ts)], diff(ts)))
}

#' Probability of event-free at time t, under piecewise constant hazard.
#'
#' \code{surv_pieces} returns the probability of being event-free at time \code{t} assuming piecewise constant hazard.
#' @param{t} time
#' @param{change_points} The breakpoints for the piecewise hazard function.
#' @param{lambdas} A vector of length one more than \code{change_points}. The hazard between consecutive breakpoints in \code{c(0, change_points, Inf)}.
#' @return The probability of being free of the event-of-interest at time \code{t}.
#' @export

surv_pieces_simple = function(t, change_points, lambdas){

  if (length(t) > 1) stop("t must be length 1")
  if (length(change_points) != length(lambdas) - 1) stop("require one event rate per time period")
  ts = c(0, change_points, Inf)

  exp(-sum(get_t_pieces(t, ts) * lambdas))

}

####################################################

########################################
## conditional probability of an event
## occuring before patient time 'min_f' and
## before calendar time 'r_period + f_period',
## given that a patient is recruited at 'r_time'.
## Note: 'r_time' is a vector here

cond_p_event = function(r_time,
                        min_f,
                        r_period,
                        f_period,
                        change_points,
                        lambdas){


  f_time = pmin(min_f, r_period + f_period - r_time)

  event_prob = numeric(length(f_time))

  for (i in seq_along(event_prob)){

    event_prob[i] = 1 - surv_pieces_simple(f_time[i], change_points, lambdas)

  }

  event_prob

}

########################################
## probability density of recruitment time

p_r = function(r_time, r_period, k) k * (r_time / r_period) ^ (k - 1) / r_period

########################################
## cond_p_event times p_r

p_event_r = function(r_time,
                     min_f,
                     r_period,
                     f_period,
                     change_points,
                     lambdas,
                     k){

  cond_p_event(r_time,
               min_f,
               r_period,
               f_period,
               change_points,
               lambdas) * p_r(r_time, r_period, k)



}
########################################
## probability of an event occuring
## prior to patient-time c(change_points, Inf) AND
## before calendar time 'r_period + f_period'

p_event = function(r_period,
                   f_period,
                   change_points,
                   lambdas,
                   k){

  min_fs = c(change_points, Inf)

  integrals = numeric(length(min_fs))

  for (i in seq_along(integrals)){

    integrals[i] = integrate(p_event_r,
                             lower = 0,
                             upper = r_period,
                             min_f = min_fs[i],
                             r_period = r_period,
                             f_period = f_period,
                             change_points = change_points,
                             lambdas = lambdas,
                             k = k)$value
  }

  integrals
}


##################################################################

#' Expected number of events
#'
#' \code{expected_events_two_arm} returns the expected number of events at a given data-cut-off time.
#' @param{dco} Time of data cut-off.
#' @param{recruitment} List of recruitment information.
#'   Containing \enumerate{
#'                 \item Sample size on control, \code{n_0}
#'                 \item Sample size on experimental, \code{n_1}
#'                 \item Recruitment period, \code{r_period}
#'                 \item Recruitment parameter for power model, \code{k}
#'               }
#' @param{model} The piecewise hazard model.
#'   A list containing the \code{change_points} and \code{lambdas}.
#' @param{total_only} By default, the function returns only the expected number of events.
#'   To see more detailed output set \code{total_only = FALSE}.
#' @return The expected number of events at \code{dco}.
#' @export

expected_events_two_arm = function(dco = 28,
                                   recruitment,
                                   model,
                                   total_only = TRUE){

  n_0 = recruitment$n_0
  n_1 = recruitment$n_1
  r_period = recruitment$r_period
  k = recruitment$k

  change_points = model$change_points
  lambdas_0 = model$lambdas_0
  lambdas_1 = model$lambdas_1

  if (length(dco) > 1) stop("dco must be length 1")
  if (length(n_0) > 1) stop("n_0 must be length 1")
  if (length(n_1) > 1) stop("n_1 must be length 1")

  # assume dco is greater than r_period.
  # what if dco < r_period

  if (dco < r_period){

    n_0 = n_0 * (dco / r_period) ^ k
    n_1 = n_1 * (dco / r_period) ^ k
    r_period = dco

  }

  f_period = dco - r_period

  expected_events_0 = n_0 * p_event(r_period,
                                    f_period,
                                    change_points,
                                    lambdas_0,
                                    k)

  expected_events_1 = n_1 * p_event(r_period,
                                    f_period,
                                    change_points,
                                    lambdas_1,
                                    k)


  total_events_0 = expected_events_0[length(expected_events_0)]
  total_events_1 = expected_events_1[length(expected_events_1)]
  total_events = total_events_0 + total_events_1

  if (total_only) return(c(total_events_0 = total_events_0,
                           total_events_1 = total_events_1,
                           total_events = total_events))


  events_per_period = cbind(events_0 = diff(c(0, expected_events_0)),
                            events_1 = diff(c(0, expected_events_1)))

  prop_events = rowSums(events_per_period) / total_events

  hrs = lambdas_1 / lambdas_0

  av_hr = exp(sum(log(hrs) * prop_events))

  list(total_events_0 = total_events_0,
       total_events_1 = total_events_1,
       total_events = total_events,
       events_per_period = cbind(c(change_points, Inf), events_per_period),
       prop_events = prop_events,
       hrs = hrs,
       av_hr = av_hr,
       n_recruited_0 = n_0,
       n_recruited_1 = n_1,
       n_recruited = n_0 + n_1)

}


########################################################################

#' Expected time of dco
#'
#' \code{expected_dco_two_arm} returns the expected calendar time of data-cut-off for a given number of events.
#' @param{total_events} Number of events.
#' @param{recruitment} List of recruitment information.
#'   Containing \enumerate{
#'                 \item Sample size on control, \code{n_0}
#'                 \item Sample size on experimental, \code{n_1}
#'                 \item Recruitment period, \code{r_period}
#'                 \item Recruitment parameter for power model, \code{k}
#'               }
#' @param{model} The piecewise hazard model.
#'   A list containing the \code{change_points} and \code{lambdas}.
#' @return The expected calendar time of reaching \code{total_events}.
#' @export

expected_dco_two_arm = function(total_events,
                                recruitment,
                                model){

  if (length(total_events) > 1) stop("total_events must be length 1")

  find_dco = function(x){


    expected_events_two_arm(dco = x,
                            recruitment,
                            model,
                            total_only = TRUE)["total_events"] - total_events

  }


  uniroot(find_dco, c(0.001, 1000))$root

}
