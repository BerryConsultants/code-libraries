
installer <- function(package){
  installed_packages <- package %in% rownames(installed.packages())
  if (any(installed_packages == FALSE)) {
    install.packages(packages[!installed_packages])
  }
}


# List of dependent packages
packages = c(
  'tibble',
  'dplyr',
  'tidyr',
  'readr',
  'stringr',
  'knitr',
  'ggplot2',
  'latex2exp',
  'RColorBrewer',
  'rstan',
  'bayesplot',
  'coda',
  'survival',
  'survminer',
  'reshape2',
  'parallel',
  'kableExtra'
)

# Install + source the packages
installer(packages)
invisible(lapply(packages, library, character.only = TRUE))


# Initialize global variables
if ('rstan' %in% packages) {
  rstan_options(auto_write = TRUE)
  options(mc.cores = parallel::detectCores()) # use multiple cores if available, speeds up simulations
}
theme_set(theme_bw(base_size = 14))
options(dplyr.summarise.inform = FALSE)


# Auxiliary ---------------------------------------------------------------


#' Extract the legend from a ggplot
#' @param plt: plot whose legend is to be extracted
#' @return legend object
ext_legend = function(plt){
  
  tmp = ggplot_gtable(ggplot_build(plt))
  leg = which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend = tmp$grobs[[leg]]
  
  return(legend)
}


#' Empty theme 
custom_theme = function() {
  
  theme_survminer() %+replace%
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      panel.grid.major = element_line(colour = "grey90"),
      panel.border = element_blank(),
      panel.background = element_blank()
    )
}


#' Rename the columns of a dataframe
#' @param df: data.frame whose column are to be renamed
#' @param names: names of the columns
#' @return data.frame with new column names
colnames_inplace = function(df, names){
  colnames(df) = names
  df
}


#' Function to calculate the integrated accrual profile
#' @param t: time
#' @param peak_rate: peak accrual rate (weeks)
#' @param ramp_up: ramp up complete at week
#' @param t0: date of last enrolled patient (in weeks, since start of trial)
#' @return data.frame with new column names
Lambda = function(t, peak_rate, ramp_up, t0 = 0) {
  if (ramp_up == 0){
    peak_rate * (t + t0)
  }
  else {
    0.5 * peak_rate/ramp_up * (t + t0)^2 * as.numeric((t + t0) < ramp_up) + 
      (peak_rate * ((t + t0) - ramp_up) + 0.5 * peak_rate * ramp_up) * as.numeric((t + t0) >= ramp_up)
  }
}


#' Calculate the hazard function for piecewise constant hazards
#' @param t: time grid where to evaluate the cumulative hazard function
#' @param s_j: cut-points of the J intervals where the hazards are constant (length J + 1)
#' @param lambda: hazard rates in each of the J intervals (if a matrix, the hazard rates are on each row)
#' @return hazard function for each possible row of lambda
ht = function(t, s_j, lambda){
  J = length(s_j) - 1
  if (any(t < 0)) { 
    stop("Cannot have negative event times!\n")
  }
  s_j[J+1] = +Inf
  if (is.vector(lambda)){
    J = length(lambda)
    out = array(NA, dim = length(t))
    for (j in 1:J){
      idx = which( (t >= s_j[j]) & (t < s_j[j+1]))
      out[idx] = lambda[j]
    }
  }
  else if (is.matrix(lambda)){
    J = ncol(lambda)
    out = array(NA, dim = c(nrow(lambda), length(t)))
    for (j in 1:J){
      idx = which( (t >= s_j[j]) & (t < s_j[j+1]))
      out[,idx] = lambda[,j]
    }
  }
  else {
    stop("The hazard rates lambda should be a vector or a matrix!\n")
  }
  
  return (out)
}


#' Calculate the cumulative hazard function for piecewise constant hazards
#' @param t: time grid where to evaluate the cumulative hazard function
#' @param s_j: cut-points of the J intervals where the hazards are constant (length J + 1)
#' @param lambda: hazard rates in each of the J intervals (if a matrix, the hazard rates are on each row)
#' @return cumulative hazard function for each possible row of lambda
Ht = function(t, s_j, lambda){
  J = length(s_j) - 1
  s_j[J+1] = +Inf
  if (is.vector(lambda)){
    J = length(lambda)
    out = array(0, dim = length(t))
  }
  else if (is.matrix(lambda)){
    J = ncol(lambda)
    out = array(0, dim = c(length(t), nrow(lambda)))
  }
  else {
    stop("The hazard rates lambda should be a vector or a matrix!\n")
  }
  
  EXP = t(outer(t, s_j[2:(J+1)], FUN = function(x, y){pmin(x, y)})) - s_j[1:J]
  EXP[which(EXP < 0)] = 0
  
  if (is.vector(lambda)){
    out = as.numeric(colSums(EXP * lambda))
  }
  else if (is.matrix(lambda)){
    out = t(lambda %*% EXP)
  }
  
  return (out)
}


#' Calculate the survival function for piecewise constant hazards
#' @param t: time grid where to evaluate the survival function
#' @param s_j: cut-points of the J intervals where the hazards are constant (length J + 1)
#' @param lambda: hazard rates in each of the J intervals (if a matrix, the hazard rates are on each row)
#' @return survival function for each possible row of lambda
St = function(t, s_j, lambda){
  out = exp(-Ht(t, s_j, lambda))
  return (out)
}



# Data generation ---------------------------------------------------------


#' Normalized integrated accrual profile
#' @param x: time
#' @param T_max: maximum follow up time
#' @param peak_rate: peak accrual rate (weeks)
#' @param ramp_up: ramp up complete at week
#' @param t0: date of last enrolled patient (in weeks, since start of trial)
#' @return data.frame with new column names
Ft = function(x, T_max, peak_rate, ramp_up, t0 = 0) {
  (Lambda(x, peak_rate, ramp_up, t0) - Lambda(0, peak_rate, ramp_up, t0)) / (Lambda(T_max, peak_rate, ramp_up, t0) - Lambda(0, peak_rate, ramp_up, t0))
}


#' Inverse of the normalized integrated accrual profile
#' @param Ft: Normalized integrated accrual profile
#' @param u: uniform draw (to use the inverse cdf method)
#' @param T_max: maximum follow up time
#' @param peak_rate: peak accrual rate (weeks)
#' @param ramp_up: ramp up complete at week
#' @param t0: date of last enrolled patient (in weeks, since start of trial)
#' @return data.frame with new column names
Ftinv = function(Ft, u, T_max, peak_rate, ramp_up, t0 = 0) {
  
  N = length(u)
  a = rep(0, N)
  b = rep(T_max, N)
  binf = rep(0, N)
  bsup = rep(0, N)
  
  for (j in 1:20) {
    idx_low = which(Ft((a + b) / 2, T_max, peak_rate, ramp_up, t0) <= u)
    binf[idx_low] = (a[idx_low] + b[idx_low])/2
    bsup[idx_low] = b[idx_low]
    
    idx_upp = which(Ft((a + b) / 2, T_max, peak_rate, ramp_up, t0) >= u)
    bsup[idx_upp] = (a[idx_upp] + b[idx_upp]) / 2
    binf[idx_upp] = a[idx_upp]
    
    a = binf
    b = bsup
  }
  
  return ( (a + b) / 2 )
}


#' Sample observations from a piecewise constant hazard model
#' @param n: sample size
#' @param X: design matrix (without the column of leading ones)
#' @param lambda: hazard rates in each of the J intervals 
#' @param theta: log HR for the different active doses (length D)
#' @param s_j: cut-points of the J intervals where the hazards are constant (length J + 1)
#' @param T_max: maximum follow up time after the last patient randomized
#' @param accr_profile: accrual profile
#' @param t0: date of last enrolled patient (in weeks, since start of trial)
#' @return (time, event, dates): event times (or censoring times, whichever is sooner), event indicators, and accrual dates
sample_pwc_data = function(n, X, lambda, theta, s_j, T_max, accr_profile, t0){
  
  # Sample times from start of accrual according to a Poisson Process
  n_max = rpois(1, Lambda(2*n/accr_profile$peak_rate, accr_profile$peak_rate, accr_profile$ramp_up, t0))
  if (n_max < n){
    stop('Cannot simulate n individuals. Consider increasing accrual!\n')
  }
  # Use the inverse cdf method to sample n_max accrual dates 
  dates_temp = Ftinv(Ft, runif(n_max), 2*n/accr_profile$peak_rate, accr_profile$peak_rate, accr_profile$ramp_up, t0)
  # Take only n accrual dates, sort them and shift them by the initial date
  dates = t0 + sort(dates_temp)[1:n]
  
  # Now we use the inverse cdf method to sample from the likelihood model
  u = runif(n)
  XtBeta = as.numeric(X %*% theta)
  time = event = dropout = array(NA, dim = n)
  J = length(lambda)
  s_j[J+1] = max(s_j[J+1], dates[length(dates)] - dates[1] + T_max + 1)
  for (i in 1:n){
    # Define the max censoring time as the gap time between ith enrollment and 
    # last enrollment plus the maximum follow up after last randomized (T_max)
    obs_window_i = dates[length(dates)] - dates[i] + T_max
    
    # Calculate the hazard rate for individual i
    lambda_i = lambda * exp(XtBeta[i])
    # Check which time interval the observation should fall into
    H_temp = c(0, cumsum(lambda_i * diff(s_j)))
    interval_idx = head(which(log(u[i]) >= - H_temp), 1) - 1
    
    if (length(interval_idx) == 0){ # censored at the observation window endpoint
      time[i] = obs_window_i
      event[i] = 0
    }
    else { # sampled from the model in that interval
      event[i] = 1
      time[i] = s_j[interval_idx] - 1/lambda_i[interval_idx] * (log(u[i]) + H_temp[interval_idx])
      if (time[i] > obs_window_i){ # it still could be censored in that interval
        time[i] = obs_window_i
        event[i] = 0
      }
    }
  }
  
  return (list('time' = time, 'event' = event, 'dates' = dates))
}


# Predictive probabilities ------------------------------------------------


#' Compute the summary statistics for running the piecewise constant hazards model
#' The methodology is described here: https://myweb.uiowa.edu/pbreheny/7210/f19/notes/10-24.pdf
#' @param y: event times or censoring times, whichever is sooner (vector of length N)
#' @param nu: event indicators (vector of length N)
#' @param s_j: cut-points of the J intervals where the hazards are constant (length J + 1)
#' @return (delta_i, Hij): interval number in which subject i failed, cumulative hazards
get_summary_stats = function(y, nu, s_j){
  
  N = length(y)
  J = length(s_j) - 1
  s_j[J+1] = max(y) + 1
  
  # Calculate the sufficient statistics for the STAN model to run
  delta_i = array(0, dim = N) # interval number in which subject i failed
  Hij = array(NA, dim = c(N, J))
  for (i in 1:N){
    if (nu[i] == 1){
      delta_i[i] = head(which(y[i] < s_j), 1) - 1
    }
    Hij[i,] = pmin(y[i], s_j[2:(J+1)]) - s_j[1:J]
    Hij[i,which(Hij[i,] < 0)] = 0
  }
  
  # Check if we have enough events in each interval
  if (!(all(table(factor(delta_i, levels = 0:J)) > 4))) {
    warning ("Some intervals have less than 5 observations. Consider collapsing into larger intervals.\n")
  }
  stopifnot(rowSums(Hij) == y) # Just an additional check
  
  return (list('delta_i' = delta_i, 'Hij' = Hij))
}


#' Followup the current patients for T_max weeks using a piecewise constant 
#' hazard model
#' @param subj_id: subject IDs for the currently enrolled patients
#' @param time: vector indicating exposure time for each patient
#' @param event: vector indicating if each patient has already had the event or not
#' @param dropout: vector indicating if each patient has already dropped out or not
#' @param lambda: hazard rates in each of the J intervals 
#' @param s_j: cut-points of the J intervals where the hazards are constant (length J + 1)
#' @param X: design matrix (without the column of leading ones)
#' @param beta: regression coefficients (length p)
#' @param T_max: maximum follow up time (censoring window)
#' @param cens_rate: weekly rate corresponding to the desired censoring rate 
#' @return imputed dataset
sample_followup_data_curr = function(subj_id, time, event, dropout, lambda, s_j, 
                                     X, beta, T_max, cens_rate){
  
  # We use the inverse cdf method to impute using the piecewise constant hazard
  # likelihood
  idx_not_to_impute = which(event == 1 | dropout == 1) # if already observed, nothing to impute
  idx_to_impute = which(event == 0 & dropout == 0) # if not observed, events could occur
  
  # Now we use the inverse cdf method to sample from the likelihood model
  n = length(time)
  J = length(lambda)
  
  s_j[J+1] = +Inf
  XtBeta = as.numeric(X %*% beta)
  
  imp_time = imp_event = array(NA, dim = n)
  imp_dropout = array(0, dim = n)
  imp_time[idx_not_to_impute] = time[idx_not_to_impute]
  imp_event[idx_not_to_impute] = event[idx_not_to_impute]
  imp_dropout[idx_not_to_impute] = dropout[idx_not_to_impute]
  u = runif(n)
  if (length(idx_to_impute) > 0){ # continue followup for T_max weeks
    for (i in idx_to_impute){
      lambda_i = lambda * exp(XtBeta[i])
      
      H_temp = c(0, cumsum(lambda_i * diff(pmax(0, s_j - time[i]))))
      interval_idx = head(which(log(u[i]) >= - H_temp), 1) - 1
      
      if (length(interval_idx) == 0){ # event is censored 
        imp_time[i] = time[i] + T_max
        imp_event[i] = 0
      }
      else { # event is sampled in one of the intervals defined by s_j
        imp_time[i] = time[i] - 1/lambda_i[interval_idx] * (log(u[i]) + H_temp[interval_idx])
        imp_event[i] = 1
        if (imp_time[i] > time[i] + T_max){ # it still could be censored in that interval
          imp_time[i] = time[i] + T_max
          imp_event[i] = 0
        }
      }
      
      # Additionally, a censoring time is sampled from an exponential distribution 
      # with a rate that results in the assumed dropout rate
      c_time = time[i] + rexp(1, cens_rate) # using the memoryless property of Exp
      
      # If this censoring time is less than the event time for the subject, then 
      # the patient will be censored at this time and will not count as an event.
      if (c_time < imp_time[i]) {
        imp_dropout[i] = 1
        imp_event[i] = 0
        imp_time[i] = c_time
      }
    }
  }
  
  imp_dat = tibble(ID = subj_id, 
                   dropout = imp_dropout, 
                   event = imp_event, 
                   time = imp_time)
  
  return (imp_dat)
}


#' Sample new observations from a piecewise constant hazard model
#' @param n: sample size to be imputed
#' @param lambda: hazard rates in each of the J intervals 
#' @param s_j: cut-points of the J intervals where the hazards are constant (length J + 1)
#' @param T_max: maximum follow up time (censoring window)
#' @param cens_rate: weekly rate corresponding to the desired censoring rate 
#' @return imputed dataset
sample_pwc_data_new = function(n, lambda, s_j, T_max, cens_rate){
  
  # Now we use the inverse cdf method to sample from the likelihood model
  u = runif(n)
  time = event = array(NA, dim = n)
  dropout = array(0, dim = n)
  J = length(lambda)
  
  for (i in 1:n){
    # Check which time interval the observation should fall into
    H_temp = c(0, cumsum(lambda * diff(s_j)))
    interval_idx = head(which(log(u[i]) >= - H_temp), 1) - 1
    
    if (length(interval_idx) == 0){ # censored at the observation window endpoint
      time[i] = T_max
      event[i] = 0
    }
    else { # sampled from the model in that interval
      event[i] = 1
      time[i] = s_j[interval_idx] - 1/lambda[interval_idx] * (log(u[i]) + H_temp[interval_idx])
      if (time[i] > T_max){ # it still could be censored in that interval
        time[i] = T_max
        event[i] = 0
      }
    }
    
    # Additionally, a censoring time is sampled from an exponential distribution 
    # with a rate that results in the assumed dropout rate at 1 year. 
    c_time = rexp(1, cens_rate) # censoring time SINCE TIME 0 (not post-blanking)
    
    # If this censoring time is less than the event time for the subject, then 
    # the patient will be censored at this time and will not count as an event.
    if (c_time < time[i]) {
      dropout[i] = 1
      event[i] = 0
      time[i] = c_time
    }
  }

  imp_dat = tibble(ID = paste0('imp_', 1:n), 
                   DPOUTFL = dropout, 
                   peefl = event, 
                   EXPTIME = time)
  
  return (imp_dat)
}



