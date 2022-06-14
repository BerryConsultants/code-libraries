

installer <- function(package){
  # install a package if it isn't installed
  installed_packages <- package %in% rownames(installed.packages())
  if (any(installed_packages == FALSE)) {
    install.packages(packages[!installed_packages])
  }
}

# List of dependent packages
packages = c(
  'MASS',
  'ggalluvial',
  'pammtools',
  'dplyr',
  'tidyr',
  'stringr',
  'purrr',
  'ggplot2',
  'latex2exp',
  'RColorBrewer',
  'rstan',
  'bayesplot',
  'coda',
  'egg',
  'parallel', 
  'forcats', 
  'MCMCpack', 
  'knitr', 
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


#' Inverse logistic function
#' @param x: argument
#' @return inverse logistic applied to x
inv_logit = function(x) {
  p = 1 / (1 + exp(-x))
  p[p < 0] = exp(x) / (1 + exp(-1))
  return(p)
}


#' Logistic function
#' @param p: argument
#' @return logistic applied to p
logit = function(p){
  log(p / (1-p))
}


#' Probability density function of the logistic distribution
#' @param x: argument
#' @return pdf of the logistic applied to x
logistic_pdf = function(x) {
  exp(-x) / (1 + exp(-x))^2
}


#' Rename the columns of a dataframe
#' @param df: data.frame whose column are to be renamed
#' @param names: names of the columns
#' @return data.frame with new column names
colnames_inplace = function(df, names){
  colnames(df) = names
  df
}


#' Calculate the class proportions in the treatment group for ordinal logistic regression
#' @param p_ctr: class proportions in the control group
#' @param OR: common odds ratio
#' @return class proportions in the treatment group
OR_transform = function(p_ctr, OR){
  
  gamma = logit(head(cumsum(p_ctr), -1)) # transform to latent parameters (cutpoints)
  theta = log(OR) # common Odds-ratio
  p_trt = round(diff(c(0, inv_logit(theta + gamma), 1)), 2)
  
  return (p_trt)
}


# Data generation ---------------------------------------------------------


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


#' @param N: sample size
#' @param TRTPN: treatment arm indicator
#' @param p_ctr_true: probabilities under the control
#' @param p_trt_true: probabilities under the treatment
#' @param P_tr: transition probability matrix from Visit t to Visit (t-1)
#' @param T_max: follow up time after the last patient randomized
#' @param accr_profile: accrual profile
#' @param t0: date of last enrolled patient (in weeks, since start of trial)
#' @return
sample_long_ordinal_data = function(N, TRTPN, p_ctr_true, p_trt_true, P_tr, T_max, accr_profile, t0){

  K = length(p_ctr_true)
  # Sample times from start of accrual according to a Poisson Process
  Ntimes = 5
  n_max = rpois(1, Lambda(Ntimes*N/accr_profile$peak_rate, accr_profile$peak_rate, accr_profile$ramp_up, t0))
  if (n_max < N){
    stop('Cannot simulate n individuals with such small accrual!\n')
  }
  # Use the inverse cdf method to sample n_max accrual dates
  dates_temp = Ftinv(Ft, runif(n_max), Ntimes*N/accr_profile$peak_rate, accr_profile$peak_rate, accr_profile$ramp_up, t0)
  # Take only n accrual dates, sort them and shift them by the initial date
  dates = t0 + sort(dates_temp)[1:N]

  OTC30FL = RESP30 = OTC90FL = RESP90 = OTC180FL = RESP180 = array(NA, dim = N)
  for (i in 1:N){
    obs_window_i = dates[length(dates)] - dates[i] + T_max

    if (obs_window_i < (180 + 21) / 7) { # if not past 180 days
      OTC180FL[i] = 0
      if (obs_window_i < (90 + 14) / 7) { # if not past 90 days
        OTC90FL[i] = 0
        if (obs_window_i < (30 + 14) / 7) { # if not past 30 days
          OTC30FL[i] = 0
        }
        else { # between 30 and 90 days
          OTC30FL[i] = 1
          if (TRTPN[i] == 'Control') {
            RESP30[i] = sample(0:(K - 1), 1, FALSE, p_ctr_true)
          }
          else if (TRTPN[i] == 'Treatment') {
            RESP30[i] = sample(0:(K - 1), 1, FALSE, p_trt_true)
          }
        }
      }
      else { # between 90 and 180 days
        OTC30FL[i] = OTC90FL[i] = 1
        if (TRTPN[i] == 'Control') {
          RESP90[i] = sample(0:(K - 1), 1, FALSE, p_ctr_true)
        }
        else if (TRTPN[i] == 'Treatment') {
          RESP90[i] = sample(0:(K - 1), 1, FALSE, p_trt_true)
        }
        RESP30[i] = sample(0:(K - 1), 1, FALSE, P_tr[RESP90[i] + 1,])
      }
    }
    else {
      OTC30FL[i] = OTC90FL[i] = OTC180FL[i] = 1
      if (TRTPN[i] == 'Control') {
        RESP180[i] = sample(0:(K - 1), 1, FALSE, p_ctr_true)
      }
      else if (TRTPN[i] == 'Treatment') {
        RESP180[i] = sample(0:(K - 1), 1, FALSE, p_trt_true)
      }
      RESP90[i] = sample(0:(K - 1), 1, FALSE, P_tr[RESP180[i] + 1,])
      RESP30[i] = sample(0:(K - 1), 1, FALSE, P_tr[RESP90[i] + 1,])
    }
  }

  data = cbind(dates, RESP30, RESP90, RESP180, OTC30FL, OTC90FL, OTC180FL) %>%
    as_tibble %>% 
    mutate(SITEID = sample(paste0('0', 101:110), N, TRUE), 
           TRTPN = factor(key_trt[dose_label + 1]),
           RANDDT = as.Date('2018-03-15') + 7 * dates,
           CUTDT = max(RANDDT), 
           .before = dates) %>% 
    group_by(SITEID) %>% 
    mutate(SUBJID = paste0(SITEID, '-', str_pad(row_number(), 3, pad = "0")), .before = SITEID) %>% 
    ungroup() %>% 
    select(-dates)

  return (data)
}


#' Function to plot the transition probability matrix
#' @param P_trans: list with three elements: p30 (vector), p90 (matrix), p180 (matrix)
#' @param digits: number of decimal digits to be used in the plot
#' @param title: title above the heatmaps
#' @param palette: color palette
#' @return plot object
plot_Ptrans = function(P_trans, digits = 1, title = NULL, palette = 'Oranges'){
  
  K = length(P_trans$p30)
  n_digits = if_else(digits == 0, "%s", paste0("%0.", digits, "f"))
  
  # Initial state probabilities
  p1 = P_trans$p30 %>%
    as_tibble() %>% 
    mutate(rowname = as.character(0:(K - 1)), 
           colname = as.character(0), .before = 1) %>% 
    mutate(p_transitions = value / sum(value)) %>% 
    mutate(label = sprintf(n_digits, value)) %>%
    ggplot(aes(colname, forcats::fct_rev(rowname), fill = p_transitions)) +
    geom_tile() +
    geom_text(aes(label = label)) +
    coord_fixed() +
    labs(x = '', y = 'Response at 30 days') + 
    scale_fill_distiller(palette = palette, direction = 1) +
    theme(legend.position = 'none',
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank()
    )
  
  # Transitions from 30 to 90
  p2 = P_trans$p90 %>%
    as_tibble() %>% 
    colnames_inplace(as.character(0:(K - 1))) %>% 
    mutate(rowname = as.character(0:(K - 1)), .before = 1) %>% 
    pivot_longer(
      cols = -rowname,
      names_to = "colname",
      values_to = "transitions"
    ) %>% 
    group_by(rowname) %>% 
    mutate(p_transitions = transitions / sum(transitions), 
           p_transitions = if_else(is.na(p_transitions), 0, p_transitions)) %>% 
    ungroup() %>% 
    mutate(label = sprintf(n_digits, transitions)) %>%
    ggplot(aes(colname, forcats::fct_rev(rowname), fill = p_transitions)) +
    geom_tile() +
    geom_text(aes(label = label)) +
    coord_fixed() +
    labs(x = 'Response at 90 days', y = 'Response at 30 days') + 
    scale_fill_distiller(palette = palette, direction = 1) +
    theme(legend.position = 'none')
  
  # Transitions from 90 to 180
  p3 = P_trans$p180 %>%
    as_tibble() %>% 
    colnames_inplace(as.character(0:(K - 1))) %>% 
    mutate(rowname = as.character(0:(K - 1)), .before = 1) %>% 
    pivot_longer(
      cols = -rowname,
      names_to = "colname",
      values_to = "transitions"
    ) %>% 
    group_by(rowname) %>% 
    mutate(p_transitions = transitions / sum(transitions), 
           p_transitions = if_else(is.na(p_transitions), 0, p_transitions)) %>% 
    ungroup() %>% 
    mutate(label = sprintf(n_digits, transitions)) %>%
    ggplot(aes(colname, forcats::fct_rev(rowname), fill = p_transitions)) +
    geom_tile() +
    geom_text(aes(label = label)) +
    coord_fixed() +
    labs(x = 'Response at 180 days', y = 'Response at 90 days') + 
    scale_fill_distiller(palette = palette, direction = 1) +
    theme(legend.position = 'none')
  
  plt = egg::ggarrange(p1, p2, p3, 
                       ncol = 3, 
                       widths = c(1, 6, 6), 
                       top = title
  )
  
  return (plt)
}


# Model fit ---------------------------------------------------------------


#' Creates a dataframe with the posterior parameters for the Dirichlet 
#' distribution used in the imputation model, for the two possible arms.
#' @param df: dataframe of study data
#' @param long_prior: prior for the Dirichlet distributions of the longitudinal model
#' @param key_trt: treatment arms labels
#' @return list of posterior hyperparameters for the Dirichlet distributions 
get_imputation_pars = function(df, long_prior, key_trt){
  
  # Number of outcome classes
  K = length(long_prior$beta)
  
  # Count how many transitions we have at 30, 90 and 180 days
  PP30 = df %>% 
    filter(OTC30FL == 1) %>% 
    mutate(RESP30 = factor(RESP30, levels = 0:(K-1))) %>% 
    count(TRTPN, RESP30, .drop = F) %>% 
    filter(!is.na(RESP30))
  PP90 = df %>% 
    filter(OTC90FL == 1) %>% 
    mutate(RESP30 = factor(RESP30, levels = 0:(K-1)), 
           RESP90 = factor(RESP90, levels = 0:(K-1))) %>% 
    count(TRTPN, RESP30, RESP90, .drop = F) %>% 
    filter(!is.na(RESP30), !is.na(RESP90))
  PP180 = df %>% 
    filter(OTC180FL == 1) %>% 
    mutate(RESP90 = factor(RESP90, levels = 0:(K-1)),
           RESP180 = factor(RESP180, levels = 0:(K-1))) %>%
    count(TRTPN, RESP90, RESP180, .drop = F) %>% 
    filter(!is.na(RESP90), !is.na(RESP180))
  
  # Calculate the updated parameters for the Dirichlet distributions
  p30_ctr = PP30 %>% filter(TRTPN == key_trt[1]) %>% pull(n) + long_prior$beta
  p30_trt = PP30 %>% filter(TRTPN == key_trt[2]) %>% pull(n) + long_prior$beta
  
  p90_ctr = PP90 %>% filter(TRTPN == key_trt[1]) %>% pull(n) %>% matrix(K, K, byrow = T) + long_prior$alpha
  p90_trt = PP90 %>% filter(TRTPN == key_trt[2]) %>% pull(n) %>% matrix(K, K, byrow = T) + long_prior$alpha
  
  p180_ctr = PP180 %>% filter(TRTPN == key_trt[1]) %>% pull(n) %>% matrix(K, K, byrow = T) + long_prior$alpha
  p180_trt = PP180 %>% filter(TRTPN == key_trt[2]) %>% pull(n) %>% matrix(K, K, byrow = T) + long_prior$alpha
  
  return(list('p30_ctr' = p30_ctr, 
              'p30_trt' = p30_trt, 
              'p90_ctr' = p90_ctr, 
              'p90_trt' = p90_trt, 
              'p180_ctr' = p180_ctr, 
              'p180_trt' = p180_trt))
}


#' Main model function. At each iteration it imputes a dataset and performs the 
#' relevant tests that declare success/futility for the clinical trial. 
#' @param df: dataframe of study data
#' @param long_prior: prior for the Dirichlet distributions of the longitudinal model
#' @param n_max: maximum sample size
#' @param key_trt: treatment arms labels
#' @param iters: number of iterations
#' @param seed: random seed
#' @return list with elements:
#'         'PPn', a vector of flags for success/failure of the trials across simulations based on the rule PPn
#'         'PPmax', a vector of flags for success/failure of the trials across simulations based on the rule PPmax
#'         'OR_n', a vector with the observed ORs across simulations based on current sample size
#'         'OR_max', a vector with the observed ORs across simulations based on maximum sample size
#'         'pvals_n', a vector with the p-values across simulations based on maximum sample size
#'         'pvals_max', a vector with the p-values across simulations based on maximum sample size
fit_model = function(df, long_prior, n_max, key_trt, iters = 1000, seed = 12345){
  
  # Prepare data structures
  K = length(long_prior$beta)
  PPn = array(NA, dim = iters)
  PPmax = array(NA, dim = iters)
  OR_n = array(NA, dim = iters)
  OR_max = array(NA, dim = iters)
  pvals_n = array(NA, dim = iters)
  pvals_max = array(NA, dim = iters)
  
  n_act = df %>%
    pull(SUBJID) %>%
    n_distinct()
  ctr_act = df %>%
    filter(TRTPN == key_trt[1]) %>%
    nrow
  trt_act = df %>%
    filter(TRTPN == key_trt[2]) %>%
    nrow
  
  ctr_max = floor(n_max/3)
  trt_max = n_max - ctr_max
  n_ctr = ctr_max - ctr_act
  n_trt = trt_max - trt_act
  
  if ( (n_ctr <= 0) | (n_trt <= 0) ) {
    stop("Number of patients in one of the arms is larger than maximum sample size\n")
  }
  
  # Get posterior parameters for the imputation model
  imp_pars = df %>%
    get_imputation_pars(long_prior, key_trt)
  
  
  set.seed(seed)
  
  # Sample all of the parameter vectors before the simulations loop
  pi30_ctr = MCMCpack::rdirichlet(iters, imp_pars$p30_ctr)
  pi30_trt = MCMCpack::rdirichlet(iters, imp_pars$p30_trt)
  
  pi90_ctr = pi90_trt = pi180_ctr = pi180_trt = array(NA, dim = c(iters, K, K))
  for (d in 1:K) {
    pi90_ctr[,d,] = MCMCpack::rdirichlet(iters, imp_pars$p90_ctr[d,])
    pi90_trt[,d,] = MCMCpack::rdirichlet(iters, imp_pars$p90_trt[d,])
    
    pi180_ctr[,d,] = MCMCpack::rdirichlet(iters, imp_pars$p180_ctr[d,])
    pi180_trt[,d,] = MCMCpack::rdirichlet(iters, imp_pars$p180_trt[d,])
  }
  
  # Calculate the relevant quantities for the current data under treatment and control
  df_ctr = df %>%
    filter(TRTPN == key_trt[1]) %>%
    select(SUBJID, TRTPN, RESP30, RESP90, RESP180)
  df_trt = df %>%
    filter(TRTPN == key_trt[2]) %>%
    select(SUBJID, TRTPN, RESP30, RESP90, RESP180)
  
  # Some things to time the implementation
  pb = txtProgressBar(min = 1, max = iters)
  start = proc.time()
  
  set.seed(seed)
  my_seeds = sample(1:100000, iters, FALSE)
  for (i in 1:iters){ # loop over simulations
    
    set.seed(my_seeds[i])
    
    ### (1) Control group missing data imputation ###
    
    pi30_temp = pi30_ctr[i,]
    pi90_temp = pi90_ctr[i,,]
    pi180_temp = pi180_ctr[i,,]
    
    # Imputes the current dataset for the treatment group
    df_imp_ctr = df_ctr %>%
      add_row(
        tibble(
          SUBJID = paste0('imp_ctr_', 1:n_ctr),
          TRTPN = key_trt[1]
        )
      ) %>%
      group_by(SUBJID) %>%
      mutate(RESP30 = if_else(is.na(RESP30), as.numeric(sample(0:(K-1), 1, TRUE, pi30_temp)), RESP30),
             RESP90 = if_else(is.na(RESP90), as.numeric(sample(0:(K-1), 1, TRUE, pi90_temp[RESP30 + 1,])), RESP90),
             RESP180 = if_else(is.na(RESP180), as.numeric(sample(0:(K-1), 1, TRUE, pi180_temp[RESP90 + 1,])), RESP180)) %>% 
      ungroup()
    
    
    ### (2) Treatment group missing data imputation ###
    
    pi30_temp = pi30_trt[i,]
    pi90_temp = pi90_trt[i,,]
    pi180_temp = pi180_trt[i,,]
    
    # Imputes the current dataset for the treatment group
    df_imp_trt = df_trt %>%
      add_row(
        tibble(
          SUBJID = paste0('imp_trt_', 1:n_trt),
          TRTPN = key_trt[2]
        )
      ) %>%
      group_by(SUBJID) %>%
      mutate(RESP30 = if_else(is.na(RESP30), as.numeric(sample(0:(K-1), 1, TRUE, pi30_temp)), RESP30),
             RESP90 = if_else(is.na(RESP90), as.numeric(sample(0:(K-1), 1, TRUE, pi90_temp[RESP30 + 1,])), RESP90),
             RESP180 = if_else(is.na(RESP180), as.numeric(sample(0:(K-1), 1, TRUE, pi180_temp[RESP90 + 1,])), RESP180)) %>% 
      ungroup()
    
    df_imp = rbind(df_imp_trt, df_imp_ctr)
    
    
    ### (3) Perform the ordinal regression test ###
    
    # Current sample size
    data_fit = df_imp %>% 
      filter(!str_detect(SUBJID, 'imp_')) %>% 
      mutate(RESP180 = factor(RESP180))

    fit_freq = MASS::polr(RESP180 ~ TRTPN, data = data_fit, Hess = TRUE)
    fit_summ = summary(fit_freq)
    
    obs_OR = exp(- fit_summ$coefficients[1,"Value"])
    OR_n[i] = obs_OR
    t_val = - fit_summ$coefficients[1,"Value"] / fit_summ$coefficients[1,"Std. Error"]
    dof = fit_freq$df.residual
    
    pval_n = 1 - pt(q = t_val, df = dof)
    pvals_n[i] = pval_n
    
    
    # Maximum sample size
    data_fit = df_imp %>% 
      mutate(RESP180 = factor(RESP180))
    
    fit_freq = MASS::polr(RESP180 ~ TRTPN, data = data_fit, Hess = TRUE)
    fit_summ = summary(fit_freq)
    
    obs_OR = exp(- fit_summ$coefficients[1,"Value"])
    OR_max[i] = obs_OR
    t_val = - fit_summ$coefficients[1,"Value"] / fit_summ$coefficients[1,"Std. Error"]
    dof = fit_freq$df.residual
    
    pval_max = 1 - pt(q = t_val, df = dof)
    pvals_max[i] = pval_max
    
    # Calculate rules to declare success/futility at the final
    PPn[i] = if_else(pval_n < 0.02, 1, 0)
    PPmax[i] = if_else(pval_max < 0.02, 1, 0)
    
    setTxtProgressBar(pb, i)
  }
  cat(sprintf('\n Imputation complete in %.1f minutes', (proc.time() - start)[3]/ 60))
  
  return(list(
    PPn = PPn,
    PPmax = PPmax, 
    OR_n = OR_n, 
    OR_max = OR_max, 
    pvals_n = pvals_n, 
    pvals_max = pvals_max
  ))
}



