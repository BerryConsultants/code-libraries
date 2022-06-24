

installer <- function(package){
  # install a package if it isn't installed
  installed_packages <- package %in% rownames(installed.packages())
  if (any(installed_packages == FALSE)) {
    install.packages(packages[!installed_packages])
  }
}

# List of dependent packages
packages = c(
  'plotly',
  'readr',
  'dplyr',
  'tidyr',
  'stringr',
  'purrr',
  'ggplot2',
  'latex2exp',
  'RColorBrewer',
  'rstan'
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
options(warn = -1)


#' Rename the columns of a dataframe
#' @param df: data.frame whose column are to be renamed
#' @param names: names of the columns
#' @return data.frame with new column names
colnames_inplace = function(df, names){
  colnames(df) = names
  df
}


#' Inverse logistic function
#' @param x: argument
#' @return inverse logistic applied to x
inv_logit = function(x) {
  p = 1 / (1 + exp(-x))
  p[p < 0] = exp(x) / (1 + exp(-1))
  return(p)
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


