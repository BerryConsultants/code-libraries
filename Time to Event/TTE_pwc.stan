functions{
  // Likelihood function
  real pc_haz(vector lambda, vector beta, int[] delta, matrix H, matrix X){
    
    vector[num_elements(delta)] loglik_i;
    real loglik;
    vector[num_elements(delta)] Lambda; 
    vector[num_elements(delta)] eta; 
    
    eta = X * beta;
    Lambda = - (H * lambda) .* exp(eta);

    loglik_i = Lambda;
    for (i in 1:num_elements(delta)) {
      if (delta[i] > 0){
        loglik_i[i] = loglik_i[i] + log(lambda[delta[i]]) + eta[i];
      }
    }
    loglik = sum(loglik_i);
    
    return loglik;
  }
}

data{
  int<lower=1> N; // number of patients
  int<lower=1> D; // number of doses
  matrix[N,D] X; // design matrix (without column of ones)
  int<lower=1> J; // number of intervals
  
  int delta[N]; // vector indicating in which interval the event happens for each patient (0 if it does never happen)
  matrix[N,J] H; // matrix with the sufficient statistics to calculate the cumulative hazard
  
  vector<lower=0>[J] a_j; // Weight hyperparameters for each component of the gamma prior on lambda
  vector<lower=0>[J] b_j; // Mean hyperparameters for each component of the gamma prior on lambda
  real<lower=0> alpha; // Weight hyperparameter for the inv-gamma prior on tau^2
  real<lower=0> beta; // Central value hyperparameter for the inv-gamma prior on tau^2
}

parameters{
  vector<lower=0>[J] lambda; // piecewise constant baseline hazards
  vector[D] theta; // covariate effects
  real<lower=0> tau2;
}

transformed parameters {
  real<lower=0> tau = sqrt(tau2);
}

model{
  // Priors
  lambda ~ gamma(a_j, b_j); // prior for the piecewise constant baseline hazards
  tau2 ~ scaled_inv_chi_square(alpha, beta);
  theta[1] ~ normal(0, 1);
  for (d in 2:D){
    theta[d] ~ normal(theta[d-1], tau);
  }
  
  // Likelihood
  target += pc_haz(lambda, theta, delta, H, X);
}
