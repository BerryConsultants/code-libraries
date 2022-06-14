data {
  int<lower=1> N;             // Number of observations
  int<lower=1> D;             // Max of ordinal categories
  
  int<lower=0, upper=1> X[N];
  int<lower=1, upper=D> y[N]; // Observed ordinals {0, ..., K}
}

parameters {
  real theta;   // Latent effect
  ordered[D - 1] c; // (Internal) cut points
}

model {
  // Prior model
  // Implicit uniform priors
  
  // Observational model
  for (i in 1:N) {
    y[i] ~ ordered_logistic(0.0, c + X[i] * theta);
  }
}
