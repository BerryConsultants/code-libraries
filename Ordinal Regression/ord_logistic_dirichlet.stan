functions {
  real induced_dirichlet_lpdf(vector c, vector alpha) {
    int D = num_elements(c) + 1;
    vector[D - 1] sigma = inv_logit(- c);
    vector[D] p;
    matrix[D, D] J = rep_matrix(0, D, D);
    
    // Induced ordinal probabilities
    p[1] = 1 - sigma[1];
    for (k in 2:(D - 1)) {
      p[k] = sigma[k - 1] - sigma[k];
    }
    p[D] = sigma[D - 1];
    
    // Baseline column of Jacobian
    for (k in 1:D) {
      J[k, 1] = 1;
    }
    
    // Diagonal entries of Jacobian
    for (k in 2:D) {
      J[k, k] = - sigma[k - 1] * (1 - sigma[k - 1]);
      J[k - 1, k] = sigma[k - 1] * (1 - sigma[k - 1]);
    }
    
    return (dirichlet_lpdf(p | alpha) + log_determinant(J));
  }
}

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
  theta ~ normal(0, 2); // prior on the OR parameter 
  c ~ induced_dirichlet(rep_vector(1, D));
  
  // Observational model
  for (i in 1:N) {
    y[i] ~ ordered_logistic(0.0, c + X[i] * theta);
  }
}
