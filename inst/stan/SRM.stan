data {
  int<lower=0>     m;        // Number of small areas
  int<lower=0>     p;        // Number of auxiliary variables
  real           y[m];       // Direct Estimate
  matrix[m,p]   X;           // Auxiliary variable
  real<lower=0> sDi[m];      // Sampling Error
  matrix[m,m]   R;           // Spatial matrix
  matrix[m,m]   I;           // Identity matrix
}
parameters {
  real<lower=0.0, upper=0.99999> rho;   // Spatial parameter 0 <= rho <= 1
  real<lower=0>  sigma_sq;                   // Scale parameter
  vector[p]          beta;                   // Regression parameter
  vector[m]         theta;                   // Area characteristics
}
transformed parameters {
  vector[m] mu;
  mu = X*beta;
}
model {
  //target += tau_sq;
  theta    ~ multi_normal_prec(mu, (1/sigma_sq)*( rho*(R) + (1-rho)*I ));
  for( i in 1:m)
  y[i]     ~ normal(theta[i], sDi[i]);
}
