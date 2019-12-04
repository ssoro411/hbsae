data {
  int<lower=0>     m;        // Number of small areas
  int<lower=0>     p;        // Number of auxiliary variables
  real           y[m];       // Direct Estimate
  matrix[m,p]   X;           // Auxiliary variable
  real<lower=0> sDi[m];      // Sampling Error
  matrix[m,m]   Dw;           // Diagonal matrix with rowSums(W)
  matrix[m,m]   W;           // Spatial matrix
  matrix[m,m]   I;           // Identity matrix
}
parameters {
  real<lower=-1, upper=1> rho;  // Spatial parameter
  real<lower=0>  sigma_sq;       // Scale parameter
  vector[p]          beta;       // Regression parameter
  vector[m]         theta;       // Characteristic of areas
}
transformed parameters {
  vector[m] mu;
  mu = X*beta;
}
model {
  theta    ~ multi_normal_prec(mu, (1/sigma_sq)*( (Dw-rho*(W)) ) );
  for( i in 1:m)
  y[i]     ~ normal(theta[i], sDi[i]);
}
