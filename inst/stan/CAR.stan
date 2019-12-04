data {
  int<lower=0>     m;        // Number of small areas
  int<lower=0>     p;        // Number of auxiliary variables
  real           y[m];       // Direct Estimate
  real<lower=0> sDi[m];      // Sampling Error
  matrix[m,m]   W;           // Spatial matrix
  matrix[m,p]   X;           // Auxiliary variable
  matrix[m,m]   I;           // Identity matrix
  real uu;                   // Upper bound of rho
  real ll;                   // Lower bound of rho
}
parameters {
  real<lower=ll, upper=uu> rho;  // Spatial parameter
  real<lower=0>  sigma_sq;       // Scale parameter
  vector[p]          beta;       // Regression parameter
  vector[m]         theta;       // Characteristic of areas
}
transformed parameters {
  vector[m] mu;
  mu = X*beta;
}
model {
  theta    ~ multi_normal_prec(mu, (1/sigma_sq)*( (I-rho*(W)) ) );
  for( i in 1:m)
  y[i]     ~ normal(theta[i], sDi[i]);
}
