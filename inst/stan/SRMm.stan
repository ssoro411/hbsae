data {
  int<lower=0>        m2;       // Number of sampled small areas
  int<lower=0>        m1;       // Number of non-sampled small areas
  int<lower=0>        p;        // Number of auxiliary variables
  real                y[m2];    // Direct Estimate
  real<lower=0>       sDi[m2];  // Sampling error
  matrix[m1+m2,p]     X;        // Auxiliary variable
  matrix[m2,m1+m2] M;
  matrix[m1+m2,m1+m2]   R;           // Spatial matrix
  matrix[m1+m2,m1+m2]   I;           // Identity matrix
}
parameters {
  real<lower=0.0, upper=0.99999> rho;   // Spatial parameter 0 <= rho <= 1
  real<lower=0>  sigma_sq;                   // Scale parameter
  //real             tau_sq;                   // Scale parameter
  vector[p]          beta;                   // Regression parameter
  vector[m1+m2]         theta;                   // Area characteristics
}
transformed parameters {
  vector[m2]       Mtheta;
  vector[m1+m2] mu;
  mu = X*beta;
  Mtheta = M*theta;
}
model {
  //target += tau_sq;
  theta    ~ multi_normal_prec(mu, (1/sigma_sq)*( rho*(R) + (1-rho)*I ));
  //theta    ~ multi_normal_prec(mu, (exp(-tau_sq))*( rho*(R) + (1-rho)*I ));
  for( i in 1:m2)
  y[i]     ~ normal( Mtheta[i], sDi[i]);
}
