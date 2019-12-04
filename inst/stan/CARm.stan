data {
  int<lower=0>        m2;       // Number of sampled small areas
  int<lower=0>        m1;       // Number of non-sampled small areas
  int<lower=0>        p;        // Number of auxiliary variables
  real                y[m2];    // Direct Estimate
  real<lower=0>       sDi[m2];  // Sampling error
  matrix[m1+m2,p]     X;        // Auxiliary variable
  matrix[m2,m1+m2] M;
  matrix[m1+m2,m1+m2]   W;           // Spatial matrix
  matrix[m1+m2,m1+m2]   I;           // Identity matrix
  real uu;                   // Upper bound of rho
  real ll;                   // Lower bound of rho
}
parameters {
  real<lower=ll, upper=uu> rho;  // Spatial parameter
  real<lower=0>      sigma_sq;   // Scale parameter
  //real                 tau_sq;   // Scale parameter
  vector[p]              beta;   // Regression parameter
  vector[m1+m2]         theta;   // Characteristic of areas
}
transformed parameters {
  vector[m2]       Mtheta;
  vector[m1+m2] mu;
  mu = X*beta;
  Mtheta = M*theta;
}
model {
  //target += tau_sq;
  theta    ~ multi_normal_prec(mu, (1/sigma_sq)*( (I-rho*(W)) ) );
  //theta    ~ multi_normal_prec(mu, (exp(-tau_sq))*( (I-rho*(W)) ) );
  for( i in 1:m2)
  y[i]     ~ normal( Mtheta[i], sDi[i]);
}
