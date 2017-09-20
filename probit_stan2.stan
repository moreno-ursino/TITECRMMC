data {
  real a; // intercept
  int<lower=0> N; // number of previous patients
  int<lower=0> K; // number of doses
  int<lower=1,upper=K> x[N]; // dose given to patients
  int<lower=1,upper=3> y[N]; // tox = 1,2,3
  real<lower=0,upper=1> w1[N]; // weights y=2
  real<lower=0,upper=1> w2[N]; // weights y=3
  real d[K]; // transformed doses for CRMMC
}
parameters {
  real beta;
  real gam;
}
transformed parameters {
  real theta[2];
  theta[1] <- (0 - a)/exp(beta);//  0.5^(1/exp(beta));
  theta[2] <- (exp(gam) -0.675 - a)/exp(beta); //     0.25^(1/(exp(beta)+exp(gam)));
}
model {
  real p1p; // probabilities first constraint
  real p2p; // probabilities second constraint
  real p1; // probabilities first constraint
  real p2; // probabilities second constraint
  real p3; // probabilities first constraint
  //real x1[N]; // logistic transformation T1
  real loglik;
  
  for (n in 1:N){
    p1p <- w1[n]*normal_cdf(a+exp(beta)*d[x[n]],0,1);
    //p1[n] <- if_else(p1[n]==0, 2^(-900), p1[n]);
    p2p <- w2[n]*normal_cdf(a+exp(beta)*d[x[n]]-exp(gam),0,1); 
    
    p1 <- 1-p1p;
    p1 <- if_else(p1==0, 2^(-900), p1);
    
    p2 <- p1p-p2p;
    p2 <- if_else(p2==0, 2^(-900), p2);
    
    p3 <-if_else(p2p==0, 2^(-900), p2p);
    
    loglik <- (y[n]==1)*log(p1) + (y[n]==2)*log(p2) + (y[n]==3)*log(p3);
    increment_log_prob(loglik);
  }
  
  beta ~ normal(0, sqrt(1.34));
  gam ~ normal(0, sqrt(1.34));
  //lambda ~ gamma(1, 0.1);
}