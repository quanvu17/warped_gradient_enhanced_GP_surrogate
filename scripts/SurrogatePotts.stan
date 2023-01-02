functions {
  vector ft(vector t, real tC, real e0, real ecrit, real v0, real vmaxLo, real vmaxHi, real phi1, real phi2) {
    vector[num_elements(t)] mu;
    real sqrtBcritPhi = sqrt(tC)*phi1;
    for (i in 1:num_elements(t)) {
      if (t[i] <= tC) {
        real sqrtBdiffPhi = sqrt(tC - t[i])*phi1;
        mu[i] = e0 + t[i]*v0 - ((2*(vmaxLo-v0))/(phi1^2))*((sqrtBcritPhi + 1)/exp(sqrtBcritPhi) - (sqrtBdiffPhi + 1)/exp(sqrtBdiffPhi));
      } else {
        real sqrtBdiff = sqrt(t[i] - tC);
        mu[i] = ecrit - ((2*vmaxHi)/phi2)*(sqrtBdiff/exp(phi2*sqrtBdiff) + (exp(-phi2*sqrtBdiff) - 1)/phi2);
      }
    }
    return mu;
  }

  vector dfdt(vector t, real tC, real v0, real vmaxLo, real vmaxHi, real r1, real r2) {
    vector[num_elements(t)] dmu;
    for (i in 1:num_elements(t)) {
      if (t[i] <= tC) {
        dmu[i] = v0 + (vmaxLo-v0)*exp(-r1*sqrt(tC - t[i]));
      } else {
        dmu[i] = vmaxHi*exp(-r2*sqrt(t[i] - tC));
      }
    }
    return dmu;
  }
}
data {
  int<lower = 1> M;
  int<lower = 1> N;
  real<lower = 1> maxY;
  real<lower = 1> Vlim;
  real<lower = 0> e0;
  real<lower = 0> v0;
  real tcrit;
  matrix<lower=0, upper=maxY>[M,N] y;
  vector[M] t;
}
parameters {
  real<lower = 0> a;
  real<lower = 0> b;
  real<lower = e0, upper=maxY> ecrit;
  real<lower = 0, upper=Vlim> vmaxLo;
  real<lower = 0, upper=Vlim> vmaxHi;
}
transformed parameters {
  vector[M] curr_mu;
  vector[M] curr_var;
  curr_mu = ft(t, tcrit, e0, ecrit, v0, vmaxLo, vmaxHi, a, b);
  curr_var = dfdt(t, tcrit, v0, vmaxLo, vmaxHi, a, b);
}
model {
  for (i in 1:M) {
    y[i,] ~ normal(curr_mu[i], sqrt(curr_var[i]));
  }
}
