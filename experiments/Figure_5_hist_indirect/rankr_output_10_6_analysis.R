tau_tildesall = NULL
cases = NULL

seed = 1

for (seed in 1:40){
    filename = sprintf("output10_6/output_rankr_106_25_%i.Rdata",seed)
    if (file.exists(filename)) {
      load(filename)
      tau_tildesall = c(tau_tildesall,tau_tildes)
    }
}
tau_tildes = tau_tildesall

n= 1000000
B = 200000 ##number of simulations to estimate sigma_IND
rhon = n^(-2/5)
propensity= 0.4
r = 3


withinp = 0.8
betweenp = 0.2
##stochatic block model
outcome_fun = function(W, theta, prop_treated_ngb){
  Y = 1/2*(W + theta*prop_treated_ngb)^2
  return(Y)
}
Psi_fun = function(x,y){
  same_group = (x < y/3) & (x > (y-1)/3)
  return(same_group)
}
lambdas = c(withinp - betweenp, withinp - betweenp, withinp - betweenp, betweenp)

outcome_deriv = function(W, theta, prop_treated_ngb){
  Y = (W + theta*prop_treated_ngb)*theta
  return(Y)
}


##compute tau_IND and sigma_IND
tau_INDs = NULL
sigma_INDs = NULL
for (l in 1:B){
  u1 = runif(1)
  u2 = runif(1)
  base_deriv = propensity*outcome_deriv(1, u1, propensity) + (1-propensity)*outcome_deriv(0, u1, propensity)
  tau_INDs = c(tau_INDs,base_deriv)
  ##normal errors get canceled, assume Y(0) and Y(1) have the same error
  c_u1u2 = 1/4*(outcome_fun(1, u1, propensity)
                -outcome_fun(0, u1, propensity)
                +outcome_fun(1, u2, propensity)
                -outcome_fun(0, u2, propensity)
  )^2
  Psi_u1u2 = cbind(outer(c(u1,u2), 1:r, Psi_fun),c(1,1))
  G_u1u2 = (Psi_u1u2 %*% diag(lambdas) %*% t(Psi_u1u2))[1,2]
  sigma_INDs = c(sigma_INDs,sqrt(c_u1u2*G_u1u2))
}
tau_IND = mean(tau_INDs)
sigma_IND = sqrt(mean(sigma_INDs^2))
tau_IND
sigma_IND
#hist((tau_tildes - tau_IND)/(sqrt(rhon)*sigma_IND))

library(latex2exp)
library(ggplot2)

dat_his = data.frame(tau_tildes = tau_tildes)
range_abs = max(-min(tau_tildes-tau_IND), max(tau_tildes-tau_IND))
x = seq(tau_IND-range_abs, tau_IND+range_abs, length.out=100)
dat_dens = data.frame(x = x, den = dnorm(x, tau_IND, sqrt(rhon)*sigma_IND))
dfgg = with(dat_his,dat_dens)


colors <- c("CLT predictions" = "red")
fills <- c("Estimators" = "grey")
p = ggplot(dfgg)
p = p + geom_histogram(data = dat_his, aes(x = tau_tildes, y = ..density..,fill = "Estimators"), bins = 20, color = "black")
p = p + xlim(tau_IND-range_abs, tau_IND+range_abs)
p = p + geom_line(data = dat_dens, aes(x = x, y = den, color = "CLT predictions"))
p = p + scale_color_manual(name="Normal Approximations (CLT)", values = colors)
p = p + scale_fill_manual(name="PC balancing Estimator", values = fills)
p = p + ggtitle(TeX('Histogram of $\\tilde{\\tau}_{IND}$')) + xlab(TeX('$\\tilde{\\tau}_{IND}$'))
print(p)

