


##compute tau_IND and sigma_IND

B = 10000
propensity = 0.4
tau_IND_cases = NULL
sigma_IND_cases = NULL
for (case in 1:10){
  
  withinp = 0.8
  betweenp = 0.2
  
  
  if (case == 1){
    r = 3
    Psi_fun = function(x,y){
      same_group = (x < y/r) & (x > (y-1)/r)
      return(same_group)
    }
    lambdas = c(1,1,1)
    outcome_fun = function(W, theta, prop_treated_ngb){
      Y = 1/2*(W + theta*prop_treated_ngb)^2
      return(Y)
    }
    outcome_deriv = function(W, theta, prop_treated_ngb){
      Y = (W + theta*prop_treated_ngb)*theta
      return(Y)
    }
  } else if (case == 2){
    r = 3
    Psi_fun = function(x,y){
      return(x^y)
    }
    lambdas = c(1,-2,1)*27/4
    outcome_fun = function(W, theta, prop_treated_ngb){
      Y = cos(3*W*prop_treated_ngb)
      return(Y)
    }
    outcome_deriv = function(W, theta, prop_treated_ngb){
      Y = -3*W*sin(3*W*prop_treated_ngb)
      return(Y)
    }
  }else if (case == 3){
    r = 3
    Psi_fun = function(x,y){
      return(x^y)
    }
    lambdas = c(1,-2,1)*27/4
    outcome_fun = function(W, theta, prop_treated_ngb){
      Y = -exp(theta)*cos(3*W*prop_treated_ngb)
      return(Y)
    }
    outcome_deriv = function(W, theta, prop_treated_ngb){
      Y = exp(theta)*3*W*sin(3*W*prop_treated_ngb)
      return(Y)
    }
  }else if (case == 4){
    r = 3
    Psi_fun = function(x,y){
      return(x >= (y-1)/r)
    }
    lambdas = c(1,1,1)/3
    outcome_fun = function(W, theta, prop_treated_ngb){
      Y = (1 + W)*exp(prop_treated_ngb)
      return(Y)
    }
    outcome_deriv = function(W, theta, prop_treated_ngb){
      Y = (1 + W)*exp(prop_treated_ngb)
      return(Y)
    }
  } else if (case == 5){
    r = 3
    Psi_fun = function(x,y){
      return(x >= (y-1)/r)
    }
    lambdas = c(1,1,1)/3
    outcome_fun = function(W, theta, prop_treated_ngb){
      Y = 1/5*(1+theta)^2*(1 + W)*exp(prop_treated_ngb)
      return(Y)
    }
    outcome_deriv = function(W, theta, prop_treated_ngb){
      Y = 1/5*(1+theta)^2*(1 + W)*exp(prop_treated_ngb)
      return(Y)
    }
  } else if (case == 6){
    r = 1
    Psi_fun = function(x,y){
      same_group = 0.3 + 0.6*(x>0.5)
      return(same_group)
    }
    lambdas = 1
    outcome_fun = function(W, theta, prop_treated_ngb){
      Y = 1/2*(W + theta*prop_treated_ngb)^2
      return(Y)
    }
    outcome_deriv = function(W, theta, prop_treated_ngb){
      Y = (W + theta*prop_treated_ngb)*theta
      return(Y)
    }
  } else if (case == 7){
    r = 1
    Psi_fun = function(x,y){
      return(sin(2*pi*x)*0.3 + 0.5)
    }
    lambdas = 1
    outcome_fun = function(W, theta, prop_treated_ngb){
      Y = cos(3*W*prop_treated_ngb)
      return(Y)
    }
    outcome_deriv = function(W, theta, prop_treated_ngb){
      Y = -3*W*sin(3*W*prop_treated_ngb)
      return(Y)
    }
  }else if (case == 8){
    r = 1
    Psi_fun = function(x,y){
      return(sin(2*pi*x)*0.3 + 0.5)
    }
    lambdas = 1
    outcome_fun = function(W, theta, prop_treated_ngb){
      Y = -exp(theta)*cos(3*W*prop_treated_ngb)
      return(Y)
    }
    outcome_deriv = function(W, theta, prop_treated_ngb){
      Y = exp(theta)*3*W*sin(3*W*prop_treated_ngb)
      return(Y)
    }
  }else if (case == 9){
    r = 1
    Psi_fun = function(x,y){
      return((x+1)^4/20 + 0.1)
    }
    lambdas = 1
    outcome_fun = function(W, theta, prop_treated_ngb){
      Y = (1 + W)*exp(prop_treated_ngb)
      return(Y)
    }
    outcome_deriv = function(W, theta, prop_treated_ngb){
      Y = (1 + W)*exp(prop_treated_ngb)
      return(Y)
    }
  } else if (case == 10){
    r = 1
    Psi_fun = function(x,y){
      return((x+1)^4/20 + 0.1)
    }
    lambdas = 1
    outcome_fun = function(W, theta, prop_treated_ngb){
      Y = 1/5*(1+theta)^2*(1 + W)*exp(prop_treated_ngb)
      return(Y)
    }
    outcome_deriv = function(W, theta, prop_treated_ngb){
      Y = 1/5*(1+theta)^2*(1 + W)*exp(prop_treated_ngb)
      return(Y)
    }
  }

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
    Psi_u1u2 = outer(c(u1,u2), 1:r, Psi_fun)
    G_u1u2 = (Psi_u1u2 %*% diag(lambdas) %*% t(Psi_u1u2))[1,2]
    if (case == 1){
      G_u1u2 = (Psi_u1u2 %*% diag(lambdas) %*% t(Psi_u1u2))[1,2]*(withinp-betweenp) + betweenp
    }
    sigma_INDs = c(sigma_INDs,sqrt(c_u1u2*G_u1u2))
  }
  tau_IND = mean(tau_INDs)
  sigma_IND = sqrt(mean(sigma_INDs^2))
  tau_IND_cases = c(tau_IND_cases,tau_IND)
  sigma_IND_cases = c(sigma_IND_cases,sigma_IND)
}

##########load the tau_tildes
tau_tilde_vars = NULL
tau_tilde_means = NULL
cases = NULL
ns = NULL
rhons = NULL

for (seed in 1:200){
  filename = sprintf("output/output_var_vs_n%i.Rdata",seed)
  if (file.exists(filename)) {
    load(filename)
    if(case >= 6){
      tau_tilde_vars = c(tau_tilde_vars,var(tau_tilde1s))
      tau_tilde_means = c(tau_tilde_means,mean(tau_tilde1s))
    } else{
      tau_tilde_vars = c(tau_tilde_vars,var(tau_tildes))
      tau_tilde_means = c(tau_tilde_means,mean(tau_tildes))
    }
    cases = c(cases, case)
    ns = c(ns, n)
    rhons = c(rhons, rhon)
  }
}

dat1 = data.frame(tau_tilde_vars = tau_tilde_vars, tau_tilde_means = tau_tilde_means,
                 cases = cases, ns = ns, rhons = rhons,
                 true_means = tau_IND_cases[cases],
                 true_vars = rhons*(sigma_IND_cases[cases])^2
                 )
library(dplyr)

for (rate_of_converge in 1:2){
  dat2 = data.frame(vars = tau_tilde_vars, means = tau_tilde_means,
                    Setting = as.factor(cases), ns = ns, rhons = rhons,
                    biases = tau_tilde_means - tau_IND_cases[cases],
                    Rank = as.factor(3-2*(cases>=6)))
  dat3 = dat2%>% mutate(logmse = log(vars + biases^2)/log(10),
                         logn = log(ns)/log(10),
                        rates = -5*log(rhons)/log(ns)
                        ) %>%
    filter(abs(rates - rate_of_converge) < 0.5)
  
  
  #library(jcolors)
  #library(gridExtra)
  library(latex2exp)
  library(ggplot2)
  #library(directlabels)
  library(tidyverse)
  p = ggplot(data = dat3, aes(x = logn, y = logmse,
                              group = interaction(Setting,Rank),
                              color = Setting,
                              linetype = Rank
                              ))
  p = p + theme_bw(base_size = 15)
  p = p + theme(line = element_blank(), legend.key.width=unit(1.5,"cm"))
  p = p + geom_line(size = 1) #+ scale_color_brewer(palette="Set1")
  my_cols = c("black",  "darkgray", "burlywood3",  "purple2",
    "blue2", "green2",  "yellow2", "orange1", "red3", "maroon2")
  p = p + scale_color_manual(values = my_cols)
  p = p + scale_linetype_manual(values=c("solid", "longdash"))
  #p = p + scale_color_jcolors(palette = "pal12")
  #title = sprintf('PC balancing Estimator: $\\log_{10}( MSE(\\tilde{\\tau}))$ when $\\rho_n = n^{-%i/5}$', rate_of_converge)
  #p = p + ggtitle(TeX(title))
  p = p + xlab(TeX('$\\log_{10}(n)$'))
  p = p + ylab(TeX('$\\log_{10}( MSE $)$'))
  subtitle = sprintf('Slope of the reference lines  $= %i/5$', -rate_of_converge)
  p = p + labs(caption = TeX(subtitle))
  for (i in (-10:10)/5){
    p = p + geom_abline(intercept = i, slope = -rate_of_converge/5, color = "grey", linetype = 3)
  }
  #p = p + geom_dl(aes(label = Setting), method = list(dl.trans(x = x + .2), "last.points"))
  print(p)
}

