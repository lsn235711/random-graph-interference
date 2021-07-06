args = commandArgs(TRUE)
seed = as.numeric(args[1])
set.seed(seed)

rate = (seed-1)%/%100 + 1

seed1 = seed %% 100 + 1
case = seed1%%10 + 1
mult = (seed1-1)%/%10

library(RSpectra)
library(Matrix)
ptm <- proc.time()
n = floor(1000*(1.3^mult))
N = 500
rhon = n^(-rate/5)


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



tau_hats = NULL
tau_tildes = NULL
tau_tilde1s = NULL
for (i in 1:N){
  U = runif(n)
  propensity = 0.4
  ## rank r case
  Psi = outer(U, 1:r, Psi_fun)
  if(case == 1){
    Gij = (withinp-betweenp)*Psi %*% diag(lambdas) %*% t(Psi) + betweenp
  }else{
    Gij = Psi %*% diag(lambdas) %*% t(Psi)
  }
  Eij = (matrix(runif(n^2), ncol = n) < rhon*Gij) + 0
  ##symmetrize adjacency matrix
  Eij.diag = diag(Eij)
  Eij[lower.tri(Eij,diag=T)] = 0
  Eij = Eij + t(Eij) + diag(rep(0,n))
  
  theta = U
  W = (runif(n) < propensity) + 0
  Ni = apply(Eij,2,sum)
  Mi = apply(Eij*W,2,sum)
  prop_treated_ngb = Mi/pmax(Ni,1)
  Y = outcome_fun(W, theta, prop_treated_ngb) + 0.2*rnorm(n)
  
  ## unbiased estimator
  tau_hat = mean(Y*(Mi/propensity - (Ni - Mi)/(1-propensity)))
  tau_hats = c(tau_hats, tau_hat)
  ## PC balancing estimator
  res = eigs_sym(Eij,r)
  #res = eigs_sym(Gij,r)
  Psi_hat = res$vectors
  gamma_original = Mi/propensity - (Ni - Mi)/(1-propensity)
  beta = -t(Psi_hat) %*% gamma_original
  tau_tilde = mean(Y*(Mi/propensity - (Ni - Mi)/(1-propensity) + Psi_hat %*% beta))
  tau_tildes = c(tau_tildes, tau_tilde)
  
  ##assume rank1
  beta1 = -sum(Ni*(Mi - propensity*Ni))/(propensity*(1-propensity)*sum(Ni^2))
  tau_tilde1 = mean(Y*(Mi/propensity - (Ni - Mi)/(1-propensity) + beta1*Ni))
  tau_tilde1s = c(tau_tilde1s, tau_tilde1)
  rm(list=setdiff(ls(), c("tau_tildes","tau_tilde1s","tau_hats","n","N","r","rhon","propensity","withinp","betweenp","outcome_fun","Psi_fun","lambdas","case","ptm","seed")))
  gc(verbose = TRUE)
}


var(tau_tildes)
var(tau_tilde1s)
var(tau_hats)

save(tau_tildes, tau_tilde1s,tau_hats,case,n, rhon, file = paste("output/output_var_vs_n",seed,".Rdata", sep = ""))
