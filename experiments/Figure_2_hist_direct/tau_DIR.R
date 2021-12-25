##Direct effect

set.seed(1)

N1 = 3000
propensity = 0.7
g = 0.4
n = 1000
tauhatDIRs = NULL
tauhatHajs = NULL

sigma = 0.3

for (ii in 1:N1){
  u = runif(n,0,1)
  theta = u
  W = rbinom(n,1,propensity)
  Gu1u2 = matrix(rep(g,n^2),ncol = n)
  diag(Gu1u2) = rep(0,n)
  edge = matrix(runif(n^2), ncol = n) < Gu1u2
  edge.diag = diag(edge)
  edge[lower.tri(edge,diag=T)] = 0
  edge = edge + t(edge)
  
  Ni = apply(edge,2,sum)
  Mi = apply(edge*W,2,sum)
  
  prop = Mi/pmax(1,Ni)
  y = W*prop/propensity^2 + sigma*rnorm(n)
  tauhatDIR = mean(y*W/propensity - y*(1-W)/(1-propensity))
  tauhatHaj = mean(y*W/mean(W) - y*(1-W)/(1-mean(W)))
  tauhatDIRs = c(tauhatDIRs, tauhatDIR)
  tauhatHajs = c(tauhatHajs, tauhatHaj)
}

sigma02 = propensity*(1-propensity)*(1/(propensity^4) + (sigma*(1/propensity + 1/(1-propensity)))^2)
sigma12 = propensity*(1-propensity)*(4/(propensity^4) + (sigma*(1/propensity + 1/(1-propensity)))^2)

sigmaHa02 = propensity*(1-propensity)*(sigma*(1/propensity + 1/(1-propensity)))^2
sigmaHa12 = propensity*(1-propensity)*(1/(propensity^4) + (sigma*(1/propensity + 1/(1-propensity)))^2)

tau_true = 1/propensity

hist(tauhatDIRs)
sqrt(n*mean((tauhatDIRs - tau_true)^2))
sqrt(sigma02)
sqrt(sigma12)

dat_his = data.frame(tauhat = tauhatDIRs)
range1 = max( tau_true - min(tauhatDIRs), max(tauhatDIRs) - tau_true)
x = seq(tau_true - range1, tau_true + range1, length.out=100)
dat_dens1 = data.frame(x = x, den = dnorm(x, tau_true, sqrt(sigma12/n)))
dat_dens0 = data.frame(x = x, den = dnorm(x, tau_true, sqrt(sigma02/n)))

dat_his2 = data.frame(tauhat = tauhatHajs)
range2 = max( tau_true - min(tauhatHajs), max(tauhatHajs) - tau_true)
x2 = seq(tau_true - range2, tau_true + range2, length.out=100)
dat_dens21 = data.frame(x = x2, den = dnorm(x2, tau_true, sqrt(sigmaHa12/n)))
dat_dens20 = data.frame(x = x2, den = dnorm(x2, tau_true, sqrt(sigmaHa02/n)))
dfgg = with(dat_his,dat_dens1, dat_dens0)

##################
##  Make plots  ##
##################
library(latex2exp)
library(ggplot2)

cols <- c("LINE1"="red","LINE2"="blue")
his <- c("hist"="grey")

colors <- c("Correct for intervention" = "red", "Not correct for intervention" = "blue")
fills <- c("Estimator" = "grey")
p = ggplot(dfgg)
p = p + geom_histogram(data = dat_his, aes(x = tauhat, y = ..density..,fill = "Estimator"), bins = 20, color = "black")
p = p + geom_line(data = dat_dens1, aes(x = x, y = den, color = "Correct for intervention"))
p = p + geom_line(data = dat_dens0, aes(x = x, y = den, color = "Not correct for intervention"))
p = p + scale_color_manual(name="Normal Approximations (CLT)", values = colors)
p = p + scale_fill_manual(name="Estimators", values = fills)
#p = p +  ggtitle(TeX("Histogram of $\\hat{\\tau}_{DIR}$")) + xlab(TeX("$\\hat{\\tau}_{DIR}$"))
p = p + xlab("Horvitz-Thompson estimator")
p = p + theme_bw()
#p


cols <- c("LINE1"="red","LINE2"="blue")
his <- c("hist"="grey")

colors <- c("Correct for intervention" = "red", "Not correct for intervention" = "blue")
fills <- c("Estimator" = "grey")
p2 = ggplot(dfgg)
p2 = p2 + geom_histogram(data = dat_his2, aes(x = tauhat, y = ..density..,fill = "Estimator"), bins = 20, color = "black")
p2 = p2 + geom_line(data = dat_dens21, aes(x = x, y = den, color = "Correct for intervention"))
p2 = p2 + geom_line(data = dat_dens20, aes(x = x, y = den, color = "Not correct for intervention"))
p2 = p2 + scale_color_manual(name="Normal Approximations (CLT)", values = colors)
p2 = p2 + scale_fill_manual(name="Estimators", values = fills)
#p2 = p2 +  ggtitle(TeX("Histogram of $\\hat{\\tau}_{DIR}$")) + xlab(TeX("$\\hat{\\tau}_{DIR}$"))
p2 = p2 + xlab(TeX("H\'ajek estimator"))
p2 = p2 + theme_bw()
#p2

library(ggpubr)
pp = ggarrange(
  p, p2, 
  common.legend = TRUE, legend = "none"
)

library(tidyverse)
pp %>% ggsave(file="Hist_DIR.pdf", width=8, height=4, units="in")

pp

