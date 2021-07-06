args = commandArgs(TRUE)
arg = as.numeric(args[1])
seed = arg

set.seed(seed)
library(RSpectra)
library(Matrix)
ptm <- proc.time()
n = 1000000
N = 25
r = 3
rhon = n^(-2/5)/2
#rhon = 1

#alpha0 = 0.5
#alpha1 = 0.8
#H_fun <- function(x,y){return(alpha0 + alpha1/2*sin(2*pi*x*y)/sqrt(r))}
#stochastic block model
withinp = 0.8
betweenp = 0.2
propensity = 0.4

tau_tildes = NULL

##it's the same if we reorder the indexes by U
##sample the blocks first

for (num in 1:N){
  print(num)
  n3 = rbinom(1,n,1/3)
  n1 = rbinom(1,n-n3,1/2)
  n2 = n-n3-n1
  block_sizes = c(n1,n2,n3)
  block_cum = c(0,cumsum(block_sizes)[c(1,2)])
  
  i_list = rep(0,n^2*rhon)
  j_list = rep(0,n^2*rhon)
  x_no = 0
  
  for(bi in 1:3){
    for(bj in 1:3){
      same_group = (bi == bj)
      prob = same_group*withinp + (1-same_group)*betweenp
      block_sizei = block_sizes[bi]
      block_sizej = block_sizes[bj]
      k = rbinom(1,block_sizei*block_sizej, prob*rhon)
      ijs = sample(0:(block_sizei*block_sizej-1), k)
      i_list[x_no + 1:k] = ijs %% block_sizei + 1 + block_cum[bi]
      j_list[x_no + 1:k] = ijs %/% block_sizei + 1 + block_cum[bj]
      x_no = x_no + k
    }
  }
  i_list_short = i_list[1:x_no]
  j_list_short = j_list[1:x_no]
  Eij0 = sparseMatrix(i = i_list_short, j = j_list_short, x = 1, dims = c(n,n))
  Eij = triu(Eij0, k= 1) + t(triu(Eij0, k = 1))
  #Eij = forceSymmetric(Eij)
  
  U = sort(runif(n))
  theta = U
  W = rbinom(n,1,propensity) + 0
  Ni = as.numeric(Eij %*% rep(1,n))
  Mi = as.numeric(Eij %*% W)
  prop_treated_ngb = Mi/pmax(Ni,1)
  Y = 1/2*(W + theta*prop_treated_ngb)^2 + 0.2*rnorm(n)
  # 
  # PC balancing estimator
  res = eigs_sym(Eij,r)
  #res = eigs_sym(Gij,r)
  Psi_hat = res$vectors
  gamma_original = Mi/propensity - (Ni - Mi)/(1-propensity)
  beta = -t(Psi_hat) %*% gamma_original
  tau_tilde = mean(Y*(Mi/propensity - (Ni - Mi)/(1-propensity) + Psi_hat %*% beta))
  tau_tildes = c(tau_tildes, tau_tilde)
  rm(list=setdiff(ls(), c("tau_tildes","n","N","r","rhon","propensity","withinp","betweenp","seed")))
  gc(verbose = TRUE)
}

save(tau_tildes, seed, file = paste("output_rankr_106_25_",seed,".Rdata", sep = ""))
