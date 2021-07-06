##choose which plot
which_fig = 3

n = 100
u1 = ((1:n)-0.5)/n
u2 = u1
r = 3

withinp = 0.8
betweenp = 0.2

Psi_fun = function(x,y){
  same_group = (x < y/r) & (x > (y-1)/r)
  
  if(which_fig == 1){
    return(same_group)
  }
  else{
    return(x^y)
  }
}

if(which_fig == 1){
    lambdas = c(1,1,1)
} else{
    lambdas = c(1,-2,1)*27/4
}
Psi = outer(u1, 1:r, Psi_fun)

if(which_fig == 2){
    G = Psi %*% diag(lambdas) %*% t(Psi)
}else if(which_fig == 1){
    G = (withinp-betweenp)*Psi %*% diag(lambdas) %*% t(Psi) + betweenp
}else{
    ff = function(ui,uj){
      return(0.25 + floor(3*min(ui,uj))/4)
    }
    ff_v = Vectorize(ff)
    G = outer(u1,u2,ff_v)
}

xs = NULL
ys = NULL
xs = rep(u1,n)
ys = rep(u2,each = n)
zs = as.vector(G)

library(plotly)
fig <- plot_ly(x = xs, y = ys, z = zs, size = 0.5, 
               marker = list(color = zs, line = list(color= zs, width = 0.1))
               )
fig <- fig %>% add_markers()
fig <- fig %>% layout(scene = list(
                        xaxis = list(title = "u1"),
                        yaxis = list(title = "u2"),
                        zaxis = list(title = "G(u1,u2)")
                      ))
fig
