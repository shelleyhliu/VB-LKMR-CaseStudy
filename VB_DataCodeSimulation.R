# Simulation case study

library(lars); library(SuppDists); library(MCMCpack); library(magic); library(lasso2); library(mvtnorm); library(grplasso); library(kernlab); library(MASS); library(fields)

n      		= 100
grpsize 	= 5
numtime 	= 3
p			    = n*numtime
c 			  = 2 # number of confounders
ar 			  = 0.5

autocorr.mat <- function(p, rho) {
  mat <- diag(p)
  return(rho^abs(row(mat)-col(mat)))
}
time_mat = autocorr.mat(numtime,ar)

metal_mat <- matrix( c(1, 0.19, 0.32, 0.15, 0.4, 0.19, 1, 0.11, 0.25, 0.4, 0.32, 0.11, 1, 0.17, 0.24, 0.15, 0.25, 0.17, 1, 0.35, 0.4, 0.4, 0.24, 0.35, 1), 5,5, byrow=TRUE) 

cmat = kronecker(time_mat, metal_mat)
obs <- matrix( mvrnorm(n, mu=rep(0, 5*3), Sigma=cmat), nrow = n)

Z = obs

Z.1 = Z[,1:5]
Z.2 = Z[,6:10]
Z.3 = Z[,11:15]

beta = rep(1, c)
X = cbind(scale(matrix(rnorm(n, mean = 10, sd = 1))), scale(matrix(sample(1:2, n, replace = TRUE))))

res.sd = 1 
Y = matrix(rnorm(n=n, mean=0, sd=res.sd)) + X%*%beta + 0.5*(Z.1[,1]^2 - Z.1[,2]^2 + 1/2*Z.1[,1]*Z.1[,2] + Z.1[,1] + Z.1[,2]) + 1*(Z.2[,1]^2 - Z.2[,2]^2 + 1/2*Z.2[,1]*Z.2[,2] + Z.2[,1] + Z.2[,2]) + 1.5*(Z.3[,1]^2 - Z.3[,2]^2 + 1/2*Z.3[,1]*Z.3[,2] + Z.3[,1] + Z.3[,2] + Z.3[,4])
Y = scale(Y)

offDiagonal = function(x) {
  diagBelow = diag(x)
  i = 1
  while (i <= n) {
    diagBelow=rbind(rep(0,length(x)	+i),cbind(diagBelow,rep(0,length(x) + i - 1)))
    i = i + 1
  }
  mat <- diagBelow + t(diagBelow) - diag(diag(diagBelow))
  return(mat)
}

beta = as.matrix(beta)