#############################
# Create simulation data
#############################

source("VB_DataCodeSimulation.R")

N 				= n
T 				= numtime # Number of time points
P 				= N*T
M				= grpsize # Number of metals

TauG = N 

poly = polydot(degree = 2, offset = 1)

offDiagonal = function(x) {
		diagBelow = diag(x)
		i = 1
		while (i <= TauG) { 
			diagBelow=rbind(rep(0,length(x)	+i),cbind(diagBelow,rep(0,length(x) + i - 1)))
			i = i + 1
		}
		mat <- diagBelow + t(diagBelow) - diag(diag(diagBelow))
		return(mat)
}

#############################

W 				= diag(N)
counter 		= 1
while(counter < T) {
	W 		= cbind(W, diag(N))
	counter = counter+1
}
