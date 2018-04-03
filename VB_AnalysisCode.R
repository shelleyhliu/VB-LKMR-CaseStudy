#############################
# Shelley H. Liu
# Mean Field Variational Bayes for Lagged Kernel Machine Regression - Case Study
#############################
#############################

doKernelApproach = TRUE
cofactor = 1e-7

if (doKernelApproach) source("VB_SetupCode.R")

# Set tasks:

doMFVB 			= TRUE
doMCMC 			= FALSE
doPrelimPlots 	= FALSE
plotTRUE		= FALSE
doAnalysis		= FALSE

ptm 			= proc.time()

# Set constant matrices needed for MFVB:

XTX 			= crossprod(X)
WTW 			= crossprod(W)

# Set hyperparameters
r1 				= 60 
r2 				= 45 
delta1 			= 10
delta2 			= 10

if (doMFVB)
{

   # Do iter_num iterations of MFVB updating:

   iter_num 			= 200

   E.q.recip.sigsq 		= 1
   E.q.beta 			= as.matrix(rep(1, ncol(X)))
   E.q.recip.tausq 		= rep(1, T)
   E.q.recip.omegasq 	= rep(1, T-1)
   E.q.lambda1sq 		= 1
   E.q.lambda2sq		= 1
   
   A.q.recip.tausq 		= rep(1, T)
   B.q.recip.tausq 		= rep(1, T)
   A.q.recip.omegasq 	= rep(1, T-1)
   B.q.recip.omegasq 	= rep(1, T-1)
   
   sig.shape = 85 
   sig.scale = 10
   a = sig.shape
   gamma = sig.scale
   
   list.G 				= list()
   for (g in 1:T) {
   		list.G[[g]]		= solve(kernelMatrix(poly, Z[,(M*g-(M-1)):(M*g)]) + cofactor*diag(N))
   	}
 
	for (itnum in 1:iter_num) {
	
	list.GT = list()
	for (g in 1:T) {
		list.GT[[g]] = list.G[[g]] * E.q.recip.tausq[g]
	}

	cov.1 				= do.call(adiag, list.GT)
	
	# 2. Diagonal

	list.sig.sqf 		= c()
	list.sig.sqf[1] 	= E.q.recip.omegasq[1]
	list.sig.sqf[T] 	= E.q.recip.omegasq[T-1]
	for (g in 2:(T - 1)) {
		list.sig.sqf[g] = E.q.recip.omegasq[g] + E.q.recip.omegasq[g-1]
	}

	list.mat.sig.sqf 	= list()
	for (g in 1:T) {
			list.mat.sig.sqf[[g]] = list.sig.sqf[g] * diag(TauG) 
	}

	cov.2 				= do.call(adiag, list.mat.sig.sqf)

	#3. Off-diagonals
	cov.3 				= offDiagonal((-1)*(rep(E.q.recip.omegasq, each=TauG))) 

	    E.q.Sigma.h.inverse = cov.1 + cov.2 + cov.3 
   
	Cov.q.h 			= solve( E.q.recip.sigsq * WTW + E.q.Sigma.h.inverse + cofactor*diag(P)) 
	E.q.h				= E.q.recip.sigsq * Cov.q.h %*% t(W) %*% (Y - X %*% E.q.beta)
	E.q.beta 			= solve(XTX) %*% t(X) %*% (Y - W %*% E.q.h)

        B.q.recip.sigsq = as.numeric(crossprod(Y - X %*% E.q.beta - W %*% E.q.h)/2) + sig.scale 
               
        A.q.recip.sigsq = N/2 + sig.shape 

		E.q.recip.sigsq 	= A.q.recip.sigsq/B.q.recip.sigsq
	
	list.h 				= list()
	
	for (g in 1:T) {
		list.h[[g]] 		= E.q.h[(1+N*(g-1)):(N*g)]               
		A.q.recip.tausq[g] 	= sqrt( E.q.lambda1sq / (t(list.h[[g]]) %*% list.G[[g]] %*% list.h[[g]]) ) 
		B.q.recip.tausq[g] 	= E.q.lambda1sq
		E.q.recip.tausq[g] 	= A.q.recip.tausq[g]
	}	
	
	for (g in 1:(T-1)) {
		A.q.recip.omegasq[g] = sqrt( E.q.lambda2sq / ( sum((list.h[[g+1]] - list.h[[g]])^2)) )
		B.q.recip.omegasq[g] = E.q.lambda2sq
		E.q.recip.omegasq[g] = A.q.recip.omegasq[g]
	}
		
	A.q.lambda1sq = T*(N+1)/2 + r1
    B.q.lambda1sq = 1/2*sum(sapply(E.q.recip.tausq, function(x) 1/x)) + delta1
    E.q.lambda1sq = A.q.lambda1sq/B.q.lambda1sq
   
    A.q.lambda2sq = (T - 1 + r2)
    B.q.lambda2sq = (1/2*sum(sapply(E.q.recip.omegasq, function(x) 1/x)) + delta2)
	E.q.lambda2sq = A.q.lambda2sq / B.q.lambda2sq 
	
	# log(p;y)
	# Increment
	logpy1 = -N/2*log(2*pi) - 1/2*t(E.q.h) %*% E.q.Sigma.h.inverse %*% E.q.h 
	logpy2 =  1/2*t(E.q.h - E.q.recip.sigsq * solve( E.q.recip.sigsq * WTW + E.q.Sigma.h.inverse + cofactor*diag(P)) %*% Cov.q.h %*% t(W) %*% (Y - X %*% E.q.beta)) %*% ( E.q.recip.sigsq * WTW + E.q.Sigma.h.inverse + cofactor*diag(P)) %*% (E.q.h - E.q.recip.sigsq * solve( E.q.recip.sigsq * WTW + E.q.Sigma.h.inverse + cofactor*diag(P)) %*% Cov.q.h %*% t(W) %*% (Y - X %*% E.q.beta))
	logpy3 = P/2*log(2*pi) + log(det(solve(E.q.recip.sigsq*XTX))^(1/2))
	logpy4 = 1/2*E.q.recip.sigsq * t(E.q.beta - solve(XTX) %*% t(X) %*% (Y - W %*% E.q.h)) %*% XTX %*% (E.q.beta - solve(XTX) %*% t(X) %*% (Y - W %*% E.q.h)) # this equals zero? 
	logpy5 = a * log(gamma) - log(gamma(a)) + log(N/2 + a) - (N/2 + 1)*log(  (t(Y - X %*% E.q.beta - W %*% E.q.h) %*% (Y - X %*% E.q.beta - W %*% E.q.h) )/2 + gamma)
	logpy6 = N*T/2 + T/2*log(pi) - N*sum(log(E.q.recip.tausq))
	logpy7 = -E.q.lambda1sq/2 * sum(1/E.q.recip.tausq) + E.q.lambda1sq/2 * sum( (E.q.recip.tausq - sqrt( E.q.lambda1sq / (t(list.h[[g]]) %*% list.G[[g]] %*% list.h[[g]]) ))^2 / ( sqrt( E.q.lambda1sq / (t(list.h[[g]]) %*% list.G[[g]] %*% list.h[[g]]) ) * E.q.recip.tausq) )
	logpy8 = (T-1)/2 * log( E.q.lambda2sq/2) + (T-1)/2*log(pi) - sum(E.q.recip.omegasq)/2
	logpy9 = E.q.lambda2sq/2 *sum(1/E.q.recip.omegasq) + E.q.lambda2sq/2 * sum( (E.q.recip.omegasq - sqrt( E.q.lambda2sq / ( sum((list.h[[g+1]] - list.h[[g]])^2)) ))^2 / (sqrt( E.q.lambda2sq / ( sum((list.h[[g+1]] - list.h[[g]])^2)) ) * E.q.recip.omegasq) )
	logpy10 = r2 * log(delta2)- log(gamma(r2)) + log(gamma(T - 1 + r2)) 
	logpy11 = -(T-1 + r2)*log(sum((1/E.q.recip.omegasq)/2) + delta2) - (T-1)*log(E.q.lambda2sq)
	logpy12 = sum((1/E.q.recip.omegasq)/2)*E.q.lambda2sq + r1*log(delta1) - log(gamma(r1)) 
	logpy13 = -(T*(N+1)/2 + r1)*log( sum( (1/E.q.recip.tausq)/2 ) + delta1) - T*(N+1)/2 * log(E.q.lambda1sq)
	logpy14 = sum( (1/E.q.recip.tausq)/2 ) * E.q.lambda1sq
	
	logpy = logpy1 + logpy2 + logpy3 + logpy4 + logpy5 + logpy6 + logpy7 + logpy8 + logpy9 + logpy10 + logpy11 + logpy12 + logpy13 + logpy14
	
	print(c(itnum, logpy))
	logpyvec[itnum+1] = logpy
	if (abs(logpy - logpyvec[itnum]) < 10e-8) break
	
   }

  proc.time() - ptm
}


