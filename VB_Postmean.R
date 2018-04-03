# Postmean functions for VB, all updated according to paper.
# Code last updated 1/31/2017.

#########################################
# Timepoint 1. Trimester 2 of pregnancy #
#########################################

newh.postmean.gfl.1 <- function(Time, g, numtime, Znew, sel) {
	if(is.null(dim(Znew))) Znew 			= matrix(Znew, nrow=1)
	if(class(Znew) == "data.frame") Znew 	= data.matrix(Znew)
	
	Zall 			= Z
	Z 				= Z[,Time]
	Znew 			= as.matrix(Znew)
	W 				= diag(n)
	counter 		= 1
	while(counter < numtime) {
		W 			= cbind(W, diag(n))
		counter 	= counter+1
	}
	
	Y 				= Y
	pnew 			= nrow(Znew)

	sigsq.eps 		= 1/E.q.recip.sigsq
	sigsqf.mean 	= 1/E.q.recip.omegasq # posterior mean of omegas
	tausqf.mean 	= 1/E.q.recip.tausq # posterior mean of taus
	beta.mean 		= E.q.h #posterior mean of NK betas
	conf.mean		= E.q.beta
		 
	#######################################
	# Recreate the Sigma_h inverse matrix #
	#######################################

	# Diagonals 
	vec.sigsqf 				= c()
	vec.sigsqf[1] 			= 1/sigsqf.mean[1]
	for (ind in 2:(numtime - 1)) {
		vec.sigsqf[ind] 	= 1/sigsqf.mean[ind] + 1/sigsqf.mean[ind-1]
	}
	vec.sigsqf[numtime] 	= 1/sigsqf.mean[numtime-1]
	
	# reordering the omegas so it corresponds to h_2,h_3,h_1
	list.mat.sig.sqf.diag 				= list()
	for (ind in 2:numtime) {
		list.mat.sig.sqf.diag[[ind-1]] 	= vec.sigsqf[ind]*diag(n)
	}
	list.mat.sig.sqf.diag[[numtime]] 	= vec.sigsqf[1]*diag(n)
	
	cov.omega.diag 						= do.call(adiag, list.mat.sig.sqf.diag)
   
   	# Off diagonals 
	cov.omega.off.diag 		= offDiagonal((-1)*c(rep(1/sigsqf.mean[2:2],each=n), rep(0,n)) )
	delta.diag 				= row(cov.omega.off.diag) - col(cov.omega.off.diag)
	cov.omega.off.diag[abs(delta.diag) == (n*grpsize - n)] == -1/sigsqf.mean[1]
	cov.omega 				= cov.omega.diag + cov.omega.off.diag
	
	cov.omega.new 			= "[<-"(matrix(0, (numtime*n+pnew), (numtime*n+pnew)), 1:nrow(cov.omega), 1:ncol(cov.omega), value = cov.omega) # adding rows/cols of zeroes to square matrix

	# 1. Block diagonal matrix - just T unique tau's here
    list.G 					= list()
    poly 					= polydot(degree=2, offset=1)

	for (ind in 2:(numtime)) {
		list.G[[ind-1]] 	= kernelMatrix(poly, Zall[,(grpsize*ind-(grpsize-1)):(grpsize*ind)])	
		list.G[[ind-1]] 	= solve(list.G[[ind-1]] + cofactor*diag(n))/tausqf.mean[ind]
	}   
	list.G[[numtime]] 		= solve(kernelMatrix(poly, rbind(Z, Znew)) + cofactor*diag(n+pnew))/tausqf.mean[1]

	# Contruction of the last time point, with the new z's
    cov.bf1 = do.call(adiag, list.G)   		
	cov.bf = cov.bf1 + cov.omega.new
	
	###################
   	# Sigma_h inverse #
   	###################
   	
	A = cov.bf[1:(n*numtime), 1:(n*numtime)]
	B = cov.bf[1:(n*numtime), (n*numtime + 1):(n*numtime + pnew)] #pnew is number of new subjects
	C = cov.bf[(n*numtime + 1):(n*numtime + pnew), 1:(n*numtime)]
	D = cov.bf[(n*numtime + 1):(n*numtime + pnew), (n*numtime + 1):(n*numtime + pnew)]
		
	###########	
	# Sigma_h #
	###########
	
	if (doginv) {	
		Sigma11 = ginv(A - B %*% ginv(D) %*% C)
		Sigma12 = -1*ginv(A) %*% B %*% ginv(D - C %*% ginv(A) %*% B)
		Sigma21 = -1*ginv(D) %*% C %*% ginv(A - B %*% ginv(D) %*% C)
		Sigma22 = ginv(D - C %*% ginv(A) %*% B)
	} else {
		Sigma11 = solve(A - B %*% solve(D + cofactor*diag(nrow(D))) %*% C + 					cofactor*diag(nrow(A)))
		Sigma12 = -1*solve(A + cofactor*diag(nrow(A))) %*% B %*% solve(D - C %*% 			solve(A + cofactor*diag(nrow(A))) %*% B + cofactor*diag(nrow(D)))
		Sigma21 = -1*solve(D + cofactor*diag(nrow(D))) %*% C %*% solve(A - B %*% 			ginv(D) %*% C + cofactor*diag(nrow(A)))
		Sigma22 = solve(D - C %*% solve(A + cofactor*diag(nrow(A))) %*% B + 					cofactor*diag(nrow(D)))
	}
	
	########################################
	# Approximated posterior mean of h_new #
	########################################
	
	Mu.hnew.updated = Sigma21 %*% solve(Sigma11 + cofactor*diag(nrow(Sigma11))) %*% solve(1/sigsq.eps*t(W) %*% W + solve(Sigma11 + cofactor*diag(nrow(Sigma11)))) %*% t(W) %*% (Y - X%*%conf.mean) * 1/sigsq.eps
		
	############################################
	# Approximated posterior variance of h_new #
	############################################
	
	Sigma.hnew.updated = Sigma22 - Sigma21 %*% solve(Sigma11 + cofactor*diag(nrow(Sigma11))) %*% Sigma12 + Sigma21 %*% solve(Sigma11 + cofactor*diag(nrow(Sigma11))) %*% (solve(1/sigsq.eps*t(W) %*% W + solve(Sigma11 + cofactor*diag(nrow(Sigma11))))) %*% solve(Sigma11 + cofactor*diag(nrow(Sigma11))) %*% Sigma12 
	
	##########
	# Output #
	##########
	
	list(postmean=drop(Mu.hnew.updated), postvar=drop(Sigma.hnew.updated))

}

#########################################
# Timepoint 2. Trimester 3 of pregnancy #
#########################################

newh.postmean.gfl.2 <- function(Time, g, numtime, Znew, sel) {
	
	if(is.null(dim(Znew))) Znew 			= matrix(Znew, nrow=1)
	if(class(Znew) == "data.frame") Znew 	= data.matrix(Znew)
	
	Zall 			= Z
	Z 				= Z[,Time]
	Znew 			= as.matrix(Znew)
	W 				= diag(n)
	counter 			= 1
	while(counter < numtime) {
		W 			= cbind(W, diag(n))
		counter 	= counter+1
	}
	
	Y 				= Y	
	pnew 			= nrow(Znew)
	
	sigsq.eps 		= 1/E.q.recip.sigsq
	sigsqf.mean 	= 1/E.q.recip.omegasq # posterior mean of omegas
	tausqf.mean 	= 1/E.q.recip.tausq # posterior mean of taus
	beta.mean 		= E.q.h #posterior mean of NK betas
	conf.mean		= E.q.beta #posterior mean of NK betas
	 
	#######################################
	# Recreate the Sigma_h inverse matrix #
	#######################################

	# Diagonals 
	vec.sigsqf 				= c()
	vec.sigsqf[1] 			= 1/sigsqf.mean[1]
	for (ind in 2:(numtime - 1)) {
		vec.sigsqf[ind] 	= 1/sigsqf.mean[ind] + 1/sigsqf.mean[ind-1]
	}
	vec.sigsqf[numtime] 	= 1/sigsqf.mean[numtime-1]
	
	# reordering the omegas so it corresponds to h_1,h_3,h_2
	list.mat.sig.sqf.diag 				= list()
	for (ind in 1:1) {
		list.mat.sig.sqf.diag[[ind]] 	= vec.sigsqf[ind]*diag(n)
	}
	list.mat.sig.sqf.diag[[2]] 			= vec.sigsqf[3]*diag(n)
	list.mat.sig.sqf.diag[[numtime]] 	= vec.sigsqf[2]*diag(n)
	
	cov.omega.diag = do.call(adiag, list.mat.sig.sqf.diag)
   
   	# Off diagonals 
	cov.omega.off.diag 	= offDiagonal((-1)*c(rep(1/sigsqf.mean[1], n), rep(0,n)) )
	delta.diag 			= row(cov.omega.off.diag) - col(cov.omega.off.diag) 
	cov.omega.off.diag[abs(delta.diag) == (n*grpsize - 3*n)] == c( rep(0,n), rep(-1/sigsqf.mean[2], n) )
	cov.omega 			= cov.omega.diag + cov.omega.off.diag
	
	cov.omega.new 		= "[<-"(matrix(0, (numtime*n+pnew), (numtime*n+pnew)), 1:nrow(cov.omega), 1:ncol(cov.omega), value = cov.omega) # adding rows/cols of zeroes to square matrix

	# 1. Block diagonal matrix - just T unique tau's here - reordered
    list.G 				= list()
    poly 				= polydot(degree=2, offset=1)

	for (ind in 1:1) {
		list.G[[ind]] 	= kernelMatrix(poly, Zall[,(grpsize*ind-(grpsize-1)):(grpsize*ind)])	
		list.G[[ind]] 	= solve(list.G[[ind]] + cofactor*diag(n))/tausqf.mean[ind]
	}   
	list.G[[2]] 		= kernelMatrix(poly, Zall[,(grpsize*2+1):(3*grpsize)])
    list.G[[2]] 		= solve(list.G[[2]] + cofactor*diag(n))/tausqf.mean[3]
	list.G[[numtime]] 	= solve(kernelMatrix(poly, rbind(Z, Znew)) + cofactor*diag(n+pnew))/tausqf.mean[2]
	# Contruction of the last time point, with the new z's
    cov.bf1 			= do.call(adiag, list.G)      		
	cov.bf 				= cov.bf1 + cov.omega.new
	
	###################
   	# Sigma_h inverse #
   	###################
	
	A = cov.bf[1:(n*numtime), 1:(n*numtime)]
	B = cov.bf[1:(n*numtime), (n*numtime + 1):(n*numtime + pnew)] #pnew is number of new subjects
	C = cov.bf[(n*numtime + 1):(n*numtime + pnew), 1:(n*numtime)]
	D = cov.bf[(n*numtime + 1):(n*numtime + pnew), (n*numtime + 1):(n*numtime + pnew)]
		
	###########	
	# Sigma_h #
	###########
	
	if (doginv) {	
		Sigma11 = ginv(A - B %*% ginv(D) %*% C)
		Sigma12 = -1*ginv(A) %*% B %*% ginv(D - C %*% ginv(A) %*% B)
		Sigma21 = -1*ginv(D) %*% C %*% ginv(A - B %*% ginv(D) %*% C)
		Sigma22 = ginv(D - C %*% ginv(A) %*% B)
	} else {
		Sigma11 = solve(A - B %*% solve(D + cofactor*diag(nrow(D))) %*% C + 					cofactor*diag(nrow(A)))
		Sigma12 = -1*solve(A + cofactor*diag(nrow(A))) %*% B %*% solve(D - C %*% 			solve(A + cofactor*diag(nrow(A))) %*% B + cofactor*diag(nrow(D)))
		Sigma21 = -1*solve(D + cofactor*diag(nrow(D))) %*% C %*% solve(A - B %*% 			ginv(D) %*% C + cofactor*diag(nrow(A)))
		Sigma22 = solve(D - C %*% solve(A + cofactor*diag(nrow(A))) %*% B + 					cofactor*diag(nrow(D)))
	}
	
	########################################
	# Approximated posterior mean of h_new #
	########################################
	
	Mu.hnew.updated = Sigma21 %*% solve(Sigma11 + cofactor*diag(nrow(Sigma11))) %*% solve(1/sigsq.eps*t(W) %*% W + solve(Sigma11 + cofactor*diag(nrow(Sigma11)))) %*% t(W) %*% (Y - X%*%conf.mean) * 1/sigsq.eps
		
	############################################
	# Approximated posterior variance of h_new #
	############################################
	
	Sigma.hnew.updated = Sigma22 - Sigma21 %*% solve(Sigma11 + cofactor*diag(nrow(Sigma11))) %*% Sigma12 + Sigma21 %*% solve(Sigma11 + cofactor*diag(nrow(Sigma11))) %*% (solve(1/sigsq.eps*t(W) %*% W + solve(Sigma11 + cofactor*diag(nrow(Sigma11))))) %*% solve(Sigma11 + cofactor*diag(nrow(Sigma11))) %*% Sigma12 
	
	##########
	# Output #
	##########
	
	list(postmean=drop(Mu.hnew.updated), postvar=drop(Sigma.hnew.updated))
}

#######################################
# Timepoint 3: Months 0-3 after birth #
#######################################

newh.postmean.gfl.3 = function(Time, g, numtime, Znew, sel) {
	
	if(is.null(dim(Znew))) Znew 			= matrix(Znew, nrow=1)
	if(class(Znew) == "data.frame") Znew 	= data.matrix(Znew)
	
	Zall 			= Z
	Z 				= Z[,Time]
	Znew 			= as.matrix(Znew)
	W 				= diag(n)
	counter 			= 1
	
	while(counter < numtime) {
		W 			= cbind(W, diag(n))
		counter 	= counter+1
	}
	
	Y 				= Y
	pnew 			= nrow(Znew)
	
	sigsq.eps 		= 1/E.q.recip.sigsq
	sigsqf.mean 	= 1/E.q.recip.omegasq # posterior mean of omegas
	tausqf.mean 	= 1/E.q.recip.tausq # posterior mean of taus
	beta.mean 		= E.q.h #posterior mean of NK betas
	conf.mean		= E.q.beta
	
	#######################################
	# Recreate the Sigma_h inverse matrix #
	#######################################
	
	# Diagonals
	vec.sigsqf 				= c()
	vec.sigsqf[1] 			= 1/sigsqf.mean[1]
	for (ind in 2:(numtime - 1)) {
		vec.sigsqf[ind] 	= 1/sigsqf.mean[ind] + 1/sigsqf.mean[ind-1]
	}
	vec.sigsqf[numtime] 	= 1/sigsqf.mean[numtime-1]
	
	list.mat.sig.sqf.diag 	= list()
	for (ind in 1:numtime) {
		list.mat.sig.sqf.diag[[ind]] = vec.sigsqf[ind]*diag(n)
	}
	cov.omega.diag 			= do.call(adiag, list.mat.sig.sqf.diag)
   
   	# Off diagonals
	cov.omega.off.diag 		= offDiagonal((-1)*rep(1/sigsqf.mean,each=n)) ## edited to 1/sigsqf.mean instead of sigsqf.mean
	cov.omega 				= cov.omega.off.diag + cov.omega.diag
	
	cov.omega.new 			= "[<-"(matrix(0, (numtime*n+pnew), (numtime*n+pnew)), 1:nrow(cov.omega), 1:ncol(cov.omega), value = cov.omega) # adding rows/cols of zeroes to square matrix

	# 1. Block diagonal matrix - just T unique tau's here
    list.G 					= list()
    poly 					= polydot(degree=2, offset=1)
    
	for (ind in 1:(numtime-1)) {
		list.G[[ind]] 		= solve(kernelMatrix(poly, Zall[,(grpsize*ind-(grpsize-1)):(grpsize*ind)]) + cofactor*diag(n))/tausqf.mean[ind]	
	}   
	list.G[[numtime]] 		= solve(kernelMatrix(poly, rbind(Z, Znew)) + cofactor*diag(n+pnew))/tausqf.mean[numtime]
	
	# Contruction of the last time point, with the new z's
    cov.bf1 				= do.call(adiag, list.G)
    cov.bf 					= cov.bf1 + cov.omega.new
   	
   	###################
   	# Sigma_h inverse #
   	###################
   		
	A 	= cov.bf[1:(n*numtime), 1:(n*numtime)]
	B 	= cov.bf[1:(n*numtime), (n*numtime + 1):(n*numtime + pnew)] #pnew is number of new subjects
	C 	= cov.bf[(n*numtime + 1):(n*numtime + pnew), 1:(n*numtime)]
	D 	= cov.bf[(n*numtime + 1):(n*numtime + pnew), (n*numtime + 1):(n*numtime + pnew)]
	
	###########	
	# Sigma_h #
	###########
	
	if (doginv) {	
		Sigma11 = ginv(A - B %*% ginv(D) %*% C)
		Sigma12 = -1*ginv(A) %*% B %*% ginv(D - C %*% ginv(A) %*% B)
		Sigma21 = -1*ginv(D) %*% C %*% ginv(A - B %*% ginv(D) %*% C)
		Sigma22 = ginv(D - C %*% ginv(A) %*% B)
	} else {
		Sigma11 = solve(A - B %*% solve(D + cofactor*diag(nrow(D))) %*% C + 					cofactor*diag(nrow(A)))
		Sigma12 = -1*solve(A + cofactor*diag(nrow(A))) %*% B %*% solve(D - C %*% 			solve(A + cofactor*diag(nrow(A))) %*% B + cofactor*diag(nrow(D)))
		Sigma21 = -1*solve(D + cofactor*diag(nrow(D))) %*% C %*% solve(A - B %*% solve(D + cofactor*diag(nrow(D))) %*% C + cofactor*diag(nrow(A))) 
		Sigma22 = solve(D - C %*% solve(A + cofactor*diag(nrow(A))) %*% B + 					cofactor*diag(nrow(D)))
	}
	
	########################################
	# Approximated posterior mean of h_new #
	########################################
	
	Mu.hnew.updated = Sigma21 %*% solve(Sigma11 + cofactor*diag(nrow(Sigma11))) %*% solve(1/sigsq.eps*t(W) %*% W + solve(Sigma11 + cofactor*diag(nrow(Sigma11)))) %*% t(W) %*% (Y - X%*%conf.mean) * 1/sigsq.eps
		
	############################################
	# Approximated posterior variance of h_new #
	############################################
	
	Sigma.hnew.updated = Sigma22 - Sigma21 %*% solve(Sigma11 + cofactor*diag(nrow(Sigma11))) %*% Sigma12 + Sigma21 %*% solve(Sigma11 + cofactor*diag(nrow(Sigma11))) %*% (solve(1/sigsq.eps*t(W) %*% W + solve(Sigma11 + cofactor*diag(nrow(Sigma11))))) %*% solve(Sigma11 + cofactor*diag(nrow(Sigma11))) %*% Sigma12 
	
	##########
	# Output #
	##########
	
	list(postmean=drop(Mu.hnew.updated), postvar=drop(Sigma.hnew.updated))
}




