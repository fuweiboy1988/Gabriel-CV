# demo/bench/settings.R
#
# Depends:
#   library("MASS") # mvrnorm
#
##-----------------------------------------------
Cov <- function(dim = 2, rho = 0){
	COV <- diag(dim)
	for(i in 1:dim){
		for(j in 1:dim){
			if(i!= j){COV[i,j] = rho}
		}
	}
	return (COV)
}
##-------------------------------------------------
Data_generate <- function(Dim = 10, C.num = 4, Para = 1, Rho = 0, Obs.num = 100){ 
	condition <- TRUE
	N = 0
  	while(condition){
		N=N+1
		del = 1
		nn = sample(c(Obs.num/2,Obs.num),size=C.num,replace=T)
		cl = rep(1:C.num,nn)
		 
		Centers <- mvrnorm(C.num, mu=rep(0,Dim),Sigma = (Para^2)*diag(Dim))

		for(j in 1:nrow(Centers)){
			if(j==1){
				X.data = mvrnorm(nn[j], mu = Centers[j,], Sigma = Cov(Dim,Rho))
			}else{
				DATA = mvrnorm(nn[j], mu = Centers[j,], Sigma = Cov(Dim,Rho))
				X.data = rbind(X.data, DATA)
			}
		}
		
		ss=dist(rbind(X.data,Centers))
		d=matrix(0,nrow=nrow(X.data)+C.num,ncol=nrow(X.data)+C.num)
 		d[row(d)>col(d)]=ss
 		DD=d[(nrow(X.data)+1):(nrow(X.data)+C.num),1:nrow(X.data)]
		
		for(i in 1:ncol(DD)){ DD[cl[i],i]=DD[cl[i],i]+del}
		ff=apply(DD,2,which.min)

		if(sum(ff==cl)==nrow(X.data)){
 			condition <- FALSE
 		}
	}#end while loop
	ID <- sample(1:nrow(X.data),nrow(X.data),replace = FALSE)
    	result <- list(Data = X.data[ID,], Number = N)
    	return(result)
}

##MM <- Data_generate(dim = 10, C.num = 10, Para = 0.72, Rho = 0, Obs.num = 100)
###-----------------------------------------------------
NumA <- function(P = 1, dim, c.num, rho, obs.num){
	set.seed(1)
	NN <- rep(NA,50)
	for(W in 1:50){
		Med <- Data_generate(Dim = dim, C.num = c.num, Para = P, Rho = rho, Obs.num = obs.num)
		NN[W]<- Med$Number
	}
	return(mean(NN)-2)
}

## NumA(P=3, dim=10, c.num=4, rho=0, obs.num=100)
####=----------------------------------------------------
## Auto-select the parameter Para in above function
Para.hunt <- function(Dim, C.num, Rho, Obs.num, Initial = 1){
	N.init <- NumA(P=Initial, dim=Dim, c.num=C.num, rho=Rho, obs.num=Obs.num)
	if( -0.2 <= N.init & N.init <= 0.2){
		return(Initial)
	}else{
		if(N.init > 0.2){
			return("Please use larger initial value")
		}else{
			N.old <- N.init
			P.old <- Initial
			GAP = 0.05
			Cond = TRUE
			while(Cond){
				P.new = P.old - GAP
				N.new <- NumA(P=P.new, dim=Dim, c.num=C.num, rho=Rho, obs.num=Obs.num)
				if(-0.2 <= N.new & N.new <= 0.2){
					Cond <- FALSE
					return(P.new)
				}else{
					if(N.new > 0.2){
						result <- uniroot(NumA, interval = c(P.new,P.old), tol = 0.01, dim=Dim, c.num=C.num, rho = Rho, obs.num=Obs.num)$root
						Cond <- FALSE
						return(result)
					}else{
						if(N.new != N.old & GAP > 0.02){
							GAP = GAP/2
						}
						P.old = P.new
						N.old = N.new
					}
				}
			}##end while loop
		}
	}##end if
}

## Par <- Para.hunt(Dim = 100, C.num = 10, Rho = 0, Obs.num=100, Initial = 1)
## NumA(P=Par, dim=100, c.num=10, rho=0, obs.num=100)

##==========================================================================
Data3centers <- function(Ratio = 1, Para = 1, Obs.num = 60){ 
	condition <- TRUE
	N = 0
  	while(condition){
		N=N+1
		del = 1
		nn = rep(Obs.num,3)
		cl = rep(1:3,nn)
		 
		Centers <- mvrnorm(3, mu=rep(0,20),Sigma = (Para^2)*diag(20))

		U1 <- mvrnorm(nn[1], mu = Centers[1,], Sigma = diag(20))
		U2 <- mvrnorm(nn[2], mu = Centers[2,], Sigma = (Ratio+1)/2*diag(20))
		U3 <- mvrnorm(nn[3], mu = Centers[3,], Sigma = Ratio*diag(20))
		##----------------------------------------------
		X.data = rbind(U1,U2,U3)

		ss=dist(rbind(X.data,Centers))
		d=matrix(0,nrow=nrow(X.data)+3,ncol=nrow(X.data)+3)
 		d[row(d)>col(d)]=ss
 		DD=d[(nrow(X.data)+1):(nrow(X.data)+3),1:nrow(X.data)]
		
		for(i in 1:ncol(DD)){ DD[cl[i],i]=DD[cl[i],i]+del}
		ff=apply(DD,2,which.min)

		if(sum(ff==cl)==nrow(X.data)){
 			condition <- FALSE
 		}
	}#end while loop
	ID <- sample(1:nrow(X.data),nrow(X.data),replace = FALSE)
    	result <- list(Data = X.data[ID,], Number = N)
    	return(result)
}
###-----------------------------------------------------
NumA_3centers <- function(P = 1, ratio, obs.num){
	set.seed(1)
	NN <- rep(NA,50)
	for(W in 1:50){
		Med <- Data3centers(Ratio = ratio, Para = P, Obs.num = obs.num)
		NN[W]<- Med$Number
	}
	return(mean(NN)-2)
}

## NumA_3centers(P = 3, ratio = 5, obs.num = 60)
####=----------------------------------------------------
## Auto-select the parameter Para in above function
Para_3centers <- function(Ratio = 1, Obs.num = 60, Initial = 1){
	N.init <- NumA_3centers(P = Initial, ratio = Ratio, obs.num = Obs.num)
	if( -0.2 <= N.init & N.init <= 0.2){
		return(Initial)
	}else{
		if(N.init > 0.2){
			return("Please use larger initial value")
		}else{
			N.old <- N.init
			P.old <- Initial
			GAP = 0.05
			Cond = TRUE
			while(Cond){
				P.new = P.old - GAP
				N.new <- NumA_3centers(P = P.new, ratio = Ratio, obs.num = Obs.num)
				if(-0.2 <= N.new & N.new <= 0.2){
					Cond <- FALSE
					return(P.new)
				}else{
					if(N.new > 0.2){
						result <- uniroot(NumA_3centers, interval = c(P.new,P.old), tol = 0.01, ratio = Ratio, obs.num=Obs.num)$root
						Cond <- FALSE
						return(result)
					}else{
						if(N.new != N.old & GAP > 0.02){
							GAP = GAP/2
						}
						P.old = P.new
						N.old = N.new
					}
				}
			}##end while loop
		}
	}##end
}

## PP <- Para_3centers(Ratio = 5, Obs.num = 60, Initial = 3)
## NumA_3centers(P = PP, ratio = 5, obs.num = 60)

##==========================================================================
Heavytail <- function(Df = 2, Para = 10, Obs.num = 80){
	condition <- TRUE
	N = 0
  	while(condition){
		N=N+1
		del = 1
		nn = rep(Obs.num,5)
		cl = rep(1:5,nn)
		 
		Centers <- mvrnorm(5, mu=rep(0,15),Sigma = (Para^2)*diag(15))

		U1 <- matrix(Centers[1,],nn[1], 15, byrow= T) + matrix(rt(15*nn[1],df=Df), ncol=15)
		U2 <- matrix(Centers[2,],nn[2], 15, byrow= T) + matrix(rt(15*nn[2],df=Df), ncol=15)
		U3 <- matrix(Centers[3,],nn[3], 15, byrow= T) + matrix(rt(15*nn[3],df=Df), ncol=15)
		U4 <- matrix(Centers[4,],nn[4], 15, byrow= T) + matrix(rt(15*nn[4],df=Df), ncol=15)
		U5 <- matrix(Centers[5,],nn[5], 15, byrow= T) + matrix(rt(15*nn[5],df=Df), ncol=15)
		##----------------------------------------------
		X.data = rbind(U1,U2,U3,U4,U5)

		ss=dist(rbind(X.data,Centers))
		d=matrix(0,nrow=nrow(X.data)+5,ncol=nrow(X.data)+5)
 		d[row(d)>col(d)]=ss
 		DD=d[(nrow(X.data)+1):(nrow(X.data)+5),1:nrow(X.data)]
		
		for(i in 1:ncol(DD)){ DD[cl[i],i]=DD[cl[i],i]+del}
		ff=apply(DD,2,which.min)

		if(sum(ff==cl)==nrow(X.data)){
 			condition <- FALSE
 		}
	}#end while loop
	ID <- sample(1:nrow(X.data),nrow(X.data),replace = FALSE)
    	result <- list(Data = X.data[ID,], Number = N)
    	return(result)
}
###====================================================================================
NumA_heavytail <- function(P = 10, df, obs.num){
	set.seed(1)
	NN <- rep(NA,50)
	for(W in 1:50){
		Med <- Heavytail(Df = df, Para = P, Obs.num = obs.num)
		NN[W]<- Med$Number
	}
	return(mean(NN)-2)
}
## NumA_heavytail(P = 11, df = 2, obs.num = 80)
##=================================================================
####=----------------------------------------------------
## Auto-select the parameter Para in above function
Para_heavytail <- function(Df = 2, Obs.num = 80, Initial = 12){
	N.init <- NumA_heavytail(P = Initial, df = Df, obs.num = Obs.num)
	if( -0.2 <= N.init & N.init <= 0.2){
		return(Initial)
	}else{
		if(N.init > 0.2){
			return("Please use larger initial value")
		}else{
			N.old <- N.init
			P.old <- Initial
			GAP = 0.05
			Cond = TRUE
			while(Cond){
				P.new = P.old - GAP
				N.new <- NumA_heavytail(P = P.new, df = Df, obs.num = Obs.num)
				if(-0.2 <= N.new & N.new <= 0.2){
					Cond <- FALSE
					return(P.new)
				}else{
					if(N.new > 0.2){
						result <- uniroot(NumA_heavytail, interval = c(P.new,P.old), tol = 0.01, df = Df, obs.num = Obs.num)$root
						Cond <- FALSE
						return(result)
					}else{
						if(N.new != N.old & GAP > 0.02){
							GAP = GAP/2
						}
						P.old = P.new
						N.old = N.new
					}
				}
			}##end while loop
		}
	}##end
}
## PP <- Para_heavytail(Df =11, Obs.num = 80, Initial = 2.1)
## NumA_heavytail(P = PP, df = 11, obs.num = 80)
## NumA_heavytail(P = 2.6, df = 6, obs.num = 80)
##==========================================================================

##compute the paramter for setting1 
#	Para6_10dim <- rep(NA,10)
#	for(i in 1:10){
#		Para6_10dim[i] = Para.hunt(Dim = 10, C.num = 6, Rho = i/10-0.1, Obs.num=100, Initial = 4)
#	}
##above code results following parameter set
Para6_10dim <- c(2.56,2.60,2.58,2.83,2.74,2.86,3.12,3.26,3.4,2.88)

setting1 <- function(c){
			MM<- Data_generate(Dim = 10, C.num = 6, Para = Para6_10dim[c], Rho = c/10-0.1, Obs.num = 100)
			list(x=MM$Data, ngroup = 6)
}

###### Para3_6dim = Para.hunt(Dim = 6, C.num = 3, Rho = 0, Obs.num=100, Initial = 5)
###### Para3_6dim = 2.9
## Para3_6dim = Para.hunt(Dim = 6, C.num = 3, Rho = 0, Obs.num=1000, Initial = 4)
Para3_6dim = 3.2375

setting2 <- function(c){
			MM <- Data_generate(Dim = 6, C.num = 3, Para = Para3_6dim, Rho = 0, Obs.num = 1000)
			if(c==1){
				Data.x =  MM$Data
			}else{
				Noise <- runif(nrow(MM$Data)*(c-1)*6)
			    Data.x <- cbind(MM$Data, matrix(Noise,ncol=6*(c-1)))
			}
			
			list(x=Data.x, ngroup = 3)
}

##compute the paramter for setting3 
#	Para8_centers <- rep(NA,10)
#	for(i in 1:10){
#		Para8_centers[i] = Para.hunt(Dim = i*10, C.num = 8, Rho = 0, Obs.num=100, Initial = 2.8)
#	}
##above code results following parameter set
Para8_centers <- c(2.7,1.68,1.3,1.1125,0.98,0.9,0.834,0.777,0.75,0.705)

setting3 <- function(c){
			MM<- Data_generate(Dim = 10*c, C.num = 8, Para = Para8_centers[c], Rho = 0, Obs.num = 100)
			list(x=MM$Data, ngroup = 8)
}

##compute the parameter for setting4
#	Para3_20dim <- rep(NA,10)
#	for(i in 1:10){
#		Para3_20dim[i] = Para_3centers(Ratio = (5*(i-1))^(i!=1), Obs.num = 60, Initial = 6)
#	}
Para3_20dim <- c(1.2375,2.2875,3.1250,3.7625,4.2625,4.7000,5.0875,5.4500,5.7750,6)

setting4 <- function(c){
			MM <- Data3centers(Ratio = (5*(c-1))^(c!=1), Para = Para3_20dim[c], Obs.num = 60)
			list(x=MM$Data, ngroup = 3)
}
##compute the parameter for setting5
#Para5_15dim <- rep(NA,10)
#for(i in 6:10){
#	Para5_15dim[i] = Para_heavytail(Df = i+1, Obs.num = 80, Initial = 12)
#}
# seeds taken from random.org, uniform { 1, ..., 10000000 }
Para5_15dim <-c(11.8,4.25,3.175,2.5875,2.3875,2.35,2.1875,2.1125,2.025,2.05)

setting5 <- function(c){
			MM <- Heavytail(Df = 1+c, Para = Para5_15dim[c], Obs.num = 80)
			list(x=MM$Data, ngroup = 5)
}

settings <- list("setting1" = list(simulate = setting1, seed = 4244250), 
                 "setting2" = list(simulate = setting2, seed = 5513442),
                 "setting3" = list(simulate = setting3, seed = 5685887),
                 "setting4" = list(simulate = setting4, seed = 1243061),
                 "setting5" = list(simulate = setting5, seed = 7997590))
                #  "setting6" = list(simulate = setting6, seed = 5427086),
                # "setting7" = list(simulate = setting7, seed = 1595473))


