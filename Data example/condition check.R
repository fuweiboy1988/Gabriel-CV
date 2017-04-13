
###-------------------function to retuern LHS and RHS---------------------
Cond_check <- function(data = DATA, Dim_x =1 ,Dim_y =1, K = 2){
	if( Dim_x+Dim_y != ncol(data) ){
		print("dimension doesn't match")
		return
	}
	LL <- list()
	Fit <- kmeans(data,K)
	ClusterID <- unique(Fit$cluster)
	
	for(i in ClusterID ){
		ID <- which(Fit$cluster == i)
		Data <- data[ID,]
		COV <- cov(Data)
		Sigma_xx <- COV[1:Dim_x,1:Dim_x]
		Sigma_xy <- COV[1:Dim_x,(Dim_x+1):(Dim_x+Dim_y)]
		Sigma_yx <- COV[(Dim_x+1):(Dim_x+Dim_y),1:Dim_x]
		Sigma_yy <- COV[(Dim_x+1):(Dim_x+Dim_y),(Dim_x+1):(Dim_x+Dim_y)]
		Eigen <- eigen(Sigma_yy)
		Eigen_value <- Eigen$values[1]
		Eigen_vec <- Eigen$vectors[,1]
		LHS = sqrt(Eigen_value)/2
		RHS = (t(Eigen_vec)%*%Sigma_yx%*%Sigma_xy%*%Eigen_vec)/( sqrt(t(Eigen_vec)%*%Sigma_yx%*%Sigma_xx%*%Sigma_xy%*%Eigen_vec) )
		LL <- cbind(LL,c(LHS,RHS))
	}
	return(LL)
}
##----------------------------------------------------------------------
##The Cond_check function returns a 2xK dimension matrix, where each column 
##is a pair value (LHS, RHS) of the formula in proposition 2, i.e. the first
##column shows the condition check result--(LHS, RHS) for the first cluster
##returned by k-means. If the first row is greater than the second row for 
##each column, that means the condition checks out for each cluster.
##The number of clusters using in k-means is K. 
##Dim_x is the dimension of the paritioned X and Dim_y is the dimension of 
##the paritioned Y.
##----------------------------------------------------------------------
Vote <- read.table("Voting.txt",sep=",")
Id <- vector()
for(i in 1:nrow(Vote)){
      if(sum(Vote[i,] == "?")!=0){
		Id <- c(Id,i)
	}
}

vote <- Vote[-Id,-1]

VV <- matrix(NA,232,16)

for(i in 1:232){
	for(j in 1:16){
		if(vote[i,j]=="n"){
			VV[i,j]=0
		}else if(vote[i,j]=="y"){
			VV[i,j]=1
		}
	}
}

set.seed(3)
ID <- sample(1:232, 232)
Vote <- VV[ID,]
###-----change Dim_x and Dim_y for different partition
###
Cond_check(data = Vote, Dim_x =8 ,Dim_y =8, K = 2)

###-----------------------------------------------------------------------------
Breast <- read.table("Breast cancer.txt",sep=",")
Breast <- Breast[,-1]
Id <- vector()

for(i in 1:nrow(Breast)){
      if(sum(Breast[i,] == "?")!=0){
		Id <- c(Id,i)
	}
}

Breast <- Breast[-Id,]

VV <- matrix(NA,nrow(Breast),10)

for(i in 1:nrow(Breast)){
	for(j in 1:10){
		VV[i,j] = Breast[i,j]		
	}
}

set.seed(3)
ID <- sample(1:683, 683)
Breast <- VV[ID,1:9]

###-----change Dim_x and Dim_y for different partition
###
Cond_check(data = Breast, Dim_x =4 ,Dim_y =5, K = 3)

###-----------------------------------------------------------------------------

Data <- read.table("sonar.txt",sep=",")

MM <- matrix(NA,208,61)
for(i in 1:208){
	for(j in 1:61){
		MM[i,j] = Data[i,j]
	}
}

set.seed(3)
Id <- sample(1:208,208)
Sonar <- MM[Id,1:60]

###-----change Dim_x and Dim_y for different partition
#condition fails for the second cluster (column) with equal partition
Cond_check(data = Sonar, Dim_x =30 ,Dim_y =30, K = 2)

#condition holds both clusters (columns) with unequal partition
Cond_check(data = Sonar, Dim_x =20 ,Dim_y =40, K = 2)

















