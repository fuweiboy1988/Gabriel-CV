

   library("bcv")
   library("devtools")
   library("cluster")
   library("e1071")
   library("mclust")
   library("MASS")
   library("nnet")
   load_all("../lib/fpc")
   load_all("../lib/NbClust")
   #source("../lib/NbClust.R")
   source("../code/classify.R")
   source("../code/cluster.R")
   source("../code/gabriel.R")
   source("../code/jump.R")
   source("../code/wold.R")


######------------Voting data-----------------------############
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

set.seed(1)
cv.kmeans.gabriel(Vote, 5, 2, maxcenters=10, classify.method="nearest")$centers

set.seed(1)		
gabriel_cor_correct(Vote, maxcenters=10, type = 2)

set.seed(1)		
Wold_holdout(data = Vote, CV = 5, Errortol = 0.01,  max.k =10)

set.seed(1)		
Gap <- clusGap(Vote, FUN = kmeans, K.max =10)
which(Gap[[1]][,3] == max(Gap[[1]][,3]))
			
set.seed(1)	
mcluster <- Mclust(Vote, G = 1:10)
 mcluster$G
						
set.seed(1)
Ch <- NbClust(Vote, min.nc = 2, max.nc = 10,method = "kmeans", index = "ch")
Ch$Best.nc[1]
	
set.seed(1)
Hartigan <- NbClust(Vote, min.nc = 2, max.nc = 10,method = "kmeans", index = "hartigan")
Hartigan$Best.nc[1]
			
set.seed(1)
Jump <- jump(Vote,10, plotjumps=FALSE, trace=FALSE)
Jump$maxjump
	
set.seed(1)
PS <- prediction.strength(Vote, Gmin=2, Gmax=10)
 PS$optimalk
			
set.seed(1)
SB <- nselectboot(Vote,clustermethod=kmeansCBI,classification="centroid",krange=2:10)
SB$kopt

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

set.seed(1)
cv.kmeans.gabriel(Breast, 5, 2, maxcenters=10, classify.method="nearest")$centers

set.seed(1)		
gabriel_cor_correct(Breast, maxcenters=10, type = 2)

set.seed(1)		
Wold_holdout(data = Breast, CV = 5, Errortol = 0.01,  max.k =10)

set.seed(1)		
Gap <- clusGap(Breast, FUN = kmeans, K.max =10)
which(Gap[[1]][,3] == max(Gap[[1]][,3]))
			
set.seed(1)	
mcluster <- Mclust(Breast, G = 1:10)
 mcluster$G
						
set.seed(1)
Ch <- NbClust(Breast, min.nc = 2, max.nc = 10,method = "kmeans", index = "ch")
Ch$Best.nc[1]
	
set.seed(1)
Hartigan <- NbClust(Breast, min.nc = 2, max.nc = 10,method = "kmeans", index = "hartigan")
Hartigan$Best.nc[1]
			
set.seed(1)
Jump <- jump(Breast,10, plotjumps=FALSE, trace=FALSE)
Jump$maxjump
	
set.seed(1)
PS <- prediction.strength(Breast, Gmin=2, Gmax=10)
 PS$optimalk
			
set.seed(1)
SB <- nselectboot(Breast,clustermethod=kmeansCBI,classification="centroid",krange=2:10)
SB$kopt

###-----------------------------------------------------------------------------
DATA <- read.table("promeroy.txt", header = T, sep="\t")

Data <- matrix(NA,1379,42)
for(i in 1:1379){
	for(j in 1:42){
		Data[i,j] = as.numeric(as.character(DATA[i+1,j+1]))
	}
}

DATA <- t(Data)

set.seed(1)
cv.kmeans.gabriel(DATA, 5, 2, maxcenters=10, classify.method="nearest")$centers

set.seed(1)		
gabriel_cor_correct(DATA, maxcenters=10, type = 2)

set.seed(1)		
Wold_holdout(data = DATA, CV = 5, Errortol = 0.01,  max.k =10)

set.seed(1)		
Gap <- clusGap(DATA, FUN = kmeans, K.max =10)
which(Gap[[1]][,3] == max(Gap[[1]][,3]))
			
set.seed(1)	
mcluster <- Mclust(DATA, G = 1:10)
 mcluster$G
						
set.seed(1)
Ch <- NbClust(DATA, min.nc = 2, max.nc = 10,method = "kmeans", index = "ch")
Ch$Best.nc[1]
	
set.seed(1)
Hartigan <- NbClust(DATA, min.nc = 2, max.nc = 10,method = "kmeans", index = "hartigan")
Hartigan$Best.nc[1]
			
set.seed(1)
Jump <- jump(DATA,10, plotjumps=FALSE, trace=FALSE)
Jump$maxjump
	
set.seed(1)
PS <- prediction.strength(DATA, Gmin=2, Gmax=10)
 PS$optimalk
			
set.seed(1)
SB <- nselectboot(DATA,clustermethod=kmeansCBI,classification="centroid",krange=2:10)
SB$kopt






























