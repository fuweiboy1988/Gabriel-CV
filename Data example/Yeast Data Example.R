
Data <- read.table("Yeast/cdc15.txt",header=T,sep="\t")
head(Data)
DATA = Data[,52:68]

##impute NAs by the column mean
Test <- as.matrix(DATA)
for(i in 1:ncol(Test)){
 	IDs <- which(is.na(Test[,i]));
	Impute <- median(Test[,i],na.rm=T);
	Test[,i][IDs] <- Impute;
}
##impute NAs by K-NN with k =5
##library("DMwR")
##Test <- as.matrix(DATA)
##Test <- knnImputation(Test,k=15)

##select the top 3000/2945 most variable genes(s.d/mean)
Rank <- rep(NA,nrow(Test))

for(i in 1:nrow(Test)){
	Rank[i] = sd(Test[i,])/mean(Test[i,])
}

Order <- order(Rank,decreasing = T)
Id <- Order[1:2945]
Test <- Test[Id,]
Gene <- as.character(Data$X)[Id]
Gene <- as.character(lapply(Gene,function(x)substr(x,1,7)))

##variance normalize each gene
Test <- Test[,-c(10,11)]

for(i in 1:nrow(Test)){
	Mean = mean(Test[i,]);
	SD = sd(Test[i,]);
	for(j in 1:ncol(Test)){
		Test[i,j] = (Test[i,j]-Mean)/SD
	}
}
############## abandon above ########################################
##Use processed data
Yeast <- read.table("yeast3000.txt")
Test <- as.matrix(Yeast[,2:16])
#######################################
set.seed(1)
K.original <- cv.kmeans.gabriel(Test, 5, 2, maxcenters=40, classify.method="nearest")$centers
K.original

set.seed(0)
Test2 <- Uncorrelate2(Test,K.original)
cv.kmeans.gabriel(Test2, 5, 2, maxcenters=40, classify.method="nearest")$centers

###################################################################################						
###==============================================================
Data <- Yeast
names(Data)[1] = "Gene"
Data$Gene = as.character(Data$Gene)
Data$Gene <- as.character(lapply(Data$Gene,function(x)substr(x,1,7)))

set.seed(1)
Kfit <- kmeans(Test,5,nstart = 100)
Data$Cluster <- Kfit$cluster

Names_cluster <- c("Cluster 1; size=550","Cluster 2; size=590","Cluster 3; size=654","Cluster 4; size=634","Cluster 5; size=517")

for(i in 1:5){
 	Subset <- Data[Data$Cluster==i,2:16];
	print(nrow(Subset));
	Means = colMeans(Subset);
	Sds = apply(Subset, 2, sd);
	Summary <- cbind(rep(Names_cluster[i],15),1:15,unname(Means),unname(Sds))
	if(i==1){
	  Sum_data = Summary
	}else{
	  Sum_data = rbind(Sum_data,Summary) 
	}
}
Sum_data<- as.data.frame(Sum_data)
names(Sum_data) <- c("Cluster","X","Y","Sd")
Sum_data$Y<- as.numeric(as.character(Sum_data$Y))
Sum_data$X<- as.numeric(as.character(Sum_data$X))
Sum_data$Sd<- as.numeric(as.character(Sum_data$Sd))

Sum_data$Upper = Sum_data$Y+Sum_data$Sd
Sum_data$Lower = Sum_data$Y-Sum_data$Sd

library(ggplot2)
pdf("5_clusters.pdf",width=9)
ggplot(Sum_data, aes(X, Y, group = Cluster))+labs(title = "", x = "", y = "")+geom_line(size = 1.3)+facet_wrap(~Cluster,nrow=3)+geom_errorbar(aes(ymax = Upper, ymin=Lower),colour = "black")
dev.off()
##-----------------------------------------------------------
Types <- read.table("go_slim_mapping.txt",header=F,sep="\t")
Types <- Types[Types[,4]=="P",c(1,5)] ## biological process
names(Types)<- c("Gene","Process")
Types$Gene <- as.character(Types$Gene)
Types$Gene <- as.character(lapply(Types$Gene,function(x)substr(x,1,7)))
Types$Process <- as.character(Types$Process)

Final <- merge(Data,Types,by="Gene",all.x=TRUE)
write.csv(Final, "Final_Yeast.csv")
##-----------------------------
Final2 <- read.csv("Final_Yeast.csv")
Final2$Gene <-as.character(Final2$Gene)
Final2$Process <-as.character(Final2$Process )

library(plyr)
reference_table <- ddply(Final2,.(Process),summarise, Total = length(Gene))

SIZE <- c(550,590,654,634,517)

for(k in 1:5){
 	subset = Final2[Final2$Cluster==k,]
	print(paste("cluster",k))
	Process.table = ddply(subset,.(Process),summarise, Total = length(Gene))
	Process.table$Grand = NA;
	Process.table$P = NA;
	for(j in 1:nrow(Process.table)){
		z <- Process.table[j,]$Total
		m <-  reference_table[reference_table[,1]==Process.table[j,]$Process,2][1]
		Process.table[j,]$Grand = m
		Process.table[j,]$P = 1-phyper(z-1,m,2945-m,SIZE[k])
	}
	print(Process.table[Process.table$P <= 10^(-4),])
}







