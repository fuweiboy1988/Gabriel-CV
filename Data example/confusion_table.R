Tava <- read.csv("Tava.csv", header= T)
Tava2 <- Tava[,c(1,3)]
#Tava2[,2] <- as.character(lapply(Tava2[,2],function(x)substr(x,1,7)))
Tava2[,2] <- as.character(lapply(Tava2[,2],function(x)toupper(x)))
names(Tava2)[2]<- c("Gene")

Yeast <- read.table("yeast3000.txt")
Data <- Yeast
Data[,1] <- as.character(lapply(Data[,1],function(x)toupper(x)))

set.seed(1)
Kfit <- kmeans(Data[,2:16],5,nstart = 100)
Data$Gabriel <- Kfit$cluster
names(Data)[1] <- "Gene"

##merge two data
DATA <- cbind(Data,Tava2)

##Get 30X5 confusion matrix
Confusion <- matrix(NA,31,6)

for(i in 1:30){
	Subset <- DATA[DATA$Cluster==i,];
	for(j in 1:5){
		Confusion[i,j] = sum(Subset$Gabriel==j)
	}
	Confusion[i,6] = sum(Confusion[i,1:5])
}
Confusion[31,] <- colSums(Confusion[1:30,])
Confusion <- as.data.frame(Confusion)
colnames(Confusion) <- c("Cluster 1", "Cluster 2","Cluster 3","Cluster 4","Cluster 5","Total")

for(i in 1:30){
	rownames(Confusion)[i] = paste("Cluster ", i, sep="")
}
rownames(Confusion)[31] = "Total"

##convert to latex table 
options(xtable.floating = FALSE)
xtable(Confusion,digits=0)

write.csv(Confusion,"Confusion_matrix.csv")

## reduce to 7X5
reduced <- Confusion[c(1:4,7,8,14),]
Other <- as.numeric(unname(Confusion[31,] - colSums(reduced)))
Final <- rbind(reduced,Other,Confusion[31,])
xtable(Final,digits=0)
###=============split plots==========
DATA$Tava <- NA

for(i in 1:nrow(DATA)){
	if( DATA[i,"Cluster"] == 1 ){DATA[i,"Tava"] = "1"}
	if( DATA[i,"Cluster"] == 2 ){DATA[i,"Tava"] = "2"}
	if( DATA[i,"Cluster"] == 3 ){DATA[i,"Tava"] = "3"}
	if( DATA[i,"Cluster"] == 4 ){DATA[i,"Tava"] = "4"}
	if( DATA[i,"Cluster"] == 7 ){DATA[i,"Tava"] = "7"}
	if( DATA[i,"Cluster"] == 8 ){DATA[i,"Tava"] = "8"}
	if( DATA[i,"Cluster"] == 14 ){DATA[i,"Tava"] = "14"}
	if( !(DATA[i,"Cluster"] %in% c(1:4,7,8,14)) ){DATA[i,"Tava"] = "Other"}
}

##==================for top plot================================
##function SE
SE <- function(x){sd(x)/sqrt(length(x))}

for(i in 1:5){
 	Subset <- DATA[DATA$Gabriel==i,];
	Means = colMeans(Subset[,2:16]);
	Sds = apply(Subset[,2:16], 2, sd);
	Summary <- cbind(rep(unique(Subset$Gabriel),15),1:15,unname(Means),unname(Sds))
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

#marginal plot of class - plot on gabriel clusters
top <- ggplot(Sum_data, aes(X, Y))+geom_line(size = 1)+labs(title = "", x = "", y = "")+ facet_grid(~Cluster)+geom_errorbar(aes(ymax = Upper, ymin=Lower),colour = "black")

plot(top)
##==================for right plot================================
for(i in c("1","2","3","4","7","8","14","Other")){
 	Subset <- DATA[DATA$Tava==i,];
	Means = colMeans(Subset[,2:16]);
	Sds = apply(Subset[,2:16], 2, sd);
	Summary <- cbind(rep(unique(Subset$Tava),15),1:15,unname(Means),unname(Sds))
	if(i=="1"){
	  Sum_data2 = Summary
	}else{
	  Sum_data2 = rbind(Sum_data2,Summary) 
	}
}

Sum_data2 <- as.data.frame(Sum_data2)
names(Sum_data2) <- c("Cluster","X","Y","Sd")
Sum_data2$Y<- as.numeric(as.character(Sum_data2$Y))
Sum_data2$X<- as.numeric(as.character(Sum_data2$X))
Sum_data2$Sd<- as.numeric(as.character(Sum_data2$Sd))

Sum_data2$Upper = Sum_data2$Y+Sum_data2$Sd
Sum_data2$Lower = Sum_data2$Y-Sum_data2$Sd


#marginal plot of class - plot on Tava clusters
Sum_data2$Cluster <- factor(Sum_data2$Cluster, levels=c("1","2","3","4","7","8","14","Other"))
right <- ggplot(Sum_data2, aes(X, Y))+geom_line(size = 1)+labs(title = "", x = "", y = "")+facet_grid(Cluster~.)+geom_errorbar(aes(ymax = Upper, ymin=Lower),colour = "black")
plot(right)
##==================for main plot================================
count=0

for(i in 1:5){ 
	for(j in c("1","2","3","4","7","8","14","Other")){
		main.subset = DATA[ which(DATA$Gabriel==i & DATA$Tava==j),]
		if( nrow(main.subset)<20 ) next;
		Means = colMeans(main.subset[,2:16]);
		Sds = apply(main.subset[,2:16], 2, sd);
		Summary.main <- cbind(rep(unique(main.subset$Tava),15),rep(unique(main.subset$Gabriel),15),1:15,unname(Means),unname(Sds))
		if( count!=0 ){
			Sum_data3 = rbind(Sum_data3, Summary.main)
		}else{ Sum_data3 = Summary.main }
		count=1
	}
}
Sum_data3<- as.data.frame(Sum_data3)
names(Sum_data3) <- c("Tava","Gabriel","X","Y","Sd")
Sum_data3$Y<- as.numeric(as.character(Sum_data3$Y))
Sum_data3$X<- as.numeric(as.character(Sum_data3$X))
Sum_data3$Sd<- as.numeric(as.character(Sum_data3$Sd))

Sum_data3$Upper = Sum_data3$Y+Sum_data3$Sd
Sum_data3$Lower = Sum_data3$Y-Sum_data3$Sd

##main plot
Sum_data3$Tava <- factor(Sum_data3$Tava, levels=c("1","2","3","4","7","8","14","Other"))
main <- ggplot(Sum_data3, aes(X, Y))+geom_line(size = 1)+labs(title = "", x = "Time", y = "Expression level")+facet_grid(Tava~Gabriel)+geom_errorbar(aes(ymax = Upper, ymin=Lower),colour = "black")

##=============================================================
library(ggplot2)
library(gridExtra)

#placeholder plot - prints nothing at all
empty <- ggplot()+geom_point(aes(1,1), colour="white") +
     theme(                              
       plot.background = element_blank(), 
       panel.grid.major = element_blank(), 
       panel.grid.minor = element_blank(), 
       panel.border = element_blank(), 
       panel.background = element_blank(),
       axis.title.x = element_blank(),
       axis.title.y = element_blank(),
       axis.text.x = element_blank(),
       axis.text.y = element_blank(),
       axis.ticks = element_blank()
     )

##arrange all plots together

pdf("split_plots.pdf",width = 10, height = 10)
grid.arrange( main, right, top,  empty,  ncol=2, nrow=2, widths=c(3, 1), heights=c(3, 1))
dev.off()


###




