library("PropCIs")
library(ggplot2)

nclusters <- list("setting1" = 6L,
		     "setting2" = 3L,
			 "setting3" = 8L,
			 "setting4" = 3L,
			 "setting5" = 5L)

Scale_X <-  list("setting1" = seq(0.0,0.9,0.1),
		     "setting2" = 6*seq(0,9,1),
			 "setting3" = seq(10,100,10),
			 "setting4" = c(1,seq(5,45,5)),
			 "setting5" = seq(2,11,1))

Label_X <- list("setting1" = "Rho",
		    "setting2" = "Noise Dimension",
			"setting3" = "Data Dimension",
			"setting4" = "Max Variance Ratio",
			"setting5" = "Degree of freedom")


for (s in list.dirs(full.names=FALSE, recursive=FALSE)){
	## remove the object so the if cond checkd out for next iter
	rm(Result_data)
	rm(Fit_result)

	Result.path <- file.path(s,"result.rds")

	if(!file.exists(Result.path)){
		Fit.result <- list()

		for (m in names(methods)) {
	 		Matrix_result <- matrix(NA,100,10)
			for(c in 1:10){
				f.replicates <- file.path(s, paste0("C" ,c), "method",paste0(m ,".rds"))
				if(file.exists(f.replicates)){
					Matrix_result[,c]<- readRDS(f.replicates)
				}
			}
			Fit.result[[m]] = Matrix_result
		}
		saveRDS(Fit.result, Result.path)
	}
	####----------------------Read from RDS and make plot-----------------

	Fit_result <- readRDS(Result.path)
	
	for(i in 1:length(Fit_result)){
		Proportion <- as.numeric(apply(Fit_result[[i]], 2, function(x) sum(x==nclusters[[s]],na.rm=T)))
		X_axis <- Scale_X[[s]]
		Methods <- rep(names(Fit_result)[i],10)

		if(!exists("Result_data")){
			Result_data <- as.data.frame(cbind(X_axis,Proportion,Methods))
		}else{
			Result_data <- rbind(Result_data,as.data.frame(cbind(X_axis,Proportion,Methods)))
		}	
	}
	Result_data[,1] <- as.numeric(as.character(Result_data[,1]))
	Result_data[,2] <- as.numeric(as.character(Result_data[,2]))
	
	Result_data$lower <- sapply(Result_data$Proportion,function(x) 100*scoreci(x,100,0.95)$conf.int[1])
	Result_data$upper <- sapply(Result_data$Proportion,function(x) 100*scoreci(x,100,0.95)$conf.int[2])
	###==================================================================
	if(s == "setting5"){
		ggplot(Result_data, aes(X_axis,Proportion, group = Methods, shape = Methods, colour = Methods))+labs(x = Label_X[[s]])+ geom_line(size = 1.3) + scale_x_reverse()
		ggsave(file.path(s,"Overlay.pdf"),width = 11, height = 7)

	
		ggplot(Result_data, aes(X_axis,Proportion, group = Methods, colour = Methods))+ geom_line(size = 1.3)+labs(x = Label_X[[s]])+facet_wrap(~Methods,nrow=2)+guides(colour=FALSE)+geom_errorbar(aes(ymax = upper, ymin=lower),colour = "black")+scale_x_reverse()
		ggsave(file.path(s,"Facet.pdf"),width = 11, height = 7)
	}else{
		ggplot(Result_data, aes(X_axis,Proportion, group = Methods, shape = Methods, colour = Methods))+labs(x = Label_X[[s]])+ geom_line(size = 1.3)
		ggsave(file.path(s,"Overlay.pdf"),width = 11, height = 7)

	
		ggplot(Result_data, aes(X_axis,Proportion, group = Methods, colour = Methods))+ geom_line(size = 1.3)+labs(x = Label_X[[s]])+facet_wrap(~Methods,nrow=2)+guides(colour=FALSE)+geom_errorbar(aes(ymax = upper, ymin=lower),colour = "black")
		ggsave(file.path(s,"Facet.pdf"),width = 11, height = 7)
	}
}



