#!/usr/bin/Rscript --vanilla

library("bcv")
library("cluster")
library("devtools")
library("e1071")
library("mclust")
library("MASS")
library("nnet")
library("parallel")
load_all("../../lib/fpc", export_all=FALSE)
load_all("../../lib/NbClust", export_all=FALSE)
source("../../code/classify.R")
source("../../code/cluster.R")
source("../../code/gabriel.R")
source("../../code/jump.R")
source("../../code/wold.R")
source("methods.R")


options(cl.cores=max(1, detectCores() - 1))
cl <- makeCluster(getOption("cl.cores", 2))

kmax <- 10

for (s in list.dirs(full.names=FALSE, recursive=FALSE)) {
	for(c in 1:10){
		f.replicates <- file.path(s, paste0("C" ,c), "replicates.rds")
		
		if (!file.exists(f.replicates))
			next

		replicates <- readRDS(f.replicates)
		nrep <- length(replicates)

		if (!dir.exists(file.path(s, paste0("C" ,c), "method")))
			dir.create(file.path(s, paste0("C" ,c), "method"))

		for (m in names(methods)) {
			f.method <- file.path(s, paste0("C" ,c), "method", paste0(m, ".rds"))

			if (!file.exists(f.method)) {
				cat("applying method '", m, "' to '", s, "' with Case ", c, "\n", sep="")

				method <- methods[[m]]
				nclusters <- integer(nrep)

				clusterExport(cl=cl, varlist=c("replicates", "method", "kmax",
											   "classify_lda",
											   "classify_nearest",
											   "cluster_kmeans",
											   "compute.jump",
											   "cv.kmeans.gabriel",
											   "gabriel_cor_correct", 
											   "Uncorrelate",
		 									   "Uncorrelate2",
											   "Rotate",
											   "Cov",
											   "Impute",
											   "jump",
											   "kmeans.rndstart",
											   "predict.classify_nearest",
											   "mclustBIC",
											   "WoldIter",
											   "Wold_holdout"))
				nclusters <- parSapply(cl, replicates, function(r) {
						rng <- r$rng
						set.seed(rng$seed, rng$kind, rng$normal_kind)
						tryCatch(method(r$x, kmax),
								 error=function(e) NA_integer_)
					})

				nclusters <- as.integer(nclusters)
				saveRDS(nclusters, f.method)
			}
		}
	}
}

stopCluster(cl)

