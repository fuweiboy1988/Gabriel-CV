#!/usr/bin/Rscript --vanilla

library("MASS")
library("parallel")
source("settings.R") # settings


# Use the L'Ecuyer RNG so that we can have multiple streams of
# random numbers
rng_kind <- "L'Ecuyer-CMRG"
rng_normal_kind <- "Inversion"


nrep <- 100

for (s in names(settings)) {
    if (!dir.exists(s)) {
        dir.create(s)
    }
	for(c in 1:10){
		if (!dir.exists(file.path(s, paste0("C" ,c)))){
			dir.create(file.path(s, paste0("C" ,c)))
		}
			
		filename <- file.path(s, paste0("C" ,c), "replicates.rds")
		if (!file.exists(filename)) {

			cat("generating replicates for '", s, "' with Case ", c, "\n", sep="")

			# generate the list of seeds
			set.seed(settings[[s]]$seed, rng_kind, rng_normal_kind)
			seeds <- vector("list", nrep)
			seeds[[1L]] <- .Random.seed
			for (r in seq_len(nrep - 1L))
				seeds[[r+1L]] <- nextRNGStream(seeds[[r]])

			# generate the replicates
			replicates <- vector("list", nrep)

			pb <- txtProgressBar(0, nrep, style=3)
			for (r in seq_len(nrep)) {
				# set the seed for this replicate
				.Random.seed <<- seeds[[r]]
				replicates[[r]] <- settings[[s]]$simulate(c)

				# save the state of the RNG
				replicates[[r]][["rng"]] <- list(seed = .Random.seed,
												 kind = rng_kind,
												 normal_kind = rng_normal_kind)

				setTxtProgressBar(pb, r)
			}
			close(pb)

			saveRDS(replicates, filename)
		}
	}
}
