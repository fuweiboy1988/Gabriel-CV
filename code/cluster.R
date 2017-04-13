# code/cluster.R


# default seed chosen by random.org
cluster_kmeans <- function(x, centers, seed = 2651513, iter.max = 100,
                           nstart = 100, ...)
{
    if (!is.na(seed)) {
        if ((exists0 <- exists(".Random.seed", envir=globalenv()))) {
            seed0 <- .Random.seed
        }
        set.seed(seed)
    }

    cl <- stats::kmeans(x, centers, iter.max=iter.max, nstart=nstart, ...)

    if (!is.na(seed)) {
        if (exists0) {
            assign(".Random.seed", seed0, envir=globalenv())
        } else if (exists(".Random.seed", envir=globalenv())) {
            rm(".Random.seed", envir=globalenv())
        }
    }

    cl
}
