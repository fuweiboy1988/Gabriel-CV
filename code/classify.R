# code/classify.R
#
# Depends:
#   library("MASS")
#


classify_lda <- function(x, grouping, prior=c("proportions", "equal"))
{
    g <- as.factor(grouping)
    prior <- match.arg(prior)

    n <- length(grouping)
    ng <- nlevels(g)

    if (prior == "equal") {
        prior <- rep(1/ng, ng)
    } else {
        counts <- as.vector(table(g))
        prior <- counts / n
    }

    MASS::lda(x, g, prior)
}


classify_nearest <- function(x, grouping)
{
    g <- as.factor(grouping)
    nc <- tabulate(g)

    gmat <- model.matrix( ~ g - 1)
    sums <- t(gmat) %*% x
    rownames(sums) <- levels(g)
    means <- sums / nc

    structure(list(means=means, levels=levels(g)),
              class="classify_nearest")
}


predict.classify_nearest <- function(object, x, ...)
{
    means <- object$means
    levels <- object$levels

    d <- scale(x %*% t(means), center=0.5 * rowSums(means * means),
               scale=FALSE)
    class <- factor(levels[apply(d, 1, nnet::which.is.max)],
                    levels=levels)
    list(class=class)
}
