library(qtl)
library(broman)
errors2pushbullet()

set.seed(45212815+SUB)

data(map10)

n_sim <- 1000

L_A <- sum(summary(map10)[1:19,"length"])
L_X <- ceiling(summary(map10)["X","length"])
L <- L_A + L_X
L_AA <- L_A^2/2
L_XX <- L_X^2/2
L_AX <- L_A*L_X

# number of sims for the A:A part
n_AA <- 11
# determine number of A:X and X:X parts
n_AX <- ceiling( n_AA * L_AA / L_AX ) #  17,819
n_XX <- ceiling( n_AA * L_AA / L_XX ) # 317,508

sim_results <- vector("list", 4)
names(sim_results) <- c("AA", "AX", "XX", "all")
n_sim <- c(AA=n_AA, AX=n_AX, XX=n_XX, all=n_AA)

for(i in seq_along(sim_results)) {
    sim_results[[i]] <- matrix(ncol=6, nrow=n_sim[i])
    colnames(sim_results[[i]]) <- c("full", "fv1", "int", "add", "av1", "one")
}

for(i in 1:n_XX) {
    x <- sim.cross(map10, type="bc", n.ind=250)
    sex <- x$pheno$sex <- rep(0:1, 125)
    x <- calc.genoprob(x, step=1)

    # scan the autosomes
    if(i <= n_AA) {
        out <- scantwo(x, addcovar=sex, method="hk", verbose=FALSE)

        # maximum over the whole genome (A and X)
        sim_results[["all"]][i,] <- subrousummaryscantwo(out, for.perm=TRUE)

        # maximum over just the autosomes
        sim_results[["AA"]][i,] <- subrousummaryscantwo(subset(out, chr=1:19), for.perm=TRUE)
    }

    cat("    done with AA\n")

    # do a 2d scan on the A:X part
    if(i <= n_AX) {
        out <- scantwo(x, addcovar=sex, chr=list(1:19,"X"), method="hk", verbose=FALSE)
        sim_results[["AX"]][i,] <- subrousummaryscantwo(out, for.perm=TRUE)
    }

    cat("    done with AX\n")

    # do a genome scan on the X chromosome
    out <- scantwo(x, addcovar=sex, chr="X", method="hk", verbose=FALSE)
    sim_results[["XX"]][i,] <- subrousummaryscantwo(out, for.perm=TRUE)

    saveRDS(sim_results, "backcross_sim_results_SUB.rds")
}
