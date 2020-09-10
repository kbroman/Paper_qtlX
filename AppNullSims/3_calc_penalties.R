# calculate penalties from the null simulation results

library(qtl)

# load the data
b6btbr <- readRDS("../Application/b6btbr.rds")

# find the chromosome lengths
map <- pull.map(b6btbr)
L <- sapply(map, function(a) diff(range(a)))
L_A <- sum(L[-20])
L_X <- L[20]
L_AA <- L_A^2/2
L_XX <- L_X^2/2
L_AX <- L_A*L_X

# alpha values
alpha <- 0.05
rel_pair <- c(L_AA, L_AX, L_XX)
rel_pair <- rel_pair/sum(rel_pair)
alpha_pair <- setNames(1 - (1-alpha)^rel_pair, c("AA", "AX", "XX"))

rel_one <- c(L_A, L_X)
rel_one <- rel_one/sum(rel_one)
alpha_one <- setNames(1 - (1-alpha)^rel_one, c("A", "X"))

# load the null simulation results
nullsims <- readRDS("sim_results_2019-04-25.rds")

# penalties with A and X treated separately:
# A, X for one QTL; AA_heavy, AA_light, AX_heavy, XX_heavy
pen_AneX <- c(mainA=quantile(nullsims$AA[,"one"], 1-alpha_one["A"]),
              mainX=quantile(nullsims$XX[,"one"], 1-alpha_one["X"]),
              intH=quantile(nullsims$AA[,"int"], 1-alpha_pair["AA"]),
              intL=NA,
              intAX=quantile(nullsims$AX[,"int"], 1-alpha_pair["AX"]),
              intXX=quantile(nullsims$XX[,"int"], 1-alpha_pair["XX"]))
names(pen_AneX) <- c("mainA", "mainX", "intH", "intL", "intAX", "intXX")
pen_AneX["intL"] <- quantile(nullsims$AA[,"fv1"], 1-alpha_pair["AA"]) - pen_AneX["mainA"]

# penalties with A and X treated the same
pen_AeqX <- c(main=quantile(nullsims$all[,"one"], 1-alpha),
              intH=quantile(nullsims$all[,"int"], 1-alpha),
              intL=NA)
names(pen_AeqX) <- c("main", "intH", "intL")
pen_AeqX["intL"] <- quantile(nullsims$all[,"fv1"], 1-alpha) - pen_AeqX["main"]


save(bookpen_AneX, bookpen_AeqX,
     truepen_AneX, truepen_AeqX,
     file="penalties_2019-04-25.RData")
