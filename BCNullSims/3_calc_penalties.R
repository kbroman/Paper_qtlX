# calculate penalties from the null simulation results

library(qtl)

data(map10)

# find the chromosome lengths
L_A <- sum(summary(map10)[1:19,"length"])
L_X <- ceiling(summary(map10)["X","length"])
L <- L_A + L_X
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
nullsims <- readRDS("backcross_sim_results_2020-08-10.rds")

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


save(pen_AneX, pen_AeqX,
     file="backcross_penalties_2020-08-10.RData")
