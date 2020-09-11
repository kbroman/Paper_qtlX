# look at some traits

library(qtl)

b6btbr <- readRDS("b6btbr.rds")

# find the chromosome lengths
map <- pull.map(b6btbr)
L_vec <- sapply(map, function(a) diff(range(a)))
L_A <- sum(L_vec[-20])
L_X <- L_vec[20]
L <- L_A+L_X
L_AA <- L_A^2
L_XX <- L_X^2
L_AX <- L_A*L_X


b6btbr <- calc.genoprob(b6btbr, error.prob=0.002, map.function="c-f",
                        stepwidth="max", step=1)

# load liver phenotypes
liver <- data.table::fread("Clean/liver_mlratio_clean.csv", data.table=FALSE, header=TRUE)
rownames(liver) <- liver[,1]
liver <- liver[,-1]

m <- match(rownames(liver), b6btbr$pheno[,1])
liver <- liver[!is.na(m),]
b6btbr_liver <- subset(b6btbr, ind=match(rownames(liver), b6btbr$pheno[,1]))

b6btbr_liver$pheno[,1] <- liver[,"503520"]
out <- scantwo(b6btbr_liver, chr=c("8", "X"), method="hk")
plotPXG(b6btbr_liver, marker=find.marker(b6btbr_liver, c("X", "8"), c(39.5, 33.9)))

b6btbr_liver$pheno[,1] <- liver[,"10003837305"]
out <- scantwo(b6btbr_liver, chr=c("1", "7", "X"), method="hk")
# chr X has a much bigger effect in males
# chr 7 has a much bigger effect in females
plotPXG(b6btbr_liver, marker=find.marker(b6btbr_liver, c("X", "7"), c(54.4, 2.51)))
plotPXG(b6btbr_liver, marker=find.marker(b6btbr_liver, c("X", "1"), c(54.4, 86.7)))

####################################################################################3
# start exploring

# Read in the simulation result from Karl
b6btbr_result <- readRDS("../AppNullSims/sim_results_2019-04-25.rds")

n_whole=L*(L-1)/2
n_XX=L_X*(L_X-1)/2
n_AA=L_A*(L_A-1)/2
n_AX=L_A*L_X
# n.self=sum((chrlen(map10[1:19])/10+1)*(chrlen(map10[1:19])/10)/2)
# n.pair=n.AA-n.self

# Partition alpha based on number of test in each region (area argument)
alpha=0.05
alphaAA=1-(1-alpha)^((n_AA-1)/(n_whole-1))
alphaXX=1-(1-alpha)^((n_XX-1)/(n_whole-1))
alphaAX=1-(1-alpha)^((n_AX-1)/(n_whole-1))
# alphaself=1-(1-alpha)^((n.self-1)/(n.whole-1))
# alphapair=1-(1-alpha)^((n.pair-1)/(n.whole-1))
alphaA=1-(1-alpha)^(L_A/(L_A+L_X))
alphaX=1-(1-alpha)^(L_X/(L_A+L_X))

########################
# Calculate threshold XX region
thresX <- quantile(b6btbr_result$XX[,"one"],probs=1-alphaX) #one
thresXXH <- quantile(b6btbr_result$XX[,"int"],probs=1-alphaXX) #int
thresXXHtrue <- quantile(b6btbr_result$XX[,"trueint"],probs=1-alphaXX) #int
thresXXL <- quantile(b6btbr_result$XX[,"fv1"],probs=1-alphaXX)-thresX #fv1
thresXX.full <- quantile(b6btbr_result$XX[,"full"],probs=1-alphaXX) #full
thresXX.VL <- thresXX.full-2*thresX

###########################
# Calculate threshold AA region
thresA <- quantile(b6btbr_result$AA[,"one"],probs=1-alphaA)
thresH <- quantile(b6btbr_result$AA[,"int"],probs=1-alphaAA)
#thresH.book <- quantile(max.bookint,probs=1-alphaAA)
thresHtrue <- quantile(b6btbr_result$AA[,"trueint"],probs=1-alphaAA)
thresL <- quantile(b6btbr_result$AA[,"fv1"],probs=1-alphaAA)-thresA
#thresL.book <- quantile(max.bookfv1,probs=1-alphaAA)-thresA
#thresL.new <- quantile(max.newfv1,probs=1-alphaAA)-thresA
thresAA.full <- quantile(b6btbr_result$AA[,"full"],probs=1-alphaAA)
thresVL <- thresAA.full-2*thresA

###########################
# Calculate threshold AX region
thresAXH <- quantile(b6btbr_result$AX[,"int"],probs=1-alphaAX)
thresAXHtrue <- quantile(b6btbr_result$AX[,"trueint"],probs=1-alphaAX)
thresAXL <- quantile(b6btbr_result$AX[,"fv1"],probs=1-alphaAX)-max(thresX, thresA)
thresAX.full <- quantile(b6btbr_result$AX[,"full"],probs=1-alphaAX)
thresAX.VL <- thresAX.full-thresX-thresA
penalties = c(thresA, thresX, thresH, thresL, thresAXH, thresXXH)
truepenalties = c(thresA, thresX, thresHtrue, thresL, thresAXHtrue, thresXXHtrue)
##########################
# Stepwise with 6 penalties using new penalty with sex as covariate
##############
load("../AppNullSims//penalties_2019-04-25.RData")
# 1) liver   probe 10003837305 with QTL on chr 1, 7, X on current QTL package
b6btbr_liver$pheno[,1] <- liver[,"10003837305"]
out_10003837305 <- stepwiseqtl(b6btbr_liver, method="hk", penalties = pen_AneX)
out_10003837305_curr <- stepwiseqtl(b6btbr_liver, method="hk", penalties = pen_AeqX)
#out_10003837305_true_X <- stepwiseqtl(b6btbr_liver, method="hk", penalties = truepen_AneX)
#out_10003837305_true <- stepwiseqtl(b6btbr_liver, method="hk", penalties = truepen_AeqX)
#out_10003837305_manual_pen <- stepwiseqtl(b6btbr_liver, method="hk", penalties = penalties)
#out_10003837305_manual_pen_true <- stepwiseqtl(b6btbr_liver, method="hk", penalties = truepenalties)
# Result has one qtl on 1, 9, 10, 13, 17, X, and two on 7, there is one interaction between 7 and X
# Compare to current qtl package there is no more interaction between X and 1 but the location of the
# interaction between 7 and X is the same. There is much more on autosome
summary(fitqtl(b6btbr_liver,qtl=out_10003837305_curr, formula=attr(out_10003837305_curr, "formula"), method = "hk"))
summary(fitqtl(b6btbr_liver,qtl=out_10003837305, formula=attr(out_10003837305, "formula"), method = "hk"))
plotPXG(b6btbr_liver, marker=find.marker(b6btbr_liver, c("X", "7"), c(57.4, 2.5)))
plotPXG(b6btbr_liver, marker=find.marker(b6btbr_liver, c("10", "7"), c(24.0, 18.6)))

###############################
# 2) liver   probe 503520      with QTL on chr 8, X on current QTL package
b6btbr_liver$pheno[,1] <- liver[,"503520"]
# out_503520_manual_pen <- stepwiseqtl(b6btbr_liver, method="hk", penalties = penalties)
# out_503520_manual_pen_true <- stepwiseqtl(b6btbr_liver, method="hk", penalties = truepenalties)
out_503520 <- stepwiseqtl(b6btbr_liver, method="hk", penalties = pen_AneX)
out_503520_curr <- stepwiseqtl(b6btbr_liver, method="hk", penalties = pen_AeqX)
out_503520_true_X <- stepwiseqtl(b6btbr_liver, method="hk", penalties = truepen_AneX)
# Result are similar to current qtl package with qtl on chr 8, X with location change
plotPXG(b6btbr_liver, marker=find.marker(b6btbr_liver, c("X", "8"), c(40.5, 45)))

###############################
# 3) # liver   probe 10003839327 with QTL on chr 7, X  on current QTL package
b6btbr_liver$pheno[,1] <- liver[,"10003839327"]
#out_10003839327_manual_pen <- stepwiseqtl(b6btbr_liver, method="hk", penalties = penalties)
#out_10003839327_manual_pen_true <- stepwiseqtl(b6btbr_liver, method="hk", penalties = truepenalties)
out_10003839327 <- stepwiseqtl(b6btbr_liver, method="hk", penalties = pen_AneX)
out_10003839327_curr <- stepwiseqtl(b6btbr_liver, method="hk", penalties = pen_AeqX)
out_10003839327_true_X <- stepwiseqtl(b6btbr_liver, method="hk", penalties = truepen_AneX)

# Result are significant different than current qtl package with more qtl on 3, 4, 8, 11
# There is two more interaction between qtl in same chromosome, on chr 4 and chr 11. Is that suspicious?
plotPXG(b6btbr_liver, marker=find.marker(b6btbr_liver, c("4", "4"), c(40.7, 41.6)))
plotPXG(b6btbr_liver, marker=find.marker(b6btbr_liver, c("11", "11"), c(40.6, 56.1)))


# Need to compare the different pLOD of these two models using same penalty
# There is only in 3 genes that the new truepen is different than the AneX
# However there are 524 genes that has AeqX and AneqX are different
# 5% is inline with the simulation however we not cover all the type one error since we only look at the probe that is significant in at least one qtl
# If we run this in all QTL we might discover more
save(out_10003837305, out_10003837305_curr, out_503520, out_503520_curr, out_10003839327, out_10003839327_curr,
     file="../R/analysis.RData")
