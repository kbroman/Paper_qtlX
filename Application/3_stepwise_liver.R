# stepwise QTL analysis for liver expression traits

library(qtl)
library(broman)
errors2pushbullet()

b6btbr <- readRDS("b6btbr.rds")

b6btbr <- calc.genoprob(b6btbr, error.prob=0.002, map.function="c-f",
                        stepwidth="max", step=1)


liver <- data.table::fread("Clean/liver_mlratio_clean.csv", data.table=FALSE, header=TRUE)
rownames(liver) <- liver[,1]
liver <- liver[,-1]

# convert liver traits to normal quantiles
liver <- apply(liver, 2, nqrank)

m <- match(rownames(liver), b6btbr$pheno[,1])
liver <- liver[!is.na(m),]
b6btbr_liver <- subset(b6btbr, ind=match(rownames(liver), b6btbr$pheno[,1]))

# focus on the traits with a significant QTL by scanone()
has_qtl <- readRDS("liver_hasqtl.rds") # 9,188 traits
liver <- liver[,has_qtl]

# load penalties; use pen_AeqX, pen_AneX, truepen_AneX
load("../AppNullSims/penalties_2019-04-25.RData")

# which traits to use in this analysis?
traits <- seq(SUB, ncol(liver), by=96) # SUB takes values 1, 2, ..., 96

stepwise_output <- vector("list", length(traits))
names(stepwise_output) <- colnames(liver)[traits]

b6btbr_liver$pheno <- cbind(liver[,traits,drop=FALSE], b6btbr_liver$pheno)

penalties <- list(pen_AeqX, pen_AneX, truepen_AneX)

sexcovar <- data.frame(sex=as.numeric(b6btbr_liver$pheno$Sex=="Male")*1)

for(i in seq_along(traits)) {
    message(i, " of ", length(traits))

    stepwise_output[[i]] <- vector("list", 3)
    names(stepwise_output[[i]]) <- c("AeqX", "AneX", "AneX_truepen")

    not_missing <- !is.na(b6btbr_liver$pheno[,i])
    tmpcovar <- sexcovar[not_missing,,drop=FALSE]
    tmp <- subset(b6btbr_liver, ind=not_missing)

    for(j in 1:3) {
        message("    ", names(stepwise_output[[i]])[j])
        stepwise_output[[i]][[j]] <- stepwiseqtl(tmp, pheno.col=i, method="hk",
                                                 covar=tmpcovar,
                                                 penalties=penalties[[j]], verbose=FALSE)
    }
}

saveRDS(stepwise_output, "stepwise_liver_SUB.rds")

done()
