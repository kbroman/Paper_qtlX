# run scanone with liver; just keep track of the maximum LOD on each chromosome

library(qtl)

b6btbr <- readRDS("b6btbr.rds")

b6btbr <- calc.genoprob(b6btbr, error.prob=0.002, map.function="c-f",
                        stepwidth="max", step=1)

# load liver phenotypes
liver <- data.table::fread("Clean/liver_mlratio_clean.csv", data.table=FALSE, header=TRUE)
rownames(liver) <- liver[,1]
liver <- liver[,-1]

# convert phenotypes to normal quantiles
liver <- apply(liver, 2, nqrank)

m <- match(rownames(liver), b6btbr$pheno[,1])
liver <- liver[!is.na(m),]
b6btbr_liver <- subset(b6btbr, ind=match(rownames(liver), b6btbr$pheno[,1]))

# focus on the 37,827 traits with known QTL location
annot <- read.csv("Clean/microarray_annot.csv")
known_location <- as.character(annot$a_gene_id)[!is.na(annot$chr) & !is.na(annot$pos.Mb)]
liver <- liver[,known_location]

b6btbr_liver$pheno <- cbind(liver, b6btbr_liver$pheno)
sexcovar <- (b6btbr_liver$pheno$Sex=="Male")*1
out_liver <- scanone(b6btbr_liver, pheno.col=1:ncol(liver), method="hk", addcovar=sexcovar)

# maximum LOD score on each chromosome
maxlod <- t( apply(out_liver[,-(1:2)], 2, tapply, out_liver[,1], max, na.rm=TRUE) )
saveRDS(maxlod, file="maxlod_scanone_liver.rds")

# determine the set of traits with a significant QTL
load("../AppNullSims/penalties_2019-04-25.RData")

has_qtl <- rownames(maxlod)[(apply(maxlod[,1:19], 1, max) > pen_AneX[1] | maxlod[,20] > pen_AneX[2])]
saveRDS(has_qtl, file="liver_hasqtl.rds") # 10,814 traits
