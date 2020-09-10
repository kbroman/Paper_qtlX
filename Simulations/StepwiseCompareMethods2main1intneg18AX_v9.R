#library(qtl, lib.loc="Research/qtl12514")
library(qtl)
data(map10)
#source("~/Dropbox/QTL/calculate_threshold_v6.R")
source("~/Dropbox/QTL/import_threshold.R")
#source("~/Dropbox/QTL/stepwiseqtlQuocv11.R")
#pkgEnv = getNamespace("qtl")
#attach(pkgEnv)
# ## Manually code threshold because result from scantwo is not good.
# thresX <-  3.569351
# thresA <- 2.826747
# thresXX.VL <- thresXX.full-2*thresX
# thresVL <- thresAA.full -2*thresA
# thresAX.VL <- thresAX.full-thresX-thresA

# assignInNamespace("stepwiseqtl",stepwiseqtlQuoc, ns="qtl")
# Phenotype generating model specs
qtl.loc <- rbind(c(1,25,0), c(20,25,0))
qtl.model <- "y~Q1+Q2+Q1:Q2"
# Heritability
h=0.08 # 8%
add.effect=sqrt(4*h/(1-2*h))
effect.vec <- c(0.6,0.6,-1.8)

qtl.formula <- attr(terms(as.formula(qtl.model)), "factors")[-1,,drop=FALSE]
nterm <- apply(qtl.formula, 2, sum)
qtl.int <- qtl.formula[, nterm==2, drop=FALSE] # Get second order interactions
qtl.int <- apply(qtl.int, 2, function(a) which(a==1)) # Get the location of interaction

# penalties=c(thresA,thresX,thresL.new,thresL.new,thresAX,thresXX)

stepwiseqtlcompare <- function (repid, verbose=F) {
  print(repid)
#  pkgEnv = getNamespace("qtl")
#  attach(pkgEnv)
  x <- sim.cross(map10, type="bc", n.ind=250, model=qtl.loc)
  chr.type <- sapply(map10, function(a) ifelse(class(a)=="X","X","A"))
  class(x$geno$X) <- chr.type[20]
  main.mat <- x$qtlgeno-1 # main factors model matrix
  if (!is.null(ncol(qtl.int))) {
    int.mat <- matrix(nrow=nrow(main.mat), ncol=ncol(qtl.int),0) # interactions matrix
    for (j in 1:ncol(qtl.int)){
      int.mat[,j] <- main.mat[,qtl.int[1,j]]*main.mat[,qtl.int[2,j]]
    }
  } else int.mat <- NULL
  model.mat <- cbind(main.mat, int.mat)
  x$pheno <- x$pheno+model.mat%*%effect.vec
  x <- calc.genoprob(x, step=1)
  sex <- rep(0:1,each=125) # 0 for female
  pgm <- c(rep(0,62),rep(1,63),rep(0,63),rep(1,62)) # 0 is AA, ABf, 1 is BB ABr
#   phenotype <- c(rnorm(n=62, mean=-6, sd=1), rnorm(n=63, mean=6, sd=1), rnorm(n=125, mean=0, sd=1))
  x$pheno <- data.frame(x$pheno, sex, pgm)
  qtl.make <- makeqtl(cross=x, chr=names(map10)[qtl.loc[,1]], pos=qtl.loc[,2], what="prob")
  qtl.fit <- fitqtl(cross=x, qtl=qtl.make, formula=as.formula(qtl.model), method="hk", dropone=FALSE, get.ests=TRUE)
## (1) R/qtl-1.25.14
  penalties=c(main,intH,intL)
  res.AeqX <- stepwiseqtl(x, keeplodprofile=T, method="hk", penalties=penalties,max.qtl=max.qtl, verbose=verbose, keeptrace=T) 
  qtl.plod <- calc.plod(qtl.fit$lod, countqtlterms(as.formula(qtl.model), ignore.covar=TRUE), type="bc", penalties=penalties)
  attr(res.AeqX,"trueplod") <- qtl.plod
# #   ### (2) Backward step fixed, Light penalty for all interaction, Stop forward 3*main penalties below 0
# #   penalties=c(thresA,thresA,thresL.new,thresL.new,thresL.new,thresL.new)
# #   res.L.3Mfr0 <- stepwiseqtlQuoc(x, keeplodprofile=T, method="hk", penalties=penalties,max.qtl=max.qtl, verbose=F, keeptrace=T, stop.rule=1, k_f=k_f) 
# #   qtl.plod <- calc.plodQuoc(qtl.fit$lod, type="bc", formula=as.formula(qtl.model), qtl=qtl.make, penalties=penalties)
# #   attr(res.L.3Mfr0,"truelod") <- qtl.fit$lod
# #   attr(res.L.3Mfr0,"trueplod") <- qtl.plod
# #   ### (3) Backward step fixed, Light penalty for all interaction, Stop forward 3*main penalties below best plod
# #   penalties=c(thresA,thresA,thresL.new,thresL.new,thresL.new,thresL.new)
# #   res.L.3MfrBest <- stepwiseqtlQuoc(x, keeplodprofile=T, method="hk", penalties=penalties,max.qtl=max.qtl, verbose=F, keeptrace=T, stop.rule=2, k_f=k_f) 
# #   qtl.plod <- calc.plodQuoc(qtl.fit$lod, type="bc", formula=as.formula(qtl.model), qtl=qtl.make, penalties=penalties)
# #   attr(res.L.3MfrBest,"truelod") <- qtl.fit$lod
# #   attr(res.L.3MfrBest,"trueplod") <- qtl.plod
# #   ### (4) Backward step fixed, X threshold ,Light penalty for all interaction, Stop forward 3*main penalties below 0
# #   penalties=c(thresA,thresX,thresL.new,thresL.new,thresAX,thresXX)
# #   res.AvsX.3Mfr0 <- stepwiseqtlQuoc(x, keeplodprofile=T, method="hk", penalties=penalties,max.qtl=max.qtl, verbose=F, keeptrace=T, stop.rule=1, k_f=k_f) 
# #   qtl.plod <- calc.plodQuoc(qtl.fit$lod, type="bc", formula=as.formula(qtl.model), qtl=qtl.make, penalties=penalties)
# #   attr(res.AvsX.3Mfr0,"truelod") <- qtl.fit$lod
# #   attr(res.AvsX.3Mfr0,"trueplod") <- qtl.plod
#   ### (5) Backward step fixed, X threshold ,Light penalty for all interaction, Stop forward 3*main penalties below best plod
#   penalties=c(thresA,thresX,thresL.new,thresL.new,thresAX,thresXX)
#   res.AvsX.3MfrBest <- stepwiseqtlQuoc(x, keeplodprofile=T, method="hk", penalties=penalties,max.qtl=max.qtl, verbose=F, keeptrace=T, stop.rule=2, k_f=k_f) 
# #   qtl.plod <- calc.plodQuoc(qtl.fit$lod, type="bc", formula=as.formula(qtl.model), qtl=qtl.make, penalties=penalties)
# #   attr(res.AvsX.3MfrBest,"truelod") <- qtl.fit$lod
# #   attr(res.AvsX.3MfrBest,"trueplod") <- qtl.plod
# #   if (qtl.plod > attr(res.AvsX.3MfrBest,"pLOD")) attr(res.AvsX.3MfrBest,"x.missed") <- x
  
  ### (5) Backward step fixed, X threshold ,Light penalty for all interaction, Stop forward 3*main penalties below best plod
  penalties=c(mainA,mainX,intH,intL,intAX,intXX)
  res.AneqX <- stepwiseqtl(x, keeplodprofile=T, method="hk", penalties=penalties,max.qtl=max.qtl, verbose=verbose, keeptrace=T) 
  qtl.plod <- calc.plod(qtl.fit$lod, countqtlterms(as.formula(qtl.model), ignore.covar=TRUE), type="bc", penalties=penalties)
  attr(res.AneqX,"trueplod") <- qtl.plod
  if (qtl.plod > attr(res.AneqX,"pLOD")) attr(res.AneqX,"x.missed") <- x
  
  ### (6) X threshold , Very Light penalty for all interaction, Stop forward 3*main penalties below best plod
#   penalties=c(thresA,thresX,thresVL,thresVL,thresAX.VL,thresXX.VL)
#   res.AvsX.3MfrBest.v12514.VL <- stepwiseqtlQuocv12514(x, keeplodprofile=T, method="hk", penalties=penalties,max.qtl=max.qtl, verbose=F, keeptrace=T, stop.rule=2, k_f=k_f) 
  #   qtl.plod <- calc.plodQuoc(qtl.fit$lod, type="bc", formula=as.formula(qtl.model), qtl=qtl.make, penalties=penalties)
  #   attr(res.AvsX.3MfrBest.v12514,"truelod") <- qtl.fit$lod
  #   attr(res.AvsX.3MfrBest.v12514,"trueplod") <- qtl.plod
  #   if (qtl.plod > attr(res.AvsX.3MfrBest.v12514,"pLOD")) attr(res.AvsX.3MfrBest.v12514,"x.missed") <- x
  
  ### (7) Backward step fixed, X threshold ,Heavy and Light penalty 
  
  res <- list("AeqX"=res.AeqX, "AneqX"=res.AneqX)
#   pLOD.best <- max(attr(res.R1257,"pLOD"), attr(res.L.3Mfr0,"pLOD"), attr(res.L.3MfrBest,"pLOD"))
#   attr(res,"x.data") <- x
  attr(res,"truelod") <- qtl.fit$lod
#   
  res
#   print(res)
}

library(parallel)

# Refresh memory, clean all garbage. Alter mclapply.
# rm(list = ls())
cleanMem <- function(n=2) { for (i in 1:n) gc() }
#body(mclapply)[[8]][[3]][[3]][[4]] <- substitute(cleanMem())
for (i in 1:100) gc()

n.cores=detectCores()/2
n.sim <- 1000
max.qtl <- 10
k_f <- 3
results <- mclapply(1:n.sim, stepwiseqtlcompare, mc.cores=n.cores, mc.preschedule=FALSE)
desc <- paste("QTL 1.42-8, ", n.sim," runs with max.qtl=",max.qtl, ", compare all methods")
filename <- "~/Research/QTL/result/ResultsStepwiseCompareMethods2main1intneg18AXv9.RData"
save(results, desc, n.sim, max.qtl, qtl.loc, qtl.model, effect.vec, file=filename)
