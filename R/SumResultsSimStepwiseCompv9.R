# trueplod.dist <- rep(0,n.sim)
library(qtl)
methods <- vector("list", length(results[[1]]))
methods.names <- c("AeqX", "AneqX")
# methods <- as.list(methods)
n.method = length(results[[1]])
qtl.formula <- attr(terms(as.formula(qtl.model)), "factors")[-1,,drop=FALSE]
nterm <- apply(qtl.formula, 2, sum)
qtl.int <- qtl.formula[, nterm==2, drop=FALSE] # Get second order interactions
qtl.int <- apply(qtl.int, 2, function(a) which(a==1)) # Get the location of interaction
n.int <- ncol(qtl.int)
qtl.extra.mat = matrix(0, nrow=n.sim, ncol=n.method)
qtl.missed.mat = matrix(0, nrow=n.sim, ncol=n.method)
for (method.id in 1:n.method) {
  qtl.extra <- rep(0,n.sim)
  qtl.one <- rep(0,n.sim)
  int.extra <- rep(0,n.sim)
  qtl.extra.A <- rep(0,n.sim)
  qtl.extra.X <- rep(0,n.sim)
  qtl.missed <- rep(nrow(qtl.loc),n.sim)
  fdr <- rep(1,n.sim)
  case.to.look <- NULL
  for (sim.id in 1:n.sim){
    if (length(results[[sim.id]])<length(results[[1]])) next
    result.temp <- results[[sim.id]][[method.id]]
    if (!(class(result.temp)=="try-error")) { 
      #     trueplod.dist[i] <- attr(result.temp,"trueplod")
      #     if (trueplod.dist[i]>attr(result.temp,"pLOD")) { 
      #       print(i)
      #       print(summary(result.temp))
      #       print(trueplod.dist[i])
      #       case.to.look <- c(case.to.look, i)
      #     }
      if (!is.null(attr(result.temp,"formula"))) {
        #browser()
        temp.formula <- attr(terms(as.formula(attr(result.temp,"formula"))), "factors")[-1,,drop=FALSE]
        temp.nterm <- apply(temp.formula, 2, sum)
        temp.int <- temp.formula[, temp.nterm==2, drop=FALSE] # Get second order interactions
        temp.int <- apply(temp.int, 2, function(a) which(a==1)) # Get the location of interaction
        
          
        n.qtl <- result.temp$n.qtl
        qtlinmodel <- rep(0,n.qtl)
        for (j in 1:n.qtl){
          lodint.temp <- lodint(result.temp,qtl.index=j,drop=2,expandtomarkers=TRUE)
          if (as.character(lodint.temp[1,1])=='X') qtl.temp <- which(qtl.loc[,1]==20)
          else qtl.temp <- which(qtl.loc[,1]==lodint.temp[1,1])
          for (q in qtl.temp){            
              if ((lodint.temp[1,2] <= qtl.loc[q,2]) & (qtl.loc[q,2] <= lodint.temp[3,2])){ 
                qtlinmodel[j] <- q
              } 
          }          
        }
#         browser()
        intinmodel.flag <- NULL # in case there is no interaction in model
        if (!is.null(ncol(temp.int))){
          n.int.temp <- ncol(temp.int) 
          intinmodel.flag <- rep(FALSE,n.int.temp) # initiate the vector if int is in generating model
          if (!is.null(n.int)){
            for (col.id.temp in 1:n.int.temp){
              for (col.id.model in 1:n.int){
                # switch interaction which is in generating model to TRUE 
                if (all(sort(qtlinmodel[temp.int[,1]])==qtl.int[,col.id.model])) intinmodel.flag[col.id.temp] <- TRUE
              }
            }
            #             print(intinmodel)
#             browser()
          }
        }
        qtl.extra[sim.id] <- sum(qtlinmodel==0)
        qtl.one[sim.id] <- sum(any(qtlinmodel==1))
        int.extra[sim.id] <- sum(intinmodel.flag== FALSE)
        qtl.extra.A[sim.id] <- sum(qtlinmodel==0 & result.temp$chrtype=="A")
        qtl.extra.X[sim.id] <- sum(qtlinmodel==0 & result.temp$chrtype=="X")
        qtl.missed[sim.id] <- qtl.missed[sim.id] - length(unique(qtlinmodel[qtlinmodel!=0]))
        fdr[sim.id] <- qtl.extra[sim.id]/(qtl.extra[sim.id]+nrow(qtl.loc)-qtl.missed[sim.id])
      }
    }
  }
  exact=which(qtl.extra==0&qtl.missed==0)
  fdr.ave <-mean(fdr)
  if (!is.null(names(results[[sim.id]]))) {
    methods[[method.id]] <- names(results[[sim.id]])[method.id]
  } else methods[[method.id]] <- methods.names[method.id]
  attr(methods[[method.id]],"fdr.ave") <- fdr.ave
  attr(methods[[method.id]],"qtl.extra") <- table(qtl.extra)
  attr(methods[[method.id]],"qtl.one") <- table(qtl.one)
  attr(methods[[method.id]],"int.extra") <- table(int.extra)
  attr(methods[[method.id]],"qtl.extra.A") <- table(qtl.extra.A)
  attr(methods[[method.id]],"qtl.extra.X") <- table(qtl.extra.X)
  attr(methods[[method.id]],"qtl.extra.A.vec") <- qtl.extra.A
  attr(methods[[method.id]],"qtl.extra.X.vec") <- qtl.extra.X
  attr(methods[[method.id]],"int.extra.vec") <- int.extra
  attr(methods[[method.id]],"qtl.missed") <- table(qtl.missed)
  attr(methods[[method.id]],"nqtl.exact") <- length(exact)
  # attr(methods[[method.id]],"case.to.look") <- case.to.look
  # attr(methods[[method.id]],"fdr") <- fdr
  print(methods[[method.id]])
#   browser()
}
attr(methods,"qtl.loc") <- qtl.loc
attr(methods,"qtl.model") <- qtl.model
attr(methods,"effect.vec") <- effect.vec

# stem(trueplod.dist) 

# for (i in 1:n.sim){
#   #  print(attr(result.temp,"formula"))
#   if (!(class(result.temp)=="try-error")) {
#     #    trueplod.dist[i] <- attr(result.temp,"trueplod")
#     #     if (trueplod.dist[i]>attr(result.temp,"pLOD")) 
# { 
#   print(summary(result.temp))
#   print(trueplod.dist[i])
# }
#   }
# }

# cases=which(qtl.extra!=0&qtl.missed!=0)
# count=0
# for (i in cases) {
#   print(i)
#   print(results[[i]])
#   Sys.sleep(5)
#   if (results[[i]][[3]]$chr=="X") count=count+1
# }
