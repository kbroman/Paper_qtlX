library(qtl)
data(map10)
LA <- sum(chrlen(map10[-20]))
LX <- sum(chrlen(map10[20]))
models <- c("1main", "2main", "3main", "1mainX", "2mainX")
sumup <- vector("list")
for (model.id in models){
  filename <- paste("sim_result/ResultsStepwiseCompareMethods", model.id, "8percentv9.RData", sep="")
  load(filename)
  source("SumResultsSimStepwiseCompv9.R")
  sumup[[model.id]] <- methods
}
models <- c("2main1intneg18", "2main1intneg18AX")
for (model.id in models){
  filename <- paste("sim_result/ResultsStepwiseCompareMethods", model.id, "v9.RData", sep="")
  load(filename)
  source("SumResultsSimStepwiseCompv9.R")
  sumup[[model.id]] <- methods
}
filename <- "summariesv9.RData"
save(sumup, file=filename)
