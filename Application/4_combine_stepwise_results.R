results <- NULL
for(i in 1:96) {
    j <- i
    if(i < 100) j <- paste0("0", j)
    if(i < 10) j <- paste0("0", j)
    file <- paste0("stepwise_liver_", j, ".rds")
    if(!file.exists(file)) next

    message(i)

    x <- readRDS(file)
    for(j in seq_along(x)) {
        for(k in 1:3) {
            if(is.null(x[[j]][[k]]) || length(x[[j]][[k]])==0) {
                x[[j]][[k]] <- list(chr=character(0), pos=numeric(0), formula="")
            }
            x[[j]][[k]] <- list(chr=x[[j]][[k]]$chr,
                                pos=x[[j]][[k]]$pos,
                                formula=attr(x[[j]][[k]], "formula"))
        }
    }

    if(is.null(results)) results <- x
    else results <- c(results, x)
}

saveRDS(results, "stepwise_liver.rds")
