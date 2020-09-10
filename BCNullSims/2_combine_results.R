# combine the simulation results
results <- NULL
for(i in 1:96) {
    j <- i
    if(i < 100) j <- paste0("0", j)
    if(i < 10) j <- paste0("0", j)
    file <- paste0("BackcrossSims/backcross_sim_results_", j, ".rds")
    if(!file.exists(file)) { next }
    res <- readRDS(file)

    if(is.null(results)) {
        results <- res
    } else {
        for(j in seq_along(results)) {
            results[[j]] <- rbind(results[[j]], res[[j]])
        }
    }
}

saveRDS(results, "backcross_sim_results_2020-08-10.rds")
