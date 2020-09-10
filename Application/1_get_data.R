# load the libraries
library(qtl)

# grab the B6BTBR data from the mouse phenome database
url <- "https://phenomedoc.jax.org/QTL_Archive/attie_2015/Attie_2015_eqtl_clean.zip"
local_file <- basename(url)
if(!file.exists(local_file)) {
    download.file(url, local_file)
}

# unzip (into Clean/ subdir)
files <- unzip(local_file)

# read genotype data
b6btbr <- read.cross("csv", "Clean", "genotypes_clean.csv",
                     genotypes=c("BB", "BR", "RR"), alleles=c("B", "R"))

# drop chr "un"
b6btbr <- drop.markers(b6btbr, markernames(b6btbr, chr="un"))

# jitter map
set.seed(20190130)
b6btbr <- jittermap(b6btbr)

# save to file
saveRDS(b6btbr, "b6btbr.rds")
