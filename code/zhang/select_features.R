library(Signac)

input_files <- snakemake@input

for (i in seq_along(input_files)) {
    s <- read.table(input_files[[i]], sep = "\t", header = TRUE)
    s <- s[order(s$Row), ]
    s$mean <- s$Mean
    s$variance <- s$Variance
    s <- FitMeanVar(s)
    if (i == 1) {
        ranks <- s$rank
    } else {
        ranks <- ranks + s$rank
    }
}

features <- s$Row[head(order(ranks, decreasing=FALSE), 20000)]
writeLines(as.character(features), con = snakemake@output[['features']])