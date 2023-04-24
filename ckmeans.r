#!/usr/bin/env Rscript

# Usage: ckmeans.r file.dat K
# Read in a table, run Ckmeans on each row
# for each row print out the cluster assignment for each element, and the mean (center) of each cluster

library("Ckmeans.1d.dp")

args = commandArgs(trailingOnly=TRUE)
k <- strtoi(args[2])
file <- args[1]
table <- read.table(file)
#table
for (row in 1:nrow(table)) {
 x <- table[row,]
 result <- Ckmeans.1d.dp(x, k)
 cat(result$cluster,"\n",sep="\t")
 rounded <- round(result$centers)
 cat(rounded, "\n",sep="\t")
}