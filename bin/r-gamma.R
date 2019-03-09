args <- commandArgs(trailingOnly=TRUE)
cat(rgamma(args[1],shape=2.4, scale=114), sep="\n")
