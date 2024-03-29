# This script generates the truth for the simulated data
# INPUT: Fasta file generated by simulation
# OUTPUT: Truth files

#!/usr/bin/Rscript
library(seqinr)
library(argparser, quietly=TRUE)
library(tidyverse)

p <- arg_parser("truth generator")
p <- add_argument(p, "--input", help="input file")
p <- add_argument(p, "--out", help="output file")

argv <- parse_args(p)

truth_het_homoeo <- function(input) {
	dat <- read.fasta(input)
	d <- matrix(unlist(dat), nrow = 4, byrow = TRUE)
	#het_A <- which(d[1, ] != d[2, ])
        #het_B <- which(d[3, ] != d[4, ])
	het <- which(d[1, ] != d[2, ] | d[3, ] != d[4, ]) 
	homeo <- which(d[1, ] != d[3, ] & d[2, ] != d[4, ])
	het <- as.list(het)
	name <- basename(input)
	name <- gsub(pattern = "\\.fsa$", "", name)
	#print(name)
        #write_lines(het_A, file=argv$out)
	#write_lines(het_B, file=argv$out_2)
	df <- data.frame(Genotype=rep(0,100000), Position=rep(0,100000), Truth_het=rep(0,100000), Truth_homeo=rep(0,100000))
	for (i in 1:nrow(df)) {
		df$Genotype[i] <- name
		df$Position[i] <- i
		if (i %in% het) {	
		df$Truth_het[i] <- 1	
		}
		if (i %in% homeo) {
		df$Truth_homeo[i] <- 1
		}	
	}
	print(df)
	write.table(df, file=argv$out, sep="\t", row.names = FALSE, col.names = TRUE)
} 

truth_het_homoeo(argv$input)




