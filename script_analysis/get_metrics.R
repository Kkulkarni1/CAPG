#!/usr/bin/Rscript
# @file get_metrics.R
# @author R. Kulkarni and K. S. Dorman
#
# Purpose: part of real peanut data analysis
# process info files and compute metrics
#
# Requirements: dplyr
 
library(dplyr)

capg.home <- '.'
info.dir <- paste(capg.home, '/data/peanut/info', sep='')	# directory with "info" files
write.outfiles <- T						# write the following output files
results.dir <- 'results/peanut'					# where to put the final results
het.outfile = paste(results.dir, "/CAPG_PL_peanut_het.txt", sep='')
							# output file for heterozygous genotype calling
ho.outfile = paste(results.dir, "/CAPG_PL_peanut_ho.txt", sep='')
							# output file for homoeologous and allelic SNP calling
genotype_list <- c("SRR4124062")			# there are up to 14 accessions, 1 here for testing
	#, "SRR4124066", "SRR4124068", "SRR4124074", "SRR4124078", "SRR8361734", "SRR8361735", "SRR8361736", "SRR8361737", "SRR8361738", "SRR8736998", "SRR8737008", "SRR8737061", "SRR8737062")
load.data.only <- F
do.het <- T

if (!dir.exists(results.dir))
	dir.create(results.dir)

d <- NULL
 
for (g in genotype_list) {
	cat("Reading genotype ", g, "\n")
	if (is.null(d)) {
		d <- read.table(sprintf("%s/%s_info.txt", info.dir, g), header = T, sep = "\t")
	} else {
		d.tmp <- read.table(sprintf("%s/%s_info.txt", info.dir, g), header = T, sep = "\t")
		d <- rbind(d, d.tmp)
	}
}

stopifnot(!load.data.only)

# log metrics
for (i in 11:25) {
	d[,i] <- log(d[,i])
}

# compute metric for heterozygosity
if (do.het) {
	d$PP_hetA <- d$PA.1. - apply(d[, c('PA.0.', 'PA.2.')], 1, max)
	d$PP_hetB <- d$PB.1. - apply(d[, c('PB.0.', 'PB.2.')], 1, max)
	if (write.outfiles) {
		print(paste("Writing", het.outfile, "\n"))
		write.table(d[, c('PositionA', 'PositionB', 'Genotype', 'ChromA', 'ChromB', 'PP_hetA', 'PP_hetB', 'CovA', 'CovB', 'ECTA', 'ETAA', 'ECTB', 'ETAB')], file = het.outfile, row.names = , F, sep = "\t", quote = F)
	}
}

d$Genotype_call <- gsub('/', '', d$Genotype_call)
d$max_ppG <- apply(d[, c('PP.0.0.', 'PP.0.1.', 'PP.0.2.', 'PP.1.0.', 'PP.1.1.', 'PP.1.2.', 'PP.2.0.', 'PP.2.1.', 'PP.2.2.')], 1, max)   # can also accomplish with pmax()
d$min_ppG <- apply(d[, c('PP.0.0.', 'PP.0.1.', 'PP.0.2.', 'PP.1.0.', 'PP.1.1.', 'PP.1.2.', 'PP.2.0.', 'PP.2.1.', 'PP.2.2.')], 1, min)

d$max_ppA <- apply(d[, c('PA.0.', 'PA.1.', 'PA.2.')], 1, max)
d$min_ppA <- apply(d[, c('PA.0.', 'PA.1.', 'PA.2.')], 1, min)
d$max_ppB <- apply(d[, c('PB.0.', 'PB.1.', 'PB.2.')], 1, max)
d$min_ppB <- apply(d[, c('PB.0.', 'PB.1.', 'PB.2.')], 1, min)

d.join <- d %>%
	group_by(PositionA, PositionB, ChromA, ChromB) %>%
	mutate(
		biallelic = length(unique(unlist(strsplit(Genotype_call, split="")))), # == 2,
		N1 = names(which.max(table(unlist(strsplit(Genotype_call, split=""))))), #unique(unlist(strsplit(Genotype_call, split="")))[1],
		N2 = if_else(biallelic >= 2, names(sort(table(unlist(strsplit(Genotype_call, split=""))), decreasing=T)[2]), if_else(N1 == 'A', 'C', 'A')),
		CovA = mean(CovA),
		CovB = mean(CovB),
	) %>%
	summarise(
		biallelic = unique(biallelic),
		PP_homeo = pmax(
			sum(if_else(Major.allele == N1 & Minor.allele == N2, PP.0.2.,   # N1N1/N2N2
				if_else(Major.allele == N2 & Minor.allele == N1, PP.2.0., min_ppG)) - max_ppG),
			sum(if_else(Major.allele == N2 & Minor.allele == N1, PP.0.2.,   # N2N2/N1N1
				if_else(Major.allele == N1 & Minor.allele == N2, PP.2.0., min_ppG)) - max_ppG)),
		PP_homoA = - pmax(sum(if_else(Major.allele == N1, PA.0.,        # N1N1
				if_else(Minor.allele == N1, PA.2., min_ppA)) - max_ppA),
			sum(if_else(Major.allele == N2, PA.0.,                  # N2N2
				if_else(Minor.allele == N2, PA.2., min_ppA)) - max_ppA)),
		PP_homoB = - pmax(sum(if_else(Major.allele == N1, PB.0.,        # N1N1
				if_else(Minor.allele == N1, PB.2., min_ppB)) - max_ppB),
			sum(if_else(Major.allele == N2, PB.0.,  		# N2N2
				if_else(Minor.allele == N2, PB.2., min_ppB)) - max_ppB)),
		CovA = unique(CovA),
		CovB = unique(CovB),
	)
if (write.outfiles) {
	print(paste("Writing", ho.outfile, "\n"))
	write.table(d.join, file = ho.outfile, row.names = F, sep = "\t", quote = F)
}
