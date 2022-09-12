#!/usr/bin/Rscript
# @file get_metrics_simulation.R
# @author R. Kulkarni, K. S. Dorman
#
# Purpose: Compute metrics from CAPG "info" output for SNP calling.
#

library(dplyr)

write.results <- T			# write output files
parent.dir <- '.'			# OUTPUT: CAPG repository top level directory
sim.dir <- 'data/simulation'
out.dir <-				# OUTPUT: where the metric files go
	 paste(parent.dir, '/data/simulation/results', sep='')
homoeo_rates <- c(0.005, 0.007)		# homoeologous rates (CAPG manuscript: c(0.005, 0.007, 0.010))
covg_rates <- c(5, 10)			# coverage rate (CAPG manuscript: c(5, 10, 20), only 20 for 0.005 and 0.010)
mm_rates <- c(0.000, 0.001, 0.010)	# mismatch rate (CAPG manuscript with error: c(0.000, 0.000, 0.010)
n.indiv <- 1				# number of individuals (CAPG manuscript: 50)
info.dir <- 'info'			# INPUT: where info files are stored
tail <- ".txt"				# OUTPUT: file extension

do.hets <- T
do.homoeos <- T
load.data.only <- F


for (hr in homoeo_rates) {
for (cvg in covg_rates) {
for (mm in mm_rates) {
	if (!dir.exists(out.dir))
		dir.create(out.dir)
	het.file <- sprintf("%s/CAPG_PL_%.3f_%d_mm%.3f_het%s", out.dir, hr, cvg, mm, tail)
	hom.file <- sprintf("%s/CAPG_PL_%.3f_%d_mm%.3f_hom%s", out.dir, hr, cvg, mm, tail)
	#if (cvg != 20 & (hr == 0.005 | hr == 0.010))	# CAPG manuscript
	#	next
	cat("Homoeologous rate ", hr, "; coverage level ", cvg, "; mismatch rate ", mm, ":\n")
	for (i in 1:n.indiv) {
		cat("Reading genotype ", i, "\n")
		filename <- sprintf("%s/%s/homr%.3fmm%.3f/cov%d/%s/aln%d_info.txt", parent.dir, sim.dir, hr, mm, cvg, info.dir, i-1)
		if (i == 1) {
			d <- read.table(filename, header = T, sep = "\t")
		} else {
			d.tmp <- read.table(filename, header = T, sep = "\t")
			d <- rbind(d, d.tmp)
		}
	}
	stopifnot(!load.data.only)
	stopifnot(sum(d$PositionA == d$PositionB) == nrow(d))
	for (i in 11:25) {
		d[,i] <- log(d[,i])
	}
	# PositionA	PositionB	Genotype	PP_hetA	PP_hetB	CovA	CovB
	if (do.hets) {
		d$PP_hetA <- d$PA.1. - apply(d[, c('PA.0.', 'PA.2.')], 1, max)
		d$PP_hetB <- d$PB.1. - apply(d[, c('PB.0.', 'PB.2.')], 1, max)
		if (write.results) {
			write.table(d[, c('PositionA', 'PositionB', 'Genotype', 'PP_hetA', 'PP_hetB', 'CovA', 'CovB')], file = het.file, row.names = F, sep = "\t", quote = F)
		}
	}
	if (do.homoeos) {
		d$Genotype_call <- gsub('/', '', d$Genotype_call)
		d$max_ppG <- apply(d[, c('PP.0.0.', 'PP.0.1.', 'PP.0.2.', 'PP.1.0.', 'PP.1.1.', 'PP.1.2.', 'PP.2.0.', 'PP.2.1.', 'PP.2.2.')], 1, max)	# can also accomplish with pmax()
		d$min_ppG <- apply(d[, c('PP.0.0.', 'PP.0.1.', 'PP.0.2.', 'PP.1.0.', 'PP.1.1.', 'PP.1.2.', 'PP.2.0.', 'PP.2.1.', 'PP.2.2.')], 1, min)	# can also accomplish with pmax()
		d$max_ppA <- apply(d[, c('PA.0.', 'PA.1.', 'PA.2.')], 1, max)
		d$min_ppA <- apply(d[, c('PA.0.', 'PA.1.', 'PA.2.')], 1, min)
		d$max_ppB <- apply(d[, c('PB.0.', 'PB.1.', 'PB.2.')], 1, max)
		d$min_ppB <- apply(d[, c('PB.0.', 'PB.1.', 'PB.2.')], 1, min)
		d.join <- d %>%
			group_by(PositionA, PositionB) %>%
			mutate(
				biallelic = length(unique(unlist(strsplit(Genotype_call, split="")))) == 2,
				N1 = unique(unlist(strsplit(Genotype_call, split="")))[1],
				N2 = if_else(biallelic, unique(unlist(strsplit(Genotype_call, split="")))[2], if_else(N1 == 'A', 'C', 'A')),
				CovA = mean(CovA),
				CovB = mean(CovB),
			) %>%
			summarise(
				biallelic = unique(biallelic),
				PP_homeo = pmax(
					sum(if_else(Major.allele == N1 & Minor.allele == N2, PP.0.2., 	# N1N1/N2N2
						if_else(Major.allele == N2 & Minor.allele == N1, PP.2.0., min_ppG)) - max_ppG),
					sum(if_else(Major.allele == N2 & Minor.allele == N1, PP.0.2.,	# N2N2/N1N1
						if_else(Major.allele == N1 & Minor.allele == N2, PP.2.0., min_ppG)) - max_ppG)),

				PP_homoA = - pmax(sum(if_else(Major.allele == N1, PA.0., 	# N1N1
						if_else(Minor.allele == N1, PA.2., if_else(Major.allele == N2, PA.2., PA.0.))) - max_ppA),
					sum(if_else(Major.allele == N2, PA.0.,			# N2N2
						if_else(Minor.allele == N2, PA.2., if_else(Major.allele == N1, PA.2., PA.0.))) - max_ppA)),

				PP_homoB = - pmax(sum(if_else(Major.allele == N1, PB.0., 	# N1N1
						if_else(Minor.allele == N1, PB.2., if_else(Major.allele == N2, PB.2., PB.0.))) - max_ppB),
					sum(if_else(Major.allele == N2, PB.0.,			# N2N2
						if_else(Minor.allele == N2, PB.2., if_else(Major.allele == N1, PB.2., PB.0.))) - max_ppB)),	

				CovA = unique(CovA),
				CovB = unique(CovB),
			)
		if (write.results) {
			write.table(d.join, file = hom.file, row.names = F, sep = "\t", quote = F)
		}
	}
}
}
}
