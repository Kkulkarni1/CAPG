#!/usr/bin/Rscript
# @file get_metrics_GATK.R
# @author R. Kulkarni
#
# Purpose: compute GATK metrics as per the manuscript, using input from info files and output to "PL" files
#

library(dplyr)
library(tidyr)

write.results <- T			# write output files
parent.dir <- '.'			# OUTPUT: CAPG repository top level directory
sim.dir <- 'data/simulation'		# IN/OUTPUT: simulation data input/output
out.dir <-				# OUTPUT: where the metric files go
	 paste(parent.dir, '/data/simulation/results', sep='')
homoeo_rates <- c(0.005, 0.007)		# homoeologous rates (CAPG manuscript: c(0.005, 0.007, 0.010))
covg_rates <- c(5, 10)			# coverage rate (CAPG manuscript: c(5, 10, 20), only 20 for 0.005 and 0.010)
mm_rates <- c(0.000, 0.001, 0.010)	# mismatch rate (CAPG manuscript with error: c(0.000, 0.000, 0.010)
n.indiv <- 1				# number of individuals (CAPG manuscript: 50)
info.dir <- 'info'			# INPUT: where info files are stored
info.tail <- '_GATK_info.txt'		# OUTPUT: information tail output
tail <- ".txt"				# OUTPUT: file extension

do.hets <- T
do.homoeos <- T
load.data.only <- F

if (!dir.exists(out.dir))
	dir.create(out.dir)

for (hr in homoeo_rates) {
for (cvg in covg_rates) {
for (mm in mm_rates) {
	#if (cvg != 20 & (hr == 0.005 | hr == 0.010))
	#	next
	het.file <- sprintf("%s/GATK_PL_%.3f_%d_mm%.3f_het%s", out.dir, hr, cvg, mm, tail)
	hom.file <- sprintf("%s/GATK_PL_%.3f_%d_mm%.3f_hom%s", out.dir, hr, cvg, mm, tail)
	cat("Homoeologous rate ", hr, "; coverage level ", cvg, ":\n")
	# Genotype	ChromA	ChromB	PositionA	PositionB	Genotype_Call	Ref_A_allele	Alt_A_allele	Ref_B_allele	Alt_B_allele	PLA(0)	PLA(1)	PLA(2)	PLB(0)	PLB(1)PLB(2)	CovA	CovB
        for (i in 1:n.indiv) {
		cat("Reading genotype ", i, "\n")
		filename <- sprintf("%s/%s/homr%.3fmm%.3f/cov%d/GATK/%s/aln%d%s", parent.dir, sim.dir, hr, mm, cvg, info.dir, i-1, info.tail)
		if (i == 1) {
			d <- read.table(filename, header = T, sep = "\t")
		} else {
			d.tmp <- read.table(filename, header = T, sep = "\t")
			d <- rbind(d, d.tmp)
		}
	}
	stopifnot(!load.data.only)
	stopifnot(sum(d$PositionA == d$PositionB) == nrow(d))	# true in simulation

	# PositionA	PositionB	Genotype	PP_hetA	PP_hetB	CovA	CovB
	if (do.hets) {
		d$PP_hetA <- 0.1*log(10)*(apply(d[, c('PLA.0.', 'PLA.2.')], 1, min) - d$PLA.1.)
		d$PP_hetB <- 0.1*log(10)*(apply(d[, c('PLB.0.', 'PLB.2.')], 1, min) - d$PLB.1.)
		if (write.results) {
			write.table(d[, c('PositionA', 'PositionB', 'Genotype', 'PP_hetA', 'PP_hetB', 'CovA', 'CovB')],
				file = het.file, row.names = F, sep = "\t", quote = F)
			cat("Heterozygosity metrics written to file.\n")
		}
	}
	
	# PositionA	PositionB	PP_homeo	PP_homoA	PP_homoB	CovA	CovB
	if (do.homoeos) {
		#d$Genotype_Call <- gsub('/', '', d$Genotype_Call)

		# compute PL for allotetraploid genotypes
		d$PLAB.0.0. <- d$PLA.0. + d$PLB.0.
		d$PLAB.0.1. <- d$PLA.0. + d$PLB.1.
		d$PLAB.0.2. <- d$PLA.0. + d$PLB.2.
		d$PLAB.1.0. <- d$PLA.1. + d$PLB.0.
		d$PLAB.1.1. <- d$PLA.1. + d$PLB.1.	# could be N1N3/N2N3
		d$PLAB.1.2. <- d$PLA.1. + d$PLB.2.
		d$PLAB.2.0. <- d$PLA.2. + d$PLB.0.
		d$PLAB.2.1. <- d$PLA.2. + d$PLB.1.
		d$PLAB.2.2. <- d$PLA.2. + d$PLB.2.

		# maximum PL values for allotetraploid and diploid genotypes
		d$max_ppG <- apply(d[, c('PLAB.0.0.', 'PLAB.0.1.', 'PLAB.0.2.', 'PLAB.1.0.', 'PLAB.1.1.', 'PLAB.1.2.', 'PLAB.2.0.', 'PLAB.2.1.', 'PLAB.2.2.')], 1, max)	# can also accomplish with pmax()
		d$max_ppA <- apply(d[, c('PLA.0.', 'PLA.1.', 'PLA.2.')], 1, max)
		d$max_ppB <- apply(d[, c('PLB.0.', 'PLB.1.', 'PLB.2.')], 1, max)
	        print("here_2")

		# convert genotype to AACC form
		d.join <- d %>%
			separate(col = Genotype_Call, into = c("Call_Agenome", "Call_Bgenome")) %>%
				mutate(
					Call_Agenome = if_else(Call_Agenome==0,
						paste(Ref_A_allele, Ref_A_allele, sep=""),
						if_else(Call_Agenome == 1,
							paste(Ref_A_allele, Alt_A_allele, sep=""),
							paste(Alt_A_allele, Alt_A_allele, sep=""))),
					Call_Bgenome = if_else(Call_Bgenome==0,
						paste(Ref_B_allele, Ref_B_allele, sep=""),
						if_else(Call_Bgenome == 1,
							paste(Ref_B_allele, Alt_B_allele, sep=""),
							paste(Alt_B_allele, Alt_B_allele, sep=""))),
					Genotype_Call = paste(Call_Agenome, Call_Bgenome, sep="")
				)

		cat("After computing genotypes\n")
		print(head(as.data.frame(d.join)))
	
		# extract alternate allele and then compute metrics
		d.join <- d.join %>%
			group_by(PositionA, PositionB, ChromA, ChromB) %>%
			mutate(
				biallelic = length(unique(unlist(strsplit(Genotype_Call, split="")))), # == 2,
				N1 = names(which.max(table(unlist(strsplit(Genotype_Call, split=""))))),
				N2 = if_else(biallelic>=2, names(sort(table(unlist(strsplit(Genotype_Call, split=""))), decreasing=T)[2]), if_else(N1 == 'A', 'C', 'A')),
				CovA = mean(CovA),
				CovB = mean(CovB),
			)
			#	) %>%
		cat("After computing N1, N2\n")
		print(head(as.data.frame(d.join)))
		d.join <- d.join %>%
			group_by(PositionA, PositionB, ChromA, ChromB) %>%
			summarise(
				biallelic = unique(biallelic),
				PP_homeo = -0.1 * log(10) * min(
					sum(if_else(Ref_A_allele == N1 & !is.na(Alt_B_allele) & Alt_B_allele == N2, PLAB.0.2.,
						if_else(Ref_A_allele == N1 & Ref_B_allele == N2, PLAB.0.0.,
						if_else(!is.na(Alt_A_allele) & Alt_A_allele == N1 & !is.na(Alt_B_allele) & Alt_B_allele == N2, PLAB.2.2.,
						if_else(!is.na(Alt_A_allele) & Alt_A_allele == N1 & Ref_B_allele == N2, PLAB.2.0., max_ppG))))),	# N1N1/N2N2
					sum(if_else(Ref_A_allele == N2 & Ref_B_allele == N1, PLAB.0.0.,
						if_else(Ref_A_allele == N2 & !is.na(Alt_B_allele) & Alt_B_allele == N1, PLAB.0.2.,
						if_else(!is.na(Alt_A_allele) & Alt_A_allele == N2 & Ref_B_allele == N2, PLAB.2.0.,
						if_else(!is.na(Alt_A_allele) & Alt_A_allele == N2 & !is.na(Alt_B_allele) & Alt_B_allele == N1, PLAB.2.2., max_ppG)))))),	# N2N2/N1N1
				PP_homoA = 0.1 * log(10) * min(
					sum(if_else(Ref_A_allele == N1, PLA.0.,
						if_else(!is.na(Alt_A_allele) & Alt_A_allele == N1, PLA.2., max_ppA))),	# N1N1
					sum(if_else(Ref_A_allele == N2, PLA.0.,
						if_else(!is.na(Alt_A_allele) & Alt_A_allele == N2, PLA.2., max_ppA)))),	# N2N2
				PP_homoB = 0.1 * log(10) * min(
					sum(if_else(Ref_B_allele == N1, PLB.0.,
						if_else(!is.na(Alt_B_allele) & Alt_B_allele == N1, PLB.2., max_ppB))),	# N1N1
					sum(if_else(Ref_B_allele == N2, PLB.0.,
						if_else(!is.na(Alt_B_allele) & Alt_B_allele == N2, PLB.2., max_ppB)))),	# N2N2
				CovA = unique(CovA),
				CovB = unique(CovB),
			)
		cat("After computing PPs\n")
		#print(dim(d.join))
		print(head(as.data.frame(d.join)))
		print("done")
		warnings()
		if (write.results) {
			write.table(d.join, file = hom.file, row.names = F, sep = "\t", quote = F)
		}
	}
}
}
}
