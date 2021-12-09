#!/usr/bin/Rscript
#
# Purpose: compute metrics as per the manuscript, using input from info files and output to "PL" files
#
# For heterozygote calling within individual within subgenome: Proof for Eq. (7)
# 0.1[\min_{m\in\{0,2\}}PL_{klg}(m) - PL_{klg}(1)]
#
# For allelic SNP calling within subgenome across individuals: Proof for Eq. (11)
# -0.1[\min_{\bm\in\{0,2\}} PL_{klg}(m)]
#
# For homoeologous SNP calling across individuals: Proof for Eq. (9)
# -0.1[\min_{\bm\in\{(0,2),(2,0)\}} PL_{kl}(m)]
#
# Latest run (Thu Jul 22 10:44:59 AM CDT 2021) replaces PL with worst (maximum) observed PL when the genotype does not have the PL value of a particular genotype recorded.

library(dplyr)
library(tidyr)

# homr0.005/pair/cov20/info_files/aln49_info.txt
# [1] "Genotype"      "ChromA"        "ChromB"        "PositionA"    
# [5] "PositionB"     "Genotype_call" "Call_Agenome"  "Call_Bgenome" 
# [9] "Major.allele"  "Minor.allele"  "PP.0.0."       "PP.0.1."      
#[13] "PP.0.2."       "PP.1.0."       "PP.1.1."       "PP.1.2."      
#[17] "PP.2.0."       "PP.2.1."       "PP.2.2."       "PA.0."        
#[21] "PA.1."         "PA.2."         "PB.0."         "PB.1."        
#[25] "PB.2."         "CovA"          "CovB"         

homoeo_rates <- c(0.005, 0.007, 0.010)
#homoeo_rates <- 0.007
covg_rates <- c(5, 10, 20)	# only 20 for 0.005 and 0.010
#covg_rates <- 20
n.indiv <- 50	# 50 [KSD] Only do first...for speed
parent.dir <- '/path'
out.dir <- 'results'	#'compare_RK_ksd'	# where the PL files go
info.dir <- 'info_files'	# where the info files are
tail <- "txt"

do.hets <- T
do.homoeos <- T
load.data.only <- F

for (hr in homoeo_rates) {
for (cvg in covg_rates) {
	if (cvg != 20 & (hr == 0.005 | hr == 0.010))
		next
	cat("Homoeologous rate ", hr, "; coverage level ", cvg, ":\n")
	# Genotype ChromA ChromB PositionA PositionB Genotype_Call Ref_A_allele Alt_A_allele Ref_B_allele Alt_B_allele PLA(0) PLA(1) PLA(2) PLB(0) PLB(1) PLB(2) CovA CovB
        for (i in 1:n.indiv) {
		cat("Reading genotype ", i, "\n")
		if (i == 1) {
			d <- read.table(sprintf("%s/homr%.3f/cov%d/%s/aln%d_GATK_info.txt", out.dir, hr, cvg, info.dir, i-1), header = T, sep = "\t")
		} else {
			d.tmp <- read.table(sprintf("%s/homr%.3f/cov%d/%s/aln%d_GATK_info.txt", out.dir, hr, cvg, info.dir, i-1), header = T, sep = "\t")
			d <- rbind(d, d.tmp)
		}
	}
	stopifnot(!load.data.only)
	stopifnot(sum(d$PositionA == d$PositionB) == nrow(d))
	# PositionA	PositionB	Genotype	PP_hetA	PP_hetB	CovA	CovB
	if (do.hets) {
		d$PP_hetA <- 0.1*(apply(d[, c('PLA.0.', 'PLA.2.')], 1, min) - d$PLA.1.)
		write.table(d[, c('PositionA', 'PositionB', 'Genotype', 'PP_hetA', 'PP_hetB', 'CovA', 'CovB')], file = sprintf("%s/GATK_PL_%.3f_%d_het.%s", out.dir, hr, cvg, tail), row.names = F, sep = "\t", quote = F)
		#write.table(d[, c('PositionA', 'PositionB', 'Genotype', 'PP_hetA', 'PP_hetB', 'CovA', 'CovB')], file = sprintf("%s/GATK_PL_%.3f_%d_het_filtered.%s", out.dir, hr, cvg, tail), row.names = F, sep = "\t", quote = F)
	}
	# PositionA	PositionB	PP_homeo	PP_homoA	PP_homoB	CovA	CovB
	if (do.homoeos) {
		d$Genotype_call <- gsub('/', '', d$Genotype_Call)
		# compute PL for allotetraploid genotypes
		d$PLAB.0.0. <- d$PLA.0. + d$PLB.0.
		d$PLAB.0.1. <- d$PLA.0. + d$PLB.1.
		d$PLAB.0.2. <- d$PLA.0. + d$PLB.2.
		d$PLAB.1.0. <- d$PLA.1. + d$PLB.0.
		d$PLAB.1.1. <- d$PLA.1. + d$PLB.1.
		d$PLAB.1.2. <- d$PLA.1. + d$PLB.2.
		d$PLAB.2.0. <- d$PLA.2. + d$PLB.0.
		d$PLAB.2.1. <- d$PLA.2. + d$PLB.1.
		d$PLAB.2.2. <- d$PLA.2. + d$PLB.2.

		# maximum PL values for allotetraploid and diploid genotypes
		d$max_ppG <- apply(d[, c('PLAB.0.0.', 'PLAB.0.1.', 'PLAB.0.2.', 'PLAB.1.0.', 'PLAB.1.1.', 'PLAB.1.2.', 'PLAB.2.0.', 'PLAB.2.1.', 'PLAB.2.2.')], 1, max)	# can also accomplish with pmax()
		d$max_ppA <- apply(d[, c('PLA.0.', 'PLA.1.', 'PLA.2.')], 1, max)
		d$max_ppB <- apply(d[, c('PLB.0.', 'PLB.1.', 'PLB.2.')], 1, max)

		# convert genotype to AA/CC form
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
		print(dim(d.join))
		print(head(d.join))

		# extract alternate allele and then compute metrics
		d.join <- d.join %>%
			group_by(PositionA, PositionB) %>%
			mutate(
				biallelic = length(unique(unlist(strsplit(Genotype_Call, split="")))) == 2,
				N1 = unique(unlist(strsplit(Genotype_Call, split="")))[1],
				N2 = if_else(biallelic, unique(unlist(strsplit(Genotype_Call, split="")))[2], if_else(N1 == 'A', 'C', 'A')),
				CovA = mean(CovA),
				CovB = mean(CovB),
			)
		#	) %>%
		d.join <- d.join %>%
			group_by(PositionA, PositionB) %>%
			summarise(
				biallelic = unique(biallelic),
				PP_homeo = -0.1 * min(
					sum(if_else(Ref_A_allele == N1 & Ref_B_allele == N1, PLAB.0.2.,
						if_else(Ref_A_allele == N1 & Ref_B_allele == N2, PLAB.0.0.,
						if_else(Ref_A_allele == N2 & Ref_B_allele == N2, PLAB.2.0.,
						if_else(Ref_A_allele == N2 & Ref_B_allele == N1, PLAB.2.2., max_ppG))))),	# N1N1/N2N2
					sum(if_else(Ref_A_allele == N1 & Ref_B_allele == N1, PLAB.2.0.,
						if_else(Ref_A_allele == N1 & Ref_B_allele == N2, PLAB.2.2.,
						if_else(Ref_A_allele == N2 & Ref_B_allele == N2, PLAB.0.2.,
						if_else(Ref_A_allele == N2 & Ref_B_allele == N1, PLAB.0.0., max_ppG)))))),	# N2N2/N1N1
				PP_homoA = -0.1 * min(
					sum(if_else(Ref_A_allele == N1, PLA.0.,
						if_else(!is.na(Alt_A_allele) & Alt_A_allele == N1, PLA.2., max_ppA))),	# N1N1
					sum(if_else(Ref_A_allele == N2, PLA.0.,
						if_else(!is.na(Alt_A_allele) & Alt_A_allele == N2, PLA.2., max_ppA)))),	# N2N2
				PP_homoB = -0.1 * min(
					sum(if_else(Ref_B_allele == N1, PLB.0.,
						if_else(!is.na(Alt_B_allele) & Alt_B_allele == N1, PLB.2., max_ppB))),	# N1N1
					sum(if_else(Ref_B_allele == N2, PLB.0.,
						if_else(!is.na(Alt_B_allele) & Alt_B_allele == N2, PLB.2., max_ppB)))),	# N2N2
				CovA = unique(CovA),
				CovB = unique(CovB),
			)
		write.table(d.join, file = sprintf("%s/GATK_PL_%.3f_%d_Ho.%s", out.dir, hr, cvg, tail), row.names = F, sep = "\t", quote = F)
	}
}
}
