#!/usr/bin/Rscript
#
# Purpose: Compute metrics from CAPG output for SNP calling.
#
# For heterozygous calling within individuals and subgenomes: Eq. (3)
# $\ln\Pr(M_{klg} = 1) - \max_{m\in\{0,2\}} \ln\Pr(M_{klg} = m)$, where $kinfo_files$
# indexes individual, $l$ indexes location, and $g\in\{0,1\}$ is subgenome.
#
# For allelic SNP calling across the sample within subgenomes: Eq. (5)
# Modified to use minimum posterior probability computed when the posterior
# probability of the required genotype is not computed (instead of 0).
#
# For homoeologous SNP calling across the sample: Eq. (4)
# Modified as above for allelic SNPs.
#
# Last run (Thu Jul 22 07:19:16 PM CDT 2021) uses just the *.ksd.wof.txt and
# *.ksd.wf.txt (not run) terminators in compare_RK_ksd. Previous results are
# not to be trusted. There was a bug in the allelic SNP calling.

library(dplyr)

homoeo_rates <- c(0.005, 0.007, 0.010)
covg_rates <- c(5, 10, 20)	# only 20 for 0.005 and 0.010
n.indiv <- 50	# 50 [KSD] Only do first...for speed
parent.dir <- '../CAPG_github'
out.dir <- '../analysis_simulation_real_data/simulation/CAPG/results'	# where the PL files go
wo.filter <- T			# don't use filter during info file formation
info.dir <- 'parent.dir/info_files'
# info.dir <- ifelse(wo.filter, 'info_files_new', 'info_files_new_w_filter')	# where the info files are
# tail <- paste("ksd.", ifelse(wo.filter, "wof.txt", "wf.txt"), sep = "")

do.hets <- T
do.homoeos <- T
load.data.only <- F

for (hr in homoeo_ratesinfo_files) {
for (cvg in covg_rates) {
	if (cvg != 20 & (hr == 0.005 | hr == 0.010))
		next
	cat("Homoeologous rate ", hr, "; coverage level ", cvg, ":\n")
	for (i in 1:n.indiv) {
		cat("Reading genotype ", i, "\n")
		if (i == 1) {
			d <- read.table(sprintf("homr%.3f/pair/cov%d/%s/aln%d_info.txt", hr, cvg, info.dir, i-1), header = T, sep = "\t")
		} else {
			d.tmp <- read.table(sprintf("homr%.3f/pair/cov%d/%s/aln%d_info.txt", hr, cvg, info.dir, i-1), header = T, sep = "\t")
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
		write.table(d[, c('PositionA', 'PositionB', 'Genotype', 'PP_hetA', 'PP_hetB', 'CovA', 'CovB')], file = sprintf("%s/PL_%.3f_%d_het.%s", out.dir, hr, cvg, tail), row.names = F, sep = "\t", quote = F)
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
						if_else(Minor.allele == N1, PA.2., min_ppA)) - max_ppA),
	# the intent of the following was to allow for calling invariant sites where all genotypes were N3N3, but there is nothing to require they define the same allele N3
						# if_else(Minor.allele == N1, PA.2., if_else(Major.allele == N2, PA.2., PA.0.))) - max_ppA), # N3N3 if N1 neither major nor minor allele
					sum(if_else(Major.allele == N2, PA.0.,			# N2N2
						if_else(Minor.allele == N2, PA.2., min_ppA)) - max_ppA)),
				PP_homoB = - pmax(sum(if_else(Major.allele == N1, PB.0., 	# N1N1
						if_else(Minor.allele == N1, PB.2., min_ppB)) - max_ppB),
					sum(if_else(Major.allele == N2, PB.0.,			# N2N2
						if_else(Minor.allele == N2, PB.2., min_ppB)) - max_ppB)),	
				CovA = unique(CovA),
				CovB = unique(CovB),
			)
		write.table(d.join, file = sprintf("%s/PL_%.3f_%d_Ho.%stxt", out.dir, hr, cvg, "filter"), row.names = F, sep = "\t", quote = F)
	}
}
}
