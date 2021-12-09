#!/usr/bin/Rscript

library(dplyr)

n.target <- 1000
genotype_list <- c("SRR4124062", "SRR4124066", "SRR4124068", "SRR4124074", "SRR4124078", "SRR8361734", "SRR8361735", "SRR8361736", "SRR8361737", "SRR8361738", "SRR8736998", "SRR8737008", "SRR8737061", "SRR8737062")
for (i in 1:n.target) {
 
	for (j in genotype_list) {
		cat("Reading genotype ", j, "target", i, "\n")
		if (i == 1) {
			d <- read.table(sprintf("./info_files/target_%d/%s_%d_info.txt", i,j,i), header = T, sep = "\t")
		} else {
			d.tmp <- read.table(sprintf("./info_files/target_%d/%s_%d_info.txt", i,j,i), header = T, sep = "\t")
			d <- rbind(d, d.tmp)
		}
	}
}

sum(d$PositionA == d$PositionB) == nrow(d)

for (i in 11:25) {
	d[,i] <- log(d[,i])
}

d$PP_hetA <- d$PA.1. - apply(d[, c('PA.0.', 'PA.2.')], 1, max)
d$PP_hetB <- d$PB.1. - apply(d[, c('PB.0.', 'PB.2.')], 1, max)
write.table(d[, c('PositionA', 'PositionB', 'Genotype', 'ChromA', 'ChromB', 'PP_hetA', 'PP_hetB', 'CovA', 'CovB')], file = "PL_het_peannut.txt", row.names = , F, sep = "\t", quote = F)

d$Genotype_call <- gsub('/', '', d$Genotype_call)
d$max_ppG <- apply(d[, c('PP.0.0.', 'PP.0.1.', 'PP.0.2.', 'PP.1.0.', 'PP.1.1.', 'PP.1.2.', 'PP.2.0.', 'PP.2.1.', 'PP.2.2.')], 1, max)   # can also accomplish with pmax()

d$max_ppA <- apply(d[, c('PA.0.', 'PA.1.', 'PA.2.')], 1, max)
d$max_ppB <- apply(d[, c('PB.0.', 'PB.1.', 'PB.2.')], 1, max)

d.join <- d %>%
	group_by(PositionA, PositionB, ChromA, ChromB) %>%
	mutate(
		biallelic = length(unique(unlist(strsplit(Genotype_call, split="")))) == 2,
		#print(unlist(strsplit(Genotype_call, split="/"))),
		N1 = unique(unlist(strsplit(Genotype_call, split="")))[1],
                #print(N1),
		N2 = if_else(biallelic, unique(unlist(strsplit(Genotype_call, split="")))[2], if_else(N1 == 'A', 'C', 'A')),
                #ChromA = (ChromA),
                #ChromB = (ChromB),
		CovA = mean(CovA),
		CovB = mean(CovB),
	) %>%
	summarise(
		biallelic = unique(biallelic),
		PP_homeo = pmax(
			sum(if_else(Major.allele == N1 & Minor.allele == N2, PP.0.2.,   # N1N1/N2N2
				if_else(Major.allele == N2 & Minor.allele == N1, PP.2.0., -Inf)) - max_ppG),
			sum(if_else(Major.allele == N2 & Minor.allele == N1, PP.0.2.,   # N2N2/N1N1
				if_else(Major.allele == N1 & Minor.allele == N2, PP.2.0., -Inf)) - max_ppG)),
		PP_homoA = - pmax(sum(if_else(Major.allele == N1, PA.0.,        # N1N1 or N3N3 if N1 is neither major nor minor allele

				if_else(Minor.allele == N1, PA.2., if_else(Major.allele == N2, PA.2., PA.0.))) - max_ppA),
			sum(if_else(Major.allele == N2, PA.0.,                  # N2N2
				if_else(Minor.allele == N2, PA.2., if_else(Major.allele == N1, PA.2., PA.0.))) - max_ppA)),
		PP_homoB = - pmax(sum(if_else(Major.allele == N1, PB.0.,        # N1N1
				if_else(Minor.allele == N1, PB.2., if_else(Major.allele == N2, PB.2., PB.0.))) - max_ppB),
			sum(if_else(Major.allele == N2, PB.0.,  # N2N2
				if_else(Minor.allele == N2, PB.2., if_else(Major.allele == N1, PB.2., PB.0.))) - max_ppB)),
                #ChromA = (ChromA),
                #ChromB = (ChromB),
		CovA = unique(CovA),
		CovB = unique(CovB),
		)
                print(d.join)
		write.table(d.join, file = "PL_Ho_peanut.txt", row.names = F, sep = "\t", quote = F)




