#!/usr/bin/Rscript
#
# Purpose: Plot PR curves for heterozygous SNP calls, homeologous SNPs, heterozygous SNP calls for CAPG and GATK

library(dplyr)
library(pROC)
library(PRROC)

# Load CAPG het PL files
het_filtered_CAPG_0.007_5 <- read.table("../analysis_simulation_real_data/simulation/CAPG/PL_0.007_5_het.ksd.wof.wnh.new.txt", header = T)

het_filtered_CAPG_0.007_10 <- read.table("../analysis_simulation_real_data/simulation/CAPG/PL_0.007_10_het.ksd.wof.wnh.new.txt", header = T)

het_filtered_CAPG_0.007_20 <- read.table("../analysis_simulation_real_data/simulation/CAPG/PL_0.007_20_het.ksd.wof.wnh.new.txt", header = T)

# Load CAPG homo and homeo PL files
homeo_homo_filtered_CAPG_0.007_5 <- read.table("../analysis_simulation_real_data/simulation/CAPG/PL_0.007_5_Ho.ksd.wof.wnh.new.txt", header = T)

homeo_homo_filtered_CAPG_0.007_5$PP_homeo <- ifelse(homeo_homo_filtered_CAPG_0.007_5$PP_homeo=='-Inf', min(homeo_homo_filtered_CAPG_0.007_5$PP_homeo[is.finite(homeo_homo_filtered_CAPG_0.007_5$PP_homeo)])-1, homeo_homo_filtered_CAPG_0.007_5$PP_homeo)

homeo_homo_filtered_CAPG_0.007_10 <- read.table("../analysis_simulation_real_data/simulation/CAPG/PL_0.007_10_Ho.ksd.wof.wnh.new.txt", header = T)

homeo_homo_filtered_CAPG_0.007_10$PP_homeo <- ifelse(homeo_homo_filtered_CAPG_0.007_10$PP_homeo=='-Inf', min(homeo_homo_filtered_CAPG_0.007_10$PP_homeo[is.finite(homeo_homo_filtered_CAPG_0.007_10$PP_homeo)])-1, homeo_homo_filtered_CAPG_0.007_10$PP_homeo)

homeo_homo_filtered_CAPG_0.007_20 <- read.table("../analysis_simulation_real_data/simulation/CAPG/PL_0.007_20_Ho.ksd.wof.wnh.new.txt", header = T)

homeo_homo_filtered_CAPG_0.007_20$PP_homeo <- ifelse(homeo_homo_filtered_CAPG_0.007_20$PP_homeo=='-Inf', min(homeo_homo_filtered_CAPG_0.007_20$PP_homeo[is.finite(homeo_homo_filtered_CAPG_0.007_20$PP_homeo)])-1, homeo_homo_filtered_CAPG_0.007_20$PP_homeo)

# Load truth files
truth_het_0.007 <- read.table("../analysis_simulation_real_data/simulation/CAPG/truth_het_0.007.ksd.txt", header = T)

truth_homeo_0.007 <- read.table("../analysis_simulation_real_data/simulation/CAPG/truth_homeo_0.007.ksd.txt", header = T)

truth_homo_0.007 <- read.table("../analysis_simulation_real_data/simulation/CAPG/truth_homoA_B_0.007.ksd.txt", header = T)

# Merge PL and truth files
colnames(het_filtered_CAPG_0.007_5)[1] <- 'Position'
colnames(het_filtered_CAPG_0.007_10)[1] <- 'Position'
colnames(het_filtered_CAPG_0.007_20)[1] <- 'Position'

colnames(homeo_homo_filtered_CAPG_0.007_5)[1] <- 'Position'
colnames(homeo_homo_filtered_CAPG_0.007_10)[1] <- 'Position'
colnames(homeo_homo_filtered_CAPG_0.007_20)[1] <- 'Position'

PR_file_het_filtered_CAPG_0.007_5 <- full_join(truth_het_0.007, het_filtered_CAPG_0.007_5, by=c('Genotype','Position'))
PR_file_het_filtered_CAPG_0.007_10 <- full_join(truth_het_0.007, het_filtered_CAPG_0.007_10, by=c('Genotype','Position'))
PR_file_het_filtered_CAPG_0.007_20 <- full_join(truth_het_0.007, het_filtered_CAPG_0.007_20,, by=c('Genotype','Position'))

PR_file_homeo_filtered_CAPG_0.007_5 <- full_join(truth_homeo_0.007, homeo_homo_filtered_CAPG_0.007_5, by=c('Position'))
PR_file_homeo_filtered_CAPG_0.007_10 <- full_join(truth_homeo_0.007, homeo_homo_filtered_CAPG_0.007_10, by=c('Position'))
PR_file_homeo_filtered_CAPG_0.007_20 <- full_join(truth_homeo_0.007, homeo_homo_filtered_CAPG_0.007_20, by=c('Position'))

PR_file_homo_filtered_CAPG_0.007_5 <- full_join(truth_homo_0.007, homeo_homo_filtered_CAPG_0.007_5, by=c('Position'))
PR_file_homo_filtered_CAPG_0.007_10 <- full_join(truth_homo_0.007, homeo_homo_filtered_CAPG_0.007_10, by=c('Position'))
PR_file_homo_filtered_CAPG_0.007_20 <- full_join(truth_homo_0.007, homeo_homo_filtered_CAPG_0.007_20, by=c('Position'))


PR_file_het_filtered_CAPG_0.007_5[is.na(PR_file_het_filtered_CAPG_0.007_5)] = 0
PR_file_het_filtered_CAPG_0.007_10[is.na(PR_file_het_filtered_CAPG_0.007_10)] = 0
PR_file_het_filtered_CAPG_0.007_20[is.na(PR_file_het_filtered_CAPG_0.007_20)] = 0

PR_file_homeo_filtered_CAPG_0.007_5[is.na(PR_file_homeo_filtered_CAPG_0.007_5)] = 0
PR_file_homeo_filtered_CAPG_0.007_10[is.na(PR_file_homeo_filtered_CAPG_0.007_10)]= 0
PR_file_homeo_filtered_CAPG_0.007_20[is.na(PR_file_homeo_filtered_CAPG_0.007_20)]= 0

PR_file_homo_filtered_CAPG_0.007_5[is.na(PR_file_homo_filtered_CAPG_0.007_5)] = 0
PR_file_homo_filtered_CAPG_0.007_10[is.na(PR_file_homo_filtered_CAPG_0.007_10)] = 0
PR_file_homo_filtered_CAPG_0.007_20[is.na(PR_file_homo_filtered_CAPG_0.007_20)] = 0

#Plot curves
pr.capg.hetA_filtered_0.007_5 <- pr.curve(scores.class0 = PR_file_het_filtered_CAPG_0.007_5$PP_hetA[ PR_file_het_filtered_CAPG_0.007_5$Truth_het_A==1 ], scores.class1 = PR_file_het_filtered_CAPG_0.007_5$PP_hetA[ PR_file_het_filtered_CAPG_0.007_5$Truth_het_A==0 ], curve = T)
print("hetA_0.007_5")
print(pr.capg.hetA_filtered_0.007_5)

pr.capg.hetB_filtered_0.007_5 <- pr.curve(scores.class0 = PR_file_het_filtered_CAPG_0.007_5$PP_hetB[ PR_file_het_filtered_CAPG_0.007_5$Truth_het_B==1 ], scores.class1 = PR_file_het_filtered_CAPG_0.007_5$PP_hetB[ PR_file_het_filtered_CAPG_0.007_5$Truth_het_B==0 ], curve = T)
print("hetB_0.007_5")
print(pr.capg.hetB_filtered_0.007_5)

pr.capg.hetA_filtered_0.007_10 <- pr.curve(scores.class0 = PR_file_het_filtered_CAPG_0.007_10$PP_hetA[ PR_file_het_filtered_CAPG_0.007_10$Truth_het_A==1 ], scores.class1 = PR_file_het_filtered_CAPG_0.007_10$PP_hetA[ PR_file_het_filtered_CAPG_0.007_10$Truth_het_A==0 ], curve = T)
print("hetA_0.007_10")
print(pr.capg.hetA_filtered_0.007_10)

pr.capg.hetB_filtered_0.007_10 <- pr.curve(scores.class0 = PR_file_het_filtered_CAPG_0.007_10$PP_hetB[ PR_file_het_filtered_CAPG_0.007_10$Truth_het_B==1 ], scores.class1 = PR_file_het_filtered_CAPG_0.007_10$PP_hetB[ PR_file_het_filtered_CAPG_0.007_10$Truth_het_B==0 ], curve = T)
print("hetB_0.007_10")
print(pr.capg.hetB_filtered_0.007_10)

pr.capg.hetA_filtered_0.007_20 <- pr.curve(scores.class0 = PR_file_het_filtered_CAPG_0.007_20$PP_hetA[ PR_file_het_filtered_CAPG_0.007_20$Truth_het_A==1 ], scores.class1 = PR_file_het_filtered_CAPG_0.007_20$PP_hetA[ PR_file_het_filtered_CAPG_0.007_20$Truth_het_A==0 ], curve = T)
print("hetA_0.007_20")
print(pr.capg.hetA_filtered_0.007_20)

pr.capg.hetB_filtered_0.007_20 <- pr.curve(scores.class0 = PR_file_het_filtered_CAPG_0.007_20$PP_hetB[ PR_file_het_filtered_CAPG_0.007_20$Truth_het_B==1 ], scores.class1 = PR_file_het_filtered_CAPG_0.007_20$PP_hetB[ PR_file_het_filtered_CAPG_0.007_20$Truth_het_B==0 ], curve = T)
print("hetB_0.007_20")
print(pr.capg.hetB_filtered_0.007_20)

pr.capg.homeo_filtered_0.007_5 <- pr.curve(scores.class0 = PR_file_homeo_filtered_CAPG_0.007_5$PP_homeo[ PR_file_homeo_filtered_CAPG_0.007_5$Truth_homeo == 1 ], scores.class1 = PR_file_homeo_filtered_CAPG_0.007_5$PP_homeo[ PR_file_homeo_filtered_CAPG_0.007_5$Truth_homeo == 0 ], curve = T)
print("homeo_0.007_5")
print(pr.capg.homeo_filtered_0.007_5)

pr.capg.homeo_filtered_0.007_10 <- pr.curve(scores.class0 = PR_file_homeo_filtered_CAPG_0.007_10$PP_homeo[ PR_file_homeo_filtered_CAPG_0.007_10$Truth_homeo == 1 ], scores.class1 = PR_file_homeo_filtered_CAPG_0.007_10$PP_homeo[ PR_file_homeo_filtered_CAPG_0.007_10$Truth_homeo == 0 ], curve = T)
print("homeo_0.007_10")
print(pr.capg.homeo_filtered_0.007_10)


pr.capg.homeo_filtered_0.007_20 <- pr.curve(scores.class0 = PR_file_homeo_filtered_CAPG_0.007_20$PP_homeo[ PR_file_homeo_filtered_CAPG_0.007_20$Truth_homeo == 1 ], scores.class1 = PR_file_homeo_filtered_CAPG_0.007_20$PP_homeo[ PR_file_homeo_filtered_CAPG_0.007_20$Truth_homeo == 0 ], curve = T)
print("homeo_0.007_20")
print(pr.capg.homeo_filtered_0.007_20)


pr.capg.homoA_filtered_0.007_5 <- pr.curve(scores.class0 = PR_file_homo_filtered_CAPG_0.007_5$PP_homoA[ PR_file_homo_filtered_CAPG_0.007_5$Truth_A == 1 ], scores.class1 = PR_file_homo_filtered_CAPG_0.007_5$PP_homoA[ PR_file_homo_filtered_CAPG_0.007_5$Truth_A == 0 ], curve = T)
print("homoA_0.007_5")
print(pr.capg.homoA_filtered_0.007_5)

pr.capg.homoB_filtered_0.007_5 <- pr.curve(scores.class0 = PR_file_homo_filtered_CAPG_0.007_5$PP_homoB[ PR_file_homo_filtered_CAPG_0.007_5$Truth_B == 1 ], scores.class1 = PR_file_homo_filtered_CAPG_0.007_5$PP_homoB[ PR_file_homo_filtered_CAPG_0.007_5$Truth_B == 0 ], curve = T)
print("homoB_0.007_5")
print(pr.capg.homoB_filtered_0.007_5)

pr.capg.homoA_filtered_0.007_10 <- pr.curve(scores.class0 = PR_file_homo_filtered_CAPG_0.007_10$PP_homoA[ PR_file_homo_filtered_CAPG_0.007_10$Truth_A == 1 ], scores.class1 = PR_file_homo_filtered_CAPG_0.007_10$PP_homoA[ PR_file_homo_filtered_CAPG_0.007_10$Truth_A == 0 ], curve = T)
print("homoA_0.007_10")
print(pr.capg.homoA_filtered_0.007_10)

pr.capg.homoB_filtered_0.007_10 <- pr.curve(scores.class0 = PR_file_homo_filtered_CAPG_0.007_10$PP_homoB[ PR_file_homo_filtered_CAPG_0.007_10$Truth_B == 1 ], scores.class1 = PR_file_homo_filtered_CAPG_0.007_10$PP_homoB[ PR_file_homo_filtered_CAPG_0.007_10$Truth_B == 0 ], curve = T)
print("homoB_0.007_10")
print(pr.capg.homoB_filtered_0.007_10)

pr.capg.homoA_filtered_0.007_20 <- pr.curve(scores.class0 = PR_file_homo_filtered_CAPG_0.007_20$PP_homoA[ PR_file_homo_filtered_CAPG_0.007_20$Truth_A == 1 ], scores.class1 = PR_file_homo_filtered_CAPG_0.007_20$PP_homoA[ PR_file_homo_filtered_CAPG_0.007_20$Truth_A == 0 ], curve = T)
print("homoA_0.007_20")
print(pr.capg.homoA_filtered_0.007_20)

pr.capg.homoB_filtered_0.007_20 <- pr.curve(scores.class0 = PR_file_homo_filtered_CAPG_0.007_20$PP_homoB[ PR_file_homo_filtered_CAPG_0.007_20$Truth_B == 1 ], scores.class1 = PR_file_homo_filtered_CAPG_0.007_20$PP_homoB[ PR_file_homo_filtered_CAPG_0.007_20$Truth_B == 0 ], curve = T)
print("homoB_0.007_20")
print(pr.capg.homoB_filtered_0.007_20)


curve.points_hetA_filtered_CAPG_0.007_5 <- pr.capg.hetA_filtered_0.007_5$curve
curve.points_hetA_filtered_CAPG_0.007_10 <- pr.capg.hetA_filtered_0.007_10$curve
curve.points_hetA_filtered_CAPG_0.007_20 <- pr.capg.hetA_filtered_0.007_20$curve

curve.points_hetB_filtered_CAPG_0.007_5 <- pr.capg.hetB_filtered_0.007_5$curve
curve.points_hetB_filtered_CAPG_0.007_10 <- pr.capg.hetB_filtered_0.007_10$curve
curve.points_hetB_filtered_CAPG_0.007_20 <- pr.capg.hetB_filtered_0.007_20$curve

curve.points_homeo_filtered_CAPG_0.007_5 <- pr.capg.homeo_filtered_0.007_5$curve
curve.points_homeo_filtered_CAPG_0.007_10 <- pr.capg.homeo_filtered_0.007_10$curve
curve.points_homeo_filtered_CAPG_0.007_20 <- pr.capg.homeo_filtered_0.007_20$curve

curve.points_homoA_filtered_CAPG_0.007_5 <- pr.capg.homoA_filtered_0.007_5$curve
curve.points_homoA_filtered_CAPG_0.007_10 <- pr.capg.homoA_filtered_0.007_10$curve
curve.points_homoA_filtered_CAPG_0.007_20 <- pr.capg.homoA_filtered_0.007_20$curve

curve.points_homoB_filtered_CAPG_0.007_5 <- pr.capg.homoB_filtered_0.007_5$curve
curve.points_homoB_filtered_CAPG_0.007_10 <- pr.capg.homoB_filtered_0.007_10$curve
curve.points_homoB_filtered_CAPG_0.007_20 <- pr.capg.homoB_filtered_0.007_20$curve

#png(file="Plot_PR_het_CAPG_change_cov.png", width = 4, height = 4, units = 'in', res = 300)
postscript(file="/home/peanut/peanut2/CAPG_github/Simulation/CAPG/Plot_PR_het_CAPG_change_cov.eps")
plot(curve.points_hetA_filtered_CAPG_0.007_5[,1], curve.points_hetA_filtered_CAPG_0.007_5[,2], type="l", xlab="Recall", ylab="Precision", col="orange")
legend(x= "bottomleft", inset = 0.05, legend=c("CAPG_hetA [C=5]", "CAPG_hetB [C=5]", "CAPG_hetA [C=10]", "CAPG_hetB [C=10]", "CAPG_hetA [C=20]", "CAPG_hetB [C=20]"), col=c("orange", "orange", "blue", "blue", "red", "red"), lty=c(1,2,1,2,1,2), cex=0.8, box.lty=0)
lines(curve.points_hetB_filtered_CAPG_0.007_5[,1], curve.points_hetB_filtered_CAPG_0.007_5[,2], type="l", lty=2, col="orange")
lines(curve.points_hetA_filtered_CAPG_0.007_10[,1], curve.points_hetA_filtered_CAPG_0.007_10[,2], type="l", col="blue")
lines(curve.points_hetB_filtered_CAPG_0.007_10[,1], curve.points_hetB_filtered_CAPG_0.007_10[,2], type="l", lty=2, col="blue")
lines(curve.points_hetA_filtered_CAPG_0.007_20[,1], curve.points_hetA_filtered_CAPG_0.007_20[,2], type="l", col="red")
lines(curve.points_hetB_filtered_CAPG_0.007_20[,1], curve.points_hetB_filtered_CAPG_0.007_20[,2], type="l", lty=2, col="red")
dev.off()

#png(file="Plot_PR_homeo_CAPG_change_cov.png", width = 4, height = 4, units = 'in', res = 300)
postscript("/home/peanut/peanut2/CAPG_github/Simulation/CAPG/Plot_PR_homeo_CAPG_change_cov.eps")
plot(curve.points_hetA_filtered_CAPG_0.007_5[,1], curve.points_hetA_filtered_CAPG_0.007_5[,2], type="l", xlab="Recall", ylab="Precision", col="orange")

plot(as.numeric(curve.points_homeo_filtered_CAPG_0.007_5[,1]), as.numeric(curve.points_homeo_filtered_CAPG_0.007_5[,2]), type="l", xlab="Recall", ylab="Precision", col="orange")
legend(x= "bottomleft", inset = 0.05, legend=c("CAPG_homeo [C=5]", "CAPG_homeo [C=10]", "CAPG_homeo [C=20]"), col=c("orange", "blue", "red"), lty=c(1,1,1), cex=0.8, box.lty=0)
lines(as.numeric(curve.points_homeo_filtered_CAPG_0.007_10[,1]), as.numeric(curve.points_homeo_filtered_CAPG_0.007_10[,2]), type="l", col="blue")
lines(as.numeric(curve.points_homeo_filtered_CAPG_0.007_20[,1]), as.numeric(curve.points_homeo_filtered_CAPG_0.007_20[,2]), type="l", col="red")
dev.off()

#png(file="Plot_PR_homo_CAPG_change_cov.png", width = 4, height = 4, units = 'in', res = 300)
postscript("/home/peanut/peanut2/CAPG_github/Simulation/CAPG/Plot_PR_homo_CAPG_change_cov.eps")
plot(curve.points_homoA_filtered_CAPG_0.007_5[,1], curve.points_homoA_filtered_CAPG_0.007_5[,2], type="l", xlab="Recall", ylab="Precision", col="orange")
legend(x= "bottomleft", inset = 0.05, legend=c("CAPG_homoA [C=5]", "CAPG_homoB [C=5]", "CAPG_homoA [C=10]", "CAPG_homoB [C=10]", "CAPG_homoA [C=20]", "CAPG_homoB [C=20]"), col=c("orange", "orange", "blue", "blue", "red", "red"), lty=c(1,2,1,2,1,2), cex=0.8, box.lty=0)
lines(curve.points_homoB_filtered_CAPG_0.007_5[,1], curve.points_homoB_filtered_CAPG_0.007_5[,2], type="l", lty=2, col="orange")
lines(curve.points_homoA_filtered_CAPG_0.007_10[,1], curve.points_homoA_filtered_CAPG_0.007_10[,2], type="l", col="blue")
lines(curve.points_homoB_filtered_CAPG_0.007_10[,1], curve.points_homoB_filtered_CAPG_0.007_10[,2], type="l", lty=2, col="blue")
lines(curve.points_homoA_filtered_CAPG_0.007_20[,1], curve.points_homoA_filtered_CAPG_0.007_20[,2], type="l", col="red")
lines(curve.points_homoB_filtered_CAPG_0.007_20[,1], curve.points_homoB_filtered_CAPG_0.007_20[,2], type="l", lty=2, col="red")
dev.off()
