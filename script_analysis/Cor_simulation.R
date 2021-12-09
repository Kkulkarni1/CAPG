library(dplyr)
library(ggplot2)

CAPG_homo_data <- read.table("../analysis_simulation_real_data/simulation/correlation_CAPG_GATK/PL_CAPG_homo_0.007_20_wf_subset.txt", header = TRUE)
CAPG_het_data <- read.table("../analysis_simulation_real_data/simulation/correlation_CAPG_GATK/PL_CAPG_het_0.007_20_wf_genotype1_subset.txt", header = TRUE)

GATK_homo_data <- read.table("../analysis_simulation_real_data/simulation/correlation_CAPG_GATK/GATK_PL_0.007_20_homo_filtered_subset.txt", header = TRUE)
GATK_het_data <- read.table("../analysis_simulation_real_data/simulation/correlation_CAPG_GATK/GATK_PL_0.007_20_het_filtered_genotype1_subset.txt", header = TRUE)

CAPG_homeo_data <- read.table("../analysis_simulation_real_data/simulation/correlation_CAPG_GATK/PL_CAPG_homeo_0.007_20_wf_subset.txt", header = TRUE)
GATK_homeo_data <- read.table("../analysis_simulation_real_data/simulation/correlation_CAPG_GATK/GATK_PL_0.007_20_homeo_filtered_subset.txt", header = TRUE)

merge_ho_data <- cbind(CAPG_homo_data, GATK_homo_data)
merge_het_data <- cbind(CAPG_het_data, GATK_het_data)
merge_homeo_data <- cbind(CAPG_homeo_data, GATK_homeo_data)

#print(merge_homeo_data$PP_homeo)
#merge_homeo_data$PP_homeo <- log(merge_homeo_data$PP_homeo)
#merge_homeo_data$PP_homeo_GATK <- log(merge_homeo_data$PP_homeo_GATK) 
#print(merge_homeo_data$PP_homeo)

pdf("Cor_sim_homo_CAPG_GATK_subset.pdf")
ggplot(merge_ho_data, aes(x=PP_homoA, y=PP_GATK_homoA)) + geom_point(aes(colour = factor(Truth_homoA))) + xlab("CAPG") + ylab("GATK") + ggtitle("CAPG vs. GATK (Simulation allelic SNP)")
dev.off()


pdf("Cor_sim_het_CAPG_GATK_subset.pdf")
ggplot(merge_het_data, aes(x=PP_hetA, y=PP_GATK_hetA)) + geom_point(aes(colour = factor(Truth_het_A))) + xlab("CAPG") + ylab("GATK") + ggtitle("CAPG vs. GATK (Simulation heterozygosity)")
dev.off()

pdf("Cor_sim_homeo_CAPG_GATK_subset.pdf")
ggplot(merge_homeo_data, aes(x=PP_homeo, y=PP_homeo_GATK)) + geom_point(aes(colour = factor(Truth_homeo))) + xlab("CAPG") + ylab("GATK") + ggtitle("CAPG vs. GATK (Simulation homeologous SNP)")
dev.off()
