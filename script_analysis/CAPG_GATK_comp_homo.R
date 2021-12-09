#!/usr/bin/Rscript
#
# Purpose: Plot PR curves for homoerozygous SNP calls, homeologous SNPs, homoerozygous SNP calls for CAPG and GATK
library(tidyverse)
library(dplyr)
library(pROC)
library(PRROC)


## use 0 as a threshold
add_point <- function(dat, column = "PP_homoA", truth = 'Truth_homoA') {
  tmp = (dat[column] > 0)
  idf = which(dat[truth] == 0)
  idt = which(dat[truth] == 1)
  tn = sum(tmp[idf] == dat[idf, truth])
  tp = sum(tmp[idt] == dat[idt, truth])
  inf = which(tmp == 0)
  int = which(tmp == 1)
  fn = sum(tmp[inf] != dat[inf, truth])
  fp = sum(tmp[int] != dat[int, truth])
  p = tp/(tp + fp)
  r = tp/(fn + tp)
  return(c(r, p))
}

############ write functions!###################################
# order of 'or' needs to match capg1, capg2, gatk1, gatk2
draw_plot <- function(capg1, capg2, gatk1, gatk2, file_out, type = 'PP_homoA', truth = 'Truth_homoA',
                      or = c("CAPG_homoA", "CAPG_homoA_filtered", "GATK_homoA", "GATK_homoA_filtered")) {
  
  
  pr.capg.homoA_0.007_20 <- pr.curve(scores.class0 = capg1[[type]][capg1[[truth]]==1 ], 
                                    scores.class1 = capg1[[type]][capg1[[truth]]==0 ], curve = T)
  print("CAPG")
  print(pr.capg.homoA_0.007_20)
  pr.capg.homoA_filtered_0.007_20 <- pr.curve(scores.class0 = capg2[[type]][capg2[[truth]]==1 ], 
                                             scores.class1 = capg2[[type]][capg2[[truth]]==0 ], curve = T)
  print("CAPG_filtered")
  print(pr.capg.homoA_filtered_0.007_20)
  pr.gatk.homoA_0.007_20 <- pr.curve(scores.class0 = gatk1[[type]][ gatk1[[truth]]==1 ], 
                                    scores.class1 = gatk1[[type]][ gatk1[[truth]]==0 ], curve = T)
  print("GATK")
  print(pr.gatk.homoA_0.007_20)
  pr.gatk.homoAf_0.007_20 <- pr.curve(scores.class0 = gatk2[[type]][ gatk2[[truth]]==1 ], 
                                     scores.class1 = gatk2[[type]][ gatk2[[truth]]==0 ], curve = T)
  print("GATK_filtered")
  print(pr.gatk.homoAf_0.007_20)

  a <- pr.capg.homoA_0.007_20$curve
  b <- pr.capg.homoA_filtered_0.007_20$curve
  c <- pr.gatk.homoA_0.007_20$curve
  d <- pr.gatk.homoAf_0.007_20$curve
  
  p1 = add_point(capg1, column = type, truth = truth)
  p2 = add_point(capg2, column = type, truth = truth)
  p3 = add_point(gatk1, column = type, truth = truth)
  p4 = add_point(gatk2, column = type, truth = truth)
  
  png(filename=file_out)
  col = c("magenta", "magenta", "green", "green")
  plot(a[,1], a[,2], type="l", lty=2, xlab = "Recall",ylab = "Precision", col = col[1])
  legend(x= "bottom", legend = or, col = col, lty=c(2, 1, 2, 1), cex=0.8, box.lty=0, bg="transparent")
  lines(b[,1], b[,2], type="l", col=col[2])
  lines(c[,1], c[,2], type="l", lty=2, col=col[3])
  lines(d[,1], d[,2], type="l", col=col[4])
  
  x = c(p1[1], p2[1], p3[1], p4[1])
  y = c(p1[2], p2[2], p3[2], p4[2])
  points(x, y, col = col)
  dev.off()
}

PR_file_homo_CAPG_0.007_20 <- read.table("../analysis_simulation_real_data/simulation/CAPG/merged_files/merged_homo_CAPG_0.007_20", header = T)

PR_file_homo_filtered_CAPG_0.007_20 <- read.table("../analysis_simulation_real_data/simulation/CAPG/merged_files/merged_CAPG_homo_0.007_20.txt", header = T)

PR_file_homo_GATK_0.007_20 <- read.table("../analysis_simulation_real_data/simulation/GATK/merged_files/GATK_PL_0.007_20_homo.txt", header = T)

PR_file_homo_filtered_GATK_0.007_20 <- read.table("../analysis_simulation_real_data/simulation/GATK/merged_files/GATK_PL_0.007_20_homo_filtered.txt", header = T)


draw_plot(PR_file_homo_CAPG_0.007_20, PR_file_homo_filtered_CAPG_0.007_20, 
          PR_file_homo_GATK_0.007_20, PR_file_homo_filtered_GATK_0.007_20,
          file_out = 'CAPG_GATK_comp_homo.png')
