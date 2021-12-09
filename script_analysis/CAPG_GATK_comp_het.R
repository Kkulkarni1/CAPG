#!/usr/bin/Rscript
#
# Purpose: Plot PR curves for heterozygous SNP calls, homeologous SNPs, heterozygous SNP calls for CAPG and GATK
library(tidyverse)
library(dplyr)
library(pROC)
library(PRROC)


## use 0 as a threshold
add_point <- function(dat, column = "PP_hetA", truth = 'Truth_A') {
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

############# write functions!###################################
# order of 'or' needs to match capg1, capg2, gatk1, gatk2
draw_plot <- function(capg1, capg2, gatk1, gatk2, file_out, type = 'PP_hetA', truth = 'Truth_het_A',
                      or = c("CAPG_hetA", "CAPG_hetA_filtered", "GATK_hetA", "GATK_hetA_filtered")) {


  pr.capg.hetA_0.007_20 <- pr.curve(scores.class0 = capg1[[type]][capg1[[truth]]==1 ],
                                    scores.class1 = capg1[[type]][capg1[[truth]]==0 ], curve = T)
  pr.capg.hetA_filtered_0.007_20 <- pr.curve(scores.class0 = capg2[[type]][capg2[[truth]]==1 ],
                                             scores.class1 = capg2[[type]][capg2[[truth]]==0 ], curve = T)
  pr.gatk.hetA_0.007_20 <- pr.curve(scores.class0 = gatk1[[type]][gatk1[[truth]]==1 ],
                                    scores.class1 = gatk1[[type]][gatk1[[truth]]==0 ], curve = T)
  pr.gatk.hetAf_0.007_20 <- pr.curve(scores.class0 = gatk2[[type]][gatk2[[truth]]==1 ],
                                     scores.class1 = gatk2[[type]][gatk2[[truth]]==0 ], curve = T)

  a <- pr.capg.hetA_0.007_20$curve
  b <- pr.capg.hetA_filtered_0.007_20$curve
  c <- pr.gatk.hetA_0.007_20$curve
  d <- pr.gatk.hetAf_0.007_20$curve

  p1 = add_point(capg1, column = type, truth = truth)
  p2 = add_point(capg2, column = type, truth = truth)
  p3 = add_point(gatk1, column = type, truth = truth)
  p4 = add_point(gatk2, column = type, truth = truth)

  pdf(file_out, width = 5, height = 5)
  col = c("orange", "green", "blue", "red")
  plot(a[,1], a[,2], type="l", xlab = "Recall",ylab = "Precision", col = col[1])
  legend(x= "bottom", legend = or, col = col, lty=c(1, 1, 1, 1), cex=0.8, box.lty=0, bg="transparent")
  lines(b[,1], b[,2], col=col[2])
  lines(c[,1], c[,2], col=col[3])
  lines(d[,1], d[,2], col=col[4])

  x = c(p1[1], p2[1], p3[1], p4[1])
  y = c(p1[2], p2[2], p3[2], p4[2])
  points(x, y, col = col)
  dev.off()
}

draw_plot("/home/peanut/peanut2/CAPG/WGS/simulation/analysis_CAPG_V3/merged_files/merged_het_CAPG_0.007_20", "/home/peanut/peanut2/CAPG/WGS/simulation/analysis_CAPG_V3/analysis_filtered_info_files/merged_files/merged_CAPG_het_0.007_20.txt", "/home/peanut/peanut2/CAPG/WGS/simulation/analysis_GATK_V4/GATK_PL_0.007_20_het.txt", "/home/peanut/peanut2/CAPG/WGS/simulation/analysis_GATK_V4/GATK_PL_0.007_20_Ho_filtered.txt", file_out = 'a.pdf')
