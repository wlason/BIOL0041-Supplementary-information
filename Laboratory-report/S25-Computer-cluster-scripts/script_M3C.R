#Install M3C in my directory
library(BiocManager)
BiocManager::install("M3C", ask = FALSE, lib = "~/R_packages")

#Load clustering files
load("~/Jobs_data/seq_data.RData")
load("~/Jobs_data/annot_files.RData")

library(M3C)
M3C_data <- t(na.exclude(t(rnaseq_cancer)))
M3C_annot <- rnaseq_annot[rownames(rnaseq_annot) %in% rownames(M3C_data)]
M3C_annot$ID <- rownames(M3C_annot)

M3C_cancer <- M3C(mydata = M3C_data,
                  cores = 4,
                  des = ,
                  seed = 2020,
                  removeplots = TRUE,
                  doanalysis = FALSE)

save(M3C_cancer, file = "M3C_cancer.RData")