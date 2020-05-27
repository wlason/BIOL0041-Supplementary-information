library(stringr)
library(deconstructSigs, lib.loc = "~/R_packages")
library(BSgenome.Hsapiens.UCSC.hg38, lib.loc = "~/R_packages")

load("~/R_packages/deconstructSigs/data/signatures.exome.cosmic.v3.may2019.rda")
solid_tumours <- c("blca", "brca", "cesc", "chol", "coad", "esca",  "hnsc", "kich", "kirc", "kirp", "lihc", "luad", "lusc", "prad", "stad", "thca", "ucec")
solid_tumours <- str_to_upper(solid_tumours)


for(pipe in c("muse", "mutect2", "somaticsniper", "varscan2")){
  df.maf.all <- NULL
  
  for(tumour in solid_tumours){
    assign('df.maf', get(load(paste0("~/4_pipelines_MAF_data/", tumour, "_", pipe, "_MAF_data.RData"))))
    
    # Only interested in point mutations, not small insertions/deletions:
    df.maf <- df.maf[which(df.maf$Variant_Type == "SNP"),]
    
    # Annotating sample type:
    df.maf$CodeType <- sapply(df.maf$Tumor_Sample_Barcode,
                                    function(x) substr(strsplit(x,"-")[[1]][4],1,2))
    df.maf$SampleType <- sapply(df.maf$CodeType,
                                      function(x) ifelse(x=="01","PrimaryTumour",
                                                         ifelse(x %in% c("10","11"),"Normal",
                                                                ifelse(x=="06","Metastasis","Other"))))
    
    df.maf <- df.maf[which(df.maf$SampleType %in% c("PrimaryTumour", "Normal"))]
    
    df.maf.all <- rbind(df.maf.all, df.maf)
  }
  
  assign(paste0("df.maf_all", pipe), value = df.maf.all)
  save(list = paste0("df.maf_all", pipe), file = paste0("df.maf_all", pipe, ".RData"))
  
  # Convert to deconstructSigs input:
  sigs.input <- mut.to.sigs.input(mut.ref = as.data.frame(df.maf.all), 
                                  sample.id = "Tumor_Sample_Barcode", 
                                  chr = "Chromosome", 
                                  pos = "Start_Position", 
                                  ref = "Reference_Allele", 
                                  alt = "Tumor_Seq_Allele2",
                                  bsg = BSgenome.Hsapiens.UCSC.hg38)
  
  sigs <- NULL
  for (sample in rownames(sigs.input)) {
    sigs_1 = whichSignatures(tumor.ref = sigs.input, 
                             signatures.ref = signatures.exome.cosmic.v3.may2019,
                             sample.id = sample, 
                             contexts.needed = TRUE,
                             signature.cutoff = 0,
                             tri.counts.method = 'default')
    sigs <- rbind(sigs,sigs_1$weights)
  }
  
  sigs <- data.frame(sigs)
  
  assign(paste0("sigs", pipe), value = sigs)
  save(list = paste0("sigs", pipe), file = paste0("sigs", pipe))
}
