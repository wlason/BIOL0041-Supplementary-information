library(stringr)
library(deconstructSigs, lib.loc = "~/R_packages")
library(BSgenome.Hsapiens.UCSC.hg38, lib.loc = "~/R_packages")

load("~/R_data/signatures.exome.cosmic.v3.may2019.rda")
solid_tumours <- c("blca", "brca", "cesc", "chol", "coad", "esca",  "hnsc", "kich", "kirc", "kirp", "lihc", "luad", "lusc", "prad", "stad", "thca", "ucec")
solid_tumours <- str_to_upper(solid_tumours)
  
  calculate.sigs <- function(tumours, pipe, path_to_MAFs){
    df.maf.all <- NULL
    for(tumour in tumours){
    assign('df.maf', get(load(paste0(tumour, "_", pipe, "_MAF_data.RData"))))
    
    df.maf.short <- df.maf[,c("Tumor_Sample_Barcode",
                              "Cancer",
                              "Hugo_Symbol","Chromosome",
                              "Start_Position","End_Position",
                              "Variant_Classification",
                              "Variant_Type",
                              "Reference_Allele",
                              "Tumor_Seq_Allele1",
                              "Tumor_Seq_Allele2",
                              "One_Consequence",
                              "Consequence",
                              "SIFT","PolyPhen")]
    
    # Only interested in point mutations, not small insertions/deletions:
    df.maf.short <- df.maf.short[which(df.maf.short$Variant_Type == "SNP"),]
    
    # Annotating sample type:
    df.maf.short$CodeType <- sapply(df.maf.short$Tumor_Sample_Barcode,
                                    function(x) substr(strsplit(x,"-")[[1]][4],1,2))
    df.maf.short$SampleType <- sapply(df.maf.short$CodeType,
                                      function(x) ifelse(x=="01","PrimaryTumour",
                                                         ifelse(x %in% c("10","11"),"Normal",
                                                                ifelse(x=="06","Metastasis","Other"))))
    
    df.maf.short <- df.maf.short[which(df.maf.short$SampleType %in% c("PrimaryTumour", "Normal")), ]
  
  # Convert to deconstructSigs input:
  sigs.input <- mut.to.sigs.input(mut.ref = as.data.frame(df.maf.short), 
                                  sample.id = "Tumor_Sample_Barcode", 
                                  chr = "Chromosome", 
                                  pos = "Start_Position", 
                                  ref = "Reference_Allele", 
                                  alt = "Tumor_Seq_Allele2",
                                  bsg = BSgenome.Hsapiens.UCSC.hg38)
  
  sigs <- NULL
  percent <- as.integer(seq(0.1, 1, 0.1)*nrow(sigs.input))
  percent <- rownames(sigs.input)[percent]
  
  print("Starting calculations")
  for (sample in rownames(sigs.input)) {
    sigs_1 = whichSignatures(tumor.ref = sigs.input, 
                             signatures.ref = signatures.exome.cosmic.v3.may2019,
                             sample.id = sample, 
                             contexts.needed = TRUE,
                             signature.cutoff = 0,
                             tri.counts.method = 'default')
    sigs <- rbind(sigs,sigs_1$weights)
    if(sample %in% percent){
      z <- match(sample, percent)
      print(paste0("Progress: ", z*10, "%"))
    }
  }
  
  sigs <- data.frame(sigs)
  return(sigs)
    }
  }

  
mutect2  <- calculate.sigs(tumours = solid_tumours, pipe = "mutect2", path_to_MAFs = "~/4_pipelines_MAF_data/")
saveRDS(mutect2, file = "mutect2.rds")
