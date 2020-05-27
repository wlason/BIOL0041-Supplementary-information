library(stringr)
#Load libraries from the folder on Myriad
library(TCGAbiolinks, lib.loc = "~/R_packages")

## Define cancers to download
solid_tumours <- c("blca", "brca", "cesc", "chol", "coad", "esca",  "hnsc", "kich", "kirc", "kirp", "lihc", "luad", "lusc", "prad", "stad", "thca", "ucec")
solid_tumours <- str_to_upper(solid_tumours)

## Dwonload MAFs
for (tumour in solid_tumours){
  print(paste0(tumour, ": This is ", which(tumour == solid_tumours), " out of ", length(solid_tumours)))
  
      for(pipe in c("muse", "varscan2", "somaticsniper", "mutect2")){
      
      #Move to next pipeline/tumour on error
      skip <- FALSE
      df.maf <- tryCatch(GDCquery_Maf(tumour, 
                                      pipelines = pipe, 
                                      save.csv = TRUE,
                                      directory = "~/MAF_csv"), error = function(e) skip <- TRUE)
      if(skip && pipe == "mutect2"){
        next
      } else if (skip){
        break
      }
      # Add what cancer it is to the dataframe
      df.maf$Cancer <- as.character(tumour)
      name <- paste(tumour, pipe, "MAF_data", sep = "_")
      assign(name, value = df.maf)
      rm(df.maf)
      save(list = name, file = paste0(name, ".RData"))
      }
}

