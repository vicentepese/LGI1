# Import libraries
library(jsonlite)
library(tidyverse)
library(readr)
library(data.table)
library(ggrepel)
library(viridis)
library(hrbrthemes)
library(HIBAG)
library(parallel)

########## IMPORT ##########
setwd("~/Documents/LGI1")

# Import settings
settings <- jsonlite::fromJSON('settings.json')

########## PREFIT MODELS ############

# Load pre-fit model and comvert to hlaMODEL
model.list <- get(load(settings$file$PMRA_HLA_model))
hla.id <- names(model.list)

# Import file
gname <- settings$file$filtLGI1
yourgeno <- hlaBED2Geno(bed.fn=paste(gname, ".bed", sep = ''), fam.fn=paste(gname, ".fam", sep='')
                        , bim.fn=paste(gname, ".bim", sep=''), assembly = 'hg19')
summary(yourgeno)

# Make cluster 
cl <- makeCluster(10)

# Make predictions
for (locus in hla.id){
  model.hla <- hlaModelFromObj(model.list[[locus]])
  summary(model.hla)
  pred.guess <- predict(model.hla, yourgeno, type="response+prob", nclassifiers=100, cl=cl)
  save(pred.guess, file = paste("./Data/HLA_predictions/", paste('HLA', locus, sep = ''), '.RData', sep= ''))
}


