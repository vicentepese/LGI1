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

########## COLLAPSE #########
path <- settings$directory$HLA_predictions
pred.files <- list.files(settings$directory$HLA_predictions)

# Load first file and get list of samples
init.data <- get(load(paste(path, pred.files[1], sep = '')))
sample.id <- pred.guess$value$sample.id
data <- data.frame('FID' = sample.id, 'IID' = sample.id)

# Collapse
for (file in pred.files){
  temp <- get(load(paste(path, file, sep = '')))
  name <- strsplit(file, '\\.')[[1]]
  name <- strsplit(name, '\\-')[[1]]
  name <- paste(name[1], name[2], sep = '')
  data[paste(name, 'allele1', sep = '_')] <- temp$value$allele1
  data[paste(name, 'allele2',sep = '_')] <- temp$value$allele2
}

# Write 
write.table(data, file = settings$file$HLA_total, quote = FALSE, row.names = FALSE, col.names = TRUE, sep = '\t')

############ MERGE #############

# Read eigenvectors
eigenvec <- read.table(settings$file$finalPCAeigenvectors, header = TRUE, sep = "\t", quote = '')

# Merge 
covars <- merge(eigenvec, data, by = c('FID', 'IID'))
covars.names <- colnames(covars)
covars.names <- paste(covars.names, collapse = ', ')
# Write 
write.table(covars, file = settings$file$covariates, quote = FALSE, sep = '\t', row.names = FALSE)
write(covars.names, file = 'Resources/covars_names.txt')

covars.names

