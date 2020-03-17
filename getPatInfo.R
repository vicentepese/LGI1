# Import libraries
library(xlsx)
library(readxl)
library(jsonlite)
library(dplyr)
library(stringr)

############ INITIALIZATION ############ 

# Set working directory 
setwd("~/Documents/LGI1/")

# Import settings
settings <- jsonlite::read_json('settings.json')

# Import data
paraneo.data <- read_xlsx(settings$file$GWAS_total)
table(paraneo.data$Dx)
unique(paraneo.data$Dx)
table(paraneo.data$`Dx-ling`)
unique(paraneo.data$`Dx-ling`)

########### GET PLATES CONTROLS ###############

# Get controls data and count 
ctrls.dx <- c("NMDA control", "HIMC", "F or G3", "G3", "G3 or J", "sleep study", "k", "APOE study", "control", "K", 
  "Control", "Healthy Control", "Normal", "PSG Sstudy", "PSG study", "normal")
ctrls.data <- paraneo.data[which(paraneo.data$Dx %in% ctrls.dx),]
table(ctrls.data$Dx)

# Get plates
ctrls.plates <- sapply(ctrls.data$`GWAS ID`, FUN = filterplates)
table(ctrls.plates)
ctrls.plates <- as.numeric(unique(ctrls.plates))
ctrls.plates <- ctrls.plates[which(ctrls.plates >= 77)]
ctrls.plates

# Get control files from second batch (plates 120 and 121)
kls.data <- read_xlsx(settings$file$GWAS_KLS)
kls.GWASID <- kls.data$`GWAS ID...5` %>% unique()

# Write 
write.table(as.numeric(ctrls.plates), file = settings$file$plates, row.names = FALSE, 
          col.names = FALSE, sep = ',')

# Get GWAS ID
GWAS.ID <- ctrls.data$`GWAS ID` %>% unique()
GWAS.ID <- c(GWAS.ID, kls.GWASID)

# Write
write.table(GWAS.ID, file = settings$file$GWASIDsControls, row.names =  FALSE,
            col.names = FALSE, sep = ',')

########### GET PLATES CASES ##########

# Get cases data 
paraneo.cases <- read_xlsx(settings$file$paraneoplastic)
table(paraneo.cases$Dx)
cases.data <- paraneo.cases[which(paraneo.cases$Dx == 'LGI1'),]

# Get GWAS ID
GWAS.ID.cases <- cases.data$`GWAS ID` %>% unique()

# Write 
write.table(GWAS.ID.cases, file = settings$file$GWASIDsCases, row.names = FALSE,
            col.names = FALSE, sep = ',')


########### Functions #############

filterplates <- function(GWASID){
  return(strsplit(gsub("([0-9]*)([A-Z]*)([0-9]*)", "\\1 \\2 \\3", GWASID), " ")[[1]][1])
}

