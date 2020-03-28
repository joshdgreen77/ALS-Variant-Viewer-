library(tidyverse)
library(stringr)
library(data.table)
setwd("/Users/greenjod/Documents/GitHub/ALS-Variant-Viewer-/clinvar_cache/")
args <- commandArgs(trailingOnly = FALSE)
print(args[6])

df <- fread(input = as.character(args[6]))

# replaces all analogous of this p.Thy69=) with p.Thy69Thy
df$Protein.Consequence<- gsub(pattern = "p\\.(\\w{3})(\\d*)=\\)",replacement = "p\\.\\1\\2\\1",x = df$Protein.Consequence)

# replace fs) with fs
df$Protein.Consequence<-gsub(pattern = "(.*)fs\\)",replacement = "\\1fs",x = df$Protein.Consequence)

write.csv(x = df,file = args[6])


