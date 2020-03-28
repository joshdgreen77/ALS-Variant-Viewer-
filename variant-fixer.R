library(tidyverse)
library(stringr)
library(data.table)
setwd("/Users/greenjod/Documents/GitHub/ALS-Variant-Viewer-/clinvar_cache/")


args <- commandArgs(trailingOnly = FALSE)
print(args[6])

df <- fread(input = "ALS2_clinvar.csv")
df2 <- df %>% distinct(VariationID,.keep_all = TRUE)





# replaces all analogous of this p.Thy69=) with p.Thy69Thy
df$Protein.Consequence<- gsub(pattern = "p\\.(\\w{3})(\\d*)=\\)",replacement = "p\\.\\1\\2\\1",x = df$Protein.Consequence)

# replace fs) with fs
df$Protein.Consequence<-gsub(pattern = "(.*)fs\\)",replacement = "\\1fs",x = df$Protein.Consequence)

write_csv(x = df,file = args[6])


