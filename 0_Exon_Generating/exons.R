#make pipline for comingn the exon coordinates into one single csv file
library(tidyverse)
library(data.table)
library(stringr)
library(ggplot2)
library(plotly)
setwd(dir = "~/Documents/GitHub/ALS-Variant-Viewer-/0_Exon_Generating/exons/")

# make vector with the names of all the files in the folder
files <- as.vector(system("ls -1 *.txt",intern=TRUE))

# make an empty data frame that will be appended at the end of the for loop 
exons_df <- data.frame("gene_name"=c(),"start"=c(),"stop"=c())

# for loop iterates through each element in the files vector
for (i in c(1:length(files))) {
  
  # fread the file
  file <- as.data.frame(fread(file = files[i]))  
  
  # obtain gene name from the file
  gene_name<- as.vector(str_replace(string = files[i], pattern = '_exons.txt',replacement = '')) 
  
  # retrieve the first row and only two columns
  exons <- file[1,] %>% select("exonStarts","exonEnds")
  
  # bind the gene_name vector to the exons start and stop vector--> creates a data frame
  gene_exon<- cbind(gene_name,exons)
  
  # append data frame to empty data frame 
  exons_df<- rbind(exons_df,gene_exon)
}

write_csv(exons_df,path = "~/Documents/GitHub/ALS-Variant-Viewer-/exons.csv")





