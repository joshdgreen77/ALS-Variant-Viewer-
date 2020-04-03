#make pipline for comingn the exon coordinates into one single csv file
library(tidyverse)
library(data.table)
library(stringr)
library(ggplot2)
library(plotly)
setwd(dir = "~/Documents/GitHub/ALS-Variant-Viewer-/exons/")

# make vector with the names of all the files in the folder
files <- as.vector(system("ls -1 *.txt",intern=TRUE))

# make an empty data frame that will be appended at the end of the for loop 
exons_df <- tibble("gene_name"=c(),"start"=c(),"stop"=c())

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






















## WORKING FUNCTION##
variant_graph <- function(x){
  x$ClinVar.Stars <- sapply(as.character(x$Review.Criteria), function (x) {
    switch(x,
           "practice guideline" = "****",
           "reviewed by expert panel" = "***",
           "criteria provided, multiple submitters, no conflicts" = "**",
           "criteria provided, conflicting interpretations" = "*",
           "criteria provided, single submitter" = "*",
           "no assertion for the individual variant" = "none",
           "no assertion criteria provided" = "none",
           "no assertion provided" = "none")
  })
  
  
  #graphing----------
  p<- ggplot(x, aes(x=Position, y= Clinical.Significance, text1 = Nucleotide.Consequence, text2= Protein.Consequence,text3= rsID, text4 =ClinVar.Stars))+
    
    #prevents overplotting
    geom_jitter(position = position_jitter(width = 500,height = 0,seed = 1),mapping= aes(fill= Clinical.Significance),shape = 21,size = 4,color = "black", stroke = 1,show.legend = TRUE)+
    theme(panel.background = element_blank(),
          axis.text.x = element_text(angle=70,vjust = 0.5, hjust = 0.7),
          axis.title.x.bottom = element_text(),
          axis.line.x = element_line(size = 6 ,colour = "grey",lineend = "square"),
          #y elements-----
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.title.y = element_blank(),
          #legend elements-----
          legend.background = element_rect(fill="white"),
          legend.title = element_blank()
    )+
    labs(x = "Position (GRCh38)")+
    #assign color to each value
    scale_fill_manual(limits = c("Benign","Likely benign","Likely pathogenic","Pathogenic","Uncertain significance"),values = c("#ADDDA8","#FFFFC4","#FAAE6A","#D41B25","#3483B7"))+
    #ordering the y discrete values on the y axis
    scale_y_discrete(limits = c("0","Uncertain significance","Pathogenic","Likely pathogenic","Likely benign","Benign")) +
annotate(geom="rect", xmin = exons_start, xmax = exons_stop, ymin = 0, ymax = 1, color = "pink",fill = "pink")
  
  #convert ggplot to plotly object----------
  ggplotly(p,tooltip = c("x","fill","text1","text2","text3")) %>% layout(legend = list(orientation = "h",x=0,y=-1))
}


variant_graph(sod1)




# file <- as.data.frame(fread(file = "~/Desktop/exons/sod1_exons.txt"))
# gene_name <- c("sod1")
# exons <- file[1,] %>% select("exonStarts","exonEnds")
# cbind(gene_name,exons)
# 
# 
# sod1 <- fread(file= "~/Desktop/SOD1_clinvar.csv")
# 
# # select first row 
# exons_formatted <- exons %>% select("exonStarts","exonEnds")
# 
# file <- as.data.frame(fread(x))
# gene_name<- str_replace(string = files, pattern = '_exons.txt',replacement = '')
# mod_file <- file[1] %>% select(gene_name,"exonStarts","exonEnds")
# 
# 
# exons_start<-as.numeric(strsplit(x = exons_formatted$exonStarts, split = ",")[[1]])
# exons_stop <- as.numeric(strsplit(x = exons_formatted$exonEnds, split = ",")[[1]])









