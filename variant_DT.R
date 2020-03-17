# ----------
# title: variant_DT.R
# author: Joshua D Green
# date: March 2020
# description: variant table helper function
# usage: 
# ----------
library(DT)
library(tidyverse)

tableman <- function(data2){

  #selecting the appropriate columns and then renaming them to remove the periods.
data2 <- as.data.frame(data2) %>% select("VariationID","Position","rsID","Clinical.Significance","Protein.Consequence","Nucleotide.Consequence","Allele.Frequency","Review.Criteria") %>% 
  
  #converts the allele frequency into scientific notation
  mutate("Allele.Frequency"=formatC(data2$"Allele.Frequency",format="e"))

#Links to webpages-------------------------------------------------
#iteratively replaces all the rsIDs with URL link and adds them to a new vector
rsID_link <- c()
for(i in c(1:length(data2$rsID))){
   link <-  str_replace(string = data2$rsID[i],pattern = data2$rsID[i],replacement = paste("<a href='https://www.ncbi.nlm.nih.gov/snp/",data2$rsID[i],"'>",data2$rsID[i],"</a>",sep=""))
   rsID_link[i] <- link
}

VariationID_link <- c()
for(i in c(1:length(data2$VariationID))){
   link <-  str_replace(string = as.character(data2$VariationID[i]),pattern = as.character(data2$VariationID[i]),replacement = paste("<a href='https://www.ncbi.nlm.nih.gov/clinvar/variation/",as.character(data2$VariationID[i]),"'>",as.character(data2$VariationID[i]),"</a>",sep=""))
   VariationID_link[i] <- link
}



#add the rsID_link to the data frame
data2$rsID <- rsID_link
data2$VariationID <- VariationID_link
#rename the columns and print out the data frame
DT::datatable(data = data2,escape = FALSE, options = list(pageLength = 1000,
                                                          autoWidth = TRUE
                                                          ),
              #renames the columns
              colnames = c("ClinVar ID","Position (GRCh38)","rsID","Clinical Significance","Protein Consequence","Nucleotide Consequence","Allele Frequency\n\n(gnomAD)","Review Criteria")
              # #adds buttons to the table so it can be downloaded in various file formats
               # extensions = 'Buttons', options = list(
               #   dom = 'Bfrtip',
               #   buttons = c('copy','csv')) 
)

  
}