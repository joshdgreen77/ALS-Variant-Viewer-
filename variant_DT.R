# ----------
# title: variant_DT.R
# author: Joshua D Green
# date: March 2020
# description: format_variants_dataframe helper function
# usage: 
# ----------

library(DT)
library(tidyverse)


format_variant_table <- function(variants_dataframe){
  
# Format the data frame----------
  variants_dataframe <- as.data.frame(variants_dataframe) %>% 
  # select columns to display
  select("Clinical.Significance","VariationID","Position","rsID","Protein.Consequence",
         "Nucleotide.Consequence","Allele.Frequency","Review.Criteria") %>%
    
  # convert allele frequency to scientific notation
  mutate("Allele.Frequency"=formatC(variants_dataframe$"Allele.Frequency",format="e"))

# generate the links for the rsID in the data table
rsID_link <- c()
for(i in c(1:length(variants_dataframe$rsID))){
   link <-  str_replace(string = variants_dataframe$rsID[i],pattern = variants_dataframe$rsID[i],replacement = paste("<a href='https://www.ncbi.nlm.nih.gov/snp/",variants_dataframe$rsID[i],"'>",variants_dataframe$rsID[i],"</a>",sep=""))
   rsID_link[i] <- link}

# generate the links to the variant info page on the ClinVar site
VariationID_link <- c()
for(i in c(1:length(variants_dataframe$VariationID))){
   link <-  str_replace(string = as.character(variants_dataframe$VariationID[i]),pattern = as.character(variants_dataframe$VariationID[i]),replacement = paste("<a href='https://www.ncbi.nlm.nih.gov/clinvar/variation/",as.character(variants_dataframe$VariationID[i]),"'>",as.character(variants_dataframe$VariationID[i]),"</a>",sep=""))
   VariationID_link[i] <- link
}
# add the rsID links and the variantID links to the dataframe that will be displayed
variants_dataframe$rsID <- rsID_link
variants_dataframe$VariationID <- VariationID_link

# defines the order to display variants by default
variants_dataframe<- variants_dataframe %>% arrange(factor(x = variants_dataframe$Clinical.Significance,levels =c("Pathogenic","Likely pathogenic","Likely benign","Benign","Uncertain significance")))

# Make Data Table object----------
# rename the columns and print out the data frame
DT::datatable(data = variants_dataframe,escape = FALSE,
              #renames the columns
              colnames = c("Clinical Significance","ClinVar ID","Position (GRCh38)","rsID","Protein Consequence","Nucleotide Consequence","Allele Frequency (gnomAD)","ClinVar Review Status"),
               #adds buttons to the table so it can be downloaded in various file formats
                extensions = c('Buttons','RowGroup'),  options = list(pageLength = 1000, 
                                                        autoWidth = TRUE,
                                                        #add button for downloading
                                                        dom = 'Bfrtip', 
                                                        #the button to add for downloading
                                                        buttons = c('csv'),
                                                        searchHighlight = TRUE,
                                                        #group rows by clinical significance
                                                        rowGroup = list(dataSrc = 1) 
                                                        ))

  
}