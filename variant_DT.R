# ----------
# title: variant_DT.R
# author: Joshua D Green
# date: March 2020
# description: variant table helper function
# usage: 
# ----------
library(DT)
library(tidyverse)
#source("~/Documents/GitHub/LNG_Scripts/als_app/als_variants_app/data_import.R")

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
DT::datatable(data = data2,escape = FALSE,
              #renames the columns
              colnames = c("ClinVar ID","Position (GRCh38)","rsID","Clinical Significance","Protein Consequence","Nucleotide Consequence","Allele Frequency\n\n(gnomAD)","Review Criteria"),
               #adds buttons to the table so it can be downloaded in various file formats
                extensions = 'Buttons',  options = list(pageLength = 1000,
                                                        autoWidth = TRUE,
                                                        dom = 'Bfrtip',
                                                        buttons = c('csv'))
)

  
}





##########JUNK###################
#this has been used when you add colors to the data table
#conditional statment for the clinical significance color----------
# if (input$select == "DCTN1"){ #Benign, Likely benign, pathogenic, uncertain significance
#   Clinical.Significance.Colors <<- c("#addda8","#fffdc4","#d41e25","#3483b7") 
# }else if (input$select == "VAPB"){ #Benign, Likely benign, pathogenic, uncertain significance
#   Clinical.Significance.Colors <<- c("#addda8","#fffdc4","#d41e25","#3483b7") 
# }else if (input$select == "TARDBP") { #Likely benign, Likely pathogenic, Pathogenic, Uncertain significance
#   Clinical.Significance.Colors <<- c("#fffdc4","#faae6a","#d41e25","#3483b7")
# } else { #Benign, Likely Benign, Likely pathogenic, Pathogenic, Uncertain Significance
#   Clinical.Significance.Colors <<- c("#addda8","#fffdc4","#faae6a","#d41e25","#3483b7")
# }
# 
# 
# #conditional statement for the review status color----------
# if (input$select == "FUS"){ #multiple submitters, single submitter, no assertion
#   Review.Criteria.Colors <<- c("#FFDD3C","#FFEA61","#FFFFB7")
# }else {#conflicting interpretations,multiple submitters, no conflicts,single submitter,no assertion 
#   Review.Criteria.Colors <<- c("#FFF192","#FFDD3C","#FFEA61","#FFFFB7")
# }
  
  # #adds colors to the clinical significance column
  #formatStyle(columns = "Clinical Significance",backgroundColor = styleEqual(sort(unique(data2$"Clinical Significance")),Clinical.Significance.Colors)) 
#%>%
  # #adds color to the review criteria column
  # formatStyle(columns = "Review Criteria",
  #             backgroundColor = styleEqual(sort(unique(data2$"Review Criteria")), Review.Criteria.Colors)) %>%
  # #adds a bar graph to the data table to visualize allele frequency
  # formatStyle(columns = "Allele Frequency (gnomAD)",
  #             background = styleColorBar(data2$"Allele Frequency (gnomAD)",'steelblue'),
  #             backgroundSize = 100,
  #             backgroundRepeat = 'repeat',
  #             backgroundPostion = 'center')


#rename("Position (GRCh38)"="Position", "Clinical Significance"="Clinical.Significance", "Protein Consequence" = "Protein.Consequence", "Nucleotide Consequence"="Nucleotide.Consequence","Allele Frequency (gnomAD)"="Allele.Frequency","Review Criteria"="Review.Criteria")
