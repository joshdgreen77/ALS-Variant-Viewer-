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
