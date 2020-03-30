# ----------
# title: clinvar_cleaner.R
# author: Joshua D Green
# date: October 2019
# description: extracts gene-specific information from the processed clinvar data
# usage: 
# ----------

#packages------------------
library("rvest")
library("tidyverse")
library("stringr")
library("data.table")

processed_clinvar <- fread(file ="../2_Gene_Formatting/clinvar_ALS.csv")

# Helper functions---------
# function for parsing and extracting information from the processed clinvar data----------
clinvar_parse <- function(clinvar_df){
#extracting information from the Name column in dataframe
  clinvar_df$Protein.Consequence <-str_extract(string= clinvar_df$Name, pattern = "p\\....\\d*...")
  clinvar_df$Nucleotide.Consequence <-str_extract(string= clinvar_df$Name, pattern = "c.*>[A,T,G,C]")

# data extracted from the clinvar website
review.criteria<- c()

#for loop that: 
  #1) webscrapes review criteria from cv website, 
  #2) performs text manipulation on the scrapped data
  #3) adds each scrapped element to an empty vector
  for(i in 1:length(clinvar_df$VariationID)) {
    cv.url <- paste("https://www.ncbi.nlm.nih.gov/clinvar/variation/",clinvar_df$VariationID[i],"/",sep="")
    
    #loading bar
    print(paste(clinvar_df$Gene[i],i,"of",length(clinvar_df$VariationID),clinvar_df$Variation[i],cv.url,sep=" "))
    
    #1) webscrapes review criteria from cv website,
    review<- cv.url %>% read_html() %>% html_node(css = '.no-margin-top dd:nth-child(4)') %>% html_text()
    #2) performs text manipulation on the scrapped data 
    review_f <- gsub(pattern = "\\\n || \\s",replacement = "",x = review)
    #3) adds each scrapped element to an empty vector
    review.criteria[i] <- review_f
  }
# change the Review.Criteria
review.criteria_mod <- review.criteria %>% str_replace_all(c("practice guideline" = "****",
                                                             "reviewed by expert panel" = "***",
                                                             "criteria provided, multiple submitters, no conflicts" = "**",
                                                             "criteria provided, conflicting interpretations" = "*",
                                                             "criteria provided, single submitter" = "*",
                                                             "no assertion for the individual variant" = "none",
                                                             "no assertion criteria provided"="none",
                                                             "no assertion provided" = "none"))
# appends scrapped review criteria vector to the clinvar data frame
clinvar_df$Review.Criteria <- review.criteria_mod
clinvar_df <- as.data.frame(clinvar_df) 
}

# function for joining clinvar dataset with gnomad dataset--------------------------------
gnomad_join <- function(dataframe,gene){ #read in the gnomad gene of interest

#create a file path for each gene gnomad csv.file
file <- paste("../2_Gene_Formatting/gnomad_raw/",gene,"_gnomad.csv",sep="")

#import the gnomad csv file for the appropriate gene
gnomad <- as.data.frame(fread(file))

#make a unique identifier for the gnomad dataset
gnomad$ID <- paste(gnomad$Chromosome,gnomad$Position,sep = ":")

#select only needed columns and reorder so ID is column 1 for gnomad data
gnomad <- gnomad %>% select("ID","Position","rsID","Allele Frequency")

#join data frames by the ID column in the filtered clinvar data set and the gnomad dataset
left_join(dataframe, gnomad,by = "ID")
}

# function for exracting citations for each clinvar variant--------------------------------
citation_extractor <- function(gene){
  Pubmed_ID <- c()
  #iterate through the rsID
  for (x in 1:length(gene$rsID)) {
    df <- data.frame()
    print(x)
    # if the rsID is not available put in NA in the vector and   move onto the next 
    if (is.na(gene$rsID[x])==TRUE){
      Pubmed_ID[x]<- "N/A"
      print("next")
      next
    }
    # if the rsID is available make the url then try to webscrape. Will pass error
    else
      print("else")
    url <- paste("https://www.ncbi.nlm.nih.gov/snp/",gene$rsID[x],"#publications",sep="")
    try(df<- url %>% read_html() %>% html_node(xpath = '//*[@id="publication_datatable"]') %>% html_table())
    
    try(PMID<- as.vector(df$PMID))
    
    try(PMID_string<-paste(PMID, collapse = ', '))
    Pubmed_ID[x]<- paste(gene$rsID[x],PMID_string,sep=":")
    print(paste(x,url,gene$rsID[x],PMID_string,sep="   "))
  }
  #append to the end of the clinvar data frame
  gene<-cbind(gene,Pubmed_ID)
  return(gene)
}

# applying functions to each gene---------
format_by_gene<-function(gene){
gene_filtered_clinvar <- processed_clinvar %>% filter(Gene == as.character(gene))
gene_gnomadjoin <- gnomad_join(gene_filtered_clinvar,as.character(gene))
gene_cite <- citation_extractor(gene_gnomadjoin)
gene_parsed <- clinvar_parse(gene_cite)
gene_name <- gene_parsed %>% 
  select(VariationID,Name,Gene,GRCh38Location,"rsID","Allele Frequency",Clinical.Significance,Protein.Consequence,Nucleotide.Consequence,Review.Criteria,"Pubmed_ID") %>%
  rename("Position"="GRCh38Location","Allele.Frequency"="Allele Frequency")

# fix formatting in Protein.Consequence e.g., Before:'p.Ala430=)' After:'p.Ala430Ala')
gene_name$Protein.Consequence <- gsub(pattern = "p\\.(\\w{3})(\\d*)=\\)",replacement = "p\\.\\1\\2\\1",x = gene_name$Protein.Consequence)

# fix formatting in Protein.Consequence column (e.g., Before:'p.Ala430fs)', After:'p.Ala430fs')
gene_name$Protein.Consequence <- gsub(pattern = "(.*)fs\\)",replacement = "\\1fs",x = gene_name$Protein.Consequence) 

# write the csv file for each gene
write_csv(x=gene_name,path = paste("../clinvar_cache/",gene,"_clinvar.csv",sep=""))
}
start_time <- Sys.time()
format_by_gene("NEFH")
format_by_gene("HNRNPA2B1")
format_by_gene("NEK1")
format_by_gene("KIF5A")
format_by_gene("SPG11")
format_by_gene("UBQLN2")
format_by_gene("ANG")
format_by_gene("SIGMAR1")
format_by_gene("CHCHD10")
format_by_gene("VCP")
format_by_gene("C9orf72")
format_by_gene("SLC52A3")
format_by_gene("FIG4")
format_by_gene("SQSTM1")
format_by_gene("OPTN")
format_by_gene("MATR3")
format_by_gene("VAPB")
format_by_gene("SETX")
format_by_gene("ALS2")
format_by_gene("TBK1")
format_by_gene("TARDBP")
format_by_gene("DCTN1")
format_by_gene("FUS")
format_by_gene("SOD1")
end_time <- Sys.time()


print(end_time - start_time)