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
clinvar.parse <- function(x){
#remove the date for clinical significance column
x$Clinical.Significance <- str_replace(string = x$Clinical.Significance,pattern = "\\(.*\\)",replacement = "")
  
#extracting information from the Name column in dataframe
  x$Protein.Consequence <-str_extract(string= x$Name, pattern = "p\\....\\d*...")
  x$Nucleotide.Consequence <-str_extract(string= x$Name, pattern = "c.*>[A,T,G,C]")

# data extracted from the clinvar website
review.criteria<- c()

#for loop that: 
  #1) webscrapes review criteria from cv website, 
  #2) performs text manipulation on the scrapped data
  #3) adds each scrapped element to an empty vector
  for(i in 1:length(x$VariationID)) {
    cv.url <- paste("https://www.ncbi.nlm.nih.gov/clinvar/variation/",x$VariationID[i],"/",sep="")
    
    #loading bar
    print(paste(i,"of",length(x$VariationID),x$Variation[i],cv.url,sep=" "))
    
    #1) webscrapes review criteria from cv website,
    review<- cv.url %>% read_html() %>% html_node(css = '.no-margin-top dd:nth-child(4)') %>% html_text()
    #2) performs text manipulation on the scrapped data 
    review_f <- gsub(pattern = "\\\n || \\s",replacement = "",x = review)
    #3) adds each scrapped element to an empty vector
    review.criteria[i] <- review_f
  }
  # appends scrapped review criteria vector to the clinvar data frame
  x$Review.Criteria <- review.criteria
  x <- as.data.frame(x)
}

# function for joining clinvar dataset with gnomad dataset--------------------------------
gnomad_join <- function(dataframe,gene){ #read in the gnomad gene of interest

#create a file path for each gene gnomad csv.file
file <- paste("../gnomad_raw/",gene,"_gnomad.csv",sep="")

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
gene_parsed <- clinvar.parse(gene_cite)
gene_name <- gene_parsed %>% 
  select(VariationID,Name,Gene,GRCh38Location,"rsID","Allele Frequency",Clinical.Significance,Protein.Consequence,Nucleotide.Consequence,Review.Criteria,"Pubmed_ID") %>%
  rename("Position"="GRCh38Location","Allele.Frequency"="Allele Frequency")

write_csv(x=gene_name,path = paste("../clinvar_cache/",gene,"_clinvar.csv",sep=""))
}

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
format_by_gene("TARDBP")
format_by_gene("DCT1N")
format_by_gene("FUS")
format_by_gene("SOD1")


#most recent at the top

# ##KIF51-->Generated
# KIF5A.cv <- cv %>% filter(Gene == "KIF5A")
# KIF5A.cv <- gnomad_join(KIF5A.cv,"KIF5A")
# KIF5A.cite  <- citation_extractor(KIF5A.cv)
# KIF5A.parsed <- clinvar.parse(KIF5A.cite)
# KIF5A <- KIF5A.parsed %>% 
#   select(VariationID,Name,Gene,GRCh38Location,"rsID","Allele Frequency",Clinical.Significance,Protein.Consequence,Nucleotide.Consequence,Review.Criteria,"Pubmed_ID") %>%
#   rename("Position"="GRCh38Location","Allele.Frequency"="Allele Frequency")
# write_csv(x= KIF5A,path = "~/Documents/GitHub/LNG_Scripts/als_app/als_variants_app/clinvar_cache/KIF5A_clinvar.csv")
# 
# ##SPG11-->Generated 
# SPG11.cv <- cv %>% filter(Gene == "SPG11")
# SPG11.cv <- gnomad_join(SPG11.cv,"SPG11")
# SPG11.cite  <- citation_extractor(SPG11.cv)
# SPG11.parsed <- clinvar.parse(SPG11.cite)
# SPG11 <- SPG11.parsed %>% 
#   select(VariationID,Name,Gene,GRCh38Location,"rsID","Allele Frequency",Clinical.Significance,Protein.Consequence,Nucleotide.Consequence,Review.Criteria,"Pubmed_ID") %>%
#   rename("Position"="GRCh38Location","Allele.Frequency"="Allele Frequency")
# write_csv(x= SPG11,path = "~/Documents/GitHub/LNG_Scripts/als_app/als_variants_app/clinvar_cache/SPG11_clinvar.csv")
# 
# 
# ##UBQLN2-->Generated 
# UBQLN2.cv <- cv %>% filter(Gene == "UBQLN2")
# UBQLN2.cv <- gnomad_join(UBQLN2.cv,"UBQLN2")
# UBQLN2.cite  <- citation_extractor(UBQLN2.cv)
# UBQLN2.parsed <- clinvar.parse(UBQLN2.cite)
# UBQLN2 <- UBQLN2.parsed %>% 
#   select(VariationID,Name,Gene,GRCh38Location,"rsID","Allele Frequency",Clinical.Significance,Protein.Consequence,Nucleotide.Consequence,Review.Criteria,Pubmed_ID) %>%
#   rename("Position"="GRCh38Location","Allele.Frequency"="Allele Frequency")
# write_csv(x= UBQLN2,path = "~/Documents/GitHub/LNG_Scripts/als_app/als_variants_app/clinvar_cache/UBQLN2_clinvar.csv")
# 
# 
# 
# ##ANG-->Generated 
# ANG.cv <- cv %>% filter(Gene == "ANG")
# ANG.cv <- gnomad_join(ANG.cv,"ANG")
# ANG.cite  <- citation_extractor(ANG.cv)
# ANG.parsed <- clinvar.parse(ANG.cite)
# ANG <- ANG.parsed %>% 
#   select(VariationID,Name,Gene,GRCh38Location,"rsID","Allele Frequency",Clinical.Significance,Protein.Consequence,Nucleotide.Consequence,Review.Criteria,Pubmed_ID) %>%
#   rename("Position"="GRCh38Location","Allele.Frequency"="Allele Frequency")
# write_csv(x= ANG,path = "~/Documents/GitHub/LNG_Scripts/als_app/als_variants_app/clinvar_cache/ANG_clinvar.csv")
# 
# ##SIGMAR1-->Generated 
# SIGMAR1.cv <- cv %>% filter(Gene == "SIGMAR1")
# SIGMAR1.cv <- gnomad_join(SIGMAR1.cv,"SIGMAR1")
# SIGMAR1.cite  <- citation_extractor(SIGMAR1.cv)
# SIGMAR1.parsed <- clinvar.parse(SIGMAR1.cite)
# SIGMAR1 <- SIGMAR1.parsed %>% 
#   select(VariationID,Name,Gene,GRCh38Location,"rsID","Allele Frequency",Clinical.Significance,Protein.Consequence,Nucleotide.Consequence,Review.Criteria,Pubmed_ID) %>%
#   rename("Position"="GRCh38Location","Allele.Frequency"="Allele Frequency")
# write_csv(x= SIGMAR1,path = "~/Documents/GitHub/LNG_Scripts/als_app/als_variants_app/clinvar_cache/SIGMAR1_clinvar.csv")
# 
# ##CHCHD10-->Generated 
# CHCHD10.cv <- cv %>% filter(Gene == "CHCHD10")
# CHCHD10.cv <- gnomad_join(CHCHD10.cv,"CHCHD10")
# CHCHD10.cite  <- citation_extractor(CHCHD10.cv)
# CHCHD10.parsed <- clinvar.parse(CHCHD10.cite)
# CHCHD10 <- CHCHD10.parsed %>% 
#   select(VariationID,Name,Gene,GRCh38Location,"rsID","Allele Frequency",Clinical.Significance,Protein.Consequence,Nucleotide.Consequence,Review.Criteria,Pubmed_ID) %>%
#   rename("Position"="GRCh38Location","Allele.Frequency"="Allele Frequency")
# write_csv(x= CHCHD10,path = "~/Documents/GitHub/LNG_Scripts/als_app/als_variants_app/clinvar_cache/CHCHD10_clinvar.csv")
# 
# 
# ##TBK1-->Generated 
# TBK1.cv <- cv %>% filter(Gene == "TBK1")
# TBK1.cv <- gnomad_join(TBK1.cv,"TBK1")
# TBK1.cite  <- citation_extractor(TBK1.cv)
# TBK1.parsed <- clinvar.parse(TBK1.cite)
# TBK1 <- TBK1.parsed %>% 
#   select(VariationID,Name,Gene,GRCh38Location,"rsID","Allele Frequency",Clinical.Significance,Protein.Consequence,Nucleotide.Consequence,Review.Criteria,Pubmed_ID) %>%
#   rename("Position"="GRCh38Location","Allele.Frequency"="Allele Frequency")
# write_csv(x= TBK1,path = "~/Documents/GitHub/LNG_Scripts/als_app/als_variants_app/clinvar_cache/TBK1_clinvar.csv")
# 
# ##VCP-->Generated 
# VCP.cv <- cv %>% filter(Gene == "VCP")
# VCP.cv <- gnomad_join(VCP.cv,"VCP")
# VCP.cite  <- citation_extractor(VCP.cv)
# VCP.parsed <- clinvar.parse(VCP.cite)
# VCP <- VCP.parsed %>% 
#   select(VariationID,Name,Gene,GRCh38Location,"rsID","Allele Frequency",Clinical.Significance,Protein.Consequence,Nucleotide.Consequence,Review.Criteria,Pubmed_ID) %>%
#   rename("Position"="GRCh38Location","Allele.Frequency"="Allele Frequency")
# write_csv(x= VCP,path = "~/Documents/GitHub/LNG_Scripts/als_app/als_variants_app/clinvar_cache/VCP_clinvar.csv")
# 
# 
# ##C9orf72-->Generated 
# C9orf72.cv <- cv %>% filter(Gene == "C9orf72")
# C9orf72.cv <- gnomad_join(C9orf72.cv,"C9orf72")
# C9orf72.cite  <- citation_extractor(C9orf72.cv)
# C9orf72.parsed <- clinvar.parse(C9orf72.cite)
# C9orf72 <- C9orf72.parsed %>% 
#   select(VariationID,Name,Gene,GRCh38Location,"rsID","Allele Frequency",Clinical.Significance,Protein.Consequence,Nucleotide.Consequence,Review.Criteria,Pubmed_ID) %>%
#   rename("Position"="GRCh38Location","Allele.Frequency"="Allele Frequency")
# write_csv(x= C9orf72,path = "~/Documents/GitHub/LNG_Scripts/als_app/als_variants_app/clinvar_cache/C9orf72_clinvar.csv")
# 
# ##SLC52A3-->Generated 
# SLC52A3.cv <- cv %>% filter(Gene == "SLC52A3")
# SLC52A3.cv <- gnomad_join(SLC52A3.cv,"SLC52A3")
# SLC52A3.cite  <- citation_extractor(SLC52A3.cv)
# SLC52A3.parsed <- clinvar.parse(SLC52A3.cite)
# SLC52A3 <- SLC52A3.parsed %>% 
#   select(VariationID,Name,Gene,GRCh38Location,"rsID","Allele Frequency",Clinical.Significance,Protein.Consequence,Nucleotide.Consequence,Review.Criteria,Pubmed_ID) %>%
#   rename("Position"="GRCh38Location","Allele.Frequency"="Allele Frequency")
# write_csv(x= SLC52A3,path = "~/Documents/GitHub/LNG_Scripts/als_app/als_variants_app/clinvar_cache/SLC52A3_clinvar.csv")
# 
# 
# ##FIG4-->Generated 
# FIG4.cv <- cv %>% filter(Gene == "FIG4")
# FIG4.cv <- gnomad_join(FIG4.cv,"FIG4")
# FIG4.cite  <- citation_extractor(FIG4.cv)
# FIG4.parsed <- clinvar.parse(FIG4.cite)
# FIG4 <- FIG4.parsed %>% 
#   select(VariationID,Name,Gene,GRCh38Location,"rsID","Allele Frequency",Clinical.Significance,Protein.Consequence,Nucleotide.Consequence,Review.Criteria,Pubmed_ID) %>%
#   rename("Position"="GRCh38Location","Allele.Frequency"="Allele Frequency")
# write_csv(x= FIG4,path = "~/Documents/GitHub/LNG_Scripts/als_app/als_variants_app/clinvar_cache/FIG4_clinvar.csv")
# 
# ##SQSTM1-->Generated 
# SQSTM1.cv <- cv %>% filter(Gene == "SQSTM1")
# SQSTM1.cv <- gnomad_join(SQSTM1.cv,"SQSTM1")
# SQSTM1.cite  <- citation_extractor(SQSTM1.cv)
# SQSTM1.parsed <- clinvar.parse(SQSTM1.cite)
# SQSTM1 <- SQSTM1.parsed %>% 
#   select(VariationID,Name,Gene,GRCh38Location,"rsID","Allele Frequency",Clinical.Significance,Protein.Consequence,Nucleotide.Consequence,Review.Criteria,Pubmed_ID) %>%
#   rename("Position"="GRCh38Location","Allele.Frequency"="Allele Frequency")
# write_csv(x= SQSTM1,path = "~/Documents/GitHub/LNG_Scripts/als_app/als_variants_app/clinvar_cache/SQSTM1_clinvar.csv")
# 
# ##OPTN-->Generated 
# OPTN.cv <- cv %>% filter(Gene == "OPTN")
# OPTN.cv <- gnomad_join(OPTN.cv,"OPTN")
# OPTN.cite  <- citation_extractor(OPTN.cv)
# OPTN.parsed <- clinvar.parse(OPTN.cite)
# OPTN <- OPTN.parsed %>% 
#   select(VariationID,Name,Gene,GRCh38Location,"rsID","Allele Frequency",Clinical.Significance,Protein.Consequence,Nucleotide.Consequence,Review.Criteria,Pubmed_ID) %>%
#   rename("Position"="GRCh38Location","Allele.Frequency"="Allele Frequency")
# write_csv(x= OPTN,path = "~/Documents/GitHub/LNG_Scripts/als_app/als_variants_app/clinvar_cache/OPTN_clinvar.csv")
# 
# ##MATR3-->Generated 
# MATR3.cv <- cv %>% filter(Gene == "MATR3")
# MATR3.cv <- gnomad_join(MATR3.cv,"MATR3")
# MATR3.cite  <- citation_extractor(MATR3.cv)
# MATR3.parsed <- clinvar.parse(MATR3.cite)
# MATR3 <- MATR3.parsed %>% 
#   select(VariationID,Name,Gene,GRCh38Location,"rsID","Allele Frequency",Clinical.Significance,Protein.Consequence,Nucleotide.Consequence,Review.Criteria,Pubmed_ID) %>%
#   rename("Position"="GRCh38Location","Allele.Frequency"="Allele Frequency")
# write_csv(x= MATR3,path = "~/Documents/GitHub/LNG_Scripts/als_app/als_variants_app/clinvar_cache/MATR3_clinvar.csv")
# 
# ##VAPB-->Generated 
# VAPB.cv <- cv %>% filter(Gene == "VAPB")
# VAPB.cv <- gnomad_join(VAPB.cv,"VAPB")
# VAPB.cite  <- citation_extractor(VAPB.cv)
# VAPB.parsed <- clinvar.parse(VAPB.cite)
# VAPB <- VAPB.parsed %>% 
#   select(VariationID,Name,Gene,GRCh38Location,"rsID","Allele Frequency",Clinical.Significance,Protein.Consequence,Nucleotide.Consequence,Review.Criteria,Pubmed_ID) %>%
#   rename("Position"="GRCh38Location","Allele.Frequency"="Allele Frequency")
# write_csv(x= VAPB,path = "~/Documents/GitHub/LNG_Scripts/als_app/als_variants_app/clinvar_cache/VAPB_clinvar.csv")
# 
# ##SETX-->Generated 
# SETX.cv <- cv %>% filter(Gene == "SETX")
# SETX.cv <- gnomad_join(SETX.cv,"SETX")
# SETX.cite  <- citation_extractor(SETX.cv)
# SETX.parsed <- clinvar.parse(SETX.cite)
# SETX <- SETX.parsed %>% 
#   select(VariationID,Name,Gene,GRCh38Location,"rsID","Allele Frequency",Clinical.Significance,Protein.Consequence,Nucleotide.Consequence,Review.Criteria,Pubmed_ID) %>%
#   rename("Position"="GRCh38Location","Allele.Frequency"="Allele Frequency")
# write_csv(x= SETX,path = "~/Documents/GitHub/LNG_Scripts/als_app/als_variants_app/clinvar_cache/SETX_clinvar.csv")
# 
# ##ALS2-->Generated 
# ALS2.cv <- cv %>% filter(Gene == "ALS2")
# ALS2.cv <- gnomad_join(ALS2.cv,"ALS2")
# ALS2.cite  <- citation_extractor(ALS2.cv)
# ALS2.parsed <- clinvar.parse(ALS2.cite)
# ALS2 <- ALS2.parsed %>% 
#   select(VariationID,Name,Gene,GRCh38Location,"rsID","Allele Frequency",Clinical.Significance,Protein.Consequence,Nucleotide.Consequence,Review.Criteria,Pubmed_ID) %>%
#   rename("Position"="GRCh38Location","Allele.Frequency"="Allele Frequency")
# write_csv(x= ALS2,path = "~/Documents/GitHub/LNG_Scripts/als_app/als_variants_app/clinvar_cache/ALS2_clinvar.csv")
# 
# ##TARDBP-->Generated 
# TARDBP.cv <- cv %>% filter(Gene == "TARDBP")
# TARDBP.cv <- gnomad_join(TARDBP.cv,"TARDBP")
# TARDBP.cite  <- citation_extractor(TARDBP.cv)
# TARDBP.parsed <- clinvar.parse(TARDBP.cite)
# TARDBP <- TARDBP.parsed %>% 
#   select(VariationID,Name,Gene,GRCh38Location,"rsID","Allele Frequency",Clinical.Significance,Protein.Consequence,Nucleotide.Consequence,Review.Criteria,Pubmed_ID) %>%
#   rename("Position"="GRCh38Location","Allele.Frequency"="Allele Frequency")
# write_csv(x= TARDBP,path = "~/Documents/GitHub/LNG_Scripts/als_app/als_variants_app/clinvar_cache/TARDBP_clinvar.csv")
# 
# ##DCT1N-->Generated 
# DCTN1.cv <- cv %>% filter(Gene == "DCTN1")
# DCTN1.cv <- gnomad_join(DCTN1.cv,"DCTN1")
# DCTN1.cite  <- citation_extractor(DCTN1.cv)
# DCTN1.parsed <- clinvar.parse(DCTN1.cite)
# DCTN1 <- DCTN1.parsed %>% 
#   select(VariationID,Name,Gene,GRCh38Location,"rsID","Allele Frequency",Clinical.Significance,Protein.Consequence,Nucleotide.Consequence,Review.Criteria,Pubmed_ID) %>%
#   rename("Position"="GRCh38Location","Allele.Frequency"="Allele Frequency")
# write_csv(x= DCTN1,path = "~/Documents/GitHub/LNG_Scripts/als_app/als_variants_app/clinvar_cache/DCTN1_clinvar.csv")
# 
# ##FUS-->Generated 
# FUS.cv <- cv %>% filter(Gene == "FUS")
# FUS.cv <- gnomad_join(FUS.cv,"FUS")
# FUS.cite  <- citation_extractor(FUS.cv)
# FUS.parsed <- clinvar.parse(FUS.cite)
# FUS <- FUS.parsed %>% 
#   select(VariationID,Name,Gene,GRCh38Location,"rsID","Allele Frequency",Clinical.Significance,Protein.Consequence,Nucleotide.Consequence,Review.Criteria,Pubmed_ID) %>%
#   rename("Position"="GRCh38Location","Allele.Frequency"="Allele Frequency")
# write_csv(x= FUS,path = "~/Documents/GitHub/LNG_Scripts/als_app/als_variants_app/clinvar_cache/FUS_clinvar.csv")
# 
# ##SOD1-->Generated 
# SOD1.cv <- cv %>% filter(Gene == "SOD1")
# SOD1.cv <- gnomad_join(SOD1.cv,"SOD1")
# SOD1.cite  <- citation_extractor(SOD1.cv)
# SOD1.parsed <- clinvar.parse(SOD1.cite)
# SOD1 <- SOD1.parsed %>% 
#   select(VariationID,Name,Gene,GRCh38Location,"rsID","Allele Frequency",Clinical.Significance,Protein.Consequence,Nucleotide.Consequence,Review.Criteria,Pubmed_ID) %>%
#   rename("Position"="GRCh38Location","Allele.Frequency"="Allele Frequency")
# write_csv(x= SOD1,path = "~/Documents/GitHub/LNG_Scripts/als_app/als_variants_app/clinvar_cache/SOD1_clinvar.csv")















#############
#####JUNK####
#############

# #Importing clinvar data--------------------------------------------------
# cv_raw <- as.data.frame(fread(file = "~/Documents/GitHub/LNG_Scripts/als_app/als_variants_app/clinvar_ALS.txt"))
# 
# #clinvar data clean up-----------------------------------------------------------
# cv <- cv_raw %>% 
#   select("Name","VariationID","Gene(s)","Condition(s)","Clinical significance (Last reviewed)","GRCh37Chromosome","GRCh37Location","GRCh38Chromosome","GRCh38Location") 
# #create ID for clinvar dataset
# cv$ID <- paste(cv$GRCh37Chromosome,cv$GRCh37Location,sep=":")
# 
# #reorder and rename variables
# cv <- cv %>% select(ID,VariationID,1:9) %>% rename(Gene = "Gene(s)",Condition = "Condition(s)",Clinical.Significance="Clinical significance (Last reviewed)")
# 
# #Filter out structural variants
# cv <- cv %>% 
#   filter(grepl(pattern = "-",x = cv$GRCh38Location) != TRUE)
# 
# #Filter out blank columns in Position
# cv <- cv %>% 
#   filter(grepl(pattern ="\\S",x = cv$GRCh38Location) == TRUE)
# 
# 
# #Coerce Gene and Position to appropriate class
# cv$Gene <- as.factor(cv$Gene)
# cv$GRCh38Location <- as.numeric(cv$GRCh38Location)


