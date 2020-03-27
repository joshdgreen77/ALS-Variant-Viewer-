# ----------
# title: Processing_Raw_ClinVar_Data.R
# author: Joshua D Green
# date: March 2019
# description: this workflow is used to process the raw clinvar data so that it can be read by the clinvar_cleaner.R app
# usage: 
# ----------


# import packages----------
library("tidyverse")
library("stringr")
library("data.table")


# read raw clinvar data----------
clinvar_raw <- as.data.frame(fread(file ="clinvar_ALS_raw.txt"))

# clean raw clinvar data----------
# select columns of interest
clinvar_selected <- clinvar_raw %>% 
  select("Name","VariationID","Gene(s)","Condition(s)","Clinical significance (Last reviewed)","GRCh37Chromosome","GRCh37Location","GRCh38Chromosome","GRCh38Location") 

# create a unique ID for clinvar dataset
clinvar_selected$ID <- paste(clinvar_selected$GRCh37Chromosome,clinvar_selected$GRCh37Location,sep=":")

# reorder and rename variables
clinvar_reorder_rename <- clinvar_selected %>% 
  select(ID,VariationID,1:9) %>% 
  rename(Gene = "Gene(s)",Condition = "Condition(s)",Clinical.Significance="Clinical significance (Last reviewed)")

# filter out structural variants
clinvar_no_structuralvariants <- clinvar_reorder_rename %>% 
  filter(grepl(pattern = "-",x = clinvar_reorder_rename$GRCh38Location) != TRUE)

# filter out blank columns in Position
clinvar_no_blanks <- clinvar_no_structuralvariants %>% 
  filter(grepl(pattern ="\\S",x = clinvar_no_structuralvariants$GRCh38Location) == TRUE)

# coerce "Gene" to factor and "GRCh38Location" to numeric
clinvar_no_blanks$Gene <- as.factor(clinvar_no_blanks$Gene)
clinvar_no_blanks$GRCh38Location <- as.numeric(clinvar_no_blanks$GRCh38Location)


# clean up the "Gene" column---------
# assign all variants to a single gene. clean up the gene column
clinvar_clean_gene <- clinvar_no_blanks # name new data frame


# OPTN
clinvar_clean_gene$Gene<- gsub(pattern = "LOC108903148\\|OPTN|OPTN\\|LOC108903148",replacement = "OPTN",x=clinvar_clean_gene$Gene)
# ANG
clinvar_clean_gene$Gene<- gsub(pattern = "ANG\\|RNASE4|RNASE4\\|ANG",replacement = "ANG",x=clinvar_clean_gene$Gene)
# CNTF
clinvar_clean_gene$Gene<- gsub(pattern = "CNTF\\|ZFP91-CNTF",replacement = "CNTF",x=clinvar_clean_gene$Gene)
#VCP
clinvar_clean_gene$Gene<- gsub(pattern = "FANCG\\|VCP",replacement = "VCP",x=clinvar_clean_gene$Gene)
#TARDBP
clinvar_clean_gene$Gene<- gsub(pattern = "MASP2\\|TARDBP",replacement = "TARDBP",x=clinvar_clean_gene$Gene)
#PRPH 
clinvar_clean_gene$Gene<- gsub(pattern = "PRPH\\|LOC101927267",replacement = "PRPH",x=clinvar_clean_gene$Gene)


# clean up the "Clinical.Significance" column---------
# Remove the date from the clinical significance column
clinvar_clean_gene$Clinical.Significance <- str_replace(string = clinvar_clean_gene$Clinical.Significance,pattern = "\\(.*\\)",replacement = "")

# Define column to sort clinical.significance column
clinsig_stratify <- function(gene){
  #stratify the different Clinical.Significance levels 1 of 5 categories
  gene <- gsub(pattern = "Conflicting interpretations of pathogenicity",replacement="Uncertain significance", x = gene) 
  gene <- gsub(pattern = "Pathogenic/Likely pathogenic",replacement = "Likely pathogenic",x=gene)
  gene <- gsub(pattern = "Benign/Likely benign",replacement = "Likely benign",x=gene)
  gene <- gsub(pattern = "risk factor",replacement = "Likely benign",x = gene)
}

# apply function to the data frame
clinvar_clean_gene$Clinical.Significance<-clinsig_stratify(clinvar_clean_gene$Clinical.Significance)

# write the cleaned up clinvar data to a file to be read by clinvar cleaner
write_csv(x = clinvar_clean_gene,path = "../2_Gene_Formatting/clinvar_ALS.csv")





# line of code to check that the gsub function is working
#before <- clinvar_clean_gene %>% group_by(Gene) %>% summarize(count=n())

