# ----------
# title: data_import.R
# author: Joshua D Green
# date: March 2020
# description: data import helper function for als variant viewer
# usage: 
# ----------
library(tidyverse)
library(stringr)
library(data.table)
library(DT)

#generate the gene names by listing files and extracted from listed files
gene_csv<- system("ls -1 clinvar_cache/",intern=TRUE)
gene_csv <- as.vector(gene_csv)

gene_names<- str_replace(string = gene_csv, pattern = '_clinvar.csv',replacement = '')
gene_names <- as.vector(gene_names)

# read in the csv file containing the information relevant to the gene itself
gene_info <- as.data.frame(fread("gene_info.csv"))



#Function to process the imported clinvar data----------------
clinsig_stratify <- function(gene){
  #gene <- as.character(gene) #make sure the gene column is read as a character
  
#stratify the different Clinical.Significance levels 1 of 5 categories
  gene <- gsub(pattern = "Conflicting interpretations of pathogenicity",replacement="Uncertain significance", x = gene) 
  gene <- gsub(pattern = "Pathogenic/Likely pathogenic",replacement = "Likely pathogenic",x=gene)
  gene <- gsub(pattern = "Benign/Likely benign",replacement = "Likely benign",x=gene)
  gene <- gsub(pattern = "risk factor",replacement = "Likely benign",x = gene)
}

#import datasets from the clinvar_cache--> generated from the clinvar_cleaner.R script
#reads each clinvar_csv file and reads it into a dataframe named the gene
counter = 0
for (i in gene_csv){
 counter = counter + 1 #counter for refernce which gene in the gene names list is being used
 genename<- gene_names[counter] #assign the current gene to the variable genename
 df <- as.data.frame(fread(file = paste("clinvar_cache/",i,sep="")))
 assign(x=genename,value=df)
}

#stratify the clinical significance categories into 5 
SOD1$Clinical.Significance<- clinsig_stratify(SOD1$Clinical.Significance)
FUS$Clinical.Significance<- clinsig_stratify(FUS$Clinical.Significance)
DCTN1$Clinical.Significance<- clinsig_stratify(DCTN1$Clinical.Significance)
TARDBP$Clinical.Significance<- clinsig_stratify(TARDBP$Clinical.Significance)
ALS2$Clinical.Significance<- clinsig_stratify(ALS2$Clinical.Significance)
SETX$Clinical.Significance <- clinsig_stratify(SETX$Clinical.Significance)
VAPB$Clinical.Significance <- clinsig_stratify(VAPB$Clinical.Significance)
MATR3$Clinical.Significance<- clinsig_stratify(MATR3$Clinical.Significance)
OPTN$Clinical.Significance <- clinsig_stratify(OPTN$Clinical.Significance)
SQSTM1$Clinical.Significance <- clinsig_stratify(SQSTM1$Clinical.Significance)
FIG4$Clinical.Significance <- clinsig_stratify(FIG4$Clinical.Significance)
SLC52A3$Clinical.Significance <- clinsig_stratify(SLC52A3$Clinical.Significance)
C9orf72$Clinical.Significance <- clinsig_stratify(C9orf72$Clinical.Significance)
VCP$Clinical.Significance <- clinsig_stratify(VCP$Clinical.Significance)
TBK1$Clinical.Significance <- clinsig_stratify(TBK1$Clinical.Significance)
CHCHD10$Clinical.Significance <- clinsig_stratify(CHCHD10$Clinical.Significance)
SIGMAR1$Clinical.Significance <- clinsig_stratify(SIGMAR1$Clinical.Significance)
ANG$Clinical.Significance <- clinsig_stratify(ANG$Clinical.Significance)
UBQLN2$Clinical.Significance <- clinsig_stratify(UBQLN2$Clinical.Significance)
SPG11$Clinical.Significance <- clinsig_stratify(SPG11$Clinical.Significance)
KIF5A$Clinical.Significance <- clinsig_stratify(KIF5A$Clinical.Significance)
NEK1$Clinical.Significance <- clinsig_stratify(NEK1$Clinical.Significance)