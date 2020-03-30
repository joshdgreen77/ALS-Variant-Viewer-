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


#import datasets from the clinvar_cache--> generated from the clinvar_cleaner.R script
#reads each clinvar_csv file and reads it into a dataframe named the gene
counter = 0
for (i in gene_csv){
  counter = counter + 1 #counter for refernce which gene in the gene names list is being used
  genename<- gene_names[counter] #assign the current gene to the variable genename
  gene_specific_dataframe <- as.data.frame(fread(file = paste("clinvar_cache/",i,sep="")))
  assign(x=genename,value=gene_specific_dataframe)
}