# ----------
# title: gene_scrapper.R
# author: Joshua D Green
# date: February 2020
# description: scraps data from the NCBI website for each gene. 
# usage: 
# ----------
# import the packages
library(rvest)
library(stringr)
library(tidyverse)

# gene scrape functions
generate_gene_info <- function(gene,url) {
#define the url
url <- url
# webscrape exons
exon <- url %>% read_html() %>% html_node(css = '.exon-count dd') %>% html_text()

# webscrape locus
locus <- url %>% read_html() %>% html_node(css = '.dl-chr-info span') %>% html_text()

# scrape the protein name
protein_raw<- url %>% read_html() %>% html_node(css = 'dd:nth-child(4)') %>% html_text()
# remove the unwanted text characters
protein<- str_replace(string=protein_raw, pattern="provided by HGNC",replacement = "")

# gene function
pfunction <- url %>% read_html() %>% html_node(css = 'dd:nth-child(20)') %>% html_text()

link <-paste(url) 
return(data.frame(gene,exon,locus,protein,pfunction,link))
}

# scraps the information from each gene from the NCBI  website

#NEK1 
NEK1 <- generate_gene_info("NEK1","https://www.ncbi.nlm.nih.gov/gene/4750")

#KIF5A
KIF5A <- generate_gene_info("KIF5A","https://www.ncbi.nlm.nih.gov/gene/3798")

#SPG11
SPG11 <- generate_gene_info("SPG11","https://www.ncbi.nlm.nih.gov/gene/80208")

#UBQLN2
UBQLN2 <- generate_gene_info("UBQLN2","https://www.ncbi.nlm.nih.gov/gene/29978")
#ANG
ANG <- generate_gene_info("ANG","https://www.ncbi.nlm.nih.gov/gene/283")

#SIGMAR1 
SIGMAR1 <- generate_gene_info("SIGMAR1","https://www.ncbi.nlm.nih.gov/gene/10280")

#CHCHD10 
CHCHD10 <- generate_gene_info("CHCHD10","https://www.ncbi.nlm.nih.gov/gene/400916")

# TBK1
TBK1 <- generate_gene_info("TBK1","https://www.ncbi.nlm.nih.gov/gene/29110")

# VCP 
VCP <- generate_gene_info("VCP","https://www.ncbi.nlm.nih.gov/gene/7415")

# C9orf72
C9orf72 <- generate_gene_info("C9orf72","https://www.ncbi.nlm.nih.gov/gene/203228")

#SLC52A3
SLC52A3 <- generate_gene_info("SLC52A3","https://www.ncbi.nlm.nih.gov/gene/113278")

#FIG4
FIG4 <- generate_gene_info("FIG4","https://www.ncbi.nlm.nih.gov/gene/9896")

#SQSTM1
SQSTM1 <- generate_gene_info("SQSTM1","https://www.ncbi.nlm.nih.gov/gene/8878")

#OPTN 
OPTN <- generate_gene_info("OPTN","https://www.ncbi.nlm.nih.gov/gene/10133")

#MATR3
MATR3 <- generate_gene_info("MATR3","https://www.ncbi.nlm.nih.gov/gene/9782")

#VAPB
VAPB <-  generate_gene_info("VAPB","https://www.ncbi.nlm.nih.gov/gene/9217") 

#SETX
SETX <- generate_gene_info("SETX","https://www.ncbi.nlm.nih.gov/gene/23064") 

#ALS2
ALS2 <- generate_gene_info("ALS2","https://www.ncbi.nlm.nih.gov/gene/57679") 

#TARDBP
TARDBP<- generate_gene_info("TARDBP","https://www.ncbi.nlm.nih.gov/gene/23435") 

# DCTN1 
DCTN1 <- generate_gene_info("DCTN1","https://www.ncbi.nlm.nih.gov/gene/1639")

# FUS
FUS <- generate_gene_info("FUS","https://www.ncbi.nlm.nih.gov/gene/2521")

# SOD1
SOD1 <- generate_gene_info("SOD1","https://www.ncbi.nlm.nih.gov/gene/6647/")


# for re scraping all genes-----------------
# bind the initial genes together
gene_info<-rbind(SOD1,FUS,DCTN1,TARDBP,ALS2,SETX,VAPB,MATR3,OPTN,SQSTM1,FIG4,SLC52A3,C9orf72,VCP,TBK1,CHCHD10,SIGMAR1,ANG,UBQLN2,SPG11,KIF5A,NEK1)
# write csv file to be sourced in the app.R
write.csv(x = gene_info,file = "../gene_info.csv")



# append new genes --> run this function in the console just to prevent accidentally appending same gene multiple times
# appender<-function(x){
#   josh <- read.csv(file = "../gene_info.csv")
#   josh <- josh[ -c(1)]
#   josh <- rbind(josh,x)
# 
#   write.csv(x = josh,file = "../gene_info.csv")
# }
