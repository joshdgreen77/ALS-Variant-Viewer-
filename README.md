# ALS-Variant-Viewer-
Shiny application that is used to view variants and genes reported to be associated with ALS on ClinVar
## Quick Start
### 1) Install relevant packages 
Copy and paste this into your R console.
```

install.packages(c("shiny","shinydashboard","tidyverse","stringr","data.table","DT","RColorBrewer","plotly"))

```
### 2) Clone repository to your local drive 

### 3) Open App.R and click "Run App" in top right corner of Rstudio



## Detailed Workflow 
Processing the raw ClinVar data
For the sake of transparency, below I will go through the workflow of processing the data so that it can be viewed on the ALS Variant Viewer application.

### Step 1: Run Processing_Raw_ClinVar_Data.R script
`rscript Processing_Raw_ClinVar_Data.R`
##### Load data
```R
clinvar_raw <- as.data.frame(fread(file ="clinvar_ALS_raw.txt"))
```
##### Clean raw ClinVar data
```R
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
```
##### Clean up the "Gene" column
Some variants are labeled with multiple gene designations so to make data easier to read I assigned each variant to a single gene
```R
# renamed data frame to make easier to see purpose
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
```
#### Export cleaned data as csv file
```R
# write the cleaned up clinvar data to a file to be read by clinvar cleaner
write_csv(x = clinvar_clean_gene,path = "../2_Gene_Formatting/clinvar_ALS.csv")
```

