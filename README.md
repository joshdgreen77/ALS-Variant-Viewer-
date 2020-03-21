# ALS-Variant-Viewer-
Shiny application that is used to view variants and genes reported to be associated with ALS on ClinVar
## Quick Start
### Step 1: Install relevant packages 
Copy and paste this into your R console.
```

install.packages(c("shiny","shinydashboard","tidyverse","stringr","data.table","DT","RColorBrewer","plotly"))

```
### Step 2: Clone repository to your local drive 

### Step 3: Open App.R and click "Run App" in top right corner of Rstudio



## Detailed Workflow 
Processing the raw ClinVar data
For the sake of transparency, below I will go through the workflow of processing the data so that it can be viewed on the ALS Variant Viewer application.

### Step 1: Run Processing_Raw_ClinVar_Data.R script
`rscript Processing_Raw_ClinVar_Data.R`
##### Load raw ClinVar data
```R
clinvar_raw <- as.data.frame(fread(file ="clinvar_ALS_raw.txt"))
```
##### Polish raw ClinVar data
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
Some variants are labeled with multiple gene designations so I assigned each variant to a single gene
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
#### Export data as  file
```R
# write the cleaned up clinvar data to a file to be read by clinvar cleaner
write_csv(x = clinvar_clean_gene,path = "../2_Gene_Formatting/clinvar_ALS.csv")
```
### Step 2: Run clinvar_cleaner.R script
`rscript clinvar_cleaner.R`

**Objective:** Separate ClinVar variants by gene and merge them with gnomAD data to get rsID and allele frequency.

##### Load polished ClinVar data 
```R
processed_clinvar <- fread(file ="../2_Gene_Formatting/clinvar_ALS.csv")
```
#### Helper functions
##### `clinvar.parse` function

**Purpose:** Polish, extract, and web scrape Clinical Significance, Protein Consequence, Nucleotide Consequence, and Review Criteria from the processed ClinVar data.

```R
# function for parsing and extracting information from the processed clinvar data----------
clinvar.parse <- function(x){
#remove the date for clinical significance column
x$Clinical.Significance <- str_replace(string = x$Clinical.Significance,pattern = "\\(.*\\)",replacement = "")
  
#extracting information from the Name column in dataframe
  x$Protein.Consequence <-str_extract(string= x$Name, pattern = "p\\....\\d*...")
  x$Nucleotide.Consequence <-str_extract(string= x$Name, pattern = "c.*>[A,T,G,C]")

# Web scrape variant review criteria from ClinVar website
review.criteria<- c()
  for(i in 1:length(x$VariationID)) {
  #loading bar
    print(paste(i,"of",length(x$VariationID),x$Variation[i],cv.url,sep=" "))
    
    # generate URL
    cv.url <- paste("https://www.ncbi.nlm.nih.gov/clinvar/variation/",x$VariationID[i],"/",sep="")
    # web scrap from URL
    review<- cv.url %>% read_html() %>% html_node(css = '.no-margin-top dd:nth-child(4)') %>% html_text()
    # text polishing
    review_f <- gsub(pattern = "\\\n || \\s",replacement = "",x = review)
    # append to empty vector
    review.criteria[i] <- review_f
  }
  # append vector to the data frame
  x$Review.Criteria <- review.criteria
  x <- as.data.frame(x)
}
```
##### `gnomad_join` function

**Purpose:** join the polished ClinVar dataset, that has been filtered by gene, with the gene-specific data on gnomAD to obtain the rsID and the allele frequency columns
```R
# function for joining clinvar dataset with gnomad dataset
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
```
##### `citation_extractor` function

**Purpose:** to web scrape the PubMed ID to be included in the app's variant table. (Still being developed)
```
#Not fully developed yet
```
#### Applying helper functions to each gene
**Purpose:** this function filters the processed ClinVar data by gene, applies the 3 helper functions,  selects and renames columns, and exports the generated data frame as a csv file into the "clinvar_cache" folder

```R
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
```
All genes included are below
```R
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
```
#### Web scrap information on each gene


