#the source code for the als_clinvar that is used for generating a plot. Uses ggplot. 
library(tidyverse)
library(RColorBrewer)
library(data.table)
library(plotly)

# #data import----------
variant_graph <- function(x){
 x$ClinVar.Stars <- sapply(as.character(x$Review.Criteria), function (x) {
  switch(x,
         "practice guideline" = "****",
         "reviewed by expert panel" = "***",
         "criteria provided, multiple submitters, no conflicts" = "**",
         "criteria provided, conflicting interpretations" = "*",
         "criteria provided, single submitter" = "*",
         "no assertion for the individual variant" = "none",
         "no assertion criteria provided" = "none",
         "no assertion provided" = "none")
})


#graphing----------
p<- ggplot(x, aes(x=Position,text3= rsID, y= Clinical.Significance, text1 = Nucleotide.Consequence, text2= Protein.Consequence, text4 = ClinVar.Stars))+
  
  # #adds lollipops
  # geom_linerange(position = position_jitter(width = 2000,height = 0,seed = 1),size = 1,
  #                aes(
  #                  x = Position,
  #                  ymin = 0,
  #                  ymax = Clinical.Significance
  #                ),show.legend = FALSE)+
  
  #prevents overplotting
    geom_jitter(position = position_jitter(width = 2000,height = 0,seed = 1),mapping= aes(fill= Clinical.Significance),shape = 21,size = 5,color = "black", stroke = 1,show.legend = TRUE)+
  theme(panel.background = element_blank(),
        axis.text.x = element_text(angle=70,vjust = 0.5, hjust = 0.7),
        axis.title.x.bottom = element_text(),
        axis.line.x = element_line(size = 6 ,colour = "grey",lineend = "square"),
        #y elements-----
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        #legend elements-----
        legend.background = element_rect(fill="white"),
        legend.title = element_blank()
   )+
  labs(x = "Position (GRCh38)")+
  
  #assign color to each value
  scale_fill_manual(limits = c("Benign","Likely benign","Likely pathogenic","Pathogenic","Uncertain significance"),values = c("#ADDDA8","#FFFFC4","#FAAE6A","#D41B25","#3483B7"))+
   #ordering the y discrete values on the y axis
  scale_y_discrete(limits = c("Uncertain significance","Pathogenic","Likely pathogenic","Likely benign","Benign"))

#convert ggplot to plotly object----------
ggplotly(p,tooltip = c("x","fill","text1","text2","text3","text4")) %>% layout(legend = list(orientation = "h",x=1,y=-1))
}






# SOD1<- as.data.frame(fread("~/Documents/GitHub/LNG_Scripts/als_app/clinvar_cache/SOD1_clinvar.csv"))
# 
# #data cleanup----------
# SOD1$Clinical.Significance <- gsub(pattern = "Conflicting interpretations of pathogenicity",replacement="Uncertain significance",x=SOD1$Clinical.Significance)
# SOD1$Clinical.Significance <- gsub(pattern = "Benign/Likely benign",replacement = "Likely benign",x=SOD1$Clinical.Significance)
# SOD1$Clinical.Significance <- gsub(pattern = "Pathogenic/Likely pathogenic", replacement = "Likely pathogenic", x=SOD1$Clinical.Significance)
# SOD1$Position<- as.numeric(SOD1$Position)
# SOD1$Review.Criteria <- as.factor(SOD1$Review.Criteria)
# SOD1$Clinical.Significance <- as.factor(SOD1$Clinical.Significance)
# SOD1$Review.Criteria<-as.factor(SOD1$Review.Criteria)
# 
# # #Sort SOD1 vector in the custom order
#  SOD1 <- SOD1 %>%
#    arrange(factor(Clinical.Significance, levels = c("Pathogenic","Benign","Likely benign","Likely pathogenic","Uncertain significance")))

#########EVERYTHING BELOW SHOULD BE included in the function statement when we Decide we are done revising the plot########################################## 