#the source code for the als_clinvar that is used for generating a plot. Uses ggplot. 
library(tidyverse)
library(RColorBrewer)
library(data.table)
library(plotly)

# #data import----------
variant_graph <- function(clinvar_data,exons_start,exons_stop){


#graphing----------
p<- ggplot(clinvar_data, aes(x=Position, y= Clinical.Significance, text1 = Nucleotide.Consequence, text2= Protein.Consequence,text3= rsID))+
  
  #prevents overplotting
    geom_jitter(position = position_jitter(width = 0,height = 0,seed = 1),mapping= aes(fill= Clinical.Significance),shape = 21,size = 3,color = "black", stroke = 0.5,show.legend = TRUE)+
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
  scale_y_discrete(limits = c("Uncertain significance","Pathogenic","Likely pathogenic","Likely benign","Benign"))+
  
  annotate(geom="rect", xmin = exons_start, xmax = exons_stop, ymin = 0, ymax = 0.5, color = "black",fill = "orange")

  ##WORK IN PROGRESS##
#this is for adding text to label the exons
#annotate("text",x=min(clinvar_data$Position)-7000,y=0.25,label="exons")

#convert ggplot to plotly object----------
ggplotly(p,tooltip = c("x","fill","text1","text2","text3")) %>% layout(legend = list(orientation = "h",x=0,y=-1))
}