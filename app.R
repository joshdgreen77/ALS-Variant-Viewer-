#---------
#importing the required packages----------
library(shiny)
library(shinydashboard)


# Sourcing and Importing-----------

#import the processed clinvar datasets
source("data_import.R")

#import the function that makes variant plot makes the "variant_graph" function
source("variant_graph.R")

#import the function that makes the data table "tableman" function
source("variant_DT.R")



# User interface---------------
ui <- dashboardPage(
  dashboardHeader(title ="ALS Variant Viewer v1.0"),
  dashboardSidebar(
    sidebarMenu(
      #make sidebar menu item for the variant plot
      menuItem("Variant Plot", tabName = "plot"),
      menuItem("About",tabName = "About"),
      #make selector widget in the sidebar
      selectInput(inputId = "select",
                  label = "Choose Gene to display:",
                  choices = gene_names,
                  selected = gene_names[1])
    )
  ),
  dashboardBody(
      tabItems(
        #the interface for the plot tab
        tabItem(tabName =  "plot",
                fluidRow(
                  column(width=6, box(width=12,
                                      h1(textOutput(outputId = "gene")),
                                      hr(),
                                      h2("Exons:"),h3(textOutput(outputId="exons")),
                                      h2("Locus:"),h3(textOutput(outputId="locus")),
                                      h2("Protein:"), h3(textOutput(outputId="protein")),
                                      h2("Function:"),h3(textOutput(outputId="pfunction"))
                  )),
                  column(width=6,
                         box(width=12,
                             plotlyOutput("plot")))
                ),
                fluidRow(
                  box(width=12,
                      DT::dataTableOutput("dt"))
                )
        ),
        tabItem(tabName="About",
                box(width=12,
                    h1("About"),
                    hr(),
                    h2("Authors"),
                    p(" Written and Developed by Joshua D Green"),
                    hr(),
                    h2("Data"),
                    hr(),
                  p("The ALS-linked variants displayed in this browser are from NCBI ClinVar (Last retrieved March 15 2020). The SNP identifiers (rsID) and Allele Frequencies are from gnomadAD v2.1.1. The information on each gene is from the NBCI gene search tool")
                    )
                )
                 
       
      )
    
  )
)



# Server processing----------
server <- function(input, output){
  # assign the choice in the select input to a variable
  output$gene <- renderText({
    paste(input$select)
  })
  
# generates the gene information object
  
  # number of exons
  output$exons <- renderText({
    info <- gene_info %>% filter(gene == input$select)
    paste0(info$exon)})
  
  # cytogenetic locus
  output$locus <- renderText({
    info <- gene_info %>% filter(gene == input$select)
    paste0(info$locus)})
  # protein 
  output$protein <- renderText({
    info <- gene_info %>% filter(gene == input$select)
    paste0(info$protein)})
  
  # protein function
  output$pfunction <- renderText({
    info <- gene_info %>% filter(gene == input$select)
    paste0(info$pfunction)})
  
  # generates the plotly object 
  output$plot <-renderPlotly({
    data <- switch(input$select,
                   "SOD1"=SOD1,
                   "FUS"=FUS,
                   "DCTN1"=DCTN1,
                   "TARDBP"=TARDBP,
                   "ALS2"=ALS2,
                   "SETX"=SETX,
                   "VAPB"=VAPB,
                   "MATR3"=MATR3,
                   "OPTN"=OPTN,
                   "SQSTM1"=SQSTM1,
                   "FIG4"=FIG4,
                   "SLC52A3"=SLC52A3,
                   "C9orf72"=C9orf72,
                   "VCP"=VCP,
                   "TBK1"=TBK1,
                   "CHCHD10"=CHCHD10,
                   "SIGMAR1"=SIGMAR1,
                   "ANG"=ANG,
                   "UBQLN2"=UBQLN2,
                   "SPG11"=SPG11,
                   "KIF5A"=KIF5A
                   
    )
    # the imported variant_graph function comes from the thing sourced in line 10
    variant_graph(data)
  })
  
  
  output$dt <- DT::renderDataTable({
    #switches the textual input$select data to the actual data frame variables
    data <- switch(input$select,
                   "SOD1"=SOD1,
                   "FUS"=FUS,
                   "DCTN1"=DCTN1,
                   "TARDBP"=TARDBP,
                   "ALS2"=ALS2,
                   "SETX"=SETX,
                   "VAPB"=VAPB,
                   "MATR3"=MATR3,
                   "OPTN"=OPTN,
                   "SQSTM1"=SQSTM1,
                   "FIG4"=FIG4,
                   "SLC52A3"=SLC52A3,
                   "C9orf72"=C9orf72,
                   "VCP"=VCP,
                   "TBK1"=TBK1,
                   "CHCHD10"=CHCHD10,
                   "SIGMAR1"=SIGMAR1,
                   "ANG"=ANG,
                   "UBQLN2"=UBQLN2,
                   "SPG11"=SPG11,
                   "KIF5A"=KIF5A
                   
    )
    tableman(data2 = data)
  })
  
}
shinyApp(ui = ui, server = server)






