#---------
##RUN THIS TO UPDATE SHINY APP ON THE SERVER##
# library(rsconnect)
# rsconnect::deployApp("/Users/greenjod/Documents/GitHub/ALS-Variant-Viewer-/")


#importing the required packages----------
library(shiny)
library(shinydashboard)

# Sourcing and Importing-----------

#import the processed clinvar datasets
source("data_import.R")

#import the function that makes variant plot makes the "variant_graph" function
source("variant_graph.R")

#import the function that makes the data table "format_variant_table" function
source("variant_DT.R")



# User interface---------------
ui <- dashboardPage(
  dashboardHeader(title ="ALS Variant Browser"),
  dashboardSidebar(
    sidebarMenu(
      #make sidebar menu item for the variant plot
      menuItem("Variant Browser ", tabName = "browser"),
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
        tabItem(tabName =  "browser",
                fluidRow(
                  column(width=6, box(width=12,
                                      h1(uiOutput(outputId = "gene")),
                                      hr(),
                                      h2("Exons:"),h3(textOutput(outputId="exons")),
                                      h2("Locus:"),h3(textOutput(outputId="locus")),
                                      h2("Protein:"), h3(textOutput(outputId="protein")),
                                      h2("Function:"),h3(textOutput(outputId="pfunction"))
                  )),
                  column(width=6,
                         box(width=12,
                             plotlyOutput("plot")),
                         box(width=12,
                             tableOutput("clinsig_count")))
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
                    p("Written and Developed by Joshua D Green"),
                    hr(),
                    h2("Data"),
                    hr(),
                  tags$ul(
                    tags$li("The variants displayed in this browser have been linked to Amyotrophic Lateral Sclerosis (ALS) on NCBI ClinVar (Last checked March 15 2020)."),
                    tags$li("Reference SNP cluster ID (rsID) and Allele Frequencies are retrieved from gnomadAD v2.1.1."),
                    tags$li("Information about each gene is from NCBI Gene database")),
                  hr(),
                  h2("References"),
                  tags$ul(
                    tags$li(" Karczewski, K. J., Francioli, L. C., Tiao, G., Cummings, B. B., Alföldi, J., Wang, Q., Collins, R. L., Laricchia, K. M., Ganna, A., Birnbaum, D. P., Gauthier, L. D., Brand, H., Solomonson, M., Watts, N. A., Rhodes, D., Singer-Berk, M., England, E. M., Seaby, E. G., Kosmicki, J. A., … MacArthur, D. G. (2019). Variation across 141,456 human exomes and genomes reveals the spectrum of loss-of-function intolerance across human protein-coding genes. In bioRxiv (p. 531210),",tags$a(href="https://doi.org/10.1101/531210", "https://doi.org/10.1101/531210")), 
                    br(),
                    tags$li("Landrum MJ, Lee JM, Benson M, Brown GR, Chao C, Chitipiralla S, Gu B, Hart J, Hoffman D, Jang W, Karapetyan K, Katz K, Liu C, Maddipatla Z, Malheiro A, McDaniel K, Ovetsky M, Riley G, Zhou G, Holmes JB, Kattman BL, Maglott DR. ClinVar: improving access to variant interpretations and supporting evidence. Nucleic Acids Res . 2018 Jan 4. PubMed PMID:",tags$a(href="https://www.ncbi.nlm.nih.gov/pubmed/29165669","29165669")) 
                    )
                  )
                )
        )
      )
)



# Server processing----------
server <- function(input, output){
  # assign the choice in the select input to a variable
  output$gene <- renderUI({
    info <- gene_info %>% filter(gene == input$select)
    
    url <- a(input$select, href=info$link)
    tagList(url)
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
    #retrieve the variable assigned to the input$select
    data <- get(input$select)
    variant_graph(data)
  })
  
  # generates the data table object
  output$dt <- DT::renderDataTable({
    #retrieve the object assigned to the input$select
    data <- get(input$select)
    format_variant_table(data)
  })
  
  #
  output$clinsig_count <- renderTable({
    data <- get(input$select)
    tibble("Variant Clinical Significance"=c("Pathogenic","Likely pathogenic","Likely benign","Benign","Uncertain significance"),
               "Count"=c(nrow(subset(data,data$Clinical.Significance=="Pathogenic")),
                         nrow(subset(data,data$Clinical.Significance=="Likely pathogenic")),
                         nrow(subset(data,data$Clinical.Significance=="Likely benign")),
                         nrow(subset(data,data$Clinical.Significance=="Benign")),
                         nrow(subset(data,data$Clinical.Significance=="Uncertain significance"))))
    
  })
  
}
shinyApp(ui = ui, server = server)






