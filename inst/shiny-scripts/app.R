# This script is adapted from
# https://github.com/anjalisilva/TestingPackage/blob/master/inst/shiny-scripts/app.R

library(shiny)
library(shinyalert)

# Define UI
ui <- fluidPage(

  # page title
  titlePanel("Fragmentomic Analysis of cfDNA for Cancer Detection and Subtyping:
             Finding Variation in Sample and Patient Data"),

  # Sidebar layout with input and output definitions ----
  sidebarLayout(

    # Sidebar panel for inputs ----
    sidebarPanel(

      tags$p("This is a Shiny App that is part of the CfDNAfragmentomics
              package in R. Its purpose is to perform Wilcoxon Rank Sum-test
             analysis and plotting of cfDNA fragmentome data."),
      # br() element to introduce extra vertical spacing ----
      br(),

      tags$b("Description: CfDNAfragmentomics is an R package used
          to find significance of difference in cfDNA length fragments
          between patient and control data. This Shiny App is part of the
          CfDNAfragmentomics pakage. It permits to calculate the Wilcoxon Rank Sum-test
          for the variation in cfDNA fragment sizes,
          and significance of difference between patient and control data,
          and compares the sizes of the cfDNA to predict if the patient
          data is potentially cancerous. This app also allows for plotting of
          the cfDNA fragment sizes to allow for a visual
          comparison of the cfDNA fragment lengths. For more
             details, see ?nucleosomeRatio"),

      # br() element to introduce extra vertical spacing ----
      br(),
      br(),

      # input
      tags$p("Instructions: Below, enter or select values required to perform
             the analysis. Default
             values are shown. Then press 'Run'. Navigate through
             the two tabs to the right to explore the results."),

      # br() element to introduce extra vertical spacing ----
      br(),

      # input
      shinyalert::useShinyalert(),  # Set up shinyalert
      uiOutput("tab2"),
      actionButton(inputId = "data1",
                   label = "Dataset 1 Details"),
      uiOutput("tab1"),
      actionButton(inputId = "data2",
                   label = "Dataset 2 Details"),
      fileInput(inputId = "file1",
                label = "Select a cfDNA bed file of patient data. File should
                be in .bed or .txt format rows corresponding
                to chromosome regions and tab seperated columns
                corresponding to chromosome, start, and end points.",
                accept = c(".bed", ".txt")),
      fileInput(inputId = "file2",
                label = "Select a cfDNA bed file of control data.
                File should be in .bed or .txt format with rows corresponding
                to chromosome regions and tab seperated columns corresponding
                to chromosome, start, and end points.",
                accept = c(".bed", ".txt")),
      textInput(inputId = "di_nucleosome_p_value",
                label = "Enter dinucleosome p_value.
                This should be a positive decimal value less than 1:", "0.05"),
      textInput(inputId = "mono_nucleosome_p_value",
                label = "Enter mononucleosome p_value.
                This should be a positive decimal value less than 1:", "0.05"),


      # br() element to introduce extra vertical spacing ----
      br(),

      # actionButton
      actionButton(inputId = "button2",
                   label = "Run"),

    ),

    # Main panel for displaying outputs ----
    mainPanel(

      # Output: Tabset w/ plot, summary, and table ----
      tabsetPanel(type = "tabs",
                  tabPanel("Plot of cfDNA length Dataset",
                           h3("Instructions: Enter values and click 'Run'
                              at the bottom left side."),
                           h3("Pairs Plot of patient and control length
                              densities:"),
                           br(),
                           plotOutput("density_plot_out")
                  ), #close plot panel
                  tabPanel("Result of the Wilcox-test",
                           h3("Instructions: Enter values and click 'Run' at
                              the bottom left side."),
                           h3("Result of the wilcox-test:"),
                           h5("If the result is True, then for the inputted
                              patient data set, the mononucleosome
                              and dinucleosome cfDNA  fragment sizes are
                              significantly
                              shorter than the control cfDNA fragment sizes.
                              This is a cancer biomarker and indicates that patient
                              data is potentially cancerous"),
                           h5("If the result is False, then for the inputted
                              patient data set, the monnucleosome
                              and dinucleosome cfDNA  fragment sizes are not
                              significantly
                              shorter than the control cfDNA fragment sizes.
                              This indicates that the patient
                              data does not have indication of this specific cancer
                              biomarker and is potentially non-cancerous."),
                           br(),
                           verbatimTextOutput("ttest_out_cancerous"),
                           br(),

                           h5("Calculated p-value for mononucleosome cfDNA lengths"),
                           verbatimTextOutput("ttest_out_mono"),

                           br(),

                           h5("Calculated p-value for dinucleosome cfDNA lengths"),
                           verbatimTextOutput("ttest_out_di")

                  ),#close verbatim panel
                  tabPanel("Plot of cfDNA coverage",
                           h3("Instructions: Enter values and click 'Run'
                              at the bottom left side."),
                           h3("Plot of patient genomic coverage. This function
                              might takea long time to output the graph as
                              the chromosome length could be quite long."),
                           br(),
                           plotOutput("coverage_plot_out")
                  ), #close plot panel

                  tabPanel("cfDNA coverage array",
                           h3("Instructions: Enter values and click 'Run'
                              at the bottom left side."),
                           h3("output array of coverage: this function
                              might take a long time to output an array
                              because of the long chromosome lengths"),
                           br(),
                           verbatimTextOutput("coverage_list_out")
                  ) #close plot panel



      ) #close tabset

    ) #close main

  )
) #close iu

# Define server logic for random distribution app ----
server <- function(input, output) {


  # Calculate wilcox test value
  start_ttest <- eventReactive(eventExpr = input$button2, {

    if (!is.null(input$file1) & !is.null(input$file2))
      file1_acc <- input$file1
    file2_acc <- input$file2
    sample_bed_f <- read.table(file1_acc$datapath, header = FALSE,
                               sep="\t", stringsAsFactors=FALSE, quote="")
    controls_bed_f <- read.table(file2_acc$datapath, header = FALSE, sep="\t",
                                 stringsAsFactors=FALSE, quote="")


    CfDNAfragmentomics::nucleosomeRatio(
      controls_bed = controls_bed_f,
      sample_bed = sample_bed_f,
      di_nucleosome_p_value = as.numeric(input$di_nucleosome_p_value),
      mono_nucleosome_p_value = as.numeric(input$mono_nucleosome_p_value))



  })

  start_plotting <- eventReactive(eventExpr = input$button2, {

    if (!is.null(input$file1) & !is.null(input$file2))
      file1_acc <- input$file1
    file2_acc <- input$file2
    sample_bed_f <- read.table(file1_acc$datapath, header = FALSE, sep="\t",
                               stringsAsFactors=FALSE, quote="")
    controls_bed_f <- read.table(file2_acc$datapath, header = FALSE, sep="\t",
                                 stringsAsFactors=FALSE, quote="")



    CfDNAfragmentomics::nucleosomeDensityPlot(controls_bed = controls_bed_f,
                                              sample_bed = sample_bed_f)


  })


  coverage_plotting <- eventReactive(eventExpr = input$button2, {

    if (!is.null(input$file1))
      file1_acc <- input$file1
    sample_bed_f <- read.table(file1_acc$datapath, header = FALSE, sep="\t",
                               stringsAsFactors=FALSE, quote="")



    CfDNAfragmentomics::nucleosomeCoveragePlot(sample_bed = sample_bed_f)


  })

  # Calculate coverage array
  start_coverage_array <- eventReactive(eventExpr = input$button2, {

    if (!is.null(input$file1))
      file1_acc <- input$file1
    sample_bed_f <- read.table(file1_acc$datapath, header = FALSE,
                               sep="\t", stringsAsFactors=FALSE, quote="")


    CfDNAfragmentomics::nucleosomeCoverage(
      sample_bed = sample_bed_f)



  })

  #-------------------Outputs------------------------


  # Text output for wilcox-test analysis
  output$ttest_out_cancerous <- renderPrint({
    if (! is.null(start_ttest))
      start_ttest()$cancerous
  })


  # Text output for U-test analysis
  output$ttest_out_mono <- renderPrint({
    if (! is.null(start_ttest))
      start_ttest()$ttest_mono_pvalue
  })

  # Text output for U-test analysis
  output$ttest_out_di <- renderPrint({
    if (! is.null(start_ttest))
      start_ttest()$ttest_di_pvalue
  })

  # Plotting cfDNA density values
  output$density_plot_out <- renderPlot({
    if (! is.null(start_plotting))

      start_plotting()

  })

  # Plotting cfDNA coverage values
  output$coverage_plot_out <- renderPlot({
    if (! is.null(coverage_plotting))

      coverage_plotting()

  })

  # Plotting cfDNA coverage values
  output$coverage_list_out <- renderPrint({
    if (! is.null(start_coverage_array))

      start_coverage_array()

  })


  #---------------------------------------------data set descriptions-----------

  # URLs for downloading data
  url1 <- a("Example Dataset 2 (p1.bed)",
            href="https://raw.githubusercontent.com/Yasamin-Nourijelyani/CfDNAfragmentomics/master/inst/extdata/p1.bed")
  output$tab1 <- renderUI({
    tagList("Download:", url1)
  })

  observeEvent(input$data2, {
    # Show a modal when the button is pressed
    shinyalert(title = "Example Dataset 2",
               text = "This is a control cfDNA length example dataset (not real data) BED files
               from the Illumina website. It includes the chromosome numbers,
               and fragment start and ends (columns 1-3 are used in this analysis).
               It has a size of n = 3251 observations along rows and 4 variable columns.
               To save the file, click on link, then right click and save as .txt file.
               citation: Illumina, https://support.illumina.com/downloads/nextera-flex-for-enrichment-BED-files.html\n
               More examples can be found at: https://github.com/Yasamin-Nourijelyani/CfDNAfragmentomics/tree/master/inst/extdata

               "

               ,
               type = "info")
  })

  url2 <- a("Example Dataset 1 (d1.bed)", href="https://raw.githubusercontent.com/Yasamin-Nourijelyani/CfDNAfragmentomics/master/inst/extdata/d1.bed")
  output$tab2 <- renderUI({
    tagList("Download:", url2)
  })

  observeEvent(input$data1, {
    # Show a modal when the button is pressed
    shinyalert(title = "Example Dataset 1",
               text = "This is a patient cfDNA length example (not real data) dataset BED files
               from the Illumina website. It includes the chromosome numbers,
               and fragment start and ends (columns 1-3 are used in this analysis).
               It has a size of n = 1895 observations along rows, and 4 columns.
               To save the file, click on link, then right click and save as .txt file.
               citation: Illumina, https://support.illumina.com/downloads/nextera-flex-for-enrichment-BED-files.html\n
               More examples can be found at: https://github.com/Yasamin-Nourijelyani/CfDNAfragmentomics/tree/master/inst/extdata

               ",
               type = "info")
  })



}

# Create Shiny app ----
shinyApp(ui, server)

# [END]
