#setwd("/Users/marek/JLUbox/Andreas/cTagPipe")
setwd("/home/marek/ccTop/")
library(roxygen2)
library(devtools)
install_github("mrkbrtkhn/cTagPipe")

#list.files()
#install("cTagPipeTest",reload = TRUE)

library(cTagPipeTest)

  distanceBetweenGuides <- 100
  recArmLength<-150
  removeAddGeneHits <- TRUE
  browserWindow<-200
  ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
  ### target sequences
  seqs<-readDNAStringSet(system.file('extdata','rec_sequences.fa',package='cTagPipeTest'))

  mAID <- seqs$mAID
  tripleFlag<- seqs$triple_FLAG
  t2A<- seqs$t2A
  fiveP<- seqs$`5_prime_PITCH`
  threeP<- seqs$`3_prime_PITCH`
  neoM<- seqs$NeomycinR
  hygM<- seqs$HygromycinR


  ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
  ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
  ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
  ## Only run examples in interactive R sessions
  if (interactive()) {

    ui <- shiny::fluidPage(
      shiny::pageWithSidebar(
        shiny::headerPanel('C-terminal tagging pipeline'),
        shiny::sidebarPanel(
          shiny::selectInput('myGene', 'Transcript ID', c("NM_001304504","NM_153252","NM_001273",
                                                   "NM_000237","NM_002106","NM_138635","NM_014660",
                                                   "NM_175061","NM_030665","NM_005650","NM_001141969",
                                                   "NM_003325","NM_003496","NM_014034","NM_020713",
                                                   "NM_001203258","NM_004689","NM_006565","NM_080618",
                                                   "NM_002875","NM_133487","NM_003883","NM_021975","NM_001077700",
                                                   "NM_052927","NM_012308","NM_000937","NM_001004456","NM_001278215",
                                                   "NM_001077700","NM_005349")),
          #textInput('myGene', 'Transcript ID', value = "", width = NULL, placeholder = NULL),
          shiny::selectInput('myGenome', 'Genome', c("hg38")),
          shiny::numericInput('winSize', 'Window Size', 100, min = 50, max = 150),
          shiny::selectInput('myCassette', 'Rec. template', c("PITCH/homolgy_arm/AID/tripleFLAG/T2A/homolgy_arm/PITCH")),
          shiny::actionButton('run', 'Run'),
          width = 3
        ),

        shiny::mainPanel(#width=800,
          shiny::tabsetPanel(
            shiny::tabPanel("Analysis report",shiny::verbatimTextOutput("testout")),
            shiny::tabPanel("Plot", shiny::plotOutput("plot1")),
            shiny::tabPanel("CCTop summary", DT::dataTableOutput("ccSumTable", width = 800)),
            shiny::tabPanel("CCTop results", DT::dataTableOutput("ccbindTable", width = 800)),
            shiny::tabPanel("Guide pairs", DT::dataTableOutput("pairwiseSelections", width = 800)),
            shiny::tabPanel("Alignment test", msaR::msaROutput("msa", width="100%")),
            shiny::tabPanel("Sequences", shiny::verbatimTextOutput("seqSummary")),
            shiny::tabPanel("Mutation view guide 1", DT::dataTableOutput("mutView1", width = 1000)),
            shiny::tabPanel("Mutation view guide 2", DT::dataTableOutput("mutView2", width = 1000))
            )
        )))

    shiny::shinyApp(ui, server)
  }


