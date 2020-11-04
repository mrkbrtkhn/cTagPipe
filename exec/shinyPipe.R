  Sys.info()["sysname"]->mySys
  if (mySys == "Linux") {source('/data1/COMPUTING/R_data/chipseq_functions.R')}
  if (mySys == "Darwin") {source('/Users/marek/Computing/R_stuff/functions/chipseq_functions.R')}
  #library(CRISPRseek)
  library(RCurl)
  library(DT)
  library(msa)
  library(shiny)
  library(msaR)
  library(seqinr)
  #setwd("/Users/marek/Desktop/tagging")
  setwd("/Users/marek/Documents/Projects/andreas/tagging")
  #source("Crispr_functions_4.R")
  
  
  distanceBetweenGuides <- 100
  recArmLength<-150
  removeAddGeneHits <- TRUE
  
  ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
  ### target sequences
  #mAID <- DNAString(c("AAGGAGAAGAGTGCTTGTCCTAAAGATCCAGCCAAACCTCCGGCCAAGGCACAAGTTGTGGGATGGCCACCGGTGAGATCATACCGGAAGAACGTGATGGTTTCCTGCCAAAAATCAAGCGGTGGCCCGGAGGCGGCGGCGTTCGTGAAGGTATCAATGGACGGAGCACCGTACTTGAGGAAAATCGATTTGAGGATGTATAAA"))
  
  mAID <- Biostrings::DNAString(c("AAAGAGAAGAGTGCTTGTCCTAAAGATCCAGCCAAACCTCCGGCCAAGGCACAAGTTGTGGGATGGCCACCGGTGAGATCATACCGGAAGAACGTGATGGTTTCCTGCCAAAAATCAAGCGGTGGCCCGGAGGCGGCGGCGTTCGTGAAGGTATCAATGGACGGAGCACCGTACTTGAGGAAAATCGATTTGAGGATGTATAAA"))
  tripleFlag<-Biostrings::DNAString(c("GACTACAAGGATGACGACGATAAGGACTACAAGGATGACGACGATAAGGACTACAAGGATGACGACGATAAG"))
  t2A<-Biostrings::DNAString(c("GAGGGCAGAGGAAGTCTGCTAACATGCGGTGACGTCGAGGAGAATCCTGGACCT"))
  fiveP<-Biostrings::DNAString(c("GCATCGTACGCGTACGTGTTTGG"))
  threeP<-Biostrings::DNAString(c("CCAAACACGTACGCGTACGATGCG"))
  neoM<-Biostrings::DNAString(c("GGATCGGCCATTGAACAAGATGGATTGCACGCAGGTTCTCCGGCCGCTTGGGTGGAGAGGCTATTCGGCTATGACTGGGCACAACAGACAATCGGCTGCTCTGATGCCGCCGTGTTCCGGCTGTCAGCGCAGGGGCGCCCGGTTCTTTTTGTCAAGACCGACCTGTCCGGTGCCCTGAATGAACTGCAGGACGAGGCAGCGCGGCTATCGTGGCTGGCCACGACGGGCGTTCCTTGCGCAGCTGTGCTCGACGTTGTCACTGAAGCGGGAAGGGACTGGCTGCTATTGGGCGAAGTGCCGGGGCAGGATCTCCTGTCATCTCACCTTGCTCCTGCCGAGAAAGTATCCATCATGGCTGATGCAATGCGGCGGCTGCATACGCTTGATCCGGCTACCTGCCCATTCGACCACCAAGCGAAACATCGCATCGAGCGAGCACGTACTCGGATGGAAGCCGGTCTTGTCGATCAGGATGATCTGGACGAAGAGCATCAGGGGCTCGCGCCAGCCGAACTGTTCGCCAGGCTCAAGGCGCGCATGCCCGACGGCGATGATCTCGTCGTGACCCATGGCGATGCCTGCTTGCCGAATATCATGGTGGAAAATGGCCGCTTTTCTGGATTCATCGACTGTGGCCGGCTGGGTGTGGCGGACCGCTATCAGGACATAGCGTTGGCTACCCGTGATATTGCTGAAGAGCTTGGCGGCGAATGGGCTGACCGCTTCCTCGTGCTTTACGGTATCGCCGCTCCCGATTCGCAGCGCATCGCCTTCTATCGCCTTCTTGACGAGTTCTTCTGA"))
  hygM<-Biostrings::DNAString("ATGAAAAAGCCTGAACTCACCGCGACGTCTGTCGAGAAGTTTCTGATCGAAAAGTTCGACAGCGTCTCCGACCTGATGCAGCTCTCGGAGGGCGAAGAATCTCGTGCTTTCAGCTTCGATGTAGGAGGGCGTGGATATGTCCTGCGGGTAAATAGCTGCGCCGATGGTTTCTACAAAGATCGTTATGTTTATCGGCACTTTGCATCGGCCGCGCTCCCGATTCCGGAAGTGCTTGACATTGGGGAATTCAGCGAGAGCCTGACCTATTGCATCTCCCGCCGTGCACAGGGTGTCACGTTGCAAGACCTGCCTGAAACCGAACTGCCCGCTGTTCTGCAGCCGGTCGCGGAGGCCATGGATGCGATCGCTGCGGCCGATCTTAGCCAGACGAGCGGGTTCGGCCCATTCGGACCGCAAGGAATCGGTCAATACACTACATGGCGTGATTTCATATGCGCGATTGCTGATCCCCATGTGTATCACTGGCAAACTGTGATGGACGACACCGTCAGTGCGTCCGTCGCGCAGGCTCTCGATGAGCTGATGCTTTGGGCCGAGGACTGCCCCGAAGTCCGGCACCTCGTGCACGCGGATTTCGGCTCCAACAATGTCCTGACGGACAATGGCCGCATAACAGCGGTCATTGACTGGAGCGAGGCGATGTTCGGGGATTCCCAATACGAGGTCGCCAACATCTTCTTCTGGAGGCCGTGGTTGGCTTGTATGGAGCAGCAGACGCGCTACTTCGAGCGGAGGCATCCGGAGCTTGCAGGATCGCCGCGGCTCCGGGCGTATATGCTCCGCATTGGTCTTGACCAACTCTATCAGAGCTTGGTTGACGGCAATTTCGATGATGCAGCTTGGGCGCAGGGTCGATGCGACGCAATCGTCCGATCCGGAGCCGGGACTGTCGGGCGTACACAAATCGCCCGCAGAAGCGCGGCCGTCTGGACCGATGGCTGTGTAGAAGTACTCGCCGATAGTGGAAACCGACGCCCCAGCACTCGTCCGGAGGCAAAGGAATTCGGGAGATGGGGGAGGCTAACTGAAACACGGAAGGAGACAATACCGGAAGGAACCCGCGCTATGACGGCAATAAAAAGACAGAATAAAACGCACGGGTGTTGGGTCGTTTGTTCATAA")
  
  
  ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
  ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
  ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
  #options(width = 800)
  
  server<-function(input, output, session) {
    
    options(DT.options = list(pageLength = 30, language = list(search = 'Filter:')))
    shiny::observeEvent(input$run,{
      input$goPlot # Re-run when button is clicked
      
      shiny::withProgress(message = 'loading genome', value = 0, {
        print("loading genome!")
        loadTxdb(input$myGenome)->txdb})
      
      shiny::withProgress(message = 'retrieving gene information', value = 0.15, {
        print("getting gene mappings!")
        
        getMappings(input$myGenome)->mGenIDMap
        mapRefSeqToSymbol(input$myGene,mGenIDMap)->mySymbol
        
      })
      
      shiny::withProgress(message = 'extracting ranges and sequences', value = 0.25, {
        extractSeqAroundStop(txdb,input$myGene,input$winSize,input$myGenome)->genRanges
        isPlus<-as.character(strand(genRanges$cds))[1]=="+"
        createSubGtf(genRanges,input$myGenome)->winTx
      })
      
      shiny::withProgress(message = 'running CCTop', value = 0.4, {
        if (!file.exists(paste(input$myGene,".xls",sep=""))) {
          CCTop(name = input$myGene, radQ = "single", sequence = as.character(genRanges$sequence), pamType = "NGG", targetLength = "20", sgRNA5 = "NN",
                sgRNA3 = "NN", inVitroTx = "SP6", totalMismatches = "4", useCore = "on", coreLength = "12", coreMismatches = "2", 
                species = input$myGenome, downloadDir = ".")
        } else {print("file exists!")}
        parseCCTopData(paste(input$myGene,".xls",sep=""))->v
        correctCCTopGeneAnno(v,genRanges,mySymbol) -> v
        
        firstPassFilterCCResults(v)->v      
        evaluateUpCCTop(v,mySymbol,genRanges,input$winSize,mGenIDMap)->v
        
      })
      
      shiny::withProgress(message = 'calculating mutations', value = 0.5, {
        collectMutations(genRanges,v,input$myGenome)->v
        summarizeReport(v)->sumReport
        ### clean up
        sumReport<-sumReport[as.logical(sumReport$hasHitinTarget) &
                               (as.logical(unlist(sumReport$pamMutated)) | as.logical(unlist(sumReport$mutated)) | as.logical(unlist(sumReport$pamOverStop))) 
                             &ifelse(!is.na(sumReport$mintotalMMsExonic),as.numeric(as.character(sumReport$mintotalMMsExonic)),TRUE) 
                             ,]
        do.call(rbind,lapply(v$data,bindCCTop))->bindSumReport
        
      })
      
      shiny::withProgress(message = 'combining guides', value = 0.6, {
        combinedGuideAnalysis(sumReport,v)->gDat
        print(1)
        gSelect <- gDat[!gDat$guidesOverlap,][1,]
        extractArmIntervals(gSelect,v,genRanges,input$myGenome)->armIntervals
        mutateSeqs(gSelect,v,genRanges,armIntervals)->armIntervals
      })
      
      shiny::withProgress(message = 'compiling sequences', value = 0.75, {
        if (as.character(strand(genRanges$cds))[[1]]=="+") {
          recTempA <- xscat(fiveP,armIntervals$seqs[[1]],mAID,tripleFlag,t2A,neoM,armIntervals$seqs[[2]],threeP )
          recTempB <- xscat(fiveP,armIntervals$seqs[[1]],mAID,tripleFlag,t2A,hygM,armIntervals$seqs[[2]],threeP )
        } else {
          recTempA <- xscat(fiveP,reverseComplement(armIntervals$seqs[[2]]),mAID,tripleFlag,t2A,neoM,reverseComplement(armIntervals$seqs[[1]]),threeP )
          recTempB <- xscat(fiveP,reverseComplement(armIntervals$seqs[[2]]),mAID,tripleFlag,t2A,hygM,reverseComplement(armIntervals$seqs[[1]]),threeP )
        }
      })
      
      withProgress(message = 'alignment control', value = 0.85, {
        omegaAl <- alignmentControl(v,gSelect,recTempA,mAID,fiveP,threeP,genRanges,neoM,hygroM,t2A,isPlus)
        
      })
      
     # withProgress(message = 'creating genome view', value = 0.9, {
        output$plot1 <- renderPlot({  
          plotBrowserView(v,genRanges,winTx,gSelect)
        })
    #  })
      
      output$testout<- shiny::renderText(paste("your gene:",mySymbol))
        
      output$ccSumTable <- DT::renderDataTable(
        DT::datatable(sumReport, options = list(scrollX='600px', scrollCollapse=TRUE), filter = 'top'))
      
      output$ccbindTable <- DT::renderDataTable(
        DT::datatable(bindSumReport, options = list(scrollX='600px', scrollCollapse=TRUE), filter = 'top'))  
      
      output$pairwiseSelections <- DT::renderDataTable(data.table(gDat), 
                                                       options = list(scrollX='600px', scrollCollapse=TRUE), filter = 'top')
      output$seqSummary <- renderText(paste(c("NEO:\n\n",as.character(recTempA),
                                              "\n\nHYGRO:\n\n",as.character(recTempB),
                                              "\n\nGuide 1:\n\n",
                                              v$data[[as.character(gSelect[1,1])]]$seq,
                                              "\n\nGuide 2:\n\n",
                                              v$data[[as.character(gSelect[1,2])]]$seq,sep="")))
      
      output$msa <- renderMsaR(msaR("xx.aln", menu=F, overviewbox = F))
    })
    #output$testout<- renderText(c(is.character(input$myGene),is.character(input$myGenome),is.numeric(input$winSize)))
  }
  
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
                                                   "NM_052927","NM_012308","NM_000937","NM_001004456","NM_005349")),
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
            shiny::tabPanel("Sequences", shiny::verbatimTextOutput("seqSummary"))
            )
        )))
    
    shiny::shinyApp(ui, server)
  }
  
