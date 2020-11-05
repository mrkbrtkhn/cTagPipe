server<-function(input, output, session) {

  options(DT.options = list(pageLength = 30, language = list(search = 'Filter:')))
  shiny::observeEvent(input$run,{
    input$goPlot # Re-run when button is clicked

    shiny::withProgress(message = 'loading genome', value = 0, {
      print("loading genome!")
      loadTxdb2(input$myGenome)->txdb})

    shiny::withProgress(message = 'retrieving gene information', value = 0.15, {
      print("getting gene mappings!")

      getMappings(txdb)->mGenIDMap
      #print(head(mGenIDMap))
      #print(input$myGene)
      mapRefSeqToSymbol(input$myGene,mGenIDMap)->mySymbol
      #print(mySymbol)
    })

    shiny::withProgress(message = 'extracting ranges and sequences', value = 0.25, {
      extractSeqAroundStop(txdb,input$myGene,input$winSize,input$myGenome)->genRanges
      isPlus<-as.character(strand(genRanges$cds))[1]=="+"
      print(isPlus)
      createSubGtf(genRanges,input$myGenome,browserWindow)->winTx
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
      save(v,file="/home/marek/ccTop/v.RData")

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
      plotBrowserView(v,genRanges,winTx,gSelect,browserWindow)
    })
    #  })

    output$testout<- shiny::renderText(paste("your gene:",mySymbol))

    output$ccSumTable <- DT::renderDataTable(
      DT::datatable(sumReport, options = list(scrollX='600px', scrollCollapse=TRUE), filter = 'top'))

    output$ccbindTable <- DT::renderDataTable(
      DT::datatable(bindSumReport, options = list(scrollX='600px', scrollCollapse=TRUE), filter = 'top'))

    output$pairwiseSelections <- DT::renderDataTable(data.table::data.table(gDat),
                                                     options = list(scrollX='600px', scrollCollapse=TRUE), filter = 'top')
    output$seqSummary <- renderText(paste(c("NEO:\n\n",as.character(recTempA),
                                            "\n\nHYGRO:\n\n",as.character(recTempB),
                                            "\n\nGuide 1:\n\n",
                                            v$data[[as.character(gSelect[1,1])]]$seq,
                                            "\n\nGuide 2:\n\n",
                                            v$data[[as.character(gSelect[1,2])]]$seq,sep="")))

    output$msa <- renderMsaR(msaR("xx.aln", menu=F, overviewbox = F))

    output$mutView1 <- DT::renderDataTable(DT::datatable(annotMuts(genRanges,input$myGenome,v$data[[1]])),
                                          options = list(scrollX='600px', scrollCollapse=TRUE), filter = 'top')

    output$mutView2 <- DT::renderDataTable(DT::datatable(annotMuts(genRanges,input$myGenome,v$data[[2]])),
                                          options = list(scrollX='600px', scrollCollapse=TRUE), filter = 'top')
  })
  #output$testout<- renderText(c(is.character(input$myGene),is.character(input$myGenome),is.numeric(input$winSize)))
}
