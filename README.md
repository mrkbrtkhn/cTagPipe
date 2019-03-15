cTagePipe\_readme
================

## Installation

Please install the latest version of R from <https://www.r-project.org>.
On Ubuntu you will have to install a number of system packages in order
to make the subsequent installation of R packages possible.

``` sh
sudo apt-get install libssl-dev
sudo apt-get install libcurl4-openssl-dev
sudo apt-get install libxml2-dev
```

In order to install this package, you require R’s devtools:

``` r
listOfPackages <- c("devtools")
newPackages <- listOfPackages[!(listOfPackages %in% installed.packages()[,"Package"])]
if (length(newPackages)>0) {
  source("https://bioconductor.org/biocLite.R")
BiocInstaller::biocLite(newPackages)
}
```

This enables you to install directly from GitHub repository:

``` r
install_github("mrkbrtkhn/cTagPipe")
```

Most likely this will produce a lot of errors. Dependencies are still a
mess and I currently have no time to work on them. Please install them
manually.

## Running cTagPipe functionality:

Assuming you have managed to install cTagPipe, the following code will
start a Shiny server with the functionality as it is yet implemented:

``` r
library(cTagPipeTest)

distanceBetweenGuides <- 100
  recArmLength<-150
  removeAddGeneHits <- TRUE
  browserWindow<-200
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
                                                   "NM_001077700")),
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
```

## Session Information

``` r
sessionInfo()
```

    ## R version 3.4.0 (2017-04-21)
    ## Platform: x86_64-apple-darwin15.6.0 (64-bit)
    ## Running under: macOS  10.13.6
    ## 
    ## Matrix products: default
    ## BLAS: /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib
    ## LAPACK: /Library/Frameworks/R.framework/Versions/3.4/Resources/lib/libRlapack.dylib
    ## 
    ## locale:
    ## [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
    ## 
    ## attached base packages:
    ##  [1] grid      stats4    parallel  stats     graphics  grDevices utils    
    ##  [8] datasets  methods   base     
    ## 
    ## other attached packages:
    ##  [1] cTagPipeTest_0.1.0                knitr_1.21                       
    ##  [3] rmarkdown_1.11                    BSgenome.Hsapiens.UCSC.hg38_1.4.1
    ##  [5] BSgenome_1.44.0                   Gviz_1.20.0                      
    ##  [7] org.Hs.eg.db_3.4.1                shiny_1.2.0                      
    ##  [9] GenomicFeatures_1.28.0            AnnotationDbi_1.38.2             
    ## [11] Biobase_2.36.2                    rtracklayer_1.36.6               
    ## [13] GenomicRanges_1.28.2              GenomeInfoDb_1.12.0              
    ## [15] seqinr_3.4-5                      msaR_0.3.0                       
    ## [17] msa_1.8.0                         Biostrings_2.44.0                
    ## [19] XVector_0.16.0                    IRanges_2.10.5                   
    ## [21] S4Vectors_0.14.7                  BiocGenerics_0.22.1              
    ## [23] DT_0.5                            RCurl_1.95-4.12                  
    ## [25] bitops_1.0-6                      devtools_1.13.1                  
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] nlme_3.1-131                  ProtGenerics_1.8.0           
    ##  [3] matrixStats_0.52.2            bit64_0.9-7                  
    ##  [5] RColorBrewer_1.1-2            progress_1.2.0               
    ##  [7] httr_1.3.1                    backports_1.1.3              
    ##  [9] tools_3.4.0                   R6_2.4.0                     
    ## [11] rpart_4.1-11                  Hmisc_4.0-3                  
    ## [13] DBI_1.0.0                     lazyeval_0.2.0               
    ## [15] colorspace_1.3-2              ade4_1.7-13                  
    ## [17] nnet_7.3-12                   withr_2.1.2                  
    ## [19] gridExtra_2.2.1               prettyunits_1.0.2            
    ## [21] bit_1.1-14                    curl_2.6                     
    ## [23] compiler_3.4.0                git2r_0.18.0                 
    ## [25] htmlTable_1.9                 DelayedArray_0.2.3           
    ## [27] checkmate_1.8.2               scales_0.4.1                 
    ## [29] stringr_1.4.0                 digest_0.6.18                
    ## [31] Rsamtools_1.28.0              foreign_0.8-68               
    ## [33] dichromat_2.0-0               base64enc_0.1-3              
    ## [35] pkgconfig_2.0.2               htmltools_0.3.6              
    ## [37] ensembldb_2.0.1               htmlwidgets_1.3              
    ## [39] rlang_0.3.1                   rstudioapi_0.9.0             
    ## [41] RSQLite_2.1.1                 BiocInstaller_1.26.1         
    ## [43] crosstalk_1.0.0               BiocParallel_1.10.1          
    ## [45] acepack_1.4.1                 VariantAnnotation_1.22.0     
    ## [47] magrittr_1.5                  GenomeInfoDbData_0.99.0      
    ## [49] Formula_1.2-1                 Matrix_1.2-10                
    ## [51] Rcpp_1.0.0                    munsell_0.4.3                
    ## [53] ape_5.1                       yaml_2.2.0                   
    ## [55] stringi_1.3.1                 MASS_7.3-47                  
    ## [57] SummarizedExperiment_1.6.1    zlibbioc_1.22.0              
    ## [59] AnnotationHub_2.8.1           plyr_1.8.4                   
    ## [61] blob_1.1.1                    promises_1.0.1               
    ## [63] crayon_1.3.4                  lattice_0.20-35              
    ## [65] splines_3.4.0                 hms_0.3                      
    ## [67] pillar_1.1.0                  biomaRt_2.39.2               
    ## [69] XML_3.98-1.16                 evaluate_0.10                
    ## [71] biovizBase_1.24.0             latticeExtra_0.6-28          
    ## [73] data.table_1.12.0             httpuv_1.4.5.1               
    ## [75] gtable_0.2.0                  assertthat_0.2.0             
    ## [77] ggplot2_2.2.1                 xfun_0.4                     
    ## [79] mime_0.5                      xtable_1.8-2                 
    ## [81] AnnotationFilter_1.0.0        later_0.8.0                  
    ## [83] survival_2.41-3               tibble_1.4.2                 
    ## [85] GenomicAlignments_1.12.1      memoise_1.1.0                
    ## [87] cluster_2.0.6                 interactiveDisplayBase_1.14.0

<br><br> <br><br>
