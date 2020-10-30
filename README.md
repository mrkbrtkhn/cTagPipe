cTagePipe\_readme
================

## About cTagPipe

This application aims at the identification of DNA sequences that allow
the C-terminal tagging of endogenous loci in order to express fusion
proteins. To achieve this, a combination of CRISPR/Cas9-mediated
induction of DNA breaks and homologous recombination mediated repair
with engineered repair templates is applied.

The application mainly works as a Shiny app but can be run from the R
command line, too. The user selects the gene of interest and the
pipeline calculates CRISPR guide RNAs, selects recombination templates
and modifies/mutates the templates if required as such, that the guide
RNAs will not be able to allow cutting of the recombined locus.

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
install.packages(newPackages)
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
                                                        "NM_002875","NM_133487","NM_003883","NM_021975",
                                                        "NM_001077700","NM_052927","NM_012308","NM_000937",
                                                        "NM_001004456","NM_001278215",
                                                        "NM_001077700")),
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
    ##  [1] cTagPipeTest_0.1.0                roxygen2_6.1.1.9000              
    ##  [3] knitr_1.21                        rmarkdown_1.11                   
    ##  [5] BSgenome.Hsapiens.UCSC.hg38_1.4.1 BSgenome_1.44.0                  
    ##  [7] Gviz_1.20.0                       org.Hs.eg.db_3.4.1               
    ##  [9] shiny_1.2.0                       GenomicFeatures_1.28.0           
    ## [11] AnnotationDbi_1.38.2              Biobase_2.36.2                   
    ## [13] rtracklayer_1.36.6                GenomicRanges_1.28.2             
    ## [15] GenomeInfoDb_1.12.0               seqinr_3.4-5                     
    ## [17] msaR_0.3.0                        msa_1.8.0                        
    ## [19] Biostrings_2.44.0                 XVector_0.16.0                   
    ## [21] IRanges_2.10.5                    S4Vectors_0.14.7                 
    ## [23] BiocGenerics_0.22.1               DT_0.5                           
    ## [25] RCurl_1.95-4.12                   bitops_1.0-6                     
    ## [27] devtools_1.13.1                  
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] colorspace_1.3-2              biovizBase_1.24.0            
    ##  [3] htmlTable_1.9                 base64enc_0.1-3              
    ##  [5] dichromat_2.0-0               rstudioapi_0.9.0             
    ##  [7] bit64_0.9-7                   interactiveDisplayBase_1.14.0
    ##  [9] xml2_1.2.0                    splines_3.4.0                
    ## [11] ade4_1.7-13                   jsonlite_1.6                 
    ## [13] Formula_1.2-1                 Rsamtools_1.28.0             
    ## [15] cluster_2.0.6                 compiler_3.4.0               
    ## [17] httr_1.3.1                    backports_1.1.3              
    ## [19] assertthat_0.2.0              Matrix_1.2-10                
    ## [21] lazyeval_0.2.0                later_0.8.0                  
    ## [23] acepack_1.4.1                 htmltools_0.3.6              
    ## [25] prettyunits_1.0.2             tools_3.4.0                  
    ## [27] gtable_0.2.0                  GenomeInfoDbData_0.99.0      
    ## [29] Rcpp_1.0.0                    ape_5.1                      
    ## [31] nlme_3.1-131                  crosstalk_1.0.0              
    ## [33] xfun_0.4                      stringr_1.4.0                
    ## [35] mime_0.5                      ensembldb_2.0.1              
    ## [37] XML_3.98-1.16                 AnnotationHub_2.8.1          
    ## [39] zlibbioc_1.22.0               MASS_7.3-47                  
    ## [41] scales_0.4.1                  VariantAnnotation_1.22.0     
    ## [43] BiocInstaller_1.26.1          hms_0.3                      
    ## [45] promises_1.0.1                ProtGenerics_1.8.0           
    ## [47] SummarizedExperiment_1.6.1    AnnotationFilter_1.0.0       
    ## [49] RColorBrewer_1.1-2            yaml_2.2.0                   
    ## [51] curl_2.6                      memoise_1.1.0                
    ## [53] gridExtra_2.2.1               ggplot2_2.2.1                
    ## [55] biomaRt_2.39.2                rpart_4.1-11                 
    ## [57] latticeExtra_0.6-28           stringi_1.3.1                
    ## [59] RSQLite_2.1.1                 checkmate_1.8.2              
    ## [61] BiocParallel_1.10.1           commonmark_1.7               
    ## [63] rlang_0.3.1                   pkgconfig_2.0.2              
    ## [65] matrixStats_0.52.2            evaluate_0.10                
    ## [67] lattice_0.20-35               GenomicAlignments_1.12.1     
    ## [69] htmlwidgets_1.3               bit_1.1-14                   
    ## [71] plyr_1.8.4                    magrittr_1.5                 
    ## [73] R6_2.4.0                      Hmisc_4.0-3                  
    ## [75] DelayedArray_0.2.3            DBI_1.0.0                    
    ## [77] pillar_1.1.0                  foreign_0.8-68               
    ## [79] withr_2.1.2                   survival_2.41-3              
    ## [81] nnet_7.3-12                   tibble_1.4.2                 
    ## [83] crayon_1.3.4                  progress_1.2.0               
    ## [85] data.table_1.12.0             blob_1.1.1                   
    ## [87] git2r_0.18.0                  digest_0.6.18                
    ## [89] xtable_1.8-2                  httpuv_1.4.5.1               
    ## [91] munsell_0.4.3

<br><br> <br><br>
