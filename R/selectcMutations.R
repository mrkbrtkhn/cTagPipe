selectcMutations<-function(genRanges,exonContext,deBug=FALSE) {
  
  # GENETIC_CODE
  pamMutated<-FALSE
  mutated<-FALSE
  pamOverStop<-FALSE
  
  mutList<-list()
  
  ### first try to mutate the PAM
  # is the PAM not across the coding region, or STOP or splice acceptor site: mutate it
  if (all((exonContext[(min(which(!is.na(exonContext$pamSeq)))):(max(which(!is.na(exonContext$pamSeq)))),]$aa)=="none") & 
      any(is.na(exonContext[!is.na(exonContext$pamSeq),]$splice)) & !any(exonContext[!is.na(exonContext$pamSeq),]$aa == "*") ) {
    exonContext[which(!is.na(exonContext$pamSeq)),]->pamContext
    list(type="non-coding",coodinate=as.numeric(pamContext$coordinate[1]),base=as.character(reverseComplement(DNAString(pamContext$genomeSequence[1]))),
      mutation=paste("mutates",as.character(pamContext$genomeSequence[1]),"at position",pamContext$coordinate[1],"to",as.character(reverseComplement(DNAString(pamContext$genomeSequence[1])))))->mutList[[as.character(pamContext$coordinate[1])]]
    pamMutated<-TRUE
  }
  if (deBug) {print("made it here")}
  
  ### if PAM is over stop codon
  ### don't touch it as the stop gets deleted anyway
  if (!pamMutated & !mutated &!pamOverStop) {
    stopContext<-exonContext[exonContext$aa=="*",]
    if (nrow(stopContext[!is.na(stopContext$pamSeq),]) ) {
      mutList[["pamOverStop"]]<-list(type="pamOverStop")
      pamOverStop<-TRUE
    }
  }
  
  ### try to mutate PAM  when in CDS
  if (!pamMutated & !mutated & !pamOverStop) {
    pamContext<-exonContext[!is.na(exonContext$pamSeq),]
    pamContext<-exonContext[exonContext$aaPos %in% pamContext$aaPos[!is.na(pamContext$aaPos)],]
    ### can we mutate triplet as such that PAM is defunct?
    if (nrow(pamContext)>0) {
      i<-pamContext$aaPos[!duplicated(pamContext$aa)][1]
      ### iterate over triplets
      while ((!pamMutated | !mutated |!pamOverStop) & i %in% pamContext$aaPos[!duplicated(pamContext$aa)] ) {
        pamContext[pamContext$aaPos == i,]->sPamContext
        myAA<-sPamContext$aa[1]
        pamOverTriplet<-sPamContext$tripletPos[!is.na(sPamContext$pamSeq)]
        actualTriplet <- (sPamContext$base)
        ### find triplet that differs at any pamOverTriplet Position
        ### e.g. mutate any G that keeps aa constant
        names(GENETIC_CODE[GENETIC_CODE==myAA])->myTriplets
        j<-1

        ### iterate over candidate synonymous triplets
        while (j<= length(myTriplets) & (!pamMutated)) {
          candTriplet<-strsplit(myTriplets[j],"")[[1]]
          if (any((actualTriplet != candTriplet)[c(pamOverTriplet)]) & sum(actualTriplet != candTriplet)==1 ) {
            ### found mutation of GG keeping amino acid!!!
            which(!(actualTriplet == candTriplet))->gToChange
            tripletContext<-subset(exonContext,aaPos==i)
            
            ### DEFINE MUTATION!!!!!            ###             ###             ### 
            mutCoord<-tripletContext[gToChange,]$coordinate
            mutBase<-ifelse(tripletContext[gToChange,]$base==tripletContext[gToChange,]$genomeSequence,
                            candTriplet[gToChange],
                            as.character(reverseComplement(DNAString(candTriplet[gToChange])))
            )
            mutList[[as.character(mutCoord)]]<-list(type=c("PAM in CDS"),coodinate=as.numeric(mutCoord),
                                                 base=mutBase,
                                                 mutation=paste("mutates PAM sequence",tripletContext[gToChange,]$genomeSequence,"at position",as.numeric(mutCoord),"to",mutBase,
                                                                "keeping aa",myAA,"(",paste(actualTriplet,collapse=""),"to",paste(candTriplet,collapse=""),")"))
            pamMutated<-TRUE
            
          }
          j<-j+1
        }
        i<-i+1
      }
    }
    
  }
  ### otherwise mutate coding sequence accordingly
  ### start from PAM onwards (not yet implemented)
  exonContext[!is.na(exonContext$aaPos) & !is.na(exonContext$guideSeq),]->gExonContext
  if ((!pamMutated & !mutated & !pamOverStop) & nrow(gExonContext)>18) {
    
    ### can we mutate triplet maximally?
    head(exonContext)
    gAaPos <- exonContext[!is.na(exonContext$guideSeq),]$aaPos
    
    #gExonContext <- exonContext[exonContext$aaPos %in% gAaPos, ]
    if (nrow(gExonContext)>0) {
      aaPos<-gExonContext$aaPos[!duplicated(gExonContext$aa)]
      aaPos<-aaPos[!is.na(aaPos)]
      i<- min(aaPos)
      #print(i)
      nOfMuts <- 0
      while (nOfMuts < 6 & i <= max(aaPos)) {
        gExonContext[gExonContext$aaPos == i,]->sExonContext
          myAA<-sExonContext$aa[1]
          guideOverTriplet<-sExonContext$tripletPos[!is.na(sExonContext$guideSeq)]
          actualTriplet <- (sExonContext$base)
          ### find triplet that differs maximally at guideOverTriplet Positions 
          names(GENETIC_CODE[GENETIC_CODE==myAA])->myTriplets
          mutScan<-list()
          for (j in 1:length(myTriplets)) {
            candTriplet<-strsplit(myTriplets[j],"")[[1]]
            
            if (sum(!(actualTriplet == candTriplet))==1) {
              which(!(actualTriplet == candTriplet))[1]->baseToChange
              tripletContext<-subset(exonContext,aaPos==i)
              ### DEFINE MUTATION!!!!!            ###             ###             ### 
              mutCoord<-tripletContext[baseToChange,]$coordinate[1]
              
              mutBase<-ifelse((tripletContext[baseToChange,]$base==tripletContext[baseToChange,]$genomeSequence),
                              candTriplet[baseToChange],
                              as.character(reverseComplement(DNAString(candTriplet[baseToChange])))
              )
              
              
              
              mutScan[[j]]<-list(type=c("alternative triplet"),coodinate=as.numeric(mutCoord),base=mutBase,
                                 mutation=paste("mutates",paste(actualTriplet,collapse=""),"to",
                                                paste(candTriplet,collapse="") ,"at position",as.numeric(mutCoord)),
              mutCount=sum((actualTriplet != candTriplet)[c(guideOverTriplet)]))
            } else {
              mutScan[[j]]<-list(mutCount=0,
                                 mutations=list())
            }
            
          }
          ### evaluate mutations
          if (any(lapply(mutScan,function(x) return (x$mutCount))>0)) {
            mutSelect<-which(lapply(mutScan,function(x) return (x$mutCount))==max(unlist(lapply(mutScan,function(x) return (x$mutCount)))))[1]
            #print(mutScan[[mutSelect]])
            id<-paste(mutScan[[mutSelect]]$coodinate,subset(exonContext,coordinate==mutScan[[mutSelect]]$coodinate)$genomeSequence,"to",
                      mutScan[[mutSelect]]$base,sep="_")
            mutList[[id]]<-mutScan[[mutSelect]]
            nOfMuts<-nOfMuts+mutScan[[mutSelect]]$mutCount
          }
          i<-i+1
          
      }
      
    }
    if (length(mutList)>5) {mutated<-TRUE}
  }
  ### special cases?
  ### PAM not mutable due to to overlap with splice site and majority of guide is over intron!
  if ((!pamMutated & !mutated & !pamOverStop)) {
    iExonContext<-exonContext[is.na(exonContext$pamSeq) & is.na(exonContext$splice) & !is.na(exonContext$coordinate) & (exonContext$aa=="none") & !is.na(exonContext$guideSeq),]
    if (nrow(iExonContext)>5) {
      ### randomly mutate 6 bases
      sample(rownames(iExonContext),6)->ids
      mutScan<-list()
      for (i in ids) {
        print(i)
        mutBase<-as.character(complement(DNAString(iExonContext[i,]$genomeSequence)))
        id<-paste(as.numeric(iExonContext[i,]$coordinate),iExonContext[i,]$genomeSequence,"to",mutBase,sep="_")
        mutList[[id]]<-list(type=c("intronic exchange"),coodinate=as.numeric(iExonContext[i,]$coordinate),base=mutBase,
                           mutation=paste(iExonContext[i,]$genomeSequence,"to",mutBase),
                           mutCount=6)
      }
    }
    mutated<-TRUE
  } 
    
  
  
  return(list(pamMutated=pamMutated,mutated=mutated,pamOverStop=pamOverStop,mutList=mutList))
}
