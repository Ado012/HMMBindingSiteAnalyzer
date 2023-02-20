library(dplyr)


#chromosome picking function
chromPicker <- function(targetGeneChromo)
{
  #Load Proper Chromosome File
  
  if (targetGeneChromo == 'M')
  {
    retrievedSeq <-
      readDNAStringSet(filepath = "Rdatafiles/Arabidopsis_thaliana.TAIR10.dna.chromosome.MT.fa", format =
                         "fasta")
    targetGeneChromoName <- 'ChrM'
  }
  else if (targetGeneChromo == 1)
  {
    retrievedSeq <-
      readDNAStringSet(filepath = "Rdatafiles/Arabidopsis_thaliana.TAIR10.dna.chromosome.1.fa", format =
                         "fasta")
    targetGeneChromoName <- 'Chr1'
  }
  else if (targetGeneChromo == 2)
  {
    retrievedSeq <-
      readDNAStringSet(filepath = "Rdatafiles/Arabidopsis_thaliana.TAIR10.dna.chromosome.2.fa", format =
                         "fasta")
    targetGeneChromoName <- 'Chr2'
  }
  else if (targetGeneChromo == 3)
  {
    retrievedSeq <-
      readDNAStringSet(filepath = "Rdatafiles/Arabidopsis_thaliana.TAIR10.dna.chromosome.3.fa", format =
                         "fasta")
    targetGeneChromoName <- 'Chr3'
  }
  else if (targetGeneChromo == 4)
  {
    retrievedSeq <-
      readDNAStringSet(filepath = "Rdatafiles/Arabidopsis_thaliana.TAIR10.dna.chromosome.4.fa", format =
                         "fasta")
    targetGeneChromoName <- 'Chr4'
  }
  else if (targetGeneChromo == 5)
  {
    retrievedSeq <-
      readDNAStringSet(filepath = "Rdatafiles/Arabidopsis_thaliana.TAIR10.dna.chromosome.5.fa", format =
                         "fasta")
    targetGeneChromoName <- 'Chr5'
  }
  else if (targetGeneChromo == 'C')
  {
    retrievedSeq <-
      readDNAStringSet(filepath = "Rdatafiles/Arabidopsis_thaliana.TAIR10.dna.chromosome.Pt.fa", format =
                         "fasta")
    targetGeneChromoName <- 'ChrC'
  }
  else
    print("NO recognized chromosome")
  
  
  #may fail if no recognized chromasome.
  retrievedSeqInfo <- list("sequence"= retrievedSeq, "chrname"=targetGeneChromoName)
  
  
}








#main function
targetAnnotatorHMM <- function(targetfile="downregulatedwusgenes_cyclo.csv", resultsOutput="coreMotifOutput")
{
  fileStream <- file(resultsOutput, 'w')
  
  #read in annotation file: Main source of info on the genes we'll be scanning
  geneAnnotations <-
    read.csv(
      "Rdatafiles/TAIR10_GFF3_genes_transposons.csv",
      header = TRUE,
      sep = ','
    )
  #Break up gene annotation data into relevant variables.
  geneChromList <- unlist(geneAnnotations[1]) #chromosome of gene feature
  geneStartList <- unlist(geneAnnotations[4]) #start position of gene features
  geneEndList <- unlist(geneAnnotations[5]) #end position of gene feature
  geneNamesData <- geneAnnotations[9] #a column listingt gene names the feature is associated with
  
  
  
  targetsAnnotated<- data.frame(
    chromosome = c(),
    gene = c(),
    chainstart = c(),
    chainend = c(),
    seqCenter = c(),
    seqLeft = c(),
    seqRight = c()
  )
  
  #read in file of gene targets that will be examined
  targetdata <-
    read.csv(targetfile,
             header = TRUE,
             sep = ',')
  
  targets <- targetdata[2]
  targetsul <- unlist(targets)
  targetEntry <- as.character(targetsul)
  
  targetMeta<- unlist(targetdata[3])
  targetMetachar<-as.character(targetMeta)
  
  
  geneNamesDataUL<- unlist(geneNamesData[1])
  
  #get the names of the gene transcripts.
  geneNamesList <- lapply(geneNamesDataUL, function(geneEntry)
  {
    geneEntryChar<-as.character(geneEntry)
    geneNameDataSplit <- strsplit(geneEntryChar, "[[:punct:][:space:]]+")[[1]]
    
    #print("Gene Data")
    #print(geneEntryChar)
    
    geneName <- geneNameDataSplit[2]
    
    
  })
  
  
  #extract targets one by one
  for(j in 1:length(targetEntry))
  {
    
    
    target = targetEntry[j]
    
    
    #find target gene in Arabidopsis gene list
    targetLocs <- grep(target, geneNamesList)
    
    #take the first hit, possibly modify later
    targetLoc <- targetLocs[1]
    
    
    #load gene end and beginning and search range
    seqStart <- geneStartList[targetLoc]
    seqEnd  <- geneEndList[targetLoc]
    
    seqStart <- as.integer(seqStart)
    seqEnd <- as.integer(seqEnd)
    
    chrom <- substr(target,3,3)
    chrom <- as.integer(chrom)
    
  
    #load proper chromasome file and find chrom location.
    arabSeqInfo <- chromPicker(targetGeneChromo = chrom)
    
    
    #extract sequence
    arabSeq <- arabSeqInfo$sequence
    
    arabSeqUnlist <- unlist(arabSeq[1])
    
    targetArabSeqCenter <- arabSeqUnlist[seqStart:seqEnd]
    
    if (seqStart - 1000 > 1)
    targetArabSeqLeft <-arabSeqUnlist[(seqStart-1000):seqStart]
    
    else 
      targetArabSeqLeft <-arabSeqUnlist[1:seqStart]
    
    #if (seqEnd + 1000 < length(arabSeqUnlist))
    targetArabSeqRight <- arabSeqUnlist[seqEnd:(seqEnd+1000)]
    
 #   else 
  #    targetArabSeqRight <-arabSeqUnlist[seqEnd:length(arabSeqUnlist)]
    
    print("locations")
    print(targetArabSeqRight)
    
    
    #+3 adjust end to account for motif size
    #Bundle info to write to file
    annotationEntry <-
      list(
        "chromosome"=chrom,
        "gene"=target,
        "seqStart"=seqStart,
        "seqEnd"=seqEnd,
        "seqCenter" = toString(targetArabSeqCenter),
        "seqLeft" = toString(targetArabSeqLeft),
        "seqRight" = toString(targetArabSeqRight)
      )
    
    print("motifCluster")
    print(annotationEntry)
    
    #arrange info into data frame  
    annotationEntry<-as.data.frame(do.call(cbind,annotationEntry))
    
    print("motifClusterdf")
    print(annotationEntry)
    
    targetsAnnotated<-rbind(targetsAnnotated,annotationEntry )
    
    
    #write info to file
    #write.table(motifClusterdf,append=TRUE, col.names=FALSE, row.names=FALSE, quote=FALSE, sep=' ', fileStream)
    print("check")
    
    
    
  }
  
  
  targetsAnnotated
}



