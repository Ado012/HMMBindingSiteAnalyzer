#coveragefinder
#currently only set up to handle one chromosome file fix later

library(stringdist)
library(Biostrings)

setwd("/home/ace/Documents/Software_Projects/HMModeler")

#shreds genome files
coverageFinder <- function()
{
  #set up thresholds for enriched and nonenriched regions
  clearThreshold=7
  enrichedThreshold=30
  resumeMark=0
  hits=0
  currentChrom=1
  j=1
  
  #main setup

  
  #find chromosome files
  chromFiles <- as.character(list.files(path = "Rdatafiles/chrom/"))
  
  #for each chromosome file
  for (i in 1:length(chromFiles))
  {
    
    #create dataframes for low and highly enriched coverage segments
    #can move this outside loop for single results file
    clearFrame <- data.frame(chromosome=integer(),startBase=integer(),endBase=integer(),stringsAsFactors=FALSE)
    
    enrichedFrame <- data.frame(chromosome=integer(),startBase=integer(),endBase=integer(),stringsAsFactors=FALSE)
    
    
    #create chromosome name
    fileNameSplit <- strsplit(chromFiles[i], split = '[.]')
    
    print(fileNameSplit)
    
    #extract BED formatted chromosome name
    targetChromoName <- paste("Chr", fileNameSplit[[1]][5], sep = "")
    
    #extract chromosome number
    chromNum<-fileNameSplit[[1]][5]
    currentChrom<-as.numeric(fileNameSplit[[1]][5])
    
    chromFileName <-
      paste("Rdatafiles/chrom/", chromFiles[i], sep = "")
    
    #extract chromosome
    retrievedInfo <-
      readDNAStringSet(filepath = chromFileName, format = "fasta")
    
    retrievedSeq <- retrievedInfo[[1]]
    
    #set starting position and size of target fragment
    targetStart = 1
    targetEnd = 1000
    k = 1
    l = 1
    #open file of WUS binding peaks
    peakfile <- read.csv("Data/sig1.bed", sep = '\t', header = T)
    
    
    
    #take a chunk of dna and compare to fragments in bam file
    #1:nrow(peakfile)
      while (j < 40000000 && targetEnd < width(retrievedInfo)) 
      {
        #see if the peak's chromosome still matches up with the target sequences chromosome
        peakChrom=as.numeric(peakfile[j,][1]) 
        
    
        
        
        if (peakChrom != currentChrom)
          break
        
        #track fragment overlaps
        if (peakfile[j,][2] > targetStart &&
            peakfile[j,][2] < targetEnd ||
            peakfile[j,][3] > targetStart && peakfile[i,][3] < targetEnd || peakfile[j,][2] < targetStart && peakfile[j,][3] > targetEnd)
        {
          #Note hit on stream
          hits<-hits+1
        }
        
        
        #once the peak end goes beyond the target fragment mark the location to resume the search for overlaps once we move on to the next fragment. Beware may cause permanent loop
        if (peakfile[j,][3] > targetEnd && resumeMark == 0)
        {
          resumeMark = j
        }
        
        #once the beginning if the fragments passes the target fragment. It is safe to switch to aother one
        if (peakfile[j,][2] > targetEnd)
        {
          if (hits > enrichedThreshold)
          {
            enrichedFrameRow <- data.frame(targetChromoName,targetStart, targetEnd)
            names(enrichedFrameRow)<-c("chromosome","startBase","endBase")
            
            enrichedFrame <- rbind(enrichedFrame, enrichedFrameRow)
            
            print("enriched")
            print(k)
            k<-k+1
          }
          
          if (hits < clearThreshold)
          {
            clearFrameRow <- data.frame(targetChromoName,targetStart, targetEnd)
            names(clearFrameRow)<-c("chromosome","startBase","endBase")
            
            clearFrame <- rbind(clearFrame, clearFrameRow)
            l<-l+1
          }
          
          hits = 0
          j = resumeMark
          resumeMark = 0
          #set range to search around gene
          targetStart = targetStart + 1000
          targetEnd = targetEnd + 1000
          
          
          if(targetEnd>width(retrievedInfo))
            targetEnd=width(retrievedInfo)
          
         
          
        }
          
        if (peakChrom !=1)
        {
          print("chromosome matching")
          print(currentChrom)
          print(peakChrom)
        }
          
      j<-j+1
      }
      
      
    
    
    clearResultsFileName <-
      paste("clearCoordinates",targetChromoName, sep = "")
    
    enrichedResultsFileName <-
      paste("enrichedCoordinates",targetChromoName, sep = "")
    
    
    write.table(enrichedFrame,enrichedResultsFileName, row.names = FALSE, sep = "\t", quote=FALSE)
    
    write.table(clearFrame,clearResultsFileName, row.names = FALSE, sep = "\t", quote=FALSE)
    
    
  }
  
}




    
    
    
    
    
  
  
  


