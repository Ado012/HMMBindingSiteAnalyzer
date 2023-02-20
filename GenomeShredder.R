#genomeshredder

discard<-function()
{
  discardBorderFragmentStream <- file("discardedBorderFragments", 'w')
  trimmedFragmentStream <- file("trimmedGenomeFragments", 'w')
  
  
  untrimmedGenomeFragments = read.table("untrimmedGenomeFragments.bed", sep="\t")
  
  
  for(i in untrimmedGenomeFragments)
  {
    for(i in peakResult)
    {
      if (untrimmedGenomeFragments[i][,start]>=peakResults[i][,end] || untrimmedGenomeFragments[i][,end]<=peakResults[i][,end] )
        discard=1
      
      if (discard==1)
      {
        break
      }
    }
    
    if (discard==1)
    {write.table(data.frame(chromMW, targetStart ,targetEnd),append=TRUE, col.names= FALSE, row.names=FALSE, quote=FALSE, sep="\t",discardBorderFragmentStream)
      discard=0}
    
    else {
      write.table(data.frame(chromMW, targetStart ,targetEnd),append=TRUE, col.names= FALSE, row.names=FALSE, quote=FALSE, sep="\t",trimmedFragmentStream)
      
    }
    
  }
  
  
  
}


library(stringdist)
library(Biostrings)

setwd("/home/ace/Documents/Software_Projects/HMModeler")
#shreds genome files
GenomeShredderForSortedRegions<-function(regionCoordFile,rawSeqFile, observedCoordSeqSampleFile)
{
  switchChrom=0
  segmentSize=1000
  
  rawSeqStream <- file(rawSeqFile, 'w')
  
  #sample validation option
  observedCoordSeqSampleStream <- file(observedCoordSeqSampleFile, 'w')
  
  testFragmentTitle<-paste('fastatitle','targetSeq', 'transitionSeq', sep='\t')
  
  sampleTitleLine<-paste('targetSeq','startBase', 'endBase', sep='\t')
  writeLines(sampleTitleLine,observedCoordSeqSampleStream)
  writeLines(testFragmentTitle,rawSeqStream)
  
  #read in your coordinates from coverageFinder()
  sortedFrame <- read.csv(regionCoordFile, sep = '\t', header = T)


#read in chromosome data
  chromFiles <- as.character(list.files(path="Rdatafiles/chrom/"))
  
#set reference chromosome
  currentChrom=1
  
  #creating the chromasome name
#split chromsome file title
  fileNameSplit<-strsplit(chromFiles[1], split='[.]' )
  
  targetChromoName <- paste("Chr",fileNameSplit[[1]][5], sep="")
  

#load chromosome seq files
  chromFileName <- paste("Rdatafiles/chrom/",chromFiles[currentChrom], sep="")
  
  fileNum<-currentChrom
  
  print(nrow(sortedFrame))
  #for each region
   for (i in 1:nrow(sortedFrame))
   {
     #get chromosome of protoresult file
     
     regionChrom <- as.integer(substr(sortedFrame[i,][1],4,4)) #changed from chrom 021123
     #regionChrom<-as.integer(sortedFrame[i,][1]) won't work for some reason
     
     
     if (currentChrom != regionChrom) #if reference and protoresult chrom don't match
     {
       switchChrom=1
       
      currentChrom<-regionChrom
       
     fileNameSplit<-strsplit(chromFiles[regionChrom], split='[.]' ) #switch to chrom structure file indicated by region from result file
     
     print(fileNameSplit)
     fileNum<-regionChrom
     fileChrom<-fileNameSplit[[1]][5]
     
     
     if (fileChrom != currentChrom) #if new chrom doesn't match ref chrom
     {
       #if 
       while (fileChrom != currentChrom && i <= length(chromFiles))
       {
         
         fileNameSplit<-strsplit(chromFiles[i], split='[.]' )
         fileChrom<-fileNameSplit[[1]][5]
         fileNum<-i
         i<-i+1
       }
       
     }
     
     if (fileChrom==currentChrom)#if ref and new chrom match: create name variables
     {
       switchChrom=0
     targetChromoName <- paste("Chr",fileChrom, sep="")
     chromFileName <- paste("Rdatafiles/chrom/",chromFiles[fileNum], sep="")
     }
  
     }
     
     if (switchChrom==0) #if chroms match
     {#get start and end of seq from proto result line
     targetStart<-as.numeric(sortedFrame[i,]['startBase'])
     targetEnd<-as.numeric(sortedFrame[i,]['endBase'])
     
#retrieve sequence
     retrievedInfo <-
       readDNAStringSet(filepath = chromFileName, format ="fasta")
     
     retrievedSeq<-retrievedInfo[[1]]
     

#if target end is greater than sequence. Realign it.
     if(targetEnd>width(retrievedInfo))
     {
       targetEnd=width(retrievedInfo)
     transitionSeqLength=0
     }
     
     else 
       transitionSeqLength=segmentSize
     
     #get target sequence plus another 1kb to measure transitions

     
     #ERROR: does not check bounds so fails when end of chromosome is read, fix later
     targetSeq <-
       as.character(retrievedSeq[targetStart:targetEnd])
     transitionSeq<-as.character(retrievedSeq[targetEnd:(targetEnd+transitionSeqLength)])
     

     
     #write a new entry to the fasta consisting of a name and the sequence. 
     titleLine<-paste(">",targetChromoName,"_",targetStart,"_",targetEnd,sep="")
     

#write fetched target and transition sequence and write to results file
     testFragmentEntry<-paste(titleLine,targetSeq, transitionSeq, sep='\t')
     writeLines(testFragmentEntry,rawSeqStream)
     
     #writeLines(titleLine,rawSeqStream)
     #writeLines(targetSeq,rawSeqStream)
     
     #sample validation option
     sampleLine<-paste(targetSeq,targetStart, targetEnd, sep='\t')
     writeLines(sampleLine,observedCoordSeqSampleStream)
     
     }
     
     else 
     {
       print("No matching chromosome sorry")
       
     }
     
     

     
     
   }  
  
  
  
  close(rawSeqStream)
  close(observedCoordSeqSampleStream)
  
  
  
}




#bedtobamseextractor

  #function to extract sequences from peak files 
  BedSeqExtractor<-function()
  {
    
    peakSeqStream <- file("peakSeq", 'w')
    
    bedfile = read.table("Rdatafiles/NA_peaks.bed", sep=c("\t"," "))
    
    
    
    chromFiles <- as.character(list.files(path="Rdatafiles/chrom"))
    
    
#for each file in input 
    
    for (i in 1:length(chromFiles))
    {
      
      chromFileName <- paste("Rdatafiles/chrom/",chromFiles[i], sep="")
      
      retrievedInfo <-
        readDNAStringSet(filepath = chromFileName, format ="fasta")
      

      retrievedSeq <-retrievedInfo[[1]]
      
      
      fileNameSplit<-strsplit(chromFiles[i], split='[.]')
      targetChromoName <- paste("Chr",fileNameSplit[[1]][5], sep="")#check
      
      fileChrom=fileNameSplit[[1]][5]
      
#for each row in the bedfile 
      for (i in 1:nrow(bedfile))
      {
        bedlineChrom=as.character(bedfile[i,][1])
        
    
        print(bedlineChrom)
        print(fileChrom)
        chrommatch= amatch(bedlineChrom,fileChrom)
        
        if (is.na(chrommatch)==1)#amatch does not have its NA evaluated for some reason
          chrommatch=0
        
        if( chrommatch ==1 )
        {
          peakStart=as.integer(bedfile[i,][2])
          peakEnd=as.integer(bedfile[i,][3])
          
          print("peakCoord")
          print(peakStart)
          print(peakEnd)
          targetSeq=as.character(retrievedSeq[peakStart:peakEnd])
          print(unlist(targetSeq))
          #BAM write
          #write.table(data.frame(targetChromoName, peakStart ,peakEnd,width(targetSeq),targetSeq),append=TRUE, col.names= FALSE, row.names=FALSE, quote=FALSE, sep="\t",peakSeqStream)
          #FASTA write
          titleLine<-paste(">",targetChromoName,"_",peakStart,"_",peakEnd,sep="")
          writeLines(titleLine,peakSeqStream)
          writeLines(targetSeq,peakSeqStream)
        }
      }
      
      
    }
    
    close(peakSeqStream)
  }



#create model from peaks file with hmmrbuilder
#use model and fragment file in nhmmer
#convert back to bed



#print("comparison")
#print(bedfile[,i][1])
#print(fileNameSplit[[1]][5])
#print(width(retrievedSeq))


#print("peakcoord")
#print(peakStart)
#print(peakEnd)

