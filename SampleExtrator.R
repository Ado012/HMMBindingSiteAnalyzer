
SampleExtractor('enrichedSegmentsResults','clearSegmentsResults','enrichedSegmentsResultsAbridged','clearSegmentsResultsAbridged','observedSegmentRaw')





SampleExtractor<-function(inputFile,inputFile2, abridgedFileName, abridgedFileName2, sampleFileName)
{
  
  inputFrame <- read.csv(inputFile, sep = '\t', header = T)
  inputFrame2 <- read.csv(inputFile2, sep = '\t', header = T)
  sample<-c(5,15,26,33,45,59,61,77,82,91)
  
  
  #convertedStream <- file(saFileName, 'w')
  
  sampleFrame <- data.frame( "state"= character(0), "start" = integer(0), "end" = integer(0))

  for(i in 1:10)
  {
    sampleFrame<-rbind(sampleFrame,inputFrame[(sample[i]),])
    
    
  }
  
  
  for(i in 1:10)
  {
    sampleFrame<-rbind(sampleFrame,inputFrame2[(sample[i]),])
    
    
  }
  
  
  
  abridgedFrame <- inputFrame[-sample, ]
  abridgedFrame2 <- inputFrame2[-sample, ]
 
   colnames(sampleFrame)<-c('fastatitle','targetSeq','transitionSeq')
  
  
  write.table(abridgedFrame,abridgedFileName, col.names=TRUE, row.names=FALSE, append=FALSE, quote=FALSE, sep='\t')
  write.table(abridgedFrame2,abridgedFileName2, col.names=TRUE, row.names=FALSE, append=FALSE, quote=FALSE, sep='\t')
  write.table(sampleFrame,sampleFileName, col.names=TRUE, row.names=FALSE, append=FALSE, quote=FALSE, sep='\t')
  
  #close(convertedStream)
  
}