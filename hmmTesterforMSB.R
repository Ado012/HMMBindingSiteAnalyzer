

#umbrella function for testing sequences


HMMTest<-function(controlHMM,peakHMM,observedSeqFile, rawResultsFileName)
{
  
  
  resultsDataFrame <- data.frame(
    seqNameFull = c(),
    seqNameBase = c(),
    rawSeq = c(),
    crmModelProb = c(),
    crmModelLengthPredict = c(),
    nonModelProb = c(),
    nonModelLengthPredict = c(),
    enrichedStrength1 = c(),
    enrichedStrength2 = c(),
    predict = c()
  )
  
  
  
  #nonvalidation option
  #titleLine<-paste("translatedSeq","rawSeq","crmModelProb","crmModelLengthPredict","nonModelProb","nonModelLengthPredict", "prediction", sep="\t")
  j=0
  #validation option
  #colnames(resultsDataFrame)<-c("translatedSeq","rawSeq","crmModelProb","crmModelLengthPredict","nonModelProb","nonModelLengthPredict", "prediction")
  
  #writeLines(titleLine,rawResultsStream)
 
   retrievedSeqList<-readDNAStringSet(observedSeqFile, format="fasta")
  
  
  classifiedSeqList<-c()
  #for each sequence in input
  for (i in 1:length(retrievedSeqList))
  {
    print("test start")
    #translate sequence
    translatedSeq<-Translator(retrievedSeqList[i])
    print(translatedSeq)
    translatedSeq[translatedSeq=='miscoregroup']<-'miscore' #adding a result<- screws it up
    print("bob")
   # print(translatedSeq)
    print("bob2")
    #print(translatedSeq[1])
    #format sequence for hmm calculations
    translatedSeq<-toString(translatedSeq[,1])
    translatedSeq<-strsplit(translatedSeq, ', ')
    translatedSeqFormatted<-unlist(translatedSeq)
   #pull out names of sequences
    fullname <- names(retrievedSeqList[i])
    basename <- strsplit(fullname,"_")
    basename <- unlist(basename)
    basename <- basename[1]
    #calculate probabilities    
    seqRating<-HMMCalculator(translatedSeqFormatted,controlHMM,peakHMM,retrievedSeqList[i])
    
    j=j+1
    
    resultslist <- strsplit(seqRating, split = "\t" )
    resultsunlist <- unlist(resultslist)
  
    resultsEntry <-
      list(
        "fullName" = fullname,
        "baseName" = basename,
        "rawSeq" = resultsunlist[2],
        "crmModelProb" = resultsunlist[3],
        "crmModelLengthPredict" = resultsunlist[4],
        "nonModelProb" = resultsunlist[5],
        "nonModelLengthPredict" = resultsunlist[6],
        "enrichedStrength1" = as.numeric(resultsunlist[7]),
        "enrichedStrength2" = as.numeric(resultsunlist[8]),
        "predict" = resultsunlist[9]
      )
    
    
    #arrange info into data frame  
    resultsEntry<-as.data.frame(do.call(cbind,resultsEntry))
    
    print("motifClusterdf")
    print(resultsEntry)
    
    resultsDataFrame<-rbind(resultsDataFrame,resultsEntry )
  }
  
  
  write.table(resultsDataFrame, rawResultsFileName, append = FALSE, sep = " ", dec = ".",
              row.names = TRUE, col.names = TRUE)

  
resultsDataFrame    
}


