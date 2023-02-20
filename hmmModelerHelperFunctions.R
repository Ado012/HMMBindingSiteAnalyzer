library(stringr)
library(Biostrings)
library(dplyr)
library(HMM)




DataFrame2FastaConverter<-function(dataframecol1,dataframecol2,convertedFileName)
{
  convertedStream <- file(convertedFileName, 'w')
  
  for(i in 1:nrow(dataframecol1))
  {
    if (i==1)
      write.table(unlist(as.character(dataframecol2[i,])),file=convertedStream, col.names=FALSE, row.names=FALSE, append=FALSE, quote=FALSE)
    else 
      write.table(unlist(as.character(dataframecol2[i,])),file=convertedStream, col.names=FALSE, row.names=FALSE, append=TRUE, quote=FALSE)
    
    write.table(unlist(as.character(dataframecol1[i,])),file=convertedStream, col.names=FALSE, row.names=FALSE, append=TRUE, quote=FALSE)
    
  }
  
  close(convertedStream)
  
}


DataFrame2FastaConverterFile<-function(inputFile,convertedFileName)
{
  
  inputFrame <- read.csv(inputFile, sep = '\t', header = T)
  
  convertedStream <- file(convertedFileName, 'w')
  
  targetcol<-as.data.frame(inputFrame[,'targetSeq'])
  #transitioncol<-as.data.frame(inputFrame[,'transitionSeq'])
  titlecol<-as.data.frame(inputFrame[,'fastatitle'])
  
  for(i in 1:nrow(inputFrame))
  {
    if (i==1)
      write.table(unlist(as.character(titlecol[i,])),file=convertedStream, col.names=FALSE, row.names=FALSE, append=FALSE, quote=FALSE)
    else 
      write.table(unlist(as.character(titlecol[i,])),file=convertedStream, col.names=FALSE, row.names=FALSE, append=TRUE, quote=FALSE)
    
    write.table(unlist(as.character(targetcol[i,])),file=convertedStream, col.names=FALSE, row.names=FALSE, append=TRUE, quote=FALSE)
    
  }
  
  close(convertedStream)
  
}




DataFrame2FastaConverterFileMSB<-function(inputFrame,combinedFileName,centerFileName,leftFileName,rightFileName)
{
  
  
  
  combinedStream <- file(combinedFileName, 'w')
  centerStream <- file(centerFileName, 'w')
  leftStream <- file(leftFileName, 'w')
  rightStream <- file(rightFileName, 'w')
  
  
  for(i in 1:nrow(inputFrame))
  {
    
    seqProtoName <- paste('>',targetlist[i,][2], sep="")
    seqName  <- paste(seqProtoName,"_C", sep="")
      write.table(seqName,file=combinedStream, col.names=FALSE, row.names=FALSE, append=FALSE, quote=FALSE)
      write.table(unlist(targetlist[i,][5]),file=combinedStream, col.names=FALSE, row.names=FALSE, append=FALSE, quote=FALSE)
    
      write.table(seqName,file=centerStream, col.names=FALSE, row.names=FALSE, append=FALSE, quote=FALSE)
      write.table(unlist(targetlist[i,][5]),file=centerStream, col.names=FALSE, row.names=FALSE, append=FALSE, quote=FALSE)
      
      
      seqProtoName <- paste('>',targetlist[i,][2], sep="")
      seqName  <- paste(seqProtoName, "_L", sep="")
      write.table(seqName,file=combinedStream, col.names=FALSE, row.names=FALSE, append=FALSE, quote=FALSE)
      write.table(unlist(targetlist[i,][6]),file=combinedStream, col.names=FALSE, row.names=FALSE, append=FALSE, quote=FALSE)
      
      write.table(seqName,file=leftStream, col.names=FALSE, row.names=FALSE, append=FALSE, quote=FALSE)
    write.table(unlist(targetlist[i,][6]),file=leftStream, col.names=FALSE, row.names=FALSE, append=FALSE, quote=FALSE)
    
    
    seqProtoName <- paste('>',targetlist[i,][2], sep="")
    seqName  <- paste(seqProtoName, "_R", sep="")
    write.table(seqName,file=combinedStream, col.names=FALSE, row.names=FALSE, append=FALSE, quote=FALSE)
    write.table(unlist(targetlist[i,][7]),file=combinedStream, col.names=FALSE, row.names=FALSE, append=FALSE, quote=FALSE)
    
   
    write.table(seqName,file=rightStream, col.names=FALSE, row.names=FALSE, append=FALSE, quote=FALSE)
    write.table(unlist(targetlist[i,][7]),file=rightStream, col.names=FALSE, row.names=FALSE, append=FALSE, quote=FALSE)
    
    
  }
  
  close(combinedStream)
  close(centerStream)
  close(leftStream)
  close(rightStream)
}







#hmmbuilder
HMMEmissionCounter<-function(corematrix)
{
  #instantiate hmm dataframe
  state<-c('Core','Miscore','ComplexCore','Spacer','Noise')
  stateInstances<-c(0,0,0,0,0)
  CoreEmission<-c(0,0,0,0,0)
  MisCoreEmission<-c(0,0,0,0,0)
  ComplexCoreEmission<-c(0,0,0,0,0)
  SpacerEmission<-c(0,0,0,0,0)
  NoiseEmission<-c(0,0,0,0,0)
  
  
  hmm <- data.frame(state,stateInstances,CoreEmission,MisCoreEmission,ComplexCoreEmission,SpacerEmission,NoiseEmission)
  
  
  
  #loop through and for each look at the succeding entry and fill in hmm dataframe
  
  for (i in 1:nrow(corematrix)-1)#find a more consistent way to handle edges other than deleting the last emission which maybe is the reason for the -1
  {
    
    #in the future take into account coresize
    coresize<-(corematrix[i,][,'start']-corematrix[i,][,'end'])/4
    
    corestate<-toString(corematrix[i,][,'state'])
    emittedstate<-toString(corematrix[i+1,][,'state'])
    
    
    
    
    if (corestate=='core')
    {
      hmm[1,][2]=hmm[1,][2]+1
      
      
      if (emittedstate=='core')
        hmm[1,][3]=hmm[1,][3]+1
      else if (emittedstate=='miscore' || emittedstate=='miscoregroup')
        hmm[1,][4]=hmm[1,][4]+1
      else if (emittedstate=='complexcore')
        hmm[1,][5]=hmm[1,][5]+1
      else if (emittedstate=='spacer')
        hmm[1,][6]=hmm[1,][6]+1
      else if (emittedstate=='noise')
        hmm[1,][7]=hmm[1,][7]+1
    }
    
    else if(corestate=='miscore' || corestate=='miscoregroup')
    {
      hmm[2,][2]=hmm[2,][2]+1
      
      if (emittedstate=='core')
        hmm[2,][3]=hmm[2,][3]+1
      else if (emittedstate=='miscore' || emittedstate=='miscoregroup')
        hmm[2,][4]=hmm[2,][4]+1
      else if (emittedstate=='complexcore')
        hmm[2,][5]=hmm[2,][5]+1
      else if (emittedstate=='spacer')
        hmm[2,][6]=hmm[2,][6]+1
      else if (emittedstate=='noise')
        hmm[2,][7]=hmm[2,][7]+1
    }
    
    
    else if (corestate=='complexcore')
    {
      hmm[3,][2]=hmm[3,][2]+1
      
      if (emittedstate=='core')
        hmm[3,][3]=hmm[3,][3]+1
      else if (emittedstate=='miscore' || emittedstate=='miscoregroup')
        hmm[3,][4]=hmm[3,][4]+1
      else if (emittedstate=='complexcore')
        hmm[3,][5]=hmm[3,][5]+1
      else if (emittedstate=='spacer')
        hmm[3,][6]=hmm[3,][6]+1
      else if (emittedstate=='noise')
        hmm[3,][7]=hmm[3,][7]+1
    }
    
    else if (corestate=='spacer')
    {
      hmm[4,][2]=hmm[4,][2]+1
      
      if (emittedstate=='core')
        hmm[4,][3]=hmm[4,][3]+1
      else if (emittedstate=='miscore' || emittedstate=='miscoregroup')
        hmm[4,][4]=hmm[4,][4]+1
      else if (emittedstate=='complexcore')
        hmm[4,][5]=hmm[4,][5]+1
      else if (emittedstate=='spacer')
        hmm[4,][6]=hmm[4,][6]+1
      else if (emittedstate=='noise')
        hmm[4,][7]=hmm[4,][7]+1
    }
    
    
    
    else if (corestate=='noise')
    {
      hmm[5,][2]=hmm[5,][2]+1
      
      if (emittedstate=='core')
        hmm[5,][3]=hmm[5,][3]+1
      else if (emittedstate=='miscore' || emittedstate=='miscoregroup')
        hmm[5,][4]=hmm[5,][4]+1
      else if (emittedstate=='complexcore')
        hmm[5,][5]=hmm[5,][5]+1
      else if (emittedstate=='spacer')
        hmm[5,][6]=hmm[5,][6]+1
      else if (emittedstate=='noise')
        hmm[5,][7]=hmm[5,][7]+1
    }
    
    
    
  }
  
  hmm
  
}










#Counts up state transitions
HMMStateSorter<-function(hmmCount, hiddenState, hiddenStateParent, highCoreState, medCoreState, lowCoreState)
{
  print(hiddenState)
  print(hiddenStateParent)
  
  
  #if previous state was highcore increment what states it transitions to. 
  if (hiddenState=='highCore')
  {
    if (hiddenStateParent=='highCore')
      highCoreState[,2][6]=highCoreState[,2][6]+1
    else if (hiddenStateParent=='medCore')
      highCoreState[,2][7]=highCoreState[,2][7]+1
    else if (hiddenStateParent=='lowCore')
      highCoreState[,2][8]=highCoreState[,2][8]+1
  }
  
  
  #if previous state was medcore increment what states it transitions to. 
  if (hiddenState=='medCore')
  {
    if (hiddenStateParent=='highCore')
      medCoreState[,2][6]=medCoreState[,2][6]+1
    else if (hiddenStateParent=='medCore')
      medCoreState[,2][7]=medCoreState[,2][7]+1
    else if (hiddenStateParent=='lowCore')
      medCoreState[,2][8]=medCoreState[,2][8]+1
    
  }
  
  #if previous state was lowcore increment what states it transitions to. 
  if (hiddenState=='lowCore')
  {
    if (hiddenStateParent=='highCore')
      lowCoreState[,2][6]=lowCoreState[,2][6]+1
    else if (hiddenStateParent=='medCore')
      lowCoreState[,2][7]=lowCoreState[,2][7]+1
    else if (hiddenStateParent=='lowCore')
      lowCoreState[,2][8]=lowCoreState[,2][8]+1
    
  }
  
  
  #return updated counts
  sortedCoreStateList<-list(highCoreState,medCoreState,lowCoreState)
  
  
}





Percenthmm<-function(stateEmissionCounts)#converts emission counts to percentages
{
  print(stateEmissionCounts)
  
  
  totalEmission<-sum(stateEmissionCounts[,2][1:5])#retrieve total emissions
  
  #For each entry. Divide by total 
  for (i in 1:5)
  {
    if (totalEmission==0)
      stateEmissionCounts[,2][i]=0
    
    else
      stateEmissionCounts[,2][i]<-as.numeric(stateEmissionCounts[,2][i]/totalEmission)
  }
  
  totalTransmission<-sum(stateEmissionCounts[,2][6:8])
  
  for (i in 6:8)
  {
    if (totalTransmission==0)
      stateEmissionCounts[,2][i]=0
    
    else
      stateEmissionCounts[,2][i]<-as.numeric(stateEmissionCounts[,2][i]/totalTransmission)
  }
  
  
  
  stateEmissionCounts
  
}


