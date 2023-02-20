#copyright of Albert Do, do not copy or distribute or use without permission




#Umbrella function for constructing an HMM
HMMModeler <- function(inputFile,
                       fastaFileTitle,
                       transitionFileTitle)
{
  setwd("/home/ace/Documents/Software_Projects/HMModeler")
  
  inputFrame <- read.csv(inputFile, sep = '\t', header = T)
  segSize=250
  
  #create generic dataframe to store emissions and state transitions.
  state <-
    c(
      'CoreEmission',
      'MisCoreEmission',
      'ComplexCoreEmission',
      'SpacerEmission',
      'NoiseEmission',
      'toHighCoreTransit',
      'toMedCoreTransit',
      'toLowCoreTransit'
    )
  probabilities <- c(0, 0, 0, 0, 0, 0, 0, 0)
  
  targetcol <- as.data.frame(inputFrame[, 'targetSeq'])
  transitioncol <- as.data.frame(inputFrame[, 'transitionSeq'])
  titlecol <- as.data.frame(inputFrame[, 'fastatitle'])
  
  #convert target and transition columns from input file 
  dataRowNum <- nrow(inputFrame)
  DataFrame2FastaConverter(targetcol, titlecol, fastaFileTitle)
  DataFrame2FastaConverter(transitioncol, titlecol, transitionFileTitle)
  
  
  
  #read in input sequence
  retrievedSeq <- readDNAStringSet(fastaFileTitle, format = "fasta")
  transitionSeq <-
    readDNAStringSet(transitionFileTitle, format = "fasta")
  
  #create lists to store emission and state transition counts for each state
  highCoreState <- data.frame(state, probabilities)
  medCoreState <- data.frame(state, probabilities)
  lowCoreState <- data.frame(state, probabilities)
  
  #for each sequence in the input...
  for (i in 1:length(retrievedSeq))
  {
    #translate and count emissions
    write(i, file = "looptracker",
          append = TRUE)
    write(length(retrievedSeq), file = "looptracker",
          append = TRUE)
    
    retrievedSegment <- retrievedSeq[[i]] #take a chunk of the sequence (already divided but segSize provides an additional layer of control)
    segmentBegin = 1
    segmentEnd = segSize
    #while you are still in the segment. Currently should only loop once with 1000 segSize
    while (segmentEnd <= length(retrievedSegment))
      #introduce to measure state changes
    {
      retrievedFragment <- retrievedSegment[segmentBegin:segmentEnd]
      
      translatedSeq <-
        Translator(retrievedFragment) #translate raw sequences to cores, miscores, complexcores, noise, and spacer
      hmmCount <- HMMEmissionCounter(translatedSeq)
      print("print 1")

      
      
      if (segmentBegin > 1)#first sequence doesn't have an emitter of its own
        hiddenStateParent = hiddenState
      
      #calculations for core enrichment
      allEmissions = sum(hmmCount[, 2])
      coreEmissions = hmmCount[, 2][1] + hmmCount[, 2][3] #add regular and complex cores
      coreEnrichment = as.numeric(coreEmissions / allEmissions)
      
      #print(highCoreState[,2][1])
      #print(hmmCount[,2][1])
      
      #update emission counts for each state.
      if (coreEnrichment > 0.05)
      {
        #use j instead of i so nested loops don't clash
        for (j in 1:nrow(hmmCount))
          highCoreState[2][j, ] = highCoreState[2][j, ] + hmmCount[, 2][j] #fill out the emissions for an individual highcore into the highcore basket
        
        hiddenState <- 'highCore'
      }
      
      else if (coreEnrichment > 0.01)
      {
        for (j in 1:nrow(hmmCount))
          medCoreState[2][j, ] = medCoreState[2][j, ] + hmmCount[, 2][j]
        
        hiddenState <- 'medCore'
      }
      
      else
      {
        for (j in 1:nrow(hmmCount))
          lowCoreState[2][j, ] = lowCoreState[2][j, ] + hmmCount[, 2][j]
        
        hiddenState <- 'lowCore'
      }
      
      
      
      
      #update state transitions if beyond first sequence, the last three entries in the emission probabilities dataframes
      if (segmentBegin > 1)
      {
        updatedCoreStateList <-
          HMMStateSorter(
            hmmCount,
            hiddenState,
            hiddenStateParent,
            highCoreState,
            medCoreState,
            lowCoreState
          )
        
        print(updatedCoreStateList)
        print("updatedCoreStateListPart")
        print(updatedCoreStateList[1])
        
        highCoreState <- updatedCoreStateList[[1]]
        medCoreState <- updatedCoreStateList[[2]]
        lowCoreState <- updatedCoreStateList[[3]]
        
      }
      
      segmentBegin = segmentBegin + segSize
      segmentEnd = segmentEnd + segSize
      
      
    }
    
  }
  
  print(highCoreState)
  #Convert Emission and State transition counts into percentages
  
  highCoreState <- Percenthmm(stateEmissionCounts = highCoreState)
  
  medCoreState <- Percenthmm(medCoreState)
  
  lowCoreState <- Percenthmm(lowCoreState)
  
  #bind together to create overall HMM
  hmmodel <- bind_rows(highCoreState, medCoreState, lowCoreState)
  
  hmmodel
}

