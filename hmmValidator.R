#validates predictions with processed peak enrichment data
HMMValidator<-function(rawResultsFile, enrichedRegionsFile, clearRegionsFile, refinedResultsFile, rawResultsCoordFile)
{
  #set up loop tracking variables
  enriched=0
  clear=0
  correct=0
  incorrect=0
  
  #open input files and start writing file for refined results
  
  rawResultsFrame <- read.csv(rawResultsFile, sep = '\t', header = T)
  
  rawResultsCoordFrame <- read.csv(rawResultsCoordFile, sep = '\t', header = T)
  
  
  enrichedRegionsFrame <-
    read.csv(enrichedRegionsFile, sep = '\t', header = T)
  
  clearRegionsFrame <- read.csv(clearRegionsFile, sep = '\t', header = T)
  
  refinedResultsStream <- file(refinedResultsFile, 'w')
  
  
  #write title line to refined results file
  titleLine <-
    paste(
      'raw sequence',
      'translated sequence',
      'startBase',
      'endBase',
      "crmModelProb",
      "crmModelLengthPredict",
      "nonModelProb",
      "nonModelLengthPredict",
      'prediction',
      'identity',
      'correctPrediction',
      sep = "\t"
    )
  writeLines(titleLine, refinedResultsStream)
  
  #for each line in prediction file
  for (i in 1:nrow(rawResultsFrame))
  {#examine to see if tested segment lines up with an enriched segement
    for (j in 1:nrow:(enrichedRegionsFrame))
    {
      
      if (rawResultCoordFrame['startBase'][i, ] >= enrichedRegionsFrame['startBase'][j, ] &&
          rawResultsCoordFrame['startBase'][i, ] <= enrichedRegionsFrame['endBase'][j, ] ||
          rawResultCoordFrame['endBase'][i, ] >= enrichedRegionsFrame['startBase'][j, ] &&
          rawResultsCoordFrame['endBase'][i, ] <= enrichedRegionsFrame['endBase'][j, ])
      {
        enriched = 1
        status='enriched'
        break
      }
      
      if (rawResultCoordFrame['endBase'][i, ] < enrichedRegionsFrame['startBase'][j, ])
      {
        status='indeterminate'
        break
      }
      
    }
    #examine to see if tested segment lines up with clear segement
    for (k in 1:nrow:(clearRegionsFrame))
    {
      if (coordRawResultFrame['startBase'][i, ] >= clearRegionsFrame['startBase'][k, ] &&
          coordRawResultsFrame['startBase'][i, ] <= clearRegionsFrame['endBase'][k, ] ||
          coordRawResultFrame['endBase'][i, ] >= clearRegionsFrame['startBase'][k, ] &&
          coordRawResultsFrame['endBase'][i, ] <= clearRegionsFrame['endBase'][k, ])
      {
        clear = 1
        status='clear'
        break
      }
      
      if (coordRawResultFrame['endBase'][i, ] < clearRegionsFrame['startBase'][j, ])
      {
        status='indeterminate'
        break
      }
      
    }
    #if the tested segment overlaps with a segment that is the same as its prediction it counts as a correct prediction
    if (enriched==1 && rawResultsFrame['prediction'][i,]=='enriched' || clear==1 && rawResultsFrame['prediction'][i,]=='clear')
    {
      correct=correct+1
      correctPredict=1
    }
    #otherwise its an incorrect prediction
    else if (enriched==1 && rawResultsFrame['prediction'][i,]=='clear' || clear==1 && rawResultsFrame['prediction'][i,]=='enriched')
    {
      incorrect=incorrect+1
      correctPredict=0
    }
    else
    {
      unclassified=unclassified=+1
      correctPredict=0
    }
    
    enriched=0
    clear=0
    
    
    #write prediction results to refined results file
    refinedResultLine <-
      paste(
        rawResultsFrame['rawSeq'][i,],
        rawResultsFrame['translatedSeq'][i,],
        coordRawResultsFrame['startBase'][i,],
        coordRawResultsFrame['endBase'][i,],
        
        rawResultsFrame['crmModelProb'][i,],
        rawResultsFrame['crmModelLengthPredict'][i,],
        rawResultsFrame['nonModelProb'][i,],
        rawResultsFrame['nonModelLengthPredict'][i,],
        rawResultsFrame['prediction'][i,],
        status,
        correctPredict,
        sep = "\t"
      )
    writeLines(refinedResultLine, refinedResultsStream)
    
  }
  #close refined results file and calculate total accuracy
  close(refinedResultsStream)
  totalPredictions=unclassified+correct+incorrect
  accuracy=correct/totalPredictions
}