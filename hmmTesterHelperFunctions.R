library(stringr)
library(Biostrings)
library(dplyr)
library(HMM)

#initializes and calculates probabilies for sequences being generated 
#no validation option
#HMMCalculator<-function(observationSequence,hmmModelCRM,hmmModelNon,rawSequence)

#validation option
HMMCalculator<-function(observationSequence,hmmModelCRM,hmmModelNon,rawSequence)
{
  startingOptionsCRM=3
  startingOptionsNon=3
  m=1
  n=1
  crmProbFinalSum=0
  nonProbFinalSum=0
  finalCRMColumn=0
  finalNonColumn=0
  
  
  # Initialise HMM
  hmmModelCRMTest = initHMM(c("HighCore","MedCore","LowCore"), c("core","miscore","complexcore","spacer","noise"), transProbs=matrix(c(as.numeric(hmmModelCRM[2][6,]),as.numeric(hmmModelCRM[2][7,]),as.numeric(hmmModelCRM[2][8,]),
                                                                                                                                       as.numeric(hmmModelCRM[2][14,]),as.numeric(hmmModelCRM[2][15,]),as.numeric(hmmModelCRM[2][16,]),
                                                                                                                                       as.numeric(hmmModelCRM[2][22,]),as.numeric(hmmModelCRM[2][23,]),as.numeric(hmmModelCRM[2][24,])),nrow=3,ncol=3,byrow=TRUE),
                            emissionProbs=matrix(c(as.numeric(hmmModelCRM[2][1,]),as.numeric(hmmModelCRM[2][2,]),as.numeric(hmmModelCRM[2][3,]),as.numeric(hmmModelCRM[2][4,]),as.numeric(hmmModelCRM[2][5,]),
                                                   as.numeric(hmmModelCRM[2][9,]),as.numeric(hmmModelCRM[2][10,]),as.numeric(hmmModelCRM[2][11,]),as.numeric(hmmModelCRM[2][12,]),as.numeric(hmmModelCRM[2][13,]),
                                                   as.numeric(hmmModelCRM[2][17,]),as.numeric(hmmModelCRM[2][18,]),as.numeric(hmmModelCRM[2][19,]),as.numeric(hmmModelCRM[2][20,]),as.numeric(hmmModelCRM[2][21,])),nrow=3,ncol=5,byrow=TRUE))
  
  #if state doesn't show up then starting probabilities should be zeroed out in model
  if (hmmModelCRM[2][1,]==0 && hmmModelCRM[2][2,]==0 && hmmModelCRM[2][3,]==0 && hmmModelCRM[2][4,]==0 && hmmModelCRM[2][5,]==0 )
  { hmmModelCRMTest$startProbs[1]=0
  startingOptionsCRM = startingOptionsCRM-1
  }
  
  if (hmmModelCRM[2][9,]==0 && hmmModelCRM[2][10,]==0 && hmmModelCRM[2][11,]==0 && hmmModelCRM[2][12,]==0 && hmmModelCRM[2][13,]==0 )
  {
    hmmModelCRMTest$startProbs[2]==0
    startingOptionsCRM = startingOptionsCRM-1
  }
  
  if (hmmModelCRM[2][17,]==0 && hmmModelCRM[2][18,]==0 && hmmModelCRM[2][19,]==0 && hmmModelCRM[2][20,]==0 && hmmModelCRM[2][21,]==0 )
  {
    hmmModelCRMTest$startProbs[3]==0
    startingOptionsCRM = startingOptionsCRM-1
  }
  
  if (startingOptionsCRM==2)
  {
    for (m in 1:3)
    {
      if (hmmModelCRMTest$startProbs[m] != 0)
        hmmModelCRMTest$startProbs[m]=0.5
    }
    
  }
  
  if (startingOptionsCRM==1)
  {
    for (m in 1:3)
    {
      if (hmmModelCRMTest$startProbs[m] != 0)
        hmmModelCRMTest[m]=1
    }
    
  }
  
  hmmModelNonTest = initHMM(c("HighCore","MedCore","LowCore"), c("core","miscore","complexcore","spacer","noise"), transProbs=matrix(c(as.numeric(hmmModelNon[2][6,]),as.numeric(hmmModelNon[2][7,]),as.numeric(hmmModelNon[2][8,]),
                                                                                                                                       as.numeric(hmmModelNon[2][14,]),as.numeric(hmmModelNon[2][15,]),as.numeric(hmmModelNon[2][16,]),
                                                                                                                                       as.numeric(hmmModelNon[2][22,]),as.numeric(hmmModelNon[2][23,]),as.numeric(hmmModelNon[2][24,])),nrow=3,ncol=3,byrow=TRUE),
                            emissionProbs=matrix(c(as.numeric(hmmModelNon[2][1,]),as.numeric(hmmModelNon[2][2,]),as.numeric(hmmModelNon[2][3,]),as.numeric(hmmModelNon[2][4,]),as.numeric(hmmModelNon[2][5,]),
                                                   as.numeric(hmmModelNon[2][9,]),as.numeric(hmmModelNon[2][10,]),as.numeric(hmmModelNon[2][11,]),as.numeric(hmmModelNon[2][12,]),as.numeric(hmmModelNon[2][13,]),
                                                   as.numeric(hmmModelNon[2][17,]),as.numeric(hmmModelNon[2][18,]),as.numeric(hmmModelNon[2][19,]),as.numeric(hmmModelNon[2][20,]),as.numeric(hmmModelNon[2][21,])),nrow=3,ncol=5,byrow=TRUE))
  
  #if state doesn't show up then starting probabilities should be zeroed out in model
  if (hmmModelNon[2][1,]==0 && hmmModelNon[2][2,]==0 && hmmModelNon[2][3,]==0 && hmmModelNon[2][4,]==0 && hmmModelNon[2][5,]==0 )
  {
    hmmModelNonTest$startProbs[1]=0
    startingOptionsNon = startingOptionsNon-1
  }
  
  if (hmmModelNon[2][9,]==0 && hmmModelNon[2][10,]==0 && hmmModelNon[2][11,]==0 && hmmModelNon[2][12,]==0 && hmmModelNon[2][13,]==0 )
  {
    hmmModelNonTest$startProbs[2]=0
    startingOptionsNon = startingOptionsNon-1
  }
  
  if (hmmModelNon[2][17,]==0 && hmmModelNon[2][18,]==0 && hmmModelNon[2][19,]==0 && hmmModelNon[2][20,]==0 && hmmModelNon[2][21,]==0 )
  {
    hmmModelNonTest$startProbs[3]=0
    startingOptionsNon = startingOptionsNon-1
  }
  
  if (startingOptionsNon==2)
  {
    for (n in 1:3)
    {
      if (hmmModelNonTest$startProbs[n] != 0)
        hmmModelNonTest$startProbs[n]=0.5
    }
    
  }
  
  if (startingOptionsNon==1)
  {
    for (n in 1:3)
    {
      if (hmmModelNonTest$startProbs[n] != 0)
        hmmModelNonTest$startProbs[n]=1
    }
    
  }
  
  #forward algorithm to calculate probability chart, calculates the probability of the observed sequence under the model
  CRMProb<-forward(hmmModelCRMTest, observationSequence)
  
  NonProb<-forward(hmmModelNonTest, observationSequence)
  
  print(CRMProb)
  
  CRMProb<-exp(CRMProb)
  NonProb<-exp(NonProb)
  
  #go through probability chart for model and find last nonzero column
  for (i in 1:ncol(CRMProb))
  {
    crmProbColumnSum<-sum(CRMProb[,i])
    
    #if column is zero break from loop with the last nonzero columnb
    if (crmProbColumnSum <= 0)
      break
    #otherwise set as the new nonzero column found
    else 
    {
      crmProbFinalSum<-crmProbColumnSum
      finalCRMColumn<-i
    }
    
  }
  
  for (i in 1:ncol(NonProb))
  {
    nonProbColumnSum<-sum(NonProb[,i])
    
    #if column is zero break from loop with the last nonzero columnb
    if (nonProbColumnSum <= 0)
      break
    #otherwise set as the new nonzero column found
    else 
    {
      nonProbFinalSum<-nonProbColumnSum
      finalNonColumn<-i
    }
    
  }
  

#determine enrichment experience
    if (crmProbFinalSum > nonProbFinalSum)
      prediction='enriched'
    else if (crmProbFinalSum < nonProbFinalSum)
      prediction='clear'
    else 
    {
      if (finalCRMColumn > finalNonColumn)
        prediction='enriched'
      else if (finalCRMColumn < finalNonColumn)
        prediction='clear'
      else
        prediction='indeterminate'
    
     }
  
  
  enrichedStrength1 <-finalCRMColumn - finalNonColumn
  enrichedStrength2 <- crmProbFinalSum - nonProbFinalSum
  
  #sum last column 
  #SeqProbCRM<-sum(CRMProb[,ncol(CRMProb)])
  #SeqProbNon<-sum(NonProb[,ncol(NonProb)])
  
  #predictEntry<-c(observationSequence,rawSequence,SeqProbCRM,SeqProbNon)
  #readd observed seq later
  predictEntry<-paste(observationSequence, rawSequence, crmProbFinalSum, finalCRMColumn, nonProbFinalSum, finalNonColumn, enrichedStrength1, enrichedStrength2, prediction, sep="\t")
  
  #read in sequence
  
  
}
