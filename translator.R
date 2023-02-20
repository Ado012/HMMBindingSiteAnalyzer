#takes DNA sequences and translates them into the core spacer sequences
library(stringr)
library(Biostrings)
library(dplyr)

#divides emissions up into chunks 
StateChopper<-function(chainStart, chainEnd, annotation, stateMatrix)
{#this function gets th beginning and end of complex core, miscore group, and spacer emission chains and relabels them properly 
  j=chainStart
  
  #while position is less than chain's end
  while (j<=chainEnd)
  {
  print(j)
  print(chainStart)
  print(chainEnd)
  #if chain end is less than three away from current position set chunk to end of sequence
  if ((j+3)>chainEnd)
    complexCoreFragEnd=chainEnd-j
  else
    complexCoreFragEnd=3
  
  #if chain is at its start start list to hold chunks
  if (j==chainStart)
  {
    choppedparts<-data.frame(annotation,j, j+complexCoreFragEnd)
  }
  
  else
  { #or add a new chunk to the list
    newentry <-data.frame(annotation,j, j+complexCoreFragEnd)
    choppedparts<-rbind(choppedparts,newentry)
  }
  j=j+complexCoreFragEnd+1
  }
  
  
  print('blah')#label the chunks
  colnames(choppedparts)<-c('state','start','end')
  
  #bind the chunks to the emission list
 stateMatrix <- rbind(stateMatrix,choppedparts)
}



#detects chains of cores which are intersecting or overlapping. 
ComplexCoreDetector <-function(stateFrame)
{
  
  print("begin complex core detection")
  stateFrameSize <- nrow(stateFrame)
  chain = 0


  
  chainEntry<-c()
  complexcoreseq <- data.frame( "state"= character(0), "start" = integer(0), "end" = integer(0))
  
  #not sure why it is originally (stateFrameSize>=1) > 1 seems to be better since one entry auto implies no complex cores
  if (stateFrameSize>1)
  {
    #detect complex cores
    for (i in 1:(stateFrameSize - 1))
    {
      print("enter complex core detection loop for sequence")
      
      print(stateFrame)
      
      
      #if cores starts are 4 or less spaces apart
      #subtle error where grep()==1 caused the first part to be shown as an error
      #no idea what the second part of this if statement does, look into it later
      if (( (stateFrame[i + 1, ][, 'start']- stateFrame[i, ][, 'start']) <= 4) && length(grep(stateFrame[i + 1, ][, 'state'],stateFrame[i, ][, 'state']))==1 )
      {#start chain if not already started
        print("core chain potential")
        
        if (chain == 0)
        {  
          chainStart = stateFrame[i, ][, 'start']
        chainEntry<-append(chainEntry,i)
        chain=1
        }
        
        #add to chain
        chainEntry<-append(chainEntry,i+1)
        print(chainEntry)
      }
      
      print("iterating through cores")
      
      #if chain is started and next core is more than 4 spaces apart
      #no idea what the second part of this if statement does, look into it later
      if (chain == 1 && (((stateFrame[i + 1, ][, 'start'] - stateFrame[i, ][, 'start']) > 4) || (i==stateFrameSize-1) || length(grep(stateFrame[i + 1, ][, 'state'],stateFrame[i, ][, 'state'])) ==0 ))
      {#terminate chain
        #if the chain ends because of the end of the segment. Absence was source of subtle error
        if (i==stateFrameSize-1)
          chainEnd <- stateFrame[i + 1, ][, 'end']
        #else the chain is ending simply because there is another type of emission after the current emission
        else
        chainEnd <- stateFrame[i, ][, 'end']
        
        print("chain ended")
        
        #add new complexcore
       
        
        entryState<-stateFrame[i,][,'state']
        
        if (entryState=='core')
          entryStateGroup<-c('complexcore')
        else if (entryState=='miscore')
          entryStateGroup<-c('miscoregroup')
        
        #chop up complex core into 4 nucleotide chunks and label them with the complexcore label
       stateFrame <-StateChopper(chainStart, chainEnd, entryStateGroup, stateFrame)
        
        
        
        if (nrow(stateFrame)==1)#column names are deleted if rbind to empty vector
          colnames(stateFrame)<-c('state','start','end')
        
        chain = 0
        
      }
      
      #if chain is started and the penultimate core is reached (implicitly means that the last core is part of the chain due to above if)
      if (chain>=1 && i >= (stateFrameSize-1))
      {
        print('chain started but end is reached')
        #last core is chain end
        chainEnd <- stateFrame[i+1, ][, 'end']
        
        entryState<-stateFrame[i,][,'state']
        
        if (entryState=='core')
          entryStateGroup<-c('complexcore')
        else if (entryState=='miscore')
          entryStateGroup<-c('miscoregroup')
        
        #add new complexcoremminished')
        
        chain = 0
      }
      
      
      
    }
  }
  
  print(chainEntry)
  
  #delete all cores that are part of chain
  if (is.null(chainEntry)=='FALSE')
  stateFrame <- stateFrame[-c(chainEntry), ]
  
  else 
    stateFrame=stateFrame
  
  stateFrame
  
  
}

#detects gaps among cores in list
SpaceDetector<-function(stateList,retrievedDNA)
{  
  
  stateListSize<-nrow(stateList)
  print('space detector started')
#arrange cores from beginning to end. Doing it on a 1 row df will screw it up for some reason
if (nrow(stateList)>1)
  stateList <- arrange(stateList,start)



spacerBag <- data.frame( "start" = integer(0), "end" = integer(0), "segmentEndStatus" = integer(0))
#detect intercalery spaces

if (stateListSize >1)
{
  
  for (i in 1:(stateListSize-1))
  {
    print('iterating and finding spaces')
    #if the gap between detected states is detected than it is a space
    if (stateList[i + 1, ][, 'start']-stateList[i, ][, 'end']> 1)
    {
      print("iterating and finding spaces 2")
      
      spaceStart <- stateList[i, ][, 'end'] + 1
      print("iterating and finding spaces 3")
      
      spaceEnd <- stateList[i + 1, ][, 'start'] - 1
      segmentEndStatus=0
      spacerEntry<-data.frame(spaceStart,spaceEnd,segmentEndStatus)
      spacerBag<-rbind(spacerBag,spacerEntry)
    }
    
    
  }
  
}

#detect a sequence of nothing but space
if (stateListSize==0)
{
  spaceStart=1
  spaceEnd=length(retrievedDNA)
  spacesize=spaceEnd

segmentEndStatus=1
  spacerEntry<-data.frame(spaceStart,spaceEnd,segmentEndStatus)
  spacerBag<-rbind(spacerBag,spacerEntry)
}

#detect bordering spaces
if (stateListSize>=1)
{
  
  print("size of state list is = to or greater than 1")
  #detect left border state
  if (stateList[1,][,'start']>1)
  {
    spaceStart=1
    spaceEnd <- stateList[1, ][, 'start'] - 1
    spacesize = spaceStart-spaceEnd+1
    segmentEndStatus=1
   spacerEntry<-data.frame(spaceStart,spaceEnd,segmentEndStatus)
   spacerBag<-rbind(spacerBag,spacerEntry)
  }
  
  #detect right border space
  if (stateList[stateListSize,][,'end'] < length(retrievedDNA))
  {
    print("the end of the state list is smaller than the width of the dna")
    spaceStart=stateList[stateListSize,][,'end']+1
    spaceEnd <- length(retrievedDNA)
    segmentEndStatus=1
    spacerEntry<-data.frame(spaceStart,spaceEnd,segmentEndStatus)
    spacerBag<-rbind(spacerBag,spacerEntry)
  }
  print("finished finding spaces in a statelist of one or greater")
  
}

#go through the gaps detected and classify them
for(i in 1:nrow(spacerBag))
{
  
  spaceStart=spacerBag[i,][,'spaceStart']
  spaceEnd=spacerBag[i,][,'spaceEnd']
  spacesize=spaceEnd-spaceStart+1
  spacechain <- data.frame( "state"= character(0), "start" = integer(0), "end" = integer(0))
  
#classify spaces based on size if a multiple  of 10 equaL to or below 40 (add some leeway in the future)
if (spacesize %% 10 == 0 && spacesize <= 40 && spacerBag[i,][,'segmentEndStatus']==0)
  spaceAnnotation <- ('spacer')
#otherwise its noise
else
  spaceAnnotation <- ('noise')
#chop up gaps into chunks
  stateList<-StateChopper(spaceStart, spaceEnd, spaceAnnotation, stateList)


}
print("finished space finder")
#return updated emission list
stateList
#size classifier








}

#umbrella function for translating nucleotide sequence into emission sequence
Translator <- function(retrievedDNA)
{
  #readindnafragments
  
  
  #detect TAAT cores
  cores <- str_locate_all(retrievedDNA, "(?=(TAAT|ATTA))")
  
  #record core locations
  cores <- data.frame(cores)
  
  #set end of core location properly
  cores[,'end']<-cores[,'end']+4
  
  
  #Find miscores
  #apparently str_locate_all might not work properly if not on one line.
  miscores <-
    str_locate_all(
      retrievedDNA,
      "(?=(AAAT|CAAT|GAAT|TTTA|CTTA|GTTA|TTAT|TCAT|TGAT|AATA|ACTA|AGTA|TATT|TACT|TAGT|ATAA|ATCA|ATGA|TAAA|TAAC|TAAG|ATTT|ATTC|ATTG))"
)
  
  
  print("finished miscore search")
  print(cores)
  
  miscores <- data.frame(miscores)
  miscores[,'end']<-miscores[,'end']+4
  print("entered miscore positions")
  print(miscores)
  
  #add annotations
  annotationcolumn <- rep('core', nrow(cores))
  annotationcolumn <- data.frame(annotationcolumn)
  names(annotationcolumn)<-c('state')
  
  #create dataframe statelist with states and coordinates
  cores <- cbind(annotationcolumn, cores)
  
  #create a dataframe listing out the states ie core core core 
  annotationcolumn <- rep('miscore', nrow(miscores))
  annotationcolumn <- data.frame(annotationcolumn)
  names(annotationcolumn)<-c('state')
  miscores <- cbind(annotationcolumn, miscores)
  
  #combine core and miscore list 
  stateList <- rbind(cores, miscores)
  #arrange emissions chronologically
  stateList<-arrange(stateList,start)
  #iterate through list deleting miscores which overlap with cores
  #was originally just (length(stateList)>=1 ) for some reason
  if (length(stateList)>=1 && nrow(stateList) > 1)
  {
  entriesToDelete<-c()
  for (i in 1:(nrow(stateList) - 1))
  {#will repeat entries mess up deletion? 

    if (stateList[i,]['state']=='core')
    {
      j=1
      while (j<=3)
      {#maybe first if should be seperated?
        #make sure you don't run off the edge of the state list
       if (i+j<=nrow(stateList))#checking three bases forward and backward to look for overlapping miscores
       {
        if (stateList[i+j,]['state']=='miscore' && (stateList[i,]['end'] >= stateList[i+j,]['start']) )
          entriesToDelete<-append(entriesToDelete,i+j)
       }
        #make sure you don't go beyond the first entry of state list
        if (i-j>=1)
        {
          if (stateList[i-j,]['state']=='miscore' && (stateList[i,]['start'] <= stateList[i-j,]['end']) )
            entriesToDelete<-append(entriesToDelete,i-j)
        }
          
          j=j+1
      }
      
    }
    
    
    
  }
  

  if (length(entriesToDelete)>0)#remove entries marked for deletion
  {
    entriesToDelete<-unique(entriesToDelete)
  stateList <- stateList[-c(entriesToDelete), ]
  }
  #detect core chains
  stateList<-ComplexCoreDetector(stateList)
  print(cores)

  print("finished complex core detection")
  #miscore chain finder

  #miscore chain finder
  

  print(stateList)
  
stateList<-SpaceDetector(stateList,retrievedDNA)
  print('finished space detection')
print(stateList)

  }
  
  else #chop up noise if statelist is empty
  {
    stateAnnotation <- ('noise')
    
    stateList<-StateChopper(1, length(retrievedDNA), stateAnnotation, stateList)
  }
  
  
  
  #arrange cores from beginning to end. Doing it on a 1 row df will screw it up for some reason
  if (nrow(stateList)>1)
    stateList <- arrange(stateList,start)

  
  print(stateList)
}

