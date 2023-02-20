DataFrame2FastaConverter->function(dataframecol1,dataframecol2,convertedFileName)
{
  conversionStream <- file(convertedFileName, 'w')
  
  for(i in 1:nrows(dataframecol1))
  {
        writeLines(dataframecol1[i],convertedStream)
    
    writeLines(dataframecol2[i],convertedStream)
  }
  
}