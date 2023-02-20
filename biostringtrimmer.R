retrievedSeqList<-readDNAStringSet("observedSegmentResults", format="fasta")

trimmedseqlist<-c()

i=0

while (i < 20)
{
string1<-retrievedSeqList[i]
trimmedseq<-subseq(string1, 1, 100)

print(trimmedseq)

#trimmedseqlist<-c(trimmedseqlist,trimmedseq)

trimmedseq<-as.character(trimmedseq)

header<-paste(">>",i)

write(header,                                           
      file = "trimmedseq",
      append = TRUE)

write(trimmedseq,                                           
      file = "trimmedseq",
      append = TRUE)

i=i+1
}


