setwd("/home/ace/Documents/Software_Projects/HMModeler")


library(ggplot2)
# Basic histogram
# Basic histogram

#some rows will be dropped due to no difference in the prediction in enrichedstrength2


resultsCenterEnriched <-subset(resultsCenter, predict=='enriched')


logenriched = log10(as.numeric(resultsCenterEnriched$enrichedStrength2))

logenriched <-data.frame(logenriched)

ggplot(logenriched, aes(x=logenriched)) + geom_histogram(binwidth = 10) + xlab("enriched strength") + 
  ggtitle("Enriched Strength Downregulated Wus Genes: Gene Body") + theme(plot.title = element_text(hjust = 0.5)) +
scale_x_continuous(limits = c(-400, 0)) + scale_y_continuous(limits = c(0, 100))


resultsCenterClear <-subset(resultsCenter,predict=='clear')

logclear = log10(-1* as.numeric(resultsCenterClear$enrichedStrength2))

logclear <-data.frame(logclear)

ggplot(logclear, aes(x=logclear)) + geom_histogram(binwidth = 10) + xlab("clear strength") + 
  ggtitle("Clear Strength Downregulated Wus Genes: Gene Body") + theme(plot.title = element_text(hjust = 0.5)) + 
scale_x_continuous(limits = c(-400, 0)) + scale_y_continuous(limits = c(0, 100))




ggplot(resultsCenter, aes(x=predict)) + stat_count() + xlab("binding status") + 
  ggtitle("Detected Binding vs Nonbinding Regions Downregulated: Gene Body") + theme(plot.title = element_text(hjust = 0.5)) +
  scale_y_continuous(limits = c(0, 400))



#Left


resultsLeftEnriched <-subset(resultsLeft, predict=='enriched')


logenriched = log10(as.numeric(resultsLeftEnriched$enrichedStrength2))

logenriched <-data.frame(logenriched)

ggplot(logenriched, aes(x=logenriched)) + geom_histogram(binwidth = 10) + xlab("enriched strength") + 
  ggtitle("Enriched Strength Downregulated Wus Genes: 5' End") + theme(plot.title = element_text(hjust = 0.5)) +
  scale_x_continuous(limits = c(-400, 0)) + scale_y_continuous(limits = c(0, 100))


resultsLeftClear <-subset(resultsLeft,predict=='clear')

logclear = log10(-1* as.numeric(resultsLeftClear$enrichedStrength2))

logclear <-data.frame(logclear)

ggplot(logclear, aes(x=logclear)) + geom_histogram(binwidth = 10) + xlab("clear strength") + 
  ggtitle("Clear Strength Downregulated Wus Genes: 5' End") + theme(plot.title = element_text(hjust = 0.5)) + 
  scale_x_continuous(limits = c(-400, 0)) + scale_y_continuous(limits = c(0, 100))




ggplot(resultsLeft, aes(x=predict)) + stat_count() + xlab("binding status") + 
  ggtitle("Detected Binding vs Nonbinding Regions Downregulated: 5' End") + theme(plot.title = element_text(hjust = 0.5)) +
  scale_y_continuous(limits = c(0, 400))


#Right

resultsRightEnriched <-subset(resultsRight, predict=='enriched')


logenriched = log10(as.numeric(resultsRightEnriched$enrichedStrength2))

logenriched <-data.frame(logenriched)

ggplot(logenriched, aes(x=logenriched)) + geom_histogram(binwidth = 10) + xlab("enriched strength") + 
  ggtitle("Enriched Strength Downregulated Wus Genes: 3' End") + theme(plot.title = element_text(hjust = 0.5)) +
  scale_x_continuous(limits = c(-400, 0)) + scale_y_continuous(limits = c(0, 100))


resultsRightClear <-subset(resultsRight,predict=='clear')

logclear = log10(-1* as.numeric(resultsRightClear$enrichedStrength2))

logclear <-data.frame(logclear)

ggplot(logclear, aes(x=logclear)) + geom_histogram(binwidth = 10) + xlab("clear strength") + 
  ggtitle("Clear Strength Downregulated Wus Genes: 3' End") + theme(plot.title = element_text(hjust = 0.5)) + 
  scale_x_continuous(limits = c(-400, 0)) + scale_y_continuous(limits = c(0, 100))




ggplot(resultsRight, aes(x=predict)) + stat_count() + xlab("binding status") + 
  ggtitle("Detected Binding vs Nonbinding Regions Downregulated: 3' End") + theme(plot.title = element_text(hjust = 0.5)) +
  scale_y_continuous(limits = c(0, 400))


#extract top genes down center

data <- resultsCenter1 %>%   
  arrange(desc(enrichedStrength2)) %>% 

  
  

data_new2 <- data %>%                                      # Top N highest values by group
  arrange(desc(value)) %>% 
  group_by(group) %>%
  slice(1:3)




#top enriched center
resultsCenter1 <-resultsCenter
resultsCenter1$enrichedStrength2 <-as.numeric(resultsCenter1$enrichedStrength2)

data <- resultsCenter1 %>% arrange(desc(enrichedStrength2)) 

  data <- data[1:50,]
  
write.table(data, "topgenesdownenrichedcenter", append = FALSE, sep = " ", dec = ".",
            row.names = TRUE, col.names = TRUE)


#top enriched left
resultsLeft1 <-resultsLeft
resultsLeft1$enrichedStrength2 <-as.numeric(resultsLeft1$enrichedStrength2)

data <- resultsLeft1 %>% arrange(desc(enrichedStrength2)) 

data <- data[1:50,]

write.table(data, "topgenesdownenrichedleft", append = FALSE, sep = " ", dec = ".",
            row.names = TRUE, col.names = TRUE)


#top enriched right
resultsRight1 <-resultsRight
resultsRight1$enrichedStrength2 <-as.numeric(resultsRight1$enrichedStrength2)

data <- resultsRight1 %>% arrange(desc(enrichedStrength2)) 

data <- data[1:50,]

write.table(data, "topgenesdownenrichedright", append = FALSE, sep = " ", dec = ".",
            row.names = TRUE, col.names = TRUE)






#upregulated genes





resultsCenterEnriched <-subset(resultsCenterUp, predict=='enriched')


logenriched = log10(as.numeric(resultsCenterEnriched$enrichedStrength2))

logenriched <-data.frame(logenriched)

ggplot(logenriched, aes(x=logenriched)) + geom_histogram(binwidth = 10) + xlab("enriched strength") + 
  ggtitle("Enriched Strength Upregulated Wus Genes: Gene Body") + theme(plot.title = element_text(hjust = 0.5)) +
  scale_x_continuous(limits = c(-400, 0)) + scale_y_continuous(limits = c(0, 100))


resultsCenterClear <-subset(resultsCenterUp,predict=='clear')

logclear = log10(-1* as.numeric(resultsCenterClear$enrichedStrength2))

logclear <-data.frame(logclear)

ggplot(logclear, aes(x=logclear)) + geom_histogram(binwidth = 10) + xlab("clear strength") + 
  ggtitle("Clear Strength Upregulated Wus Genes: Gene Body") + theme(plot.title = element_text(hjust = 0.5)) + 
  scale_x_continuous(limits = c(-400, 0)) + scale_y_continuous(limits = c(0, 100))




ggplot(resultsCenterUp, aes(x=predict)) + stat_count() + xlab("binding status") + 
  ggtitle("Detected Binding vs Nonbinding Regions Upregulated: Gene Body") + theme(plot.title = element_text(hjust = 0.5)) +
  scale_y_continuous(limits = c(0, 400))



#Left


resultsLeftEnriched <-subset(resultsLeftUp, predict=='enriched')


logenriched = log10(as.numeric(resultsLeftEnriched$enrichedStrength2))

logenriched <-data.frame(logenriched)

ggplot(logenriched, aes(x=logenriched)) + geom_histogram(binwidth = 10) + xlab("enriched strength") + 
  ggtitle("Enriched Strength Upregulated Wus Genes: 5' End") + theme(plot.title = element_text(hjust = 0.5)) +
  scale_x_continuous(limits = c(-400, 0)) + scale_y_continuous(limits = c(0, 100))


resultsLeftClear <-subset(resultsLeftUp,predict=='clear')

logclear = log10(-1* as.numeric(resultsLeftClear$enrichedStrength2))

logclear <-data.frame(logclear)

ggplot(logclear, aes(x=logclear)) + geom_histogram(binwidth = 10) + xlab("clear strength") + 
  ggtitle("Clear Strength Upregulated Wus Genes: 5' End") + theme(plot.title = element_text(hjust = 0.5)) + 
  scale_x_continuous(limits = c(-400, 0)) + scale_y_continuous(limits = c(0, 100))




ggplot(resultsLeftUp, aes(x=predict)) + stat_count() + xlab("binding status") + 
  ggtitle("Detected Binding vs Nonbinding Regions Upregulated: 5' End") + theme(plot.title = element_text(hjust = 0.5)) +
  scale_y_continuous(limits = c(0, 400))


#Right

resultsRightEnriched <-subset(resultsRightUp, predict=='enriched')


logenriched = log10(as.numeric(resultsRightEnriched$enrichedStrength2))

logenriched <-data.frame(logenriched)

ggplot(logenriched, aes(x=logenriched)) + geom_histogram(binwidth = 10) + xlab("enriched strength") + 
  ggtitle("Enriched Strength Upregulated Wus Genes: 3' End") + theme(plot.title = element_text(hjust = 0.5)) +
  scale_x_continuous(limits = c(-400, 0)) + scale_y_continuous(limits = c(0, 100))


resultsRightClear <-subset(resultsRightUp,predict=='clear')

logclear = log10(-1* as.numeric(resultsRightClear$enrichedStrength2))

logclear <-data.frame(logclear)

ggplot(logclear, aes(x=logclear)) + geom_histogram(binwidth = 10) + xlab("clear strength") + 
  ggtitle("Clear Strength Upregulated Wus Genes: 3' End") + theme(plot.title = element_text(hjust = 0.5)) + 
  scale_x_continuous(limits = c(-400, 0)) + scale_y_continuous(limits = c(0, 100))




ggplot(resultsRightUp, aes(x=predict)) + stat_count() + xlab("binding status") + 
  ggtitle("Detected Binding vs Nonbinding Regions Upregulated: 3' End") + theme(plot.title = element_text(hjust = 0.5)) +
  scale_y_continuous(limits = c(0, 400))


#extract top genes up center




#top enriched center
resultsCenter1Up <-resultsCenterUp
resultsCenter1Up$enrichedStrength2 <-as.numeric(resultsCenter1Up$enrichedStrength2)

data <- resultsCenter1Up %>% arrange(desc(enrichedStrength2)) 

data <- data[1:50,]

write.table(data, "topgenesupenrichedcenter", append = FALSE, sep = " ", dec = ".",
            row.names = TRUE, col.names = TRUE)


#top enriched left
resultsLeft1Up <-resultsLeftUp
resultsLeft1Up$enrichedStrength2 <-as.numeric(resultsLeft1Up$enrichedStrength2)

data <- resultsLeft1Up %>% arrange(desc(enrichedStrength2)) 

data <- data[1:50,]

write.table(data, "topgenesupenrichedleft", append = FALSE, sep = " ", dec = ".",
            row.names = TRUE, col.names = TRUE)


#top enriched right
resultsRight1Up <-resultsRightUp
resultsRight1Up$enrichedStrength2 <-as.numeric(resultsRight1Up$enrichedStrength2)

data <- resultsRight1Up %>% arrange(desc(enrichedStrength2)) 

data <- data[1:50,]

write.table(data, "topgenesupenrichedright", append = FALSE, sep = " ", dec = ".",
            row.names = TRUE, col.names = TRUE)




