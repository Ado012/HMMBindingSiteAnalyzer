




targetlist<-targetAnnotatorHMM()


DataFrame2FastaConverterFileMSB(targetlist,'targetSegmentResultsCombined', 'targetSegmentResultsCenter','targetSegmentResultsLeft','targetSegmentResultsRight' )


coverageFinder()


GenomeShredderForSortedRegions(regionCoordFile="enrichedCoordinatesChr1","enrichedSegmentsResults",'enrichedCoordinatesResults')

GenomeShredderForSortedRegions(regionCoordFile="clearCoordinatesChr1","clearSegmentsResults",'clearCoordinatesResults')



hmmThree<-HMMModeler("testSegmentsResultsAbridged","testSegmentsResultsAbridgedTargets","testSegmentsResultsAbridgedTransitions")

hmmOne<-HMMModeler("enrichedSegmentsResults","enrichedSegmentsResultsTargets","enrichedSegmentsResultsTransitions")

hmmOne<-HMMModeler("enrichedSegmentsResultsAbridged","enrichedSegmentsResultsAbridgedTargets","enrichedSegmentsResultsAbridgedTransitions")

hmmTwo<-HMMModeler("clearSegmentsResultsAbridged","clearSegmentsResultsAbridgedTargets","clearSegmentsResultsAbridgedTransitions")


#observedSegmentsResults from enrichedSegmentsResults

resultsCombined<-HMMTest(hmmOne,hmmTwo,"targetSegmentResultsCombined", 'rawPredictionResultsCombined')

resultsCenter <- HMMTest(hmmOne,hmmTwo,"targetSegmentResultsCenter", 'rawPredictionResultsCenter')

resultsLeft <- HMMTest(hmmOne,hmmTwo,"targetSegmentResultsLeft", 'rawPredictionResultsLeft')

resultsRight <- HMMTest(hmmOne,hmmTwo,"targetSegmentResultsRight", 'rawPredictionResultsRight')

HMMValidator("rawPredictionResultsCenter", "enrichedCoordinateResults", "clearCoordinateResults", "refinedResultsFile")



resultsCombinedUp<-HMMTest(hmmOne,hmmTwo,"targetSegmentResultsCombinedUp", 'rawPredictionResultsCombinedUp')

resultsCenterUp <- HMMTest(hmmOne,hmmTwo,"targetSegmentResultsCenterUp", 'rawPredictionResultsCenterUp')

resultsLeftUp <- HMMTest(hmmOne,hmmTwo,"targetSegmentResultsLeftUp", 'rawPredictionResultsLeftUp')

resultsRightUp <- HMMTest(hmmOne,hmmTwo,"targetSegmentResultsRightUp", 'rawPredictionResultsRightUp')


