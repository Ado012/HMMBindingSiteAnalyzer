# **HMM Binding Site Analyzer for Arabidopsis CRM sequences**



HMM based classifier to identify protein binding regions based on ChIP data. The model is tuned specifically for WUS CRM but can theoretically be used to identify other binding sequences for other proteins.


Example commands can be found in hmmCommandListMSB.R and are described below. Sample input files can be found in Samples directory. 



## **TARGET PREP**

Scripts which prepare sequences to be tested by the model.


------------------
File: targetAnnotatorHMM.R

targetAnnotator: prepares list of test sequences from a list of gene IDs. Test sequences are what will be tested by the model. 

Input: Arabidopsis chromosome files ex: Rdatafiles/chrom/Arabidopsis_thaliana.TAIR10.dna.chromosome.1.fa
gene list for test sequences: "downregulatedwusgenes_cyclo.csv"

output: targetlist

usage: targetlist<-targetAnnotatorHMM()


------------------
File: hmmModelerHelperFunctions.R

DataFrame2FastaConverterFileMSB: convert combined files to seperated files ready to be tested. 


input: combined target file
output: segmented target files ready for testing. 

usage: DataFrame2FastaConverterFileMSB(targetlist,'targetSegmentResultsCombinedUp', 'targetSegmentResultsCenter','targetSegmentResultsLeft','targetSegmentResultsRight' )


------------------




## **TRAINING PREP**

Scripts to prepare the training data for the model

------------------
File: coverageFindertest.R

coverageFinder: classifies sequences by coverage; get coordinates of segments from the genome which are wushel binding regions and those which are not. 
Feed in chromosome files for length
Feed in bed files for coverage 

Input: chromosome fasta files ex: Rdatafiles/chrom/Arabidopsis_thaliana.TAIR10.dna.chromosome.1.fa
Bed file: Data/sig1.bed

Output: enriched and clear coordinate files: enrichedCoordinatesChr1, clearCoordinatesChr1

Usage: coverageFinder()


------------------
File: GenomeShredder.R

genome shredder: 
feed in results from coverage finder to get sequences. 

Input: coordinate files: enrichedCoordinatesChr1
Output: sequence data: dataframe with fastatitles, sequences and transitions enrichedSegmentsResults. sequences with coordinates: enrichedCoordinatesResults

GenomeShredderForSortedRegions(regionCoordFile="enrichedCoordinatesChr1","enrichedSegmentsResults",'enrichedCoordinatesResults')
GenomeShredderForSortedRegions(regionCoordFile="clearCoordinatesChr1","clearSegmentsResults",'clearCoordinatesResults')

------------------





## **Model Construction**

Constructs HMM out of sequences. 

------------------
File: hmmModeler.R, hmmModelerHelperFunctions.R

hmmmodeler: Creates hmm model for a type of sequence through a library of fragments, specifically tuned for the CRM sequences of WUS in Arabidopsis. *Can be swapped out with alternative models for other systems.*

*Input: enriched and clear sequence fragments.
Intermediate: fasta file of sequences: enrichedSegmentResultsTargets.
fasta file of transitions: enrichedSegmentResultsTransitions
Output: HMM model


usage:
hmmOne<-HMMModeler("enrichedSegmentsResults","enrichedSegmentsResultsTargets","enrichedSegmentsResultsTransitions")
hmmTwo<-HMMModeler("clearSegmentsResults","clearSegmentsResultsTargets","clearSegmentsResultsTransitions")


-----------------


## **Model Testing and Validation
Generating predictions and calculating the accuracy of the model

-----------------
File: hmmTesterForMSB.R, hmmTesterHelperFunctions.R

HMMTest: Takes two HMM models and a set of test sequences and predicts what model those sequences are likely to belong to.
Input: two HMM models. The name of the test sequences in fasta format. 
Output: the prediction file, results data


usage: HMMTest(hmmOne,hmmTwo,"observedSegmentResults", 'rawPredictionResults')


-----------------
File: hmmValidator.R

hmmvalidator: takes as input results of hmmtest and compares them to GenomeShredder output to test accuracy of predictions if using segments from training data.
Input: prediction file, training data.
Output: refined results file

usage: HMMValidator("rawPredictionResults", "enrichedSegmentsFile", "clearSegmentsFile", "refinedResultsFile")


----------------

HMMResultsGraph: graphs results
Input: results data
Output: graphs



## **MISC
----------------

Extracts samples from sequence collections

usage: SampleExtractor('enrichedSegmentsResults','clearSegmentsResults','enrichedSegmentsResultsAbridged','clearSegmentsResultsAbridged','observedSegmentRaw')
---------------




