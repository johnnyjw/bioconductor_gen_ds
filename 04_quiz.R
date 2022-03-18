library(tidyverse)

# q1

# The yeastRNASeq experiment data package contains FASTQ files from an RNA seq experiment in yeast. 
# What fraction of reads in this file has an A nucleotide in the 5th base of the read?

library(yeastRNASeq)
library(ShortRead)
fastqFilePath <- system.file("reads", "wt_1_f.fastq.gz", package = "yeastRNASeq")

fqFile <- FastqFile(fastqFilePath)
reads <- readFastq(fqFile)

all_seq <- sread(reads)

ranges <- IRanges(start = 5,
                  width = 1)
all_seq[,ranges]
all_seq
unlist(extractAt(all_seq, ranges))
letter_five <- as.character(unlist(extractAt(all_seq, ranges)))
letter_five_a <- letter_five[letter_five == 'A']
length(letter_five_a)/length(letter_five)
# 0.363841
# Correct

# q2
# This is a continuation of Question 1
# Question: What is the average numeric quality value of the 5th base of these reads?
all_qual <- as(quality(reads), "matrix")[,5]
mean(all_qual)
# 28.93
# Correct

# q3
# The leeBamViews experiment data package contains aligned BAM files from an RNA seq experiment in yeast 
# (the same experiment as in Questions 1 and 2, but that is not pertinent to the question)

# These reads are short reads (36bp) and have been aligned to the genome using a standard aligner, 
# ie. potential junctions have been ignored (this makes some sense as yeast has very few junctions and the reads are very short).

# A read duplicated by position is a read where at least one more read shares the same position.

# We will focus on the interval from 800,000 to 801,000 on yeast chromosome 13.

# Question: In this interval, how many reads are duplicated by position?

library(leeBamViews)
bamFilePath <- system.file("bam", "isowt5_13e.bam", package="leeBamViews")
bamFile <- bamFilePath
bamFile

aln <- scanBam(bamFile)
aln[[1]]
input_region <- GRanges(seqnames='Scchr13', ranges=IRanges(start = 800000, end=801000))
aln[[1]]$qname
params <- ScanBamParam(which=input_region, what = scanBamWhat())
aln2 <- scanBam(bamFile, param=params)

multi <- tibble(pos = aln2[[1]]$pos) %>% 
  group_by(pos) %>% 
  mutate(count = n()) %>%
  ungroup() %>% 
  filter(count>1)
# 129
# Correct

# q4
# The package contains 8 BAM files in total, representing 8 different samples from 4 groups. 
# An objective of the original paper was the discovery of novel transcribed regions in yeast.
# One such region is Scchr13:807762-808068.
# What is the average number of reads across the 8 samples falling in this interval?
bpaths <- list.files(system.file("bam", package="leeBamViews"), pattern = "bam$", full=TRUE)
bamView <- BamViews(bpaths)
bamView
input_region2 <- GRanges(seqnames='Scchr13', ranges=IRanges(start = 807762, end=808068))
bamRanges(bamView) <- input_region2
#params2 <- ScanBamParam(which=input_region2, what = scanBamWhat())
aln2 <- scanBam(bamView)
aln2$seq
seq(aln2)
aln2[,]

get_lengths <- function(x){
  subname <- names(x)[1]
  return(length(x[[subname]]$seq))
}

lapply(aln2, get_lengths)
lengths <- unlist(lapply(aln2, get_lengths))
mean(lengths)
# 90.25
# Correct

# q5
# In the lecture on the oligo package an ExpressionSet with 18 samples is constructed, 
# representing normalized data from an Affymetrix gene expression microarray. 
# The samples are divided into two groups given by the \verb|group|group variable.

# What is the average expression across samples in the control group for the “8149273” probeset 
# (this is a character identifier, not a row number).
library(oligo)
library(GEOquery)

# control 1-8, treated 1-10
celfiles <- list.files("GSE38792/CEL", full=TRUE)
celfiles

rawData <- read.celfiles(celfiles)
rawData
assayData(rawData)

sampleNames(rawData)

# clean up
filename <- sampleNames(rawData)
filename
#annotate this into pData
pData(rawData)$filename <- filename

# clean up using regex
sampleNames <- sub(".*_", "", filename)
sampleNames <- sub(".CEL.gz$", "", sampleNames)
sampleNames
sampleNames(rawData) <- sampleNames
pData(rawData)$group <- ifelse(grepl("^OSA", sampleNames(rawData)),
                               "OSA", "Control")
normData <- rma(rawData)

boxplot(normData)
mean(exprs(normData)["8149273", 1:8])
# 7.02183
# Correct

# q6
# This is a continuation of Question 5.

# Use the limma package to fit a two group comparison between the control group and the OSA group, 
# and borrow strength across the genes using eBayes(). Include all 18 samples in the model fit.

# What is the absolute value of the log foldchange (logFC) of the gene with the lowest P.value.
library(limma)

# clean away levels that no longer exist in the data (would would mess up models)
normData$group <- factor(normData$group)
normData$group

# setting up design matrix
design <- model.matrix(~ normData$group)
head(design)

# fit basic model
fit <- lmFit(normData, design)
fit <- eBayes(fit)
topTable(fit)

# 0.7126
# Correct

# q7
# How many genes are differentially expressed between the two groups at an adj.P.value cutoff of 0.05?
goulet <- topTable(fit, number=10000)
goulet2 <- goulet %>% 
  filter(adj.P.Val <= 0.05)

# 0
# Correct

# q8
# An example 450k dataset is contained in the minfiData package. This dataset contains 6 samples; 3 cancer and 3 normals. 
# Cancer has been shown to be globally hypo-methylated (less methylated) compared to normal tissue of the same kind.

# Take the RGsetEx dataset in this package and preprocess it with the preprocessFunnorm function. 
# For each sample, compute the average Beta value (percent methylation) across so-called OpenSea loci.

# What is the mean difference in beta values between the 3 normal samples and the 3 cancer samples, across OpenSea CpGs?
library(minfiData)
data(RGsetEx)
RGsetEx
pData(RGsetEx)
head(sampleNames(RGsetEx))

# normalise data
grSet = preprocessFunnorm(RGsetEx)
# genomic ratio set
grSet

granges(grSet)
# cpg island status
island_statii <- getIslandStatus(grSet)

granges_set <- granges(grSet)
granges_set$islands <-  getIslandStatus(grSet)
granges_set_opensea <- subset(granges_set, islands=='OpenSea')
opensea_names <- names(granges_set_opensea)

grSet_open <- grSet[opensea_names,]

#analyse
#getBeta(grSet_open)[1:3, 1:3]
betie <- getBeta(grSet_open)
beta_df <- tibble(data.frame(betie))
beta_df2 <- beta_df %>%   
  mutate(mean_control = (`X5723646052_R02C02` + `X5723646052_R04C01` + `X5723646053_R05C02`)/3) %>% 
  mutate(mean_cancer = (`X5723646052_R05C02` + `X5723646053_R04C02` + `X5723646053_R06C02`)/3) %>% 
  mutate(mean_difference = mean_control - mean_cancer)

mean(beta_df2$mean_difference)

#  0.08863657
# correct
