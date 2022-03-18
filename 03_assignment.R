#q1 
# What is the mean expression across all features for sample 5 in the ALL dataset (from the ALL package)?
library(tidyverse)
library(ALL)
data(ALL)
ALL

mean(exprs(ALL)[,5])
#5.629627
# CORRECT

#q2
# We will use the biomaRt package to annotate an Affymetrix microarray. 
# We want our results in the hg19 build of the human genome 
# and we therefore need to connect to Ensembl 75 which is the latest release on this genome version. 
# How to connect to older versions of Ensembl is described in the biomaRt package vignette; 
# it can be achived with the command 
# mart <- useMart(host=’feb2014.archive.ensembl.org’, biomart = "ENSEMBL_MART_ENSEMBL").

# Question: Using this version of Ensembl, annotate each feature of the ALL dataset with the Ensembl gene id. 
# How many probesets (features) are annotated with more than one Ensembl gene id?
library(biomaRt)
mart <- useMart(host='feb2014.archive.ensembl.org', biomart = "ENSEMBL_MART_ENSEMBL")
mart
ensembl <- useDataset("hsapiens_gene_ensembl", mart)
attributes <- listAttributes(ensembl, page = "feature_page")

head(pData(ALL))

all.features <-  featureNames(ALL)
crosswalk <- getBM(attributes = c("ensembl_gene_id", "affy_hg_u95av2"),
      filters = "affy_hg_u95av2", values =all.features, mart = ensembl)
class(crosswalk)
summary <- crosswalk %>% 
  group_by(affy_hg_u95av2) %>% 
  summarize(n = n()) %>%
  filter(n>1)

# 1045
# CORRECT

#q3
# How many probesets (Affymetrix IDs) are annotated with one or more genes on the autosomes (chromosomes 1 to 22).
affy_chr <- getBM(attributes = c("affy_hg_u95av2", "chromosome_name"),
                   filters = "affy_hg_u95av2", values =all.features, mart = ensembl)

autosomes <- c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", 
               "14", "15", "16", "17", "18", "19", "20", "21", "22")


auto_only <- affy_chr %>% 
  filter(chromosome_name %in% autosomes)

summary <- auto_only %>% 
  group_by(affy_hg_u95av2) %>% 
  summarize(n = n()) %>%
  filter(n>=1)

# 11448
# WRONG
# 11016
# CORRECT

#q4
# Use the MsetEx dataset from the minfiData package. 
# Part of this question is to use the help system to figure out how to address the question.

# What is the mean value of the Methylation channel across the features for sample “5723646052_R04C01”?
# BiocManager::install("minfiData")
library(minfiData)
data(MsetEx)

mean(getMeth(MsetEx)[, "5723646052_R04C01"])
# 7228.277
# CORRECT

#q5
# Access the processed data from NCBI GEO Accession number GSE788. 
# What is the mean expression level of sample GSM9024?
# geo = geneexpression ominbus
# BiocManager::install("GEOquery")
library(GEOquery)

eList = getGEO("GSE788")
length(eList)
names(eList)
eData = eList[[1]]
eData

mean(exprs(eData)[,"GSM9024"])
# 756.432
# CORRECT

#q6
# We are using the airway dataset from the airway package.
# What is the average of the average length across the samples in the expriment?
# BiocManager::install("airway")
library(airway)
data("airway")

mean(airway$avgLength)
# 113.75
# CORRECT

#q7
# We are using the airway dataset from the airway package. 
# The features in this dataset are Ensembl genes.
# What is the number of Ensembl genes which have a count of 1 read or more in sample SRR1039512?
colData(airway)
assay(airway)
srr_assay <- assay(airway)[,"SRR1039512"]
srr_assay_plus <- srr_assay[srr_assay >= 1]
# 25699
# CORRECT

#q8
#  The airway dataset contains more than 64k features. 
# How many of these features overlaps with transcripts on the autosomes (chromosomes 1-22) 
# as represented by the TxDb.Hsapiens.UCSC.hg19.knownGene package?
pData(airway)
airwayFeatures <- rowRanges(airway)

library("TxDb.Hsapiens.UCSC.hg19.knownGene")
txdb = TxDb.Hsapiens.UCSC.hg19.knownGene
txdb

tx_cds <- exonsBy(txdb, 'tx')
#tx_cds <- cds(txdb)

autosomes <- paste0("chr", 1:22)
autosomes
tx_cds_auto <- keepSeqlevels(tx_cds, autosomes, pruning.mode = "coarse")
seqlevelsStyle(airwayFeatures) <- "UCSC"
airway_auto <- subsetByOverlaps(airwayFeatures, tx_cds_auto)
# 26276
# CORRECT

#q9
# For sample SRR1039508, how big a percentage (expressed as a number between 0 and 1) of the total reads 
# in the airway dataset for that sample, are part of a feature which overlaps an autosomal 
# TxDb.Hsapiens.UCSC.hg19.knownGene transcript?


srr_assay <- assay(airway)[,"SRR1039508"]
srr_assay
airway_auto_genes <- names(airway_auto)
sum(srr_assay)
srr_assay_subset <- srr_assay[names(srr_assay) %in% airway_auto_genes]

sum(srr_assay_subset) / sum(srr_assay)
#  0.9004193
# CORRECT

