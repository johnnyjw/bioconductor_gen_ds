#q1 
# What is the mean expression across all features for sample 5 in the ALL dataset (from the ALL package)?
library(tidyverse)
library(ALL)
data(ALL)
ALL

mean(exprs(ALL)[,5])
#5.629627

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

#q3
# How many probesets (Affymetrix IDs) are annotated with one or more genes on the autosomes (chromosomes 1 to 22).
affy_chr <- getBM(attributes = c("affy_hg_u95av2", "chromosome_name"),
                   filters = "affy_hg_u95av2", values =all.features, mart = ensembl)

autosomes <- c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", 
               "14", "15", "16", "17", "18", "19", "20", "21", "22")


auto_only <- affy_chr %>% 
  filter(chromosome_name %in% autosomes)

# 11448

#q4
# Use the MsetEx dataset from the minfiData package. 
# Part of this question is to use the help system to figure out how to address the question.

# What is the mean value of the Methylation channel across the features for sample “5723646052_R04C01”?
BiocManager::install("minfiData")
library(minfiData)
data(minfiData)

nphenoData(ALL)
featureData(ALL)
feature_df <- featureData(ALL)
feature_df$data
colnames(feature_df$data)
rownames(feature_df$data)
