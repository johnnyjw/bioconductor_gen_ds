

library(oligo)
library(GEOquery)

getGEOSuppFiles("GSE38792")
list.files("GSE38792")
untar("GSE38792/GSE38792_RAW.tar", exdir = "GSE38792/CEL")
list.files("GSE38792/CEL")
# control 1-8, treated 1-10
celfiles <- list.files("GSE38792/CEL", full=TRUE)
celfiles

rawData <- read.celfiles(celfiles)
rawData
#more than 1 million features here
# feature=probe to measure RNA transcript

getClass("GeneFeatureSet")

exprs(rawData)[1:4, 1:3]
max(exprs(rawData))

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
pData(rawData)

oligo::boxplot(rawData, target="core")

# normalise data
normData <- rma(rawData)
normData
boxplot(normData)
