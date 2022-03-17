library(minfi)
library(GEOquery)

#example - relationship between methylation and acute mania
getGEOSuppFiles("GSE68777")
untar("GSE68777/GSE68777_RAW.tar", exdir = "GSE68777/idat")
head(list.files("GSE68777/idat", pattern = "idat"))

idatFiles <- list.files("GSE68777/idat") 
idatFiles_path <- paste0("GSE68777/idat/", idatFiles)
sapply(idatFiles_path, gunzip, overwrite=TRUE)

rgSet <- read.metharray.exp("GSE68777/idat")
rgSet

#these two not very useful
pData(rgSet)
head(sampleNames(rgSet))

#so get source info
geomat <- getGEO("GSE68777")
pD.all <- pData(geomat[[1]])
pD <- pD.all[, c("title", "geo_accession", "characteristics_ch1.1", "characteristics_ch1.2")]
head(pD)

# clean up 
names(pD)[c(3, 4)] <- c("group", "sex")
pD$group <- sub("^diagnosis: ", "", pD$group)
pD$sex <- sub("^Sex: ", "", pD$sex)
head(pD)
sampleNames(rgSet) <- sub(".*_5", "5", sampleNames(rgSet))
rownames(pD) <- pD$title
pD <- pD[sampleNames(rgSet),]
head(sampleNames(rgSet))
head(pD)

pData(rgSet) = DataFrame(pD)

rgSet

# normalise data
grSet = preprocessQuantile(rgSet)
# genomic ranges set
grSet

granges(grSet)
# cpg island status
head(getIslandStatus(grSet))

#analyse
getBeta(grSet)[1:3, 1:3]
