library(ALL)

data(ALL)
ALL

experimentData(ALL)

?ALL

exprs(ALL)[1:4, 1:4]
head(sampleNames(ALL))
head(featureNames(ALL))

head(pData(ALL))

pData(ALL)$sex
ALL$sex

#subsetting
ALL[, 1:5]

ALL[1:10, ]

ALL[1:10, 1:5]

#often empty
featureData(ALL)

ids = featureNames(ALL)[1:5]
ids

# annotating
# from details below - annotation = hgu95av2
library(hgu95av2.db)
# BiocManager::install('hgu95av2.db')

# map ids from array to entrez
as.list(hgu95av2ENTREZID[ids])

phenoData(ALL)
names(pData(ALL))

head(pData(ALL))

pData(phenoData(ALL))

# summarised experiment
# BiocManager::install('airway')
library(airway)
data(airway)
airway
colData(airway)
airway$cell

colnames(airway)
head(rownames(airway))

#assay accesser
# step one find assay names
assayNames(airway)
# access
assay(airway, "counts")[1:4, 1:4]

#rowranges
length(rowRanges(airway))
rowRanges(airway)

sum(elementNROWS(rowRanges(airway)))

start(airway)

# subsetting
gr = GRanges("1", range = IRanges(start =1 , end = 10^7))
gr
subsetByOverlaps(airway, gr)



