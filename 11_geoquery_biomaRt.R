# geo = geneexpression ominbus
# BiocManager::install("GEOquery")
library(GEOquery)

eList = getGEO("GSE11675")
length(eList)
names(eList)
eData = eList[[1]]
eData

names(pData(eData))

elist2 = getGEOSuppFiles("GSE11675")
# raw data stored as tar files.  Can see them below
elist2

# biomaRt
# BiocManager::install('biomaRt')
library(biomaRt)
head(listMarts())

mart <- useMart("ensembl")
mart

head(listDatasets(mart))

ensembl <- useDataset("hsapiens_gene_ensembl", mart)

# affymetrix probe id
values <- c("202763_at", "209310_s_at", "207500_at")
getBM(attributes = c("ensembl_gene_id", "affy_hg_u133_plus_2"),
      filters = "affy_hg_u133_plus_2", values =values, mart = ensembl)

attributes <- listAttributes(ensembl)
nrow(attributes)
head(attributes)
tail(attributes, n=500)

# lots of filters to search with
filters <- listFilters(ensembl)
head(filters)
nrow(filters)
attributePages(ensembl)
attributes <- listAttributes(ensembl, page = "feature_page")
attributes

# watch out for attributes that span more than one page.
# pull one attribute at a time and merge back
