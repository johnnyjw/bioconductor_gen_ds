library(AnnotationHub)

# create a local annotation hub
ah <- AnnotationHub()
ah

ah[1]
ah[[1]]
unique(ah$dataprovider)
unique(ah$species)

#searching / narrowing down data
ah <- subset(ah, species == "Homo sapiens")
ah
query(ah, "H3K4me3")

query(ah, c("H3K4me3", "Gm12878"))

# to display on shiny app
# selected item goes to object
ah2 <- display(ah)
ah2

##########################################################################
# use case
# histone mark - active promotors?
ahub <- AnnotationHub()
ahub <- subset(ahub, species == "Homo sapiens")

#query for a certain type of histone methylation on Gm12878 cell line
qhs <- query(ahub, c("H3K4me3", "Gm12878"))
qhs

gr1 <- qhs[[2]]
gr2 <- qhs[[4]]
gr1
summary(width(gr1))

summary(width(gr2))
table(width(gr2))

peaks = gr2
qhs[4]

# looking for a refseq item
qhs <- query(ahub, "RefSeq")
qhs
qhs$genome
genes <- qhs[[1]]
genes

table(table(genes$name))
prom <- promoters(genes)
table(width(prom))
args(promoters)

prom
peaks

# is this histone modification (particular kind in 'peaks') enriched in promotors?
findOverlaps(prom, peaks)
# we are interested in how big a % do promotors overlap
ov <- findOverlaps(prom, peaks)
length(unique(queryHits(ov)))
length(unique(subjectHits(ov)))

length(subsetByOverlaps(peaks, prom, ignore.strand = TRUE))

length(subsetByOverlaps(peaks, prom, ignore.strand = TRUE)) / length(peaks)
# 30% -- is that good?
length(subsetByOverlaps(prom, peaks, ignore.strand = TRUE)) / length(prom)
# 50 % of all promotors  - pretty good

# how many bases do the peaks cover
sum(width(reduce(peaks, ignore.strand = TRUE))) / 10^6

# how much length do the promotors cover
sum(width(reduce(prom, ignore.strand = TRUE))) / 10^6

# how big is the overlap
sum(width(intersect(peaks, prom, ignore.strand =TRUE))) / 10^6

# 3billion bases in human genome
inOut = matrix(0, ncol=2, nrow=2)
colnames(inOut) = c("in", "out")
rownames(inOut) = c("in", "out")
inOut
inOut[1,1] = sum(width(intersect(peaks, prom, ignore.strand = TRUE)))
inOut[1,2] = sum(width(setdiff(peaks, prom, ignore.strand = TRUE)))
inOut[2,1] = sum(width(setdiff(prom, peaks, ignore.strand =TRUE)))
inOut
colSums(inOut)
rowSums(inOut)

# last column - total genome (3bill)
inOut[2,2] = 3*10^9 - sum(inOut)
inOut

# with this 2x2 table can calculate odds ratio
fisher.test(inOut)$statistic
# but this doesnt work as the integer number too high...so need to do it by hand
oddsRatio = inOut[1,1] * inOut[2,2] / (inOut[2,1] * inOut[1,2])
oddsRatio

# is 3 billion the right number? Maybe the mappable part of human genome more important
# lets guess 1/2 human genome  1.5 billion
inOut[2,2] = 0
inOut[2,2] = 1.5 * 10^9 - sum(inOut)
inOut
oddsRatio = inOut[1,1] * inOut[2,2] / (inOut[2,1] * inOut[1,2])
oddsRatio
# 8.9 suggests some enrichment