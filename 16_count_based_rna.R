library(DESeq2)
library(edgeR)

library(airway)
data(airway)
# 8 samples, 64000 features
airway
assay(airway, "counts")[1:3, 1:3]

# dex covariant
airway$dex

# the reference level in R tends to be the +ve.  So to make more sense relevel reference
airway$dex <- relevel(airway$dex, "untrt")
airway$dex

# gene details of the features
granges(airway)

# getting gene by sample/count matrix not easy
# Get an annotation, select features and count overlap
# What do we do with genes with multiple transcripts?
# How you select has impact on the statistical results!
library(edgeR)

# take data from expt and convert to limma class
dge <- DGEList(counts = assay(airway, "counts"),
               group = airway$dex)

dge$samples

# by=0 means merging on the basis of rownames
dge$samples <- merge(dge$samples,
                     as.data.frame(colData(airway)),
                     by = 0)
dge$samples
# nothing in here yet
head(dge$genes)
# annotate
dge$genes <- data.frame(name = names(rowRanges(airway)),
                        stringsAsFactors = FALSE)
head(dge$genes)

# calculate normalisation factor
dge <- calcNormFactors(dge)
# now has norm.factors
dge$samples

#estimate dispersal/variablility
dge <- estimateGLMCommonDisp(dge)
dge <- estimateGLMTagwiseDisp(dge)

#design matrix
design <- model.matrix(~dge$samples$group)
head(design)

#fit model
fit <- glmFit(dge, design)

# figure out for coffs - which genes are most signficant
# likelihood ratio test
# coef=2 - testing the second coefficient in design matrix (here=trt/nontrt)
lrt <- glmLRT(fit, coef = 2)
#toptags - features most differentially expressed across treatment groups
topTags(lrt)


##### USING DESEQ2
library(DESeq2)
# bring in data
# this needs to include the design matrix
# when you are analysing multiple variables in design, deseq focusses on last variable
dds <- DESeqDataSet(airway, design = ~ dex)

# fit
dds <- DESeq(dds)

# results
res <- results(dds)
res
# order by adjusted p value
res <- res[order(res$padj),]
res
# compare with edgeR
#NOTE: There is no overlaps in the top genes found by the two methods!!!
topTags(lrt)
