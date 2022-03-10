#q1
library(AnnotationHub)

# create a local annotation hub
ah <- AnnotationHub()

ah <- subset(ah, species == "Homo sapiens")
ah
qhs <- query(ah, "CpG Islands")
qhs$genome
#item 1 which is hg19 is our favourite
cpg <- qhs[[1]]
cpg
seqlevels(cpg)

chromosome_list <- paste("chr", 1:22, sep = "")
chromosome_list

cpg_auto <- cpg
seqlevels(cpg_auto, pruning.mode = "coarse") = chromosome_list

seqlevels(cpg_auto)
cpg_auto
#26641

#q2
cpg_c4 <- cpg_auto
seqlevels(cpg_c4, pruning.mode = "coarse") = c("chr4")
cpg_c4
# 1031

#q3 
qhs <- query(ah, c("H3K4me3", "H1"))
qhs
qhs2 <- query(ah, c("H3K4me3", "E003"))
qhs2
# use narrow peak
gr1 <- qhs2[[2]]
gr1
seqlevels(gr1, pruning.mode = "coarse") = chromosome_list
gr1
sum(width(reduce(gr1, ignore.strand = TRUE)))
# 41135164

#q4
qhs3 <- query(ah, c("H3K27me3", "E003"))
qhs3
gr2 <- qhs3[[2]]
gr2
seqlevels(gr2, pruning.mode = "coarse") = chromosome_list
gr2
mean(gr2$signalValue)

#  4.770728

# q5
length(subsetByOverlaps(gr1, gr2, ignore.strand = TRUE))
findOverlaps(gr1, gr2)
subsetByOverlaps(gr1, gr2, ignore.strand = TRUE)
sum(width(intersect(gr1, gr2, ignore.strand =TRUE)))
# 10289096

# q6
bivalent <- intersect(gr1, gr2, ignore.strand =TRUE)
bivalent
length(subsetByOverlaps(bivalent, cpg_auto, ignore.strand = TRUE)) / length(bivalent)
#  0.5383644

# q7
length(subsetByOverlaps(cpg_auto, bivalent, ignore.strand = TRUE)) / length(cpg_auto)
# 0.2915806
#WRONG
sum(width(intersect(bivalent, cpg_auto, ignore.strand = TRUE))) / sum(width(cpg_auto))
# 0.241688


# q8

bivalent_resiz <- flank(bivalent, 10000)
#length(subsetByOverlaps(cpg_auto, bivalent,  ignore.strand = TRUE))
#length(subsetByOverlaps(bivalent_resiz, cpg_auto,   ignore.strand = TRUE))
sum(width(intersect(bivalent_resiz, cpg_auto, ignore.strand =TRUE)))
# 7397073
#WRONG

# q9
sum(width(cpg_auto)) / sum(seqlengths(cpg_auto))
# 0.007047481

# q10
inOut = matrix(0, ncol=2, nrow=2)
colnames(inOut) = c("in", "out")
rownames(inOut) = c("in", "out")
inOut[1,1] = sum(width(intersect(bivalent, cpg_auto, ignore.strand = TRUE)))
inOut[1,2] = sum(width(setdiff(bivalent, cpg_auto, ignore.strand = TRUE)))
inOut[2,1] = sum(width(setdiff(cpg_auto, bivalent, ignore.strand =TRUE)))
inOut
colSums(inOut)
rowSums(inOut)
inOut[2,2] = sum(seqlengths(cpg_auto)) - sum(inOut)
inOut
oddsRatio = inOut[1,1] * inOut[2,2] / (inOut[2,1] * inOut[1,2])
oddsRatio
# 169.0962