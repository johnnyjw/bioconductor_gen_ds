library(GenomicRanges)

gr = GRanges(seqnames = c("chr1"), strand = c("+", "-", "+"), ranges = IRanges(start = c(1,3,5), width =3))
gr

flank(gr, 5)
promoters(gr)
seqinfo(gr)
seqlengths(gr) = c("chr1" = 10)
seqinfo(gr)
seqlevels(gr)
gaps(gr)

# assign new chromosomes
seqlevels(gr) = c("chr1", "chr2")
seqnames(gr) = c("chr1", "chr2", "chr1")
gr

genome(gr) = "hg19"
seqinfo(gr)

gr2=gr
genome(gr2) = "hg18"
findOverlaps(gr, gr2)

library(BiocGenerics)

ir = IRanges(start = 1:3, width=2)
df = DataFrame(ir = ir, score = rnorm(3))
df
df$ir


gr <- GRanges(seqnames  ="chr1", strand = c("+", "-", "+"),
              ranges = IRanges(start = c(1,3,5), width = 3))
gr
values(gr) = DataFrame(score = rnorm(3))
gr
values(gr)
mcols(gr)
gr$score
gr$score2 <- gr$score/3
gr

gr2 <- GRanges(seqnames = c("chr1", "chr2", "chr1"), strand = "*",
               ranges = IRanges(start = c(1,3,5), width = 3))
gr2
findOverlaps(gr, gr2)
findOverlaps(gr, gr2, ignore.strand = TRUE)

subsetByOverlaps(gr2, gr)

df = data.frame(chr = "chr", start = 1:3, end = 4:6, score = rnorm(3))
makeGRangesFromDataFrame(df)
makeGRangesFromDataFrame(df, keep.extra.columns = TRUE)


gr = GRanges(seqnames = c("chr1", "chr2"), ranges = IRanges(start = 1:2, end = 4:5))
gr
seqlevels(gr, pruning.mode="coarse") = "chr1"
gr

gr = GRanges(seqnames = c("chr1", "chr2"), ranges = IRanges(start = 1:2, end = 4:5))
dropSeqlevels(gr, "chr1", pruning.mode = "coarse")
gr

gr = GRanges(seqnames = c("chr1", "chr2"), ranges = IRanges(start = 1:2, end = 4:5))
newStyle <- mapSeqlevels(seqlevels(gr), "NCBI")
newStyle

gr = renameSeqlevels(gr, newStyle)
gr
