# Rle

library(GenomicRanges)
# construct rle
rl = Rle(c(1,1,1,1,1,2,2,2,2,2,4,4,2))
rl
runLength(rl)
runValue(rl)
as.numeric(rl)

ir = IRanges(start = c(2,8), width = 4)
aggregate(rl, ir, FUN = mean)

vec = as.numeric(rl)
mean(vec[2:5])
mean(vec[8:11])

ir = IRanges(start = 1:5, width=3)
ir

# coverage of the sequence
coverage(ir)

# slice
# vector greater or equal to 2
slice(rl, 2)

slice(rl, 3)

vi = Views(rl, IRanges(2,8))
vi

vi = Views(rl, IRanges(c(2,8), width = 2))
vi

gr <- GRanges(seqnames = "chr1", ranges = IRanges(start = 1:10, width=3))
rl <- coverage(gr)
rl

vi = Views(rl, GRanges("chr1", ranges = IRanges(3,7)))
vi

vi = Views(rl, as(GRanges("chr1", ranges = IRanges(3,7)), "Ranges"))
vi
vi$chr1

gr


# granges lists
gr1 <- GRanges(seqnames = "chr1", ranges = IRanges(start = 1:4, width = 3))
gr2 <- GRanges(seqnames = "chr2", ranges = IRanges(start = 1:4, width = 3))
gL = GRangesList(gr1=gr1, gr2=gr2)
gL
gL[[1]]
gL$gr1

start(gL)
seqnames(gL)

#much faster than sapply
elementNROWS(gL)

sapply(gL, length)

shift(gL, 10)

findOverlaps(gL, gr2)
