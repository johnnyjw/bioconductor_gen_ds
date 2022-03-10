library(rtracklayer)

#?import

library(AnnotationHub)
ahub <- AnnotationHub()

# shows the type of files available
#table(ahub$rdataclass)
ahub.bw = subset(ahub, rdataclass == "BigWigFile" & species == "Homo sapiens")
#example loading first one
bw = ahub.bw[[1]]
#bw

gr.chr22 = import(bw, which = GRanges("chr22", ranges = IRanges(1, 10^8)))                                                                                                                         
#gr.chr22

rle.chr22 = import(bw, which = GRanges("chr22", ranges = IRanges(1, 10^8)), as = "Rle")
#wrong returns all chromosomes
rle.chr22
# right
rle.chr22$chr22

#chainfiles
ahub.chain = subset(ahub, rdataclass == "ChainFile")
ahub.chain
ahub.chain = subset(ahub, species == "Homo sapiens")
ahub.chain
query(ahub.chain, c("hg18", "hg19"))
chain = query(ahub.chain, c("hg18", "hg19"))[[3]]
gr.chr22 = import(bw, which = GRanges("chr22", ranges = IRanges(1, 10^8)))

#convert
gr.hg18 = liftOver(gr.chr22, chain)
class(gr.hg18)
length(gr.hg18)
length(gr.chr22)
# the same length

table(elementNROWS(gr.hg18))
# shows value=0 no map, value=1 1:1 map value=2orabove -> fragments got broken up

