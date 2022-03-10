library(GenomicFeatures)
BiocManager::install("TxDb.Hsapiens.UCSC.hg19.knownGene")
library(TxDb.Hsapiens.UCSC.hg19.knownGene)

txdb = TxDb.Hsapiens.UCSC.hg19.knownGene
txdb

gr <- GRanges(seqnames = "chr1", strand = "+", ranges = IRanges(start = 11874, end = 14409))

genes(txdb)

subsetByOverlaps(genes(txdb), gr)

subsetByOverlaps(genes(txdb), gr, ignore.strand = TRUE)

subsetByOverlaps(transcripts(txdb), gr)

subsetByOverlaps(exons(txdb), gr)

c

# coding sequences in different splice transcripts. Not very useful
subsetByOverlaps(cds(txdb), gr)
# better
subsetByOverlaps(cdsBy(txdb, by = "tx"), gr)
# outputs exons of the second transcript

subsetByOverlaps(exonsBy(txdb, by = "tx"), gr)["2"]

subset(transcriptLengths(txdb, with.cds_len=TRUE), gene_id == "100287102")

sum(width(subsetByOverlaps(cdsBy(txdb, by = "tx"), gr) [["2"]]))
