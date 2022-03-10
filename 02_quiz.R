# q1
# What is the GC content of “chr22” in the “hg19” build of the human genome?
# Tip: The reference genome includes “N” bases; you will need to exclude those.
library(BSgenome)
library(Biostrings)

available.genomes()
# BiocManager::install("BSgenome.Hsapiens.UCSC.hg19")
library("BSgenome.Hsapiens.UCSC.hg19")
Hsapiens
letterFrequency(Hsapiens$chr22, "GC")
length(Hsapiens$chr22)
letterFrequency(Hsapiens$chr22, "N")

ans1 = letterFrequency(Hsapiens$chr22, "GC") / (length(Hsapiens$chr22) - letterFrequency(Hsapiens$chr22, "N"))
ans1
# 0.4798807
# CORRECT

#q2
# Background: In the previous assessment we studied H3K27me3 “narrowPeak” regions from the H1 cell line (recall that the Roadmap ID for this cell line is “E003”). 
# We want to examine whether the GC content of the regions influence the signal; 
# in other words wether the reported results appear biased by GC content.

# Question: What is mean GC content of H3K27me3 “narrowPeak” regions from Epigenomics Roadmap from the H1 stem cell line on chr 22.

# Clarification: Compute the GC content for each peak region as a percentage and then average those percentages to compute a number between 0 and 1.
library(AnnotationHub)
ah <- AnnotationHub()
ah <- subset(ah, species == "Homo sapiens")
qhs3 <- query(ah, c("H3K27me3", "E003"))
qhs3
gr2 <- qhs3[[2]]
gr2.chr22 <- subset(gr2, seqnames == "chr22")
gr2.chr22
peakSeq  = getSeq(Hsapiens, gr2.chr22)
peakSeq
#param = new("BSParams", X = peakSeq, FUN=letterFrequency)
#bsapply(param, "GC")
mean(letterFrequency(peakSeq, "GC") / unlist(seqlengths(peakSeq)))
# 0.528866
# CORRECT

#q3
# The “narrowPeak” regions includes information on a value they call “signalValue”.
# Question: What is the correlation between GC content and “signalValue” of these regions (on chr22)?
cor(unlist(letterFrequency(peakSeq, "GC")) / unlist(seqlengths(peakSeq)), gr2.chr22$signalValue)
# 0.004467924
# CORRECT

#q4
# 
qhs3
fc_sig = qhs3[[4]]
rle.obj = import(fc_sig,
                 which=gr2.chr22, as="Rle")

fc.signal.views = Views(rle.obj$chr22,
                        start=start(gr2.chr22),
                        end=end(gr2.chr22))

# fc.signal.mean = aggregate(fc.signal.views, gr2.chr22, FUN = mean)

fc.signal.mean = mean(fc.signal.views)
cor(fc.signal.mean, gr2.chr22$signalValue)
#  0.9149614
# CORRECT

# q5 
# Question: How many bases on chr22 have an fc.signal greater than or equal to 1?
gs.fc = import(fc_sig, which = GRanges("chr22", ranges = IRanges(1, 51304566)))
rle.obj2 = import(fc_sig, which = GRanges("chr22", ranges = IRanges(1, 51304566)), as="Rle")
rle.obj2 = rle.obj2[rle.obj2>=1]
length(rle.obj2$chr22)

# 10914671
# CORRECT

# q6
# fc signal for E055
# Question: Identify the regions of the genome where the signal in E003 is 0.5 or lower and the signal in E055 is 2 or higher.
qhs4 <- query(ah, c("H3K27me3", "E055"))
qhs4
fc.signal.e055 = qhs4[[4]]
gr.obj.e055 = import(fc.signal.e055, which = GRanges("chr22", ranges = IRanges(1, 51304566)))
gr.obj.e003 = import(fc_sig, which = GRanges("chr22", ranges = IRanges(1, 51304566)))
gr.obj.e055.filt = gr.obj.e055[gr.obj.e055$score >= 2]
gr.obj.e003.filt = gr.obj.e003[gr.obj.e003$score <= 0.5]
intie = intersect(gr.obj.e055.filt, gr.obj.e003.filt)
sum(coverage(intie))
# 1869937
# CORRECT

# q7
# Question: What is the average observed-to-expected ratio of CpG dinucleotides for CpG Islands on chromosome 22?
# Specifically, the observed CpG frequency is just the number of “CG” dinucleotides in a region. 
# The expected CpG frequency is defined as the frequency of C multiplied by the frequency of G divided by the length of the region.
qhs <- query(ah, "CpG Islands")
qhs$genome
#item 1 which is hg19 is our favourite
cpg <- qhs[[1]]
cpg_seq <- getSeq(Hsapiens, cpg)
cpg_seq_cg_freq = as.vector(dinucleotideFrequency(cpg_seq)[,"CG"])

c_freq = letterFrequency(cpg_seq, "C")
g_freq = letterFrequency(cpg_seq, "G")
cpg_width = width(cpg_seq)

cpg_stats = data.frame(
  'cg' = as.vector(cpg_seq_cg_freq),
  c = c_freq,
  g = g_freq,
  cpg_len = cpg_width
)

observed_expected <- function(x){
  x['cg']/((x['C'] * x['G'])/x['cpg_len'])}

mean(apply(cpg_stats, 1, observed_expected))
# ? 0.8341
# CORRECT


# q8
# Question: How many TATA boxes are there on chr 22 of build hg19 of the human genome?
# Clarification: You need to remember to search both forward and reverse strands.
?oligonucleotideFrequency
of_for <- oligonucleotideFrequency(Hsapiens$chr22, 6)['TATAAA']
of_rev <- oligonucleotideFrequency(Hsapiens$chr22, 6)['TTTATA']
of_for+of_rev
# 27263
# CORRECT

# q9
# Question: How many promoters of transcripts on chromosome 22 containing a coding sequence, 
# contains a TATA box on the same strand as the transcript?
# BiocManager::install("TxDb.Hsapiens.UCSC.hg19.knownGene")

# Clarification: Use the TxDb.Hsapiens.UCSC.hg19.knownGene package to define transcripts and coding sequence. 
# Here, we defined a promoter to be 900bp upstream and 100bp downstream of the transcription start site.
library("TxDb.Hsapiens.UCSC.hg19.knownGene")
txdb = TxDb.Hsapiens.UCSC.hg19.knownGene
txdb

tx_cds <- cdsBy(txdb, 'tx')
tx_cds_ch22 <- keepSeqlevels(tx_cds, 'chr22', pruning.mode = "coarse")
tx_cds_ch22_names <- names(tx_cds_ch22)
tx_subset <- subset(transcripts(txdb), tx_id %in% tx_cds_ch22_names)
tx_subset_promoters <- promoters(tx_subset, upstream = 900, downstream = 100)
tata_in_chr22 <- subset(vmatchPattern('TATAAA', Hsapiens), seqnames = "chr22")
overlaps <- findOverlaps(tx_subset_promoters, tata_in_chr22)
subsetByOverlaps(tx_subset_promoters, tata_in_chr22)
# 193
# CORRECT

