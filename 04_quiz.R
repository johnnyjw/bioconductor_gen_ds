library(yeastRNASeq)
library(ShortRead)
fastqFilePath <- system.file("reads", "wt_1_f.fastq.gz", package = "yeastRNASeq")

fqFile <- FastqFile(fastqFilePath)
reads <- readFastq(fqFile)

all_seq <- sread(reads)

ranges <- IRanges(start = 5,
                  width = 1)
all_seq[,ranges]
all_seq
unlist(extractAt(all_seq, ranges))
