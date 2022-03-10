library(BSgenome)
library(Biostrings)

available.genomes()

library("BSgenome.Scerevisiae.UCSC.sacCer2")
Scerevisiae
seqnames(Scerevisiae)
seqlengths(Scerevisiae)
# so far sequences have not been loaded into memory

# but this loads into memory
Scerevisiae$chrI

letterFrequency(Scerevisiae$chrI, "GC")
letterFrequency(Scerevisiae$chrI, "GC", as.prob = TRUE)


# param to run across each sequence in genome in a memory controlling manner
param = new("BSParams", X = Scerevisiae, FUN=letterFrequency)
bsapply(param, "GC")
unlist(bsapply(param, "GC"))

#to get a single measure - add all gc's up and divide by the sum of the lengths of the chromosomes
sum(unlist(bsapply(param, "GC"))) / sum(seqlengths(Scerevisiae))

unlist(bsapply(param, "GC", as.prob = TRUE))

# matching
dnaseq <- DNAString("ACGTACGT")
dnaseq

matchPattern(dnaseq, Scerevisiae$chrI)

countPattern(dnaseq, Scerevisiae$chrI)

vmatchPattern(dnaseq, Scerevisiae)

dnaseq == reverseComplement(dnaseq)

# BSgenome Views

dnaseq <- DNAString("ACGTACGT")

vi = matchPattern(dnaseq, Scerevisiae$chrI)
vi
ranges(vi)
Scerevisiae$chrI[57932:57939]
alphabetFrequency(vi)
shift(vi, 10)
gr=vmatchPattern(dnaseq, Scerevisiae)
gr
vi2= Views(Scerevisiae, gr)
vi2

# loading up promotors
library(AnnotationHub)
ahub = AnnotationHub()
qh = query(ahub, c("sacCer2", "genes"))
qh
genes = ahub[["AH7048"]]
prom = promoters(genes)

prom
prom = trim(prom)
prom

promViews = Views(Scerevisiae, prom)
promViews
gcProm = letterFrequency(promViews, "GC", as.prob = TRUE)
plot(density(gcProm))
abline(v = 0.38)


