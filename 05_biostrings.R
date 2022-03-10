library(Biostrings)

dna1 <- DNAString('ACTC-G')
dna1

#dna set
dna2 <- DNAStringSet(c("ACG", "ACGT", "ACGGT"))
dna2

IUPAC_CODE_MAP

dna1[2:4]
dna2[1:2]

dna2[[1]]

names(dna2) <- paste0("Seq", 1:3)
dna2

width(dna2)


sort(dna2)
rev(dna2)

rev(dna1)

reverse(dna2)

reverseComplement(dna2)

translate(dna2)

alphabetFrequency(dna2)
letterFrequency(dna2, letters = "GC")
dinucleotideFrequency(dna2)
consensusMatrix(dna2)
