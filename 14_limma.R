library(limma)
library(leukemiasEset)

data(leukemiasEset)
leukemiasEset
table(leukemiasEset$LeukemiaType)

ourData <- leukemiasEset[, leukemiasEset$LeukemiaType %in% c("ALL", "NoL")]
ourData$LeukemiaType
# clean away levels that no longer exist in the data (would would mess up models)
ourData$LeukemiaType <- factor(ourData$LeukemiaType)
ourData$LeukemiaType

# setting up design matrix
design <- model.matrix(~ ourData$LeukemiaType)
head(design)

# fit basic model
fit <- lmFit(ourData, design)
fit <- eBayes(fit)
topTable(fit)

#look at the top gene
topTable(fit, n = 1)
# plug this gene into
genename <- rownames(topTable(fit, n=1))
genename

typeMean <- tapply(exprs(ourData)[genename,], ourData$LeukemiaType, mean)
typeMean
# log fold change
typeMean["NoL"] - typeMean["ALL"]

# design 2 : notice the minus 1 - results in different paramaeters
design2 <- model.matrix(~ ourData$LeukemiaType - 1)
head(design2)
colnames(design2) <- c("ALL", "NoL")
design2

# make a contrast
# a contrast is a hypothesis
fit2 <- lmFit(ourData, design2)
# the opposite of the ref in previous design (NoL now reference)
contrast.matrix <- makeContrasts("ALL-NoL", levels = design2)
contrast.matrix
# testing if these two parameters equal to zero

fit2C <- contrasts.fit(fit2, contrast.matrix)
fit2C <- eBayes(fit2C)
topTable(fit2C)
