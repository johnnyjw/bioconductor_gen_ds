library(IRanges)

ir1 <- IRanges(start = c(1,3,5), end = c(3,5,7))
ir1

ir2 <- IRanges(start = c(1,2,3), width=3)
ir2

names(ir1) <- paste("A", 1:3, sep = "")
ir1
#no dimensions
dim(ir1)
# but can subset
ir1[1]
ir1["A1"]

plotRanges <- function(x, xlim = x, main = deparse(substitute(x)),
                       col = "black", sep = 0.5, ...) {
  height <- 1
  if (is(xlim, "Ranges"))
    xlim <- c(min(start(xlim)), max(end(xlim)))
  bins <- disjointBins(IRanges(start(x), end(x) + 1))
  plot.new()
  plot.window(xlim, c(0, max(bins)*(height+sep)))
  ybottom <- bins * (sep + height) - height
  rect(start(x)-0.5, ybottom, end(x)+0.5, ybottom + height, col = col, ...)
  title(main)
  axis(1)
}

par(mfrow = c(2,1))
ir <- IRanges(start = c(1,3,7,9), end = c(4,4,8,10))
plotRanges(ir)

plotRanges(reduce(ir))

plotRanges(disjoin(ir))

ir

resize(ir, width = 1, fix = "start")
resize(ir, width = 1, fix = "center")

ir1 <- IRanges(start = c(1, 3, 5), width =1)
ir2 <- IRanges(start = c(4,5,6), width =1)
ir1
ir2
union(ir1, ir2)
intersect(ir1, ir2)

#test overlap function
ir1 <- IRanges(start = c(1,4,8), end = c(3,7,10))
ir2 <- IRanges(start = c(3,4), width =3)
ir1
ir2
ov <- findOverlaps(ir1, ir2)
ov
queryHits(ov)
subjectHits(ov)

countOverlaps(ir1, ir2)
ir1
ir2

nearest(ir1, ir2)
