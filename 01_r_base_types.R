x = 1:10
x
class(x)
names(x) = letters[1:10]
x[1:3]
x[c("a", "b")]

x = 1:3
names(x) = c("a", "a", "c")
x["a"]
anyDuplicated(names(x))
names(x) = letters[1:3]
anyDuplicated(names(x))

x = 1
class(x)

x=1:3
class(x)

x=1L
class(x)

.Machine$integer.max

# matrices
x = matrix(1:9, ncol = 3, nrow = 3)
x
rownames(x) = letters[1:3]
nrow(x)
ncol(x)
x
x[1:2,]
x[,1:2]
x[1:2,1:2]
x["a",]
x["a",,drop=FALSE]
x[x>5]
x = matrix(1:9, ncol = 3, nrow = 3, byrow=TRUE)
x

#lists
x = list(a = rnorm(3), b = letters[1:5], matrix)
x

x[1:2]
x[1]
x[[1]]
x["a"]
x$a

as.list(1:3)


x = list(rnorm(3), 3:9)
x
lapply(x, mean)
unlist(lapply(x, mean))
sapply(x, mean)

# data frame
x <- data.frame(
  sex = c("M", "M", "F"),
  age = c(32,34,29)
)
x
x$sex
x[["sex"]]
x[1,"sex"]
sapply(x, class)
as.matrix(x)
as.list(x)

library(methods)
as(x, "matrix")
