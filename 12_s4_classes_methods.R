library(ALL)
library(GenomicRanges)

df <- data.frame(y = rnorm(10), x = rnorm(10))

#making an lm object
lm.object <- lm(y ~ x, data = df)
lm.object
class(lm.object)
names(lm.object)

# in S3 a class is really just a list with an attribute
xx = list(a = letters[1:3], b=rnorm(4))
xx
#can make it a class
class(xx) = "lm"
xx
# but it looks nothing like a linear model!!

#S4 system has validity checking built in 
# useful for complicated data structures
data(ALL)
ALL
class(ALL)
# is it an S4 object?
isS4(ALL)

#getting help on a class either or below
class?ExpressionSet
?"ExpressionSet-class"

# constructor
xx = list(a = 1:3)
ExpressionSet()

new("ExpressionSet")

# class attributes
getClass("ExpressionSet")
#slots - either or 
ALL@annotation
slot(ALL, "annotation")
# but you dont access directly - you use an accessor function
# Often you're not supposed to access the slots
#  The accessor functions show you which ones you're meant to access
annotation(ALL)
# accessor functions found in help pages  (Methods)
?ExpressionSet

# check that object satisfies class definition
validObject(ALL)


# S4 METHODS
library(GenomicRanges)

as.data.frame
# 'standard generic' in display -> S4 method

# s3 method below doesnt have this
base::as.data.frame

### look at the classes supported
showMethods("as.data.frame")

### see the code that gets run on a particular class
getMethod("as.data.frame", "GenomicRanges")
# the above is a shortened version of ...
getMethod("as.data.frame", signature(x = "GenomicRanges"))

# help
# each class on a method might have its own help
method?"as.data.frame,DataFrame"
?"as.data.frame,GenomicRanges-method"

showMethods("findOverlaps")
# findOverlaps will run different code for combinations of query and subject
getMethod("findOverlaps", signature(query = "GenomicRanges", subject = "GRangesList"))
?"findOverlaps,GenomicRanges,GRangesList-method"
