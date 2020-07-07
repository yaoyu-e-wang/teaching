#
# This script is the final result of all the analyses we will perform in the Data Wrangling with R
# introduction section
#
############################  Introduction  ###########################################
# This is a comment that is not read by R 
# Comment is for human to read to explain the code
x=5
x<-10
# Declare a character
y="a"
# Declare a logic
z=TRUE

###########################################
#
# We can also find the data type of a variable by **class** function.
class(x)
class(y)
class(z)

###########################################
#
# For example we can declare the following vectors with different data types:
#
# declare a numeric vector
x=c(1,2,3,4)
# declare a character vector
y=c("a","b","c")
# declare a logic vector, TRUE and FALSE can be short-handed as T and F
z=c(TRUE,FALSE,T,F)

# A data type is assumed when a mixed data types are declared
a=c(1,2,3,"f")
class(a)  # a would be a character vector

b=c(5, F)   # this is a numeric vector  F becomes 0

###########################################
#
# The data type can also be changed into different types:
#
class(x)
x=as.character(x)
x
class(x)
x=as.numeric(x)
x

###########################################
#
# We can get the length of the vector as well as access values by calling for their specific positions using square brackets.
#
# length can be determined 
len=length(x)
len

x[1]
x[2]

###########################################
#
# a list variable can be declare as a vector of vectors with or without names
#
y=list(a=1, 17, b=2:5, c='a')
y[[1]]
y[[2]]
y$b

###########################################
#
# We can also get the names from the list variable and rename the elements:
#
names(y)
names(y)=c("a","b","c","d")
names(y)
y$b

###########################################
#
#**Data Frame** is a special kind of list to store rectangular data sets  Think of it as excel sheet where each column is a list.  
#
df.y=data.frame(y)
df.y
df.y$a[1]
df.y$c[1]
df.y$c

###########################################
#
# **Factor** is a data type unique to R and other statistical language, it condenses variables into unique categorical variables and store them as levels.  Makes computation of large data with repeated string/character data type more efficient
#
x=sample(letters[1:5], 10, replace = TRUE)
x
y=factor(x)
y
levels(y)
as.numeric(y)

###########################################
#
# Fundamental Operators
#
x=matrix(c(1,2,3,4), nrow=2, ncol=2)
y=matrix(c(5,6,7,8), nrow=2, ncol=2)
x+y
x/y

# or compare them
x>y
x<y

a=FALSE
isTRUE(a)

############################################################################################################################

