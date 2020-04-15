#
# getwd()
# 
# Windows format:
# setwd("c:/Documents/my/working/directory")
#
# Mac and Linux format:
# setwd('/Users/username/')
# 


#
# Open csv (comma separated values) file 
# a generic data set with categorical label with measurement value
# 999 is used to represent data not available (NA) 
#
dat <- read.csv("data/example-data.csv")

#
# perform a simple two sample t-test
# for more information type ?t.test in console
#
#  An expression of the form y ~ model is interpreted as a specification 
#  that the response y is modelled by a predictor specified symbolically 
#  by model.
#
t.test(measurement ~ type, data = dat)

# Generate a boxplot define measurement is related
boxplot(measurement ~ type, data = dat)
points(measurement ~ type, data = dat)

#
# jitter is function that add a little noise to the point
# to allow for better point visualization
#
points(measurement ~ jitter(as.numeric(type)), data = dat)

# Filter out NA values from t-test 
dat <- dat[dat$measurement!=999, ]

t.test(measurement ~ type, data = dat)

