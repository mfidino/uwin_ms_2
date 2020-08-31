
# Set seed for reproducibility
set.seed(200)

# Step 1. Fit a model with group-mean centering

# generate a vector of different number of sites
#  per city.
sites <- c(400, 30, 75, 25, 60, 200, 25, 80, 23, 100)


# the number of cities
ncity <- length(sites)

# create a vector of means (e.g., average housing density)
#  for each city and scale them.
means <- rnorm(
  length(sites),
  50, 10
)
means <- as.numeric(
  scale(
    means
  )
)

# get each city covariate. We first store this in a list
#  so that each city has it's own respective sites.
city_covar <- vector(
  "list",
  length = ncity
)
for(i in 1:ncity){
  city_covar[[i]] <- rnorm(sites[i])
}

# Parameters for the model.
#  This would be like Eq. 2 in the MS. The intercept that
#  varies by city.
b0 <- function(city_means){
  1 + 0.7 * city_means
}

# This would be like Eq. 3 in the MS. The slope term that
#  varies by city.
b1 <- function(city_means){
  -0.7 - 1*city_means
}

# This is a list object that stores our response variable for each city.
y <- vector(
  "list",
  ncity
)

# generate the response varianble for each city, making the response
#  binary in order to fit a logistic regression.
for(i in 1:ncity){
  tmp_det <- b0(means[i]) + b1(means[i]) * cl[[i]]
  y[[i]] <- rbinom(
    sites[i],
    1,
    plogis(tmp_det)
  )
}

# the data.frame that will be used to fit a logistic regression.
# y = response variable
# cm = city means = means value replicated sites times for each city
# x = site-level covariate
tmp <- data.frame(
  y = unlist(y),
  cm = rep(means, times = sites),
  x = unlist(city_covar)
)

# fit the model
example_model <- glm(y ~ cm*x, data = tmp, family = 'binomial')

ests_with_truevals <- cbind(
  c(1,0.7,-0.7, -1),
  coef(example_model),
  confint(example_model)
)

# make an object that stores the true value used to generate the data
#  the model estimate, and the 95% confidence intervals.
colnames(ests_with_truevals)[1:2] <- c(
  "true_value",
  "estimate"
  )
# inspect the object
ests_with_truevals
