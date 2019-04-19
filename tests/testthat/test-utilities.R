library(PSIlence)
context("utilities")

############################# fillMissing tests #######################################

x <- dpUnif(10)
animals <- as.factor(c('dog', 'zebra', 'bird', 'hippo'))
lowBound <- 1
upBound <- 5

scaledNumeric <- scaleValues(x, 'numeric', lower=lowBound, upper=upBound)
scaledLogical <- scaleValues(x, 'logical')
scaledInteger <- scaleValues(x, 'integer', lower=lowBound, upper=upBound)
scaledFactor <- scaleValues(x, 'factor', categories=animals)

expect_equal(length(scaledNumeric), length(x))
expect_equal(length(scaledLogical), length(x))
expect_equal(length(scaledInteger), length(x))
expect_equal(length(scaledFactor), length(x))

expect_equal(sum(scaledNumeric < lowBound), 0)
expect_equal(sum(scaledNumeric > upBound), 0)
expect_equal(sum(scaledInteger < lowBound), 0)
expect_equal(sum(scaledInteger > upBound), 0)
expect_equal(sum(scaledLogical < 0), 0)
expect_equal(sum(scaledLogical > 1), 0)
expect_equal(sum((scaledLogical > 0) & (scaledLogical < 1)), 0)

y <- rnorm(100)
y[sample(1:100, size=10)] <- NA
y_imputed <- fillMissing(x=y, var.type='numeric', impute.rng=c(-1,1))

expect_equal(sum(is.na(y_imputed)), 0)

s <- sample(animals, size=100, replace=TRUE)
s[sample(1:100, size=10)] <- NA
s_imputed <- fillMissing(x=s, var.type='factor', categories=animals)

expect_equal(sum(is.na(s_imputed)), 0)
expect_true(is.factor(s_imputed))

N <- 100
x1 <- x2 <- rnorm(N)
x1[sample(1:N, size=10)] <- NA
x2[sample(1:N, size=10)] <- NA
imp.rng <- matrix(c(-3, 3, -2, 2), ncol=2, byrow=TRUE)
df <- data.frame(x1, x2)
df_imputed <- fillMissing(x=df, var.type='numeric', impute.rng=imp.rng)

expect_equal(sum(is.na(df_imputed)), 0)
  
