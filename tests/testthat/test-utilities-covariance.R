library(PSIlence)
context('utilities-covariance')

test_that('true covariance function is correct', {
    x1 <- c(0,5,3)
    x2 <- c(1,2,3)
    data <- data.frame(x1,x2)
    
    testCovar <- covar(data, intercept=FALSE)
    x<- covarianceFormatRelease(testCovar, c("x1", "x2"))
    y <- cov(data)
    
    expect_equal(as.matrix(x), y)
})

test_that('linear regression post-processing function is correct for 2x2 covariance matrix',{
    
    columns <- c('income', 'age')
    n <- 10000
    intercept <- FALSE
    formula <- 'income~age'
    data <- PUMS5extract10000[columns]
    
    covar <- covar(data, intercept=FALSE)
    
    formattedCovar <- covarianceFormatRelease(covar, columns)
    postLnReg <- covariancePostLinearRegression(formattedCovar, n, intercept, formula)
    output <- as.numeric(postLnReg[[1]][1]) #extracts coefficient from output
    trueLinearReg <- lm(formula, data=PUMS5extract10000)
    expectedOutput <- as.numeric(trueLinearReg[[1]][2]) #extracts coefficient from output
    expect_equal(floor(output), floor(expectedOutput)) #check floor of values due to floating point errors.
})

test_that('linear regression post-processing function is correct for 3x3 covariance matrix',{
    columns <- c('income', 'age', 'educ')
    n <- 10000
    intercept <- FALSE
    formula <- 'income~educ'
    data <- PUMS5extract10000[columns]
    
    covarMatrix <- covar(data, intercept=FALSE)
    
    formattedCovar <- covarianceFormatRelease(covarMatrix, columns)
    postLnReg <- covariancePostLinearRegression(formattedCovar, n, intercept, formula)
    output <- as.numeric(postLnReg[[1]][1]) #extracts coefficient from output
    
    trueLinearReg <- lm(formula, data=PUMS5extract10000)
    expectedOutput <- as.numeric(trueLinearReg[[1]][2]) #extracts coefficient from output
    expect_equal(floor(output), floor(expectedOutput)) #check floor of values due to floating point errors.
})