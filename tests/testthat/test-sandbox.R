context('testing, do not merge')

test_that('Laplace noise generation works when scale and number of noise parameters have varying dimensions',{
    # Case 1: generate a single noisy value from a Laplace distributions
    expect_equal(qLap(p=0.4, mu=0, b=3), 3*log(2 * 0.4))
    
    # Case 2: generate multiple noisy values, with vector of probabilities and single scaling parameter
    p <- c(0.4, 0.7)
    scale <- 3
    noise <- qLap(p, mu=0, b=scale)
    
    expect_equal(noise[1], scale*log(2*p[1]))
    expect_equal(noise[2], - scale*log(2-2 * p[2]))
    
    # Case 3: generate multiple noisy values, with vector of probabilities and vector of scaling parameters
    scale <- c(1, 5)
    noise <- qLap(p, mu=0, b=scale)
    
    expect_equal(noise[1], scale[1]*log(2*p[1]))
    expect_equal(noise[2], - scale[2]*log(2-2*p[2]))
})