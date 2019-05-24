library(PSIlence)
context("mechanism")

mechanismTest <- setRefClass(
  Class = 'mechanismTest',
  contains = 'mechanism'
)

epsilon=0.1
delta=0.4
mech = mechanismTest$new(mechanism="test", epsilon=epsilon, delta=delta)

expect_equal(mech$epsilon, epsilon)
expect_equal(mech$delta, delta)
expect_equivalent(mech$getFields()['delta'], delta)