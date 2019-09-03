context('utilities-input-validation')

test_that('checkN raises proper warnings and errors',{
  expect_error(checkN(-1))
  expect_error(checkN(1.5))
  expect_error(checkN(-1.5))
  expect_error(checkN(c(1,2)))
  expect_error(checkN('foo'))
})