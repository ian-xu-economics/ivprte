test_that("bernstein_polynomial() works",{
  expect_equal(bernstein_polynomial(0.5, 0, 3), 0.125)
})
