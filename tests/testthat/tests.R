test_that("bernstein_polynomial() works",{
  expect_equal(bernstein_polynomial(0.5, 0, 3), 0.125)
})

test_that("Correct average weights: olsslope", {

  dgp <- dgp(MST2018 = TRUE)

  # From Figure 3
  expect_equal(compute_average_weights_ivlike(olsslope(dgp), dgp),
               readRDS(test_path("compute_average_weights_olsslope_test.rds")))

})

test_that("Correct average weights: ivslope", {

  dgp <- dgp(MST2018 = TRUE)

  # From Figure 1
  expect_equal(compute_average_weights_ivlike(ivslope(dgp), dgp),
               readRDS(test_path("compute_average_weights_ivslope_test.rds")))

})

test_that("Correct average weights: late", {

  dgp <- dgp(MST2018 = TRUE)

  # From Figure 1
  expect_equal(compute_average_weights(tp = "LATE",
                                       dgp = dgp,
                                       late.lb = 0.35,
                                       late.ub = 0.9),
               readRDS(test_path("compute_average_weights_late_test.rds")))

})


test_that("Correct bounds for MST2018: Figure 2",{

  dgp = dgp(MST2018 = TRUE)

  bounds2 = compute_bounds(tp = late(dgp = dgp, u1 = 0.35, u2 = 0.9),
                           bases = list(constantspline_basis(c(0, 1, 0.35, 0.9, dgp$pscoreZ)),
                                        constantspline_basis(c(0, 1, 0.35, 0.9, dgp$pscoreZ))),
                           dgp = dgp,
                           ivslope = TRUE)

  expect_equal(bounds2$upper_bound, 0.5003247, tolerance = 1e-6)
  expect_equal(bounds2$lower_bound, -0.4208874, tolerance = 1e-6)
})

test_that("Correct bounds for MST2018: Figure 3",{

  dgp = dgp(MST2018 = TRUE)

  bounds3 = compute_bounds(tp = late(dgp = dgp, u1 = 0.35, u2 = 0.9),
                           bases = list(constantspline_basis(c(0, 1, 0.35, 0.9, dgp$pscoreZ)),
                                        constantspline_basis(c(0, 1, 0.35, 0.9, dgp$pscoreZ))),
                           dgp = dgp,
                           ivslope = TRUE,
                           olsslope = TRUE)

  expect_equal(bounds3$upper_bound, 0.5003247, tolerance = 1e-6)
  expect_equal(bounds3$lower_bound, -0.4111812, tolerance = 1e-6)
})

test_that("Correct bounds for MST2018: Figure 4",{

  dgp = dgp(MST2018 = TRUE)

  bounds4 = compute_bounds(tp = late(dgp = dgp, u1 = 0.35, u2 = 0.9),
                           bases = list(constantspline_basis(c(0, 1, 0.35, 0.9, dgp$pscoreZ)),
                                        constantspline_basis(c(0, 1, 0.35, 0.9, dgp$pscoreZ))),
                           dgp = dgp,
                           ivslopeind = TRUE,
                           params = list(ivslopeind = c(1,2)))

  expect_equal(bounds4$upper_bound, 0.4074924, tolerance = 1e-6)
  expect_equal(bounds4$lower_bound, -0.3197803, tolerance = 1e-6)
})

test_that("Correct bounds for MST2018: Figure 5",{

  dgp = dgp(MST2018 = TRUE)

  bounds5 = compute_bounds(tp = late(dgp = dgp, u1 = 0.35, u2 = 0.9),
                           bases = list(constantspline_basis(c(0, 1, 0.35, 0.9, dgp$pscoreZ)),
                                        constantspline_basis(c(0, 1, 0.35, 0.9, dgp$pscoreZ))),
                           dgp = dgp,
                           saturated = TRUE)

  expect_equal(bounds5$upper_bound, 0.4074924, tolerance = 1e-6)
  expect_equal(bounds5$lower_bound, -0.1377803, tolerance = 1e-6)
})

test_that("Correct bounds for MST2018: Figure 6",{

  dgp = dgp(MST2018 = TRUE)

  bounds6 = compute_bounds(tp = late(dgp = dgp, u1 = 0.35, u2 = 0.9),
                           bases = list(constantspline_basis(c(0, 1, 0.35, 0.9, dgp$pscoreZ)),
                                        constantspline_basis(c(0, 1, 0.35, 0.9, dgp$pscoreZ))),
                           dgp = dgp,
                           saturated = TRUE,
                           decreasing.MTR = TRUE)

  expect_equal(bounds6$upper_bound, 0.07731061, tolerance = 1e-6)
  expect_equal(bounds6$lower_bound, -0.09517424, tolerance = 1e-6)
})

test_that("Correct bounds for MST2018: Figure 7",{

  dgp = dgp(MST2018 = TRUE)

  bounds7 = compute_bounds(tp = late(dgp = dgp, u1 = 0.35, u2 = 0.9),
                           bases = list(bernstein_basis(9),
                                        bernstein_basis(9)),
                           dgp = dgp,
                           saturated = TRUE,
                           decreasing.MTR = TRUE)

  expect_equal(bounds7$upper_bound, 0.06660755, tolerance = 1e-6)
  expect_equal(bounds7$lower_bound, 0.0004047045, tolerance = 1e-6)
})
