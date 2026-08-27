.duplicate_knot_cdf <- function() {
  list(
    knots = c(-2, -1, -1, 0, 0.75, 2, 2),
    values = c(0.05, 0.20, 0.18, 0.50, 0.80, 0.95, 0.90)
  )
}

test_that("step-CDF preparation resolves duplicate knots monotonically", {
  prepared <- .prepare_step_cdf(.duplicate_knot_cdf())

  expect_equal(unname(prepared$knots), c(-2, -1, 0, 0.75, 2))
  expect_equal(unname(prepared$values), c(0.05, 0.20, 0.50, 0.80, 0.95))
  expect_true(all(diff(prepared$knots) > 0))
  expect_true(all(diff(prepared$values) >= 0))
})

test_that("direct I-spline regularization is a monotone CDF with valid tails", {
  .skip_semiparametric_stack()
  smooth <- .smooth_Fhat(.duplicate_knot_cdf(), n0 = 125L, K = 8L)
  grid <- seq(-3, 3, length.out = 401L)
  probability <- smooth$Ftilde(grid)
  density <- smooth$ftilde(grid)

  expect_identical(smooth$method, "ispline")
  expect_true(all(is.finite(probability)))
  expect_true(all(diff(probability) >= -1e-10))
  expect_true(all(probability >= 0 & probability <= 1))
  expect_equal(smooth$Ftilde(c(-Inf, -3)), c(0, 0))
  expect_equal(smooth$Ftilde(c(3, Inf)), c(1, 1))
  expect_true(all(is.finite(density)))
  expect_true(all(density >= -1e-12))
  expect_equal(smooth$ftilde(c(-Inf, -3, 3, Inf)), rep(0, 4L))
  expect_true(all(smooth$knot_seq > smooth$u_range[1L]))
  expect_true(all(smooth$knot_seq < smooth$u_range[2L]))
  expect_identical(anyDuplicated(smooth$knot_seq), 0L)

  interior <- seq(
    smooth$u_range[1L] + 0.1,
    smooth$u_range[2L] - 0.1,
    length.out = 25L
  )
  h <- 1e-6
  numerical_density <-
    (smooth$Ftilde(interior + h) - smooth$Ftilde(interior - h)) / (2 * h)
  expect_equal(smooth$ftilde(interior), numerical_density, tolerance = 1e-5)
})

test_that("ridge regularization uses epsilon equals one over unit count", {
  Fhat <- .duplicate_knot_cdf()
  ridge <- .smooth_Fhat_ridge(Fhat, regularization_n = 100L)
  grid <- seq(-3, 3, length.out = 401L)

  expect_identical(ridge$method, "ridge")
  expect_null(ridge$ftilde)
  expect_equal(ridge$eps, 0.01)
  expect_equal(ridge$Ftilde(-Inf), ridge$eps)
  expect_equal(ridge$Ftilde(Inf), 1)
  expect_true(all(diff(ridge$Ftilde(grid)) >= -1e-12))
  expect_true(all(ridge$Ftilde(grid) >= ridge$eps))
  expect_true(all(ridge$Ftilde(grid) <= 1))

  prepared <- .prepare_step_cdf(Fhat)
  u <- 0.75
  step_at_u <- stats::approxfun(
    prepared$knots,
    prepared$values,
    method = "constant",
    rule = 2,
    f = 0,
    ties = max
  )(u)
  ramp_at_u <- (u - min(prepared$knots)) / diff(range(prepared$knots))
  expected <- (1 - 2 * ridge$eps) * step_at_u +
    ridge$eps + ridge$eps * ramp_at_u
  expect_equal(ridge$Ftilde(u), expected)
})

test_that("regularization dispatcher supports both simulation variants", {
  .skip_semiparametric_stack()
  Fhat <- .duplicate_knot_cdf()

  ispline <- .regularize_Fhat(Fhat, n0 = 125L, method = "ispline")
  ridge <- .regularize_Fhat(
    Fhat,
    n0 = 125L,
    method = "ridge",
    regularization_n = 250L
  )

  expect_identical(ispline$method, "ispline")
  expect_identical(ridge$method, "ridge")
  expect_equal(ridge$eps, 1 / 250)
})
