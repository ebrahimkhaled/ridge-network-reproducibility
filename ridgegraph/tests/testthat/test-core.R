test_that("ridge precision inverts S + lambda I", {
  set.seed(2)
  x <- matrix(rnorm(30 * 5), 30, 5)
  xc <- scale(x, TRUE, FALSE)
  s <- crossprod(xc) / 29
  q <- ridge_precision(x, 0.3)
  expect_equal(q %*% (s + 0.3 * diag(5)), diag(5), tolerance = 1e-10)
})

test_that("effective penalty is increasing and exceeds lambda", {
  set.seed(3)
  x <- matrix(rnorm(20 * 40), 20, 40)
  ep <- effective_penalty(x, c(0.01, 0.1, 1, 10))
  expect_true(all(diff(ep$kappa) > 0))
  expect_true(all(ep$kappa >= ep$lambda))
  fl <- effective_penalty_floor(x)
  expect_true(fl > 0 && fl < ep$kappa[1])
  l <- penalty_for_kappa(x, 2 * fl)
  expect_equal(effective_penalty(x, l)$kappa, 2 * fl, tolerance = 1e-6)
})

test_that("population extrapolation identity holds", {
  set.seed(4)
  a <- matrix(rnorm(36), 6, 6)
  sig <- crossprod(a) + diag(6)
  th <- solve(sig)
  k <- c(0.4, 1.1)
  w <- c(k[2] / (k[2] - k[1]), -k[1] / (k[2] - k[1]))
  lhs <- w[1] * solve(sig + k[1] * diag(6)) + w[2] * solve(sig + k[2] * diag(6))
  e <- eigen(th, symmetric = TRUE)
  t <- e$values
  rhs <- e$vectors %*% diag(t - prod(k) * t^3 / ((1 + k[1] * t) * (1 + k[2] * t))) %*% t(e$vectors)
  expect_equal(lhs, rhs, tolerance = 1e-10)
})

test_that("permutation selection returns a symmetric adjacency", {
  set.seed(5)
  x <- matrix(rnorm(40 * 6), 40, 6)
  sel <- mps_select(x, B = 39)
  expect_true(isSymmetric(sel$adjacency))
  expect_false(any(diag(sel$adjacency)))
  expect_length(sel$threshold, 4)
})

test_that("estimated effective penalty matches the trace formula", {
  set.seed(6)
  x <- matrix(rnorm(25 * 40), 25, 40)
  xc <- scale(x, TRUE, FALSE)
  s <- crossprod(xc) / 24
  r <- sum(diag(s %*% solve(s + 0.7 * diag(40)))) / 24
  expect_equal(effective_penalty(x, 0.7)$kappa, 0.7 / (1 - r), tolerance = 1e-10)
  ev <- eigen(s, TRUE, only.values = TRUE)$values[1:24]
  expect_equal(effective_penalty_floor(x), 24 / sum(1 / ev), tolerance = 1e-10)
})

test_that("extrapolation doubles the effective penalty", {
  set.seed(7)
  x <- matrix(rnorm(30 * 50), 30, 50)
  fit <- extrapolated_ridge(x, 0.4)
  expect_equal(fit$kappa[2], 2 * fit$kappa[1], tolerance = 1e-8)
  expect_equal(sum(fit$weights), 1, tolerance = 1e-12)
})

test_that("mps_select warns when no edge can be selected", {
  set.seed(8)
  x <- matrix(rnorm(30 * 5), 30, 5)
  expect_warning(mps_select(x, B = 19), "no edge can be selected")
})
