library(LocalPareto)
library(Pareto)

context("test functions LocalPareto")


test_that("LocalPareto_Layer_Mean", {
  alpha <- function(x) ifelse(x >= 1000, 2 - 1000 / x, 0)
  expect_equal(round(LocalPareto_Layer_Mean(8000, 2000, alpha), 1), 810.9)
  expect_equal(round(LocalPareto_Layer_Mean(8000, 2000, alpha, 2000),0), 1967)
})

test_that("LocalPareto_Layer_Var", {
  alpha <- function(x) ifelse(x >= 1000, 2 - 1000 / x, 0)
  expect_equal(round(LocalPareto_Layer_Var(8000, 2000, alpha), 0), 2966127)
  expect_equal(round(LocalPareto_Layer_Var(8000, 2000, alpha, 2000),0), 4921195)
})

test_that("LocalPareto_Layer_SM", {
  alpha <- function(x) ifelse(x >= 1000, 2 - 1000 / x, 0)
  expect_equal(round(LocalPareto_Layer_SM(8000, 2000, alpha), 0), 3623649)
  expect_equal(round(LocalPareto_Layer_SM(8000, 2000, alpha, 2000),0), 8791416)
})

test_that("pLocalPareto", {
  alpha <- function(x) ifelse(x >= 1000, 2 - 1000 / x, 0)
  expect_equal(round(pLocalPareto((1:5) * 1000, alpha), 4), c(0.0000, 0.5878, 0.7836, 0.8677, 0.9110))
  expect_equal(round(pLocalPareto((1:5) * 1000, alpha, 2000), 4), c(0.000, 0.000, 0.475, 0.679, 0.784))
})

test_that("dLocalPareto", {
  alpha <- function(x) ifelse(x >= 1000, 2 - 1000 / x, 0)
  expect_equal(round(dLocalPareto((1:12) * 500, alpha), 7), c(0, 0.0010038, 0.0005522, 0.0003094, 0.0001865, 0.0001204, 0.0000816, 0.0000578, 0.0000425, 0.0000320, 0.0000248, 0.0000195))
  expect_equal(round(dLocalPareto((4:12) * 500, alpha, 2000), 7), c(0.0007514, 0.0004528, 0.0002913, 0.0001980, 0.0001404, 0.0001030, 0.0000777, 0.0000601, 0.0000474))
})

test_that("qLocalPareto", {
  alpha <- function(x) ifelse(x >= 1000, 2 - 1000 / x, 0)
  expect_equal(round(qLocalPareto((0:5) * .2, alpha), 3), c(1000.000, 1225.962, 1537.607, 2040.263, 3144.694, Inf))
  expect_equal(round(qLocalPareto((0:5) * .2, alpha, 2000), 3), c(2000.000, 2313.002, 2767.314, 3523.213, 5217.592, Inf))
})



