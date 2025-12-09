library(LocalPareto)
library(Pareto)

context("test functions LALocPareto")

test_that("LALocPareto_Layer_Mean", {
  expect_equal(round(LALocPareto_Layer_Mean(8000, 2000, 1000, 1.5, gamma = 1), 4), 230.7318)
  expect_equal(round(LALocPareto_Layer_Mean(8000, 2000, 1000, 1.5, delta = 1), 4), 239.2577)
})

test_that("LALocPareto_Layer_SM", {
  expect_equal(round(LALocPareto_Layer_SM(8000, 2000, 1000, 1.5, gamma = 1), 0), 508709)
  expect_equal(round(LALocPareto_Layer_SM(8000, 2000, 1000, 1.5, delta = 1), 0), 541356)
})

test_that("LALocPareto_Layer_Var", {
  expect_equal(round(LALocPareto_Layer_Var(8000, 2000, 1000, 1.5, gamma = 1), 0), 455472)
  expect_equal(round(LALocPareto_Layer_Var(8000, 2000, 1000, 1.5, delta = 1), 0), 484111)
})

test_that("LALocPareto_ML_Estimator_Alpha", {
  losses <- c(1300, 1800, 3500, 2480, 1738, 2000, 1400, 1200, 4700)

  MLE <- LALocPareto_ML_Estimator(losses, 1000)
  expect_equal(round(MLE$alpha_0_hat, 7), 0.1916228)
  expect_equal(round(MLE$gamma_hat, 5), 13.24756)
  expect_equal(round(MLE$delta_hat, 6), 1.759579)
  expect_equal(round(MLE$ApproxCovMatrix[2,1], 7), -40.3779335)
  expect_equal(round(MLE$ApproxCovMatrix2[2,1], 7), -0.6845955)

  MLE <- LALocPareto_ML_Estimator(losses, 1000, gamma = 0.5)
  expect_equal(round(MLE$alpha_0_hat, 6), 1.145414)
  expect_equal(round(MLE$gamma_hat, 5), 0.5)
  expect_equal(round(MLE$delta_hat, 7), 0.3969701)
  expect_equal(round(MLE$ApproxCovMatrix[1,1], 7), 0.1457747)
  expect_equal(round(MLE$ApproxCovMatrix2[1,2], 8), 0.05052166)

  MLE <- LALocPareto_ML_Estimator(losses, 1000, delta = 0.5)
  expect_equal(round(MLE$alpha_0_hat, 7), 0.9760963)
  expect_equal(round(MLE$gamma_hat, 7), 0.7390127)
  expect_equal(round(MLE$delta_hat, 7), 0.5)
  expect_equal(round(MLE$ApproxCovMatrix[2,1], 7), -0.1642318)
  expect_equal(round(MLE$ApproxCovMatrix2[1,1], 7), 0.2169192)

  t <- c(1000, 500, 500, 1500, 800, 900, 500, 1000, 1000)
  MLE <- LALocPareto_ML_Estimator(losses, t)
  expect_equal(round(MLE$alpha_0_hat, 7), 0.2386455)
  expect_equal(round(MLE$gamma_hat, 6), 5.609853)
  expect_equal(round(MLE$delta_hat, 7), 0.9279622)
  expect_equal(round(MLE$ApproxCovMatrix[2,1], 7), -6.4981638)
  expect_equal(round(MLE$ApproxCovMatrix2[2,1], 7), -0.2127012)

})

test_that("pLALocPareto", {
  expect_equal(round(pLALocPareto(2000, 1000, 2, 0.5), 5), 0.80339)
  expect_equal(round(pLALocPareto(2000, 1000, 2, delta = 0.5), 5), 0.78978)
})

test_that("dLALocPareto", {
  expect_equal(round(dLALocPareto(2000, 1000, 2, 0.5), 8), 0.00026475)
  expect_equal(round(dLALocPareto(2000, 1000, 2, delta = 0.5), 8), 0.00026278)
})

test_that("qLALocPareto", {
  expect_equal(round(pLALocPareto(qLALocPareto((0:10)*0.1, 1000, 2, 1), 1000, 2, 1), 10), (0:10)*0.1)
  expect_equal(round(qLALocPareto(pLALocPareto((1:10)*1000, 1000, 2, 1), 1000, 2, 1), 10), (1:10)*1000)
})

test_that("HelperPPareto", {
  expect_equal(round(P_Layer_Mean_s(1000, 1000, 2, 1000), 6), round(Pareto_Layer_Mean(1000, 1000, 2, 1000), 6))
  expect_equal(round(P_Layer_Mean_s(500, 0, 2, 1000), 6), round(Pareto_Layer_Mean(500, 0, 2, 1000), 6))
  expect_equal(round(P_Layer_Mean_s(Inf, 2000, 2.5, 1000), 6), round(Pareto_Layer_Mean(Inf, 2000, 2.5, 1000), 6))
  expect_equal(round(P_Layer_Mean_s(Inf, 0, 2.5, 1000), 6), round(Pareto_Layer_Mean(Inf, 0, 2.5, 1000), 6))
  expect_equal(round(P_Layer_Mean_s(1000, 500, 1, 1000), 6), round(Pareto_Layer_Mean(1000, 500, 1, 1000), 6))

  expect_equal(round(P_Layer_Var_s(1000, 1000, 2, 1000), 6), round(Pareto_Layer_Var(1000, 1000, 2, 1000), 6))
  expect_equal(round(P_Layer_Var_s(500, 0, 2, 1000), 6), round(Pareto_Layer_Var(500, 0, 2, 1000), 6))
  expect_equal(round(P_Layer_Var_s(Inf, 2000, 2.5, 1000), 6), round(Pareto_Layer_Var(Inf, 2000, 2.5, 1000), 6))
  expect_equal(round(P_Layer_Var_s(Inf, 0, 2.5, 1000), 6), round(Pareto_Layer_Var(Inf, 0, 2.5, 1000), 6))
  expect_equal(round(P_Layer_Var_s(1000, 500, 1, 1000), 6), round(Pareto_Layer_Var(1000, 500, 1, 1000), 6))

  expect_equal(round(P_Layer_SM_s(1000, 1000, 2, 1000), 6), round(Pareto_Layer_SM(1000, 1000, 2, 1000), 6))
  expect_equal(round(P_Layer_SM_s(500, 0, 2, 1000), 6), round(Pareto_Layer_SM(500, 0, 2, 1000), 6))
  expect_equal(round(P_Layer_SM_s(Inf, 2000, 2.5, 1000), 6), round(Pareto_Layer_SM(Inf, 2000, 2.5, 1000), 6))
  expect_equal(round(P_Layer_SM_s(Inf, 0, 2.5, 1000), 6), round(Pareto_Layer_SM(Inf, 0, 2.5, 1000), 6))
  expect_equal(round(P_Layer_SM_s(1000, 500, 1, 1000), 6), round(Pareto_Layer_SM(1000, 500, 1, 1000), 6))

  Parameters_PP <- list(t = 1000, alpha = 2)
  expect_equal(round(PP_Layer_Mean(1000, 500, Parameters_PP), 6), round(PiecewisePareto_Layer_Mean(1000, 500, 1000, 2), 6))
  expect_equal(round(PP_Layer_Mean(1000, 2000, Parameters_PP), 6), round(PiecewisePareto_Layer_Mean(1000, 2000, 1000, 2), 6))
  expect_equal(round(PP_Layer_Mean(Inf, 0, Parameters_PP), 6), round(PiecewisePareto_Layer_Mean(Inf, 0, 1000, 2), 6))

  t <- (1:5) * 1000
  alpha <- 1:5
  SF <- 1 - pPiecewisePareto(t, t, alpha)
  Parameters_PP <- list(t = t, alpha = alpha, SF = SF)
  Cover <- c(500, 1500, 2000, Inf)
  AttPt <- c(0, 500, 2000, 4000)
  expect_equal(round(PP_Layer_Mean(Cover, AttPt, Parameters_PP), 6), round(PiecewisePareto_Layer_Mean(Cover, AttPt, t, alpha), 6))
  expect_equal(round(PP_Layer_SM(Cover, AttPt, Parameters_PP), 6), round(PiecewisePareto_Layer_SM(Cover, AttPt, t, alpha), 6))
  expect_equal(round(PP_Layer_Var(Cover, AttPt, Parameters_PP), 6), round(PiecewisePareto_Layer_Var(Cover, AttPt, t, alpha), 6))

  p <- c(0.01, (1:9) * 0.1, 0.99)
  expect_equal(round(qPP(p, Parameters_PP), 6), round(qPiecewisePareto(p, t, alpha), 6))
  x <- qPP(p, Parameters_PP)
  expect_equal(round(pPP(x, Parameters_PP), 5), round(p, 5))
  expect_equal(round(qPP(p, Parameters_PP), 5), round(x, 5))
  expect_equal(round(dPP(x, Parameters_PP), 6), round(dPiecewisePareto(x, t, alpha), 6))



})

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

