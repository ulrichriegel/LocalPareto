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




