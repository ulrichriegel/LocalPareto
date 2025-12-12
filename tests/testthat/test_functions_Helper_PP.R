library(LocalPareto)
library(Pareto)

context("test functions Helper Functions PP")



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



