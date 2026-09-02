test_that("simplex metrics attain the vertex-swap bounds", {
  p <- c(type_a = 1, type_b = 0)
  p_hat <- c(type_a = 0, type_b = 1)
  expect_equal(eval_TV(p, p_hat, trim_shared_zeros = FALSE), 1)
  expect_equal(eval_RMSE(p, p_hat, trim_shared_zeros = FALSE), 1)
  expect_equal(eval_MaxAE(p, p_hat, trim_shared_zeros = FALSE), 1)
  expect_equal(
    eval_angular_distance(p, p_hat, trim_shared_zeros = FALSE),
    1
  )
  expect_equal(eval_SDID(p, p_hat, trim_shared_zeros = FALSE), 1)
  expect_equal(eval_MAE(p, p_hat, trim_shared_zeros = FALSE), 1)
  expect_equal(eval_TV(p, p, trim_shared_zeros = FALSE), 0)
})

test_that("presence F1 is reduced when mass is spilled onto an absent type", {
  p <- c(type_a = 0.5, type_b = 0.5, type_c = 0)
  p_hat <- c(type_a = 0.4, type_b = 0.4, type_c = 0.2)
  expect_equal(
    eval_presence_F1(p, p_hat, trim_shared_zeros = FALSE),
    0.8,
    tolerance = 1e-12
  )
})
