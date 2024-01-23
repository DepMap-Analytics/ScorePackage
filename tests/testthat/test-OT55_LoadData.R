test_that("Remove core genes works", {
  expect_equal(removeCdata, removeCoreFitness(Qnorm,removeList))
})

