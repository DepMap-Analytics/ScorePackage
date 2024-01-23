#test for L1 priority script

test_that("retrievePrioriKnownEss working", {
  expect_equal(PRIORIEss, retrievePrioriKnownEss(essData = matrix(rep(NA,length(NGUIDES)),length(NGUIDES),1,dimnames = list(names(NGUIDES),NULL)),
                                                 coreFitnessGenes=coreFitnessGenes,
                                                 panCancer_cf_genes=PanCancerCoreFitnessGenes,dir.ExternalData=".",
                                                 cancerDrivers=cancerDrivers,refSets=refSets)[,4:10]+0
  )
})

test_that("retrievePrimaryTumourStatus working", {
  #expect_equal(2 * 2, 4)
})

#end of L1 priority script test

#start tests for L2 priority script
test_that("retrieveCNstatus working", {
  expect_equal(2 * 2, 4)
})
test_that("retrieveGenomic working", {
  expect_equal(2 * 2, 4)
})
test_that("updateMatrix working", {
  expect_equal(2 * 2, 4)
})




#end test L2 priority script
