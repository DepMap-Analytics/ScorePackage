
test_that("Input BEM features with NAs calculated correctly", {
  expect_equal(BEMDep2, Regression_createInputFeatures(BEMDepNA,oneFeatureOnly=NULL,
                                                      selectedCellLines=NULL,MASTER_LIST=NULL,featFactorPopulationTh=3)$EM

  )
})

test_that("Regression output working", {
  expect_equal(limmaout, AllRegression("./",FCs,BEMDepNA,inparallel=FALSE,geneAnnotations=geneAnnot,biomarkertype="CExpr")[,colnames(limmaout)]

  )
})

test_that("Matrix version of Glass D calc working",{
  expect_equal(rep(g1,2),Glass_Deltas(testD1,testD2)$g1)
})
