test_that("model name mapping works", {
  expect_equal(colnames(testFCs), colnames(CLnameMapping(testFCs,testFCs2,annotation = annot)$refdata))
})


test_that("Identical patterns of MoBEMs are compressed",{
  expect_equal(outputbinary,my.compress_identical_patterns(inputbinary))
})
