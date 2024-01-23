source("helper_PriorityThresholds.R")
#test for priority threshold script
test_that("Get_PriorityThreshold is working",{
  expect_equal(priority_threshold,Get_PriorityThreshold(set1,set2))
})
#do tests for saveFilesForPrioritization subfunctions
test_that("LoadBiomarkerRes working", {
  expect_equal(LoadedFile, LoadBiomarkerRes(CTYPE=ctype,ANOVA.results.dir=ANOVA.results.dir,mutationNumber="",load=FALSE,subdir="")$loadedFile)
})

test_that("GetClassIDs working", {
  expect_equal(BiomarkerRes, GetClassIDs(subdir="CExpr",TOTRESpt))
})

#do tests for priorityPlot sub functions
test_that("get_PriorityRes working", {
  expect_equal(PriorityRes, get_PriorityRes(MARKERclass,ctype,mutationNumber,TRACTABLE,PCHSYM,indications,PRIORITY_vectors=Pvec,PRIORITY_vectorsL3,th,sigOnly=TRUE))
})


test_that("combineRes is working",{
  expect_equal(testCR,combineRes(r1,r2))
})
