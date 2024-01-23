source("helper_PriorityFunctions.R")
###Starting tests for the L1_L2 combining script (Charpriority subfunctions)####

test_that("Getting L1 filename works", {
  expect_equal('TestDir/33_L1/Breast.Carcinoma_L0_L11.csv', Get_L1filename(subdir="","1",ctype="Breast.Carcinoma",omictype=NULL,dir.Results="TestDir"))
})

test_that("Getting L1 filename works, PariBD", {
  expect_equal('TestDir/33_L1/PairBD/CExpr/Breast.Carcinoma/Breast.Carcinoma_L0_L1CExpr1.csv', Get_L1filename(subdir="PairBD","1",ctype="Breast.Carcinoma",omictype="CExpr",dir.Results="TestDir"))
})

test_that("Get L1 score working",{
  expect_equal(c(200/3,100),Get_L1score(L1,componentNames=c("ANOVAfdrA","ANOVAfdrB","MutPrimTum"),componentWeights=rep(100/3,3)))
})

test_that("Error when L1 weights not sum to 100",{
  expect_error(Get_L1score(L1,componentWeights=c(20,30,10,10,5)))
})

test_that("Error when L1 weights not same length as components",{
  expect_error(Get_L1score(L1,componentWeights=c(50,50)))
})

test_that("Calculating L2 score working",{
  expect_equal(round(L2score,digits=6),round(Get_L2score(L2,"Combined"),digits=6))
})

test_that("Getting TRACTinfos working",{
  expect_equal(TRACTinfo,Get_TRACTinfo(subdir="",FINAL_priority,tractability_both)$TRACTinfos)

})

test_that("Getting TRACTinfos working with Pair data",{
  expect_equal(TRACTinfoPairs,Get_TRACTinfo(subdir="PairBD",FINAL_priorityPairs,tractability_both)$TRACTinfos)

})


test_that("Getting compounds and indications working",{
  expect_equal(CompoundIndications,Get_CompoundIndications(TRACTinfo,glist,tract_sm,tract_ab))
})

test_that("Getting Tractability buckets working",{
  expect_equal(TRACTABLE[4][[1]],suppressWarnings(Get_TractableBuckets("",tractability_both,glist,L1_NNL2,L1_NNL2L3)$TRACTABLE[4][[1]]))
})

test_that("Get warning when target not in tractability annotation",{
  expect_warning(Get_TractableBuckets("",tractability_both,c(glist,"missing"),L1_NNL2,L1_NNL2L3))
})
#test_that("Getting EFOI indications working",{
#  Get_EFOIinfo(EFOI_list,CTYPE,tractability_both,tract_ab,tract_sm)
#})


##tests for CharPriority subfunctions complete - end of testing for L1_L2 script #####


###Start tests for L3 script


test_that("Getting BM PPI results okay",{
  expect_equal(TargetOut,getBMPPIdata(BMPPIdir=NULL,CTYPE,ScoreType="Avg",GENES,biomarkertype="CN",omictype="",TargetData=T2))
})
#not using Get_BMPPIScore function at the moment. Also needs debugging if was to use.
#test_that("Getting BM PPIscore for L3 working",{
 # expect_equal(list("MYH9"=50,"MEF2D"=25),Get_BMPPIScores(L3list,GENES,subdirlist,L3names=L3names)[,"Genomic"])
#})

##End tests for L3 script###



