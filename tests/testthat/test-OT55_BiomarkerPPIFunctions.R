test_that("Pairs of dependencies split and", {
  expect_equal(splitres, splitDependencies(res)$out)
})
test_that("Pairs of dependencies split or", {
  expect_equal(splitres2, splitDependencies(res2)$out)
})


test_that("Combination biomarkers are split and", {
  expect_equal(splitresM, splitMultiFeature(resM,type="CompoundME")$out)
})

test_that("Combination biomarkers are split or", {
  expect_equal(splitresM2, splitMultiFeature(resM2,type="ME")$out)
})

test_that("RWR is correct",{
  expect_equal(rwrdr$p.inf,rwr(p0,normW,0.4,maxiter=10000))
})

test_that("Matrix is normalised correctly",{
  expect_equal(rwrdr$transition.matrix,normW)
})
test_that("Matrix is normalised correctly colSums =1",{
  expect_equal(csums,colSums(normW))
})

test_that("RWR is correct, Ring network",{
  expect_equal(rwrdr2$p.inf,rwr(p02,normW2,0.4,maxiter=10000))
})

test_that("Matrix is normalised correctly, Ring network",{
  expect_equal(rwrdr2$transition.matrix,normW2)
})
test_that("Matrix is normalised correctly colSums =1, Ring network",{
  expect_equal(csums2,colSums(normW2))
})

#getRestartProb has seed set
test_that("restart probability working",{
  expect_equal(0.35,summaryRestartProbs(rwrres)[[3]])
})

test_that("starting vectors output correctly",{
  expect_equal(startvec,RWRinputAll$startingvec)
})


test_that("normalisation vec correct from weight mat",{
  expect_equal(dlist,RWRinputAll$degreelist)
})

test_that("marker info output correctly",{
  expect_equal(MarkerInfo,RWRinputAll$AllM)
})
test_that("biomarker results annotated correctly",{
  expect_equal(test2,annotateAovPPI(TOTRES,PPInet,filename=NULL,getPPIid=FALSE,pvalCol="FEATURE_ANOVA_pval"))
})
BMS1<-ScoreRWRbiomarker(RWRthresh,MarkerInfo,PPInet=PPInetTest,pval=0.1,RWRthreshScore=RWRscoreMat,PvalthreshScore=PvalscoreMat,FDRscoreMat=FDRscoreMat,pairList=NULL,pair=FALSE,DualBM=FALSE)

test_that("ScoreRWRbiomarker is working, FDR matrix",{
  expect_equal(AnovaFdr,BMS1$AnovaFdr$mut)
})

test_that("ScoreRWRbiomarker is working, Biomarker Scores",{
  expect_equal(BiomarkerScore,BMS1$BiomarkerScore)
})
test_that("ScoreRWRbiomarker is working, Target Scores",{
  expect_equal(TargetScore,BMS1$TargetScore[[1]][1:2,])
})

test_that("Annotating biomarkers is working",{
  expect_equal(BMPriority,AnnotateBiomarkers(BMS1))
})
BMrw<-BMS1$RankWeight[[1]]
test_that("Getting individual Biomarker scores is working",{
  expect_equal(BMSindivid,getIndividScores(BMrw))
})


#start tests for the OT55_BMPriority script that caclulates final scores across cancer types


test_that("Getting weights by STRING score working",{
  #Inputs$StringScore<-apply(Inputs,1,function(x) getWeightedScore(x,PPIigraphWeighted = PPIigraphWeighted,PPIgraph = PPIigraph))

})


