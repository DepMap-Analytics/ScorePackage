
test_that("zero PC removal returns original matrix",{
  expect_equal(testPC,RemovePC(testPC,droppcanumber=0,perfCheck=FALSE))
})

test_that("FCs combined for positive correlations",{
  expect_equal(posab,getPairCombFC(testFCs,c("genea","geneb"),0.5))
})

test_that("FCs combined for negative correlations",{
  expect_equal(negac,getPairCombFC(testFCs,c("genea","genec"),-0.5))
})

test_that("Pair correlations with all cell lines",{
  expect_equal(corall_bc,subcor(sBFs,2,binarymat=allBin)[3,3])
})

test_that("Pair correlations with responsive cell lines",{
  expect_equal(corpair_ab1,subcor(sBFs,1,binarymat=binaryMat)[2,3])
})

test_that("Input List for Pair correlations responsive",{
  expect_equal(testCM,prepareInputs(sBFs,BEMDep,PPIgraphT,binaryMat2,PPInet,cmp,CancerType,maxCL=3,scaledBFMax=50,removeCore=removeCore,FCs=FCs,method="pearson",allObs=FALSE,ncores=ncores)
$subCM)
})

test_that("Input List for Pair correlations all",{
  expect_equal(testCM2,prepareInputs(sBFs,BEMDep,PPIgraphT,binaryMat2,PPInet,cmp,CancerType,maxCL=3,scaledBFMax=50,removeCore=removeCore,FCs=FCs,method="pearson",allObs=FALSE,ncores=ncores)
               $corMat)
})



test_that("Test for pairs that pass correlation thresholds/connected in PPI graph",{
 expect_equal(pairList,corPairs(testCM2,simThreshold=PaircorT,PPIgraphT,PPInet=PPInet)$CorGroups[,1:2])
})

test_that("Test for pairs that pass correlation thresholds no PPI graph requirement",{
  expect_equal(pairList,corPairs(testCM2,simThreshold=PaircorT,PPIigraph=NULL,PPInet=PPInet)$CorGroups[,1:2])
})

test_that("Combining of extra output for pair dependencies, pos",{
  expect_equal(extra_sBF,combineGeneProfiles(direction=TRUE,FC=FCs[1:2,],Binary=binaryMat,BF=sBFs)$subBF)
})

test_that("Combining of extra bdep output for pair dependencies, pos",{
  expect_equal(extra_bDep,combineGeneProfiles(direction=TRUE,FC=FCs[1:2,],Binary=binaryMat,BF=sBFs)$subbinary)
})

test_that("Combining of extra bdep output for pair dependencies, neg",{
  expect_equal(extra_bDepN,combineGeneProfiles(direction=FALSE,FC=FCs[1:2,],Binary=binaryMat,BF=sBFs)$subbinary)
})

test_that("Combining of extra sBF output for pair dependencies, neg",{
  expect_equal(extra_sBFN,combineGeneProfiles(direction=FALSE,FC=FCs[1:2,],Binary=binaryMat,BF=sBFs)$subBF)
})

