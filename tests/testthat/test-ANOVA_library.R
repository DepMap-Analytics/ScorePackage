source("helper_ANOVA_library.R")

test_that("One file for ADaM context core fitness returned",{
  expect_equal(1,length(getADMfile(admfiles,"Breast")))
})

test_that("Tissue coreFitness genes selected",{
  expect_equal("09_ADM_Breast_coreFitnessGenes.RData",getADMfile(admfiles,"Breast"))
})

test_that("binary depletions and Ess profiles same dimension",{
  expect_equal(dim(getBiomarkerInput(Qnorm,bDepletions,TOTALBEM=BEMDepNA)$bDep),dim(getBiomarkerInput(Qnorm,bDepletions,TOTALBEM=BEMDepNA)$EssProf))
})

test_that("binary depletions At least 3 depletions",{
  expect_equal(3,min(colSums(getBiomarkerInput(Qnorm,bDepletions,TOTALBEM=BEMDepNA)$bDep)))
})

test_that("All models remain",{
  expect_equal(ncol(Qnorm),nrow(getBiomarkerInput(Qnorm,bDepletions,TOTALBEM=BEMDepNA)$bDep))
})


test_that("ANOVA features inputs created correctly", {
  expect_equal(aovTbem, ANOVA_createInputFeatures(additional_features=ANOVA.additionalFeatures,
                                                               oneFeatureOnly = ANOVA.settings.oneFeatureOnly,
                                                               excludeHyperMetData = ANOVA.settings.excludeHyperMetData,
                                                                additional_features_only = ANOVA.settings.additionalFeaturesOnly,
                                                                ANOVA.settings.featFactorPopulationTh=2,
                                                               selectedCellLines = colnames(bDepletions),modelID="model_id",
                                                               MASTER_LIST=MASTER_LIST,TOTALBEM=BEMin)$BEM)
})

test_that("ANOVA results correct, check2",{
  expect_equal(nrow(InputFeatures$BEM)*nrow(Qnorm),nrow(aovRes))
})

test_that("ANOVA results correct, check3",{
  expect_equal(29,ncol(aovRes))
})

test_that("Compound ANOVA features inputs created correctly, Mutual Exclusive", {
  expect_equal(NULL,ANOVA_createInputFeaturesCompound(additional_features=ANOVA.additionalFeatures,
                                                 oneFeatureOnly = ANOVA.settings.oneFeatureOnly,
                                                 excludeHyperMetData = ANOVA.settings.excludeHyperMetData,
                                                 additional_features_only = ANOVA.settings.additionalFeaturesOnly,TOTALBEM=BEMDepC,
                                                 selectedCellLines = paste0('cl',1:20),
                                                 MASTER_LIST=MASTER_LIST,ANOVA.settings.featFactorPopulationTh=3,
                                                 ParentBEM=BEMDepNA,MutualExclusive=TRUE)$CompoundTerms)

})

test_that("Compound ANOVA features inputs created correctly, Compound", {
  #expect_true(t1)

})

test_that("Error checks for ANOVA input BEM working",{
  expect_error(ANOVA_totalANOVA(fn="",ESSprofiles=NULL,InputFeatures=InputFeaturesE,inparallel=FALSE,geneAnnotations,ANOVA.settings.analysisType="CS",ANOVA.fakeCS=FALSE,
                                ANOVA.settings.includeMSI_Factor=FALSE)
  )
})


