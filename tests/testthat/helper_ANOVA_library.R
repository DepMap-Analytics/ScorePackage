set.seed(234)
geneAnnot<-data.frame(cbind(name=paste0("gene",letters[1:20]),gene_family=paste0("group",1:20),location_sortable=NA,
                            entrez_id=NA,
                            ensembl_gene_id=NA,
                            pubmed_id=NA,string_id=paste0("SID",letters[1:20])),stringsAsFactors = FALSE)

#ANOVA.essentiality.dir<-paste0(args[3],'/_R_objs/')
ANOVA.additionalFeatures.fn<-NULL
ANOVA.additionalFeatures<-NULL
#ANOVA.BAGELth.fn<-paste0(args[6],'/_R_objs/06_BAGEL_BEST_PPV_th_TRAINING.Rdata')
ANOVA.settings.SCALEtoLIMIT<-FALSE #used by BAGELr only
ANOVA.settings.quantileNorm<-FALSE
ANOVA.settings.gaussianNorm<-FALSE
ANOVA.fakeCS<-FALSE
ANOVA.settings.CELL_LINES<-'LOCALcellLines'

ANOVA.settings.includeTISSUE_Factor<-FALSE
ANOVA.settings.includeMSI_Factor<-FALSE
ANOVA.settings.includeMEDIA_factor<-FALSE
ANOVA.settings.oneFeatureOnly<-NULL
ANOVA.settings.excludeHyperMetData<-TRUE

ANOVA.settings.additionalFeaturesOnly<-FALSE

ANOVA.settings.ESSENTIALITY_domain<-"logFCs" #BAGELr or Mageck
ANOVA.settings.cellLinesToRemove<-NULL

ANOVA.settings.featFactorPopulationTh<-3
ANOVA.minVulnerableCellLines<-2
ANOVA.settings.MSIfactorPopulationTh<-3

gdscANOVA.settings.pval_correction_method<-'qvalue'

#general QC normalised FC for testing input loading
Qnorm<-matrix(rnorm(200,0,2),nrow=20,ncol=10)
rownames(Qnorm)<-paste0("gene",letters[1:20])
colnames(Qnorm)<-paste0("cl",1:10)

#create test input BEM
BEMDep<-matrix(c(1,1,1,0,0,0,sample(c(0,1),28,replace=TRUE),c(1,1,1,0,0,0)),nrow=2,ncol=20,byrow=TRUE)
rownames(BEMDep)<-paste0("gene",letters[1:2],"_mut")
colnames(BEMDep)<-paste0("cl",1:20)

BEMDep2<-BEMDep
BEMDep2[1,7]<-NA
BEMDep2[2,8]<-NA
BEMDepNA<-rbind(BEMDep2,c(1,1,rep(NA,18)))

aovTbem<-BEMDepNA[1:2,paste0("cl",c(1:6,10:13))]
colnames(aovTbem)<-paste0("cl",1:10)
BEMin<-BEMDepNA[,paste0("cl",c(1:6,10:13))]
colnames(BEMin)<-paste0("cl",1:10)
admfiles<-c("09_ADM_Breast_coreFitnessGenes.RData","09_ADM_Breast.Carcinoma_coreFitnessGenes.RData",
            "09_ADM_Cervical.Carcinoma_coreFitnessGenes.RData","09_ADM_Cervix_coreFitnessGenes.RData",
            "09_ADM_panessProf_Breast.Carcinoma.RData","09_ADM_panessProf_Breast.RData",
            "09_ADM_panessProf_Colorectal.Carcinoma.RData","09_ADMBreast_OutputStats.RData",
            "09_ADMBreast.Carcinoma_OutputStats.RData","09_ADMEndometrium_OutputStats.RData")
ANOVA.settings.analysisType<-'CS'
ANOVA.fakeCS<-FALSE
InputFeatures<-list()
InputFeatures$BEM<-aovTbem
#add error in names of features:
aovTbemE<-aovTbem
rownames(aovTbemE)[1]<-""
InputFeaturesE<-list()
InputFeaturesE$BEM<-aovTbemE

expectFeat<-c(5,4)
names(expectFeat)<-c("N FEATURE pos","N FEATURE_neg")
aovRes<-ANOVA_totalANOVA(fn="./",ESSprofiles=t(Qnorm),InputFeatures=InputFeatures,inparallel=FALSE,
                         geneAnnotations = geneAnnot,ANOVA.settings.featFactorPopulationTh=2,
                         gdscANOVA.settings.pval_correction_method="BH",ANOVA.settings.analysisType="CS",ANOVA.fakeCS=FALSE,
                         ANOVA.settings.includeMSI_Factor=FALSE)

#compound Child BEM
BEMDepC<-matrix(c(0,0,0,1,1,1,sample(c(0,1),28,replace=TRUE),c(1,1,1,0,0,0)),nrow=2,ncol=20,byrow=TRUE)
rownames(BEMDepC)<-paste0("gene",letters[3:4],"_mut")
colnames(BEMDepC)<-paste0("cl",1:20)

intTerms<-data.frame(Child=c("genec_mut","gened_mut","gened_mut"),Parent=c("genea_mut","genea_mut","geneb_mut"),stringsAsFactors = FALSE)
#compound Parent BEM - try using BEMDep


#CT<-ANOVA_createInputFeaturesCompound(additional_features=ANOVA.additionalFeatures,
#                                      oneFeatureOnly = ANOVA.settings.oneFeatureOnly,
#                                      excludeHyperMetData = ANOVA.settings.excludeHyperMetData,
#                                      additional_features_only = ANOVA.settings.additionalFeaturesOnly,TOTALBEM=BEMDepC,
#                                      selectedCellLines = paste0('cl',1:20),
#                                      MASTER_LIST=MASTER_LIST,ANOVA.settings.featFactorPopulationTh=3,
#                                      ParentBEM=BEMDepNA,MutualExclusive=FALSE)$CompoundTerms
#t1<-all.equal(intTerms,CT)

