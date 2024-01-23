L1<-data.frame(gene=c("a","b"),ANOVAfdrA=c(1,1),ANOVAfdrB=c(1,1),MutPrimTum=c(0,1),stringsAsFactors = FALSE)

tractability_both<-data.frame( id=c("ENSG00000001561","ENSG00000004948","ENSG00000005073","ESNG1","ENSG2"), antibody=c(4,5,10,1,2), small.molecule=c(4,4,10,2,2),
                               type=c("small-molecule/antibody","small-molecule","small-molecule/antibody","small-molecule","small-molecule/antibody"), min_bucket=c(4,4,10,1,2),
                               symbol=c("ENPP4","CALCR","HOXA11","ESR1","ATF4"),stringsAsFactors = FALSE)
rownames(tractability_both)<-tractability_both$symbol

tract_ab<-data.frame(Target=c("ENPP4","CALCR","HOXA11"), Tractable.ranking.bucket=c(4,5,10), drug_name=rep("",3), max_phase=rep("",3), moa_chembl=rep("",3), indication_efo_term=rep("",3),stringsAsFactors = FALSE)
rownames(tract_ab)<-tract_ab$Target

tract_sm<-data.frame(Target=c("ENPP4","CALCR","HOXA11"), Tractable.ranking.bucket=c(4,4,10), drug_name=rep("",3), max_phase=rep("",3), moa_chembl=rep("",3), indication_efo_term=rep("",3),stringsAsFactors = FALSE)
rownames(tract_sm)<-tract_sm$Target


c1<-c("SIDM00446_BAGEL_SIG_DEP", "SIDM00446_NOT_EXPRESSED", "SIDM00446_HOMDEL", "SIDM00446_Mgk_10percFDR", "SIDM00446_Mgk_5percFDR", "SIDM00446_FOLD_1_sBF", "SIDM00446_FOLD_2_sBF", "SIDM00446_FOLD_3_sBF",
"SIDM00446_FOLD_minus2_FC", "SIDM00446_FOLD_minus3_FC", "SIDM00446_FOLD_minus5_FC", "SIDM00446_HIGHLY_EXP", "SIDM00446_MUT", "SIDM00446_COSMIC_MUT", "SIDM00446_IN_DEP_PATH")
c2<-c1
c2<-gsub("0446","1234",c2)
#currentBlock[,c(6:8,12,13,15)]*100
#column 4 is the mageck 10 fdr column
#columns 6-8 are the scaled BF columns
#column 12 is highly expressed
#column 13 is mutated patient population
#column 15 is in depleted path (slap enrich)
#  currentL2score<-( currentBlock[,4]*22+currentBlock[,6]*16+currentBlock[,7]*16+currentBlock[,8]*16+
#currentBlock[,12]*10+currentBlock[,13]*10+currentBlock[,15]*10)
#sig1 22+2*16+20 =74
sig1<-c(1,0,0,1,1,1,1,0,0,0,0,0,1,0,1)
#sig1.2 = 0
sig1.2<-c(1,0,0,1,0,1,0,0,0,0,0,0,1,0,1)
#sig2 74
sig2<-c(1,0,0,1,1,1,1,0,0,0,0,1,1,0,0)
#sig2.2 74
sig2.2<-c(1,0,0,1,1,1,1,0,0,0,0,0,1,0,1)
L2<-rbind(c(sig1,sig1.2),c(sig2,sig2.2))
colnames(L2)<-c(c1,c2)
rownames(L2)<-c("ENPP4","CALCR")
#column 5 has to be greater than zero. extra filter.
L2score<-cbind(c(74,74),c(0,74))

rownames(L2score)<-rownames(L2)
colnames(L2score)<-rep("currentL2score",ncol(L2score))

FINAL_priority<-matrix(0,nrow=2,ncol=2,dimnames=list(c('ENPP4',"CALCR"),c("dummy1","dummy2")))
FINAL_priorityPairs<-matrix(0,nrow=2,ncol=2,dimnames=list(c('ENPP4|ESR1',"CALCR&ATF4"),c("dummy1","dummy2")))
TRACTinfo<-tractability_both[1:2,]
TRACTinfoPairs<-rbind(apply(tractability_both[c(1,4),],2,function(x)paste(x[1],x[2],sep="//")),
                      apply(tractability_both[c(2,5),],2,function(x)paste(x[1],x[2],sep="//")))
rownames(TRACTinfoPairs)<-rownames(FINAL_priorityPairs)
L3scores<-NULL
#L3list
#GENES
subdirlist<-c("Cont","Expr","Genomic","Met","CN","Prot","Int","Compound","ME","CompoundME","Epi")
glist<-rownames(tractability_both)[1:2]
CompoundIndications<-list()
CompoundIndications$compounds<-c("|","")
CompoundIndications$indications<-c("|","")

L1_NNL2<-c(20,30)
L1_NNL2L3<-c(24,28)
names(L1_NNL2)<-glist
names(L1_NNL2L3)<-glist

TRACTABLE<-vector("list",10)
TRACTABLEL3<-vector("list",10)
TRACTABLE[[4]]<-sort(L1_NNL2,decreasing=T)


#Starting information for L3 script tests
GENES<-c("MYH9","MEF2D")

TargetData<-data.frame(Target=c("MYH9","MEF2D"),Score=c(37.5,37.5),FDRScore=c(15.625,15.7894736842105),
                       MaxTargetScoreFDR=c(15.625,15.7894736842105),PvalScore=c(50,100),	RWRscore=c(25,50),BMType=c("CN","CN"),
                       CancerType=c("Acute.Myeloid.Leukemia","Acute.Myeloid.Leukemia"),	ScoreType=c("Avg","Avg"),
                       BMs=c("9606.ENSP00000381891","9606.ENSP00000298139//9606.ENSP00000349275"),NumberBMs=c(1,2),	RWRclasses=c("25","25//25"),
                       PvalClasses=c("50","50//50"),	TargetID=c("9606.ENSP00000216181","9606.ENSP00000271555"),	Dependencytype=c("Single","Single"),
                       StringScore=c(92,898.25),	IDS=c("MYH9Acute.Myeloid.Leukemia","MEF2DAcute.Myeloid.Leukemia"),stringsAsFactors=FALSE)

T1<-TargetData
T2<-rbind(TargetData,TargetData)
T2[3:4,"ScoreType"]<-"Sum"
TargetOut<-TargetData[,c("Score","FDRScore","MaxTargetScoreFDR","PvalScore","RWRscore","NumberBMs",
                            "StringScore","BMs","RWRclasses","PvalClasses","BMType")]
colnames(TargetOut)[ncol(TargetOut)]<-"biomarkertype"

rownames(TargetOut)<-GENES
CTYPE<-"Acute.Myeloid.Leukemia"
otype<-rep("CN",2)
L3list<-list()
temp<-TargetData
L3list[[1]]<-TargetOut
TO2<-TargetOut
TO2[,"biomarkertype"]<-"mut"
TO2[1,"FDRScore"]<-30
L3list[[2]]<-TO2
L3scores<-c(0,0,25,0)
L3names<-c("CN","Genomic")
names(L3list)<-L3names
