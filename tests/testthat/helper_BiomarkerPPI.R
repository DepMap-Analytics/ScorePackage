res<-data.frame(assoc_id=paste0("a",1),`Depleted Gene`=c("genea&geneb"),string_id=c("PS1//PS2"),stringsAsFactors = FALSE)
colnames(res)[2]<-'Depleted Gene'
splitres<-data.frame(assoc_id=rep(paste0("a",1),2),`Depleted Gene`=c("genea","geneb"),string_id=c("PS1","PS2"),stringsAsFactors = FALSE)
colnames(splitres)[2]<-'Depleted Gene'

res2<-data.frame(assoc_id=paste0("a",2),`Depleted Gene`=c("genec|gened"),string_id=c("PS3//PS4"),stringsAsFactors = FALSE)
colnames(res2)[2]<-'Depleted Gene'
splitres2<-data.frame(assoc_id=rep(paste0("a",2),2),`Depleted Gene`=c("genec","gened"),string_id=c("PS3","PS4"),stringsAsFactors = FALSE)
colnames(splitres2)[2]<-'Depleted Gene'


resM<-data.frame(assoc_id=paste0("a",1),FEATURE=c("muta&exprb"),BMstringID=c("PS1&PS2"),stringsAsFactors = FALSE)

splitresM<-data.frame(assoc_id=rep(paste0("a",1),2),FEATURE=c("muta","exprb"),BMstringID=c("PS1","PS2"),stringsAsFactors = FALSE)


resM2<-data.frame(assoc_id=paste0("a",2),FEATURE=c("mutc|exprd"),BMstringID=c("PS3|PS4"),stringsAsFactors = FALSE)

splitresM2<-data.frame(assoc_id=rep(paste0("a",2),2),FEATURE=c("mutc","exprd"),BMstringID=c("PS3","PS4"),stringsAsFactors = FALSE)




n <- 5
# starting distribution (has to sum to one)
p0    <- as.vector(rmultinom(1, 1, prob=rep(.2, n)))
# adjacency matrix (either normalized or not)
graph <- matrix(abs(rnorm(n*n)), n, n)
diag(graph)<-0
normW<-normalizeW(graph,degreelist=colSums(graph))
csums<-rep(1,n)

# computation of stationary distribution
rwrdr<-random.walk(p0,graph,r=0.4,thresh=1e-012)
g <- make_ring(10)
V(g)$name<-paste0("PS",1:10)
p02<- as.vector(rmultinom(1, 1, prob=rep(.1, 10)))
g_adj<-as_adjacency_matrix(g)
g_adj<-as.matrix(g_adj)
diag(g_adj)<-0
rwrdr2<-random.walk(p02,g_adj,r=0.4,thresh=1e-012)
rownames(rwrdr2$p.inf)<-paste0("PS",1:10)
dimnames(rwrdr2$transition.matrix)<-list(paste0("PS",1:10),paste0("PS",1:10))
normW2<-normalizeW(g_adj,degreelist=colSums(g_adj))
csums2<-rep(1,10)
names(csums2)<-paste0("PS",1:10)

rwrres<-getRestartProb(g)
PPInetTest<-cbind(c("ABD","GFH","HFH","KJL","LO","HI","K1","K2","L3",paste0("target",1:11)),paste0("PS",1:20))
colnames(PPInetTest)<-c("symbol","STRING_id")
rownames(PPInetTest)<-PPInetTest[,1]
PPInetTest<-data.frame(PPInetTest,stringsAsFactors=FALSE)
TOTRES<-data.frame(c("a1","a2"), c("K1_mut","K2_mut, L3_mut"),	c("ABD","GFH"),	c("gene1","gene2"),	rep(NA,2),	rep(NA,2)	,rep(NA,2),	rep(NA,2), c(1,2),	c(3:4),	c("PS1","PS2"),c(6,5)	,c(7,8),	c(-1,-2),	c(0.2,0.5),c(2,3),	c(0.1,0.2),	c(0.3,0.23),c(3,4),c(3,4)	,c(2,3)	,c(0.0001,0.00003),	rep(NA,2),rep(NA,2),	rep(NA,2),c(0.0003,0.0002),c(0.0023,0.0056),c("PS7","PS8,PS9"),	c(3.4,6.7),c(4.5,10.4),stringsAsFactors = FALSE)
colnames(TOTRES)<-c("assoc_id",	"FEATURE",	"Depleted Gene","Name",	"gene family",	"location",	"entrez id",	"ensembl", "gene id,",	"pumbed id","string_id",	"N FEATURE pos",	"N FEATURE_neg",	"FEATUREpos_ESS_MEAN",	"FEATUREneg_ESS_MEAN",	"FEATURE_deltaMEAN_ESS",	"FEATUREpos_ESS_sd",	"FEATUREneg_ESS_sd",	"FEATURE_ESS_effect_size",	"FEATUREpos_Glass_delta",	"FEATUREneg_Glass_delta",	"FEATURE_ANOVA_pval",	"Tissue_ANOVA_pval",	"MEDIA_ANOVA_pval",	"MSI_ANOVA_pval",	"FEATURE_ESS_T_pval",	"FEATURE_T_pvalEV",	"BMstringID",	"ANOVA FEATURE FDR %",	"t-test (unequal var) FDR %")
RWRinputAll<-graphMatricesRWR(dir.Results="./",PPIigraph = g,PPInet = PPInetTest,single=FALSE,AnovaRes =TOTRES,subdir="")

startvec<-matrix(0,nrow=10,ncol=3)
startvec[7,1]<-1
startvec[8,2]<-1
startvec[9,3]<-1

rownames(startvec)<-paste0("PS",1:10)
colnames(startvec)<-c("PS7","PS8","PS9")
#RWRinputAll<-graphMatricesRWR(dir.Results="./",PPIigraph = g,PPInet = PPInet,single=FALSE,AnovaRes =TOTRES,subdir=subdir,pvalCol="FEATURE_compound_pval",deltaPosCol="FEATUREcompound_pos_Glass_delta",deltaNegCol="FEATUREcompound_neg_Glass_delta",Int=TRUE )
dlist<-colSums(RWRinputAll$WeightMat)

MarkerInfo<-data.frame(TOTRES[,c(3,2,28,11,1)],as.character(rep("mut",2)),TOTRES[,c(22,29,20,21)],stringsAsFactors = FALSE)
colnames(MarkerInfo)<-c("Depleted.Gene", "FEATURE", "StringID", "TstringID", "assoc_id","type", "AnovaPval", "AnovaFdr", "AnovaDeltaP", "AnovaDeltaN")
MarkerInfo<-rbind(MarkerInfo,MarkerInfo[2,])
MarkerInfo[,"FEATURE"]<-c("K1_mut","K2_mut","L3_mut")
MarkerInfo[,"StringID"]<-c("PS7","PS8","PS9")
rownames(MarkerInfo)<-NULL
test2<-TOTRES
colnames(test2)<-make.names(colnames(test2))
colnames(test2)[11]<-"TstringID"
colnames(test2)[22]<-"AnovaPval"
colnames(test2)[28]<-"StringID"
test2<-rbind(test2,test2[2,])
test2[,"FEATURE"]<-c("K1_mut","K2_mut","L3_mut")
test2[,"StringID"]<-c("PS7","PS8","PS9")


othem<-data.frame(rep("",3),c("K1","K2","L3"),c("mut","mut","mut"),stringsAsFactors = FALSE)
colnames(othem)<-c("ppi_target","markers","type")
test2<-cbind(test2,othem)

#now testing for the BM PPI scoring
classMat<-data.frame(Pval=c(0.05,1,0.001,0.001,0.1),Fdr=c(101,10,101,5,101),DeltaP=c(0,0,0,1,0),DeltaN=c(0,0,0,1,0),class=c("D","B","C","A","E"),stringsAsFactors = FALSE)
biomarkerMat<-data.frame(Score=c(0.005,0.01,0.02),class=c("C","B","A"),stringsAsFactors = FALSE)
classMatExpr<-data.frame(Pval=c(0.05,1,0.001,0.001,0.1),Fdr=c(101,10,101,5,101),Range=c(0,0,0,7,0),class=c("D","B","C","A","E"),stringsAsFactors = FALSE)
classMatCN<-data.frame(Pval=c(0.05,1,0.001,0.001,0.1),Fdr=c(101,10,101,5,101),Range=c(0,0,0,2,0),class=c("D","B","C","A","E"),stringsAsFactors = FALSE)
RWRscoreMat<-data.frame(ScoringClass=c("A","B","C","D"),RWRthresh=c(1e-3,5e-4,1e-4,5e-5),Score=c(100,75,50,25))
PvalscoreMat<-data.frame(ScoringClass=c("A","B","C","D"),Pvalthresh=c(5e-7,5e-6,5e-5,5e-4),Score=c(100,75,50,25))
FDRscoreMat<-data.frame(ScoringClass=c("A","B","C","D"),FDRthresh=c(5,10,15,20),Score=c(100,75,50,25))

RWRthresh<-matrix(0,nrow=20,ncol=5)
dimnames(RWRthresh)<-list(paste0("PS",1:20),paste0("PS",5:9))
RWRthresh[sample(1:100,20)]<-runif(10,0,0.2)
RWRthresh[5,1]<-0.4
RWRthresh[6,2]<-0.38
RWRthresh[7,3]<-0.41
RWRthresh[8,4]<-0.45
RWRthresh[9,5]<-0.423
#example thresh not accurate
AnovaFdr<-matrix(c(3.4,100,100,100,6.7,6.7),nrow=2,ncol=3,byrow = T)
dimnames(AnovaFdr)<-list(c("ABD","GFH"),c("K1_mut","K2_mut","L3_mut"))

BiomarkerScore<-c(0,0,62.5,75,0)
names(BiomarkerScore)<-c("NA","NA","K1_mut","K2_mut","L3_mut")

TargetScore<-data.frame(Score=c(62.5,75.0),FDRScore=c(100,87.5),BMtype=rep("mut",2),PvalScore=c(25,50),
                        RWRscore=rep(100,2),BMs=c("PS7","PS8"),NumberBMs=rep(1,2),RWRclasses=rep("100",2),
                        PvalClasses=c("25","50"),stringsAsFactors = FALSE)
rownames(TargetScore)<-c("ABD","GFH")

BMPriority<-list(K1_mut=c("ABD"=62.5),K2_mut=c("GFH"=75))

BMSindivid<-list(Scores=c("","","62.5","75",""),Targets=c("","","ABD","GFH",""))



