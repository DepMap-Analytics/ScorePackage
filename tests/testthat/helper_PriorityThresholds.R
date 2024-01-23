set.seed(1234)

#generate samples from two normal distributions and calculate priority threshold for them:
set1<-rnorm(50,20,5)
set2<-rnorm(40,10,2)
kess<-density(set1, kernel = "gaussian")
knon<-density(set2, kernel = "gaussian")

x <- seq(0,100,0.01)
nonfitx <- approx(knon$x,knon$y,x)$y

logratio_sample <- log2( approx(kess$x,kess$y,x)$y / approx(knon$x,knon$y,x)$y )

priority_threshold<-x[min(which(logratio_sample>=1 & x>=mean(set2)))]



TOTRESpt<-data.frame(`Depleted Gene`=paste0("gene",1:10),FDR=c(2,6,10,50,1,78,3,24,15,90),
                   `P-Value`=rep(0.000000001,10),IQR=c(6,2,4,4,1,6,7,8,2,5),logFC=c(2,1,1.2,0.55,1,7.8,3,2.4,1.5,9),
                   f2=c(2,1,1.2,0.55,1,7.8,3,2.4,1.5,9),stringsAsFactors=FALSE)
MARKERclass<-data.frame(`Depleted Gene`=paste0("gene",1:10),ANALYSIS=rep("Breast.Carcinoma",10),
                        CLASS=c("A","C","","","B","","A","","",""),stringsAsFactors = FALSE)
ctype<-"Breast.Carcinoma"
mutationNumber<-""



PRIORITY_vectors<-c(56,75,16,65,52,81,26,97,49,74)

PRIORITY_vectorsL3<-c(75,51,68,100,74,55, 94,17,43,29)
names(PRIORITY_vectors)<-paste0("gene",1:10)
names(PRIORITY_vectorsL3)<-paste0('gene',1:10)
th=30
indications<-vector("list",length=4)
names(indications)<-c("AACS","AACO","AOD","NAP")
indications[[1]]<-PRIORITY_vectors[c("gene1","gene10")]
indications[[2]]<-PRIORITY_vectors[c("gene2","gene4","gene7")]
indications[[3]]<-PRIORITY_vectors[c("gene3","gene6","gene8")]
indications[[4]]<-PRIORITY_vectors[c("gene5","gene9")]

PCHSYM<-c(23,24,22,21)
names(PCHSYM)<-names(indications)

TRACTABLE<-vector("list",length=10)
TRACTABLE[[1]]<-PRIORITY_vectors["gene10"]
TRACTABLE[[2]]<-PRIORITY_vectors[c("gene9","gene8")]
TRACTABLE[[3]]<-c()
TRACTABLE[[4]]<-PRIORITY_vectors["gene3"]
TRACTABLE[[5]]<-PRIORITY_vectors[c("gene1","gene2")]
TRACTABLE[[6]]<-c()
TRACTABLE[[7]]<-c()
TRACTABLE[[8]]<-PRIORITY_vectors[c("gene7")]
TRACTABLE[[9]]<-PRIORITY_vectors[c("gene5","gene6")]
TRACTABLE[[10]]<-PRIORITY_vectors["gene4"]
ANOVA.results.dir<-"ANOVAres"
LoadedFile<-"ANOVAres/Breast.Carcinoma/OUTPUT/ANOVA_results.rdata"

Pvec<-list()
Pvec[[1]]<-PRIORITY_vectors
Pvec2<-Pvec
Pvec2[[1]]<-PRIORITY_vectorsL3
names(Pvec)<-"Breast.Carcinoma"
names(Pvec2)<-"Breast.Carcinoma"

BiomarkerRes<-list(class_A_id=c(1,7),class_B_id=c(1,5,7),class_C_id=c(1,2,5,7),class_D_id=c(1,2,3,5,7))


#removed a final NA from the newSymbols vector 22.8.20
newSymbols<-c(8,3,4,2)
names(newSymbols)<-c('A','B','C','D')
ALLMARK<-newSymbols[MARKERclass$CLASS]
ALLBUCK<-c(5,5,4,10,9,9,8,2,2,1)
names(ALLBUCK)<-paste0("gene",1:10)
ALLPCH<-c(23,24,22,24,21,22,24,22,21,23)
PriorityRes<-data.frame(ctype=rep("Breast.Carcinoma",10),TARGET=TOTRESpt$Depleted.Gene,BUCKET=ALLBUCK,
                PRIORITY=PRIORITY_vectors,
                MARKERsymbol=ALLMARK,
                indiSymbol=ALLPCH,
                PRIORITYL3=PRIORITY_vectorsL3[names(PRIORITY_vectors)],
                stringsAsFactors = FALSE)

PriorityRes<-PriorityRes[PriorityRes$PRIORITY>=th,]
PriorityRes<-PriorityRes[c("gene10","gene8","gene9","gene2","gene1","gene6","gene5","gene4"),]

#rownames(PriorityRes)<-NULL

r1<-c(61, 40,  8, 65, 44)
r2<-c(49, 14, 51, 48,26)
names(r1)<-paste0("gene",1:5)
names(r2)<-paste0("gene",4:8)
testCR<-c(r1[c(1:5)],r2[c(3:5)])
