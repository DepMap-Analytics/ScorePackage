#create input data for testing OT55 Network Functions
testPC<-matrix(rnorm(20),5,4)

testFCs<-matrix(rnorm(20),5,4)
rownames(testFCs)<-paste0("gene",letters[1:5])

posab<-apply(testFCs[c("genea","geneb"),],2,mean)
negac<-apply(testFCs[c("genea","genec"),],2,min)

binaryMat<-matrix(sample(c(0,1),40,replace=TRUE),nrow=2,ncol=20)
binaryMat<-rbind(binaryMat,rep(1,20))
colnames(binaryMat)<-paste0("cl",1:20)
rownames(binaryMat)<-paste0("gene",letters[1:3])


responsive<-binaryMat[1,]|binaryMat[2,]
corpair_ab1<-cor(sBFs[1,responsive],sBFs[2,responsive])



allBin<-matrix(1,nrow=3,ncol=20)
dimnames(allBin)<-dimnames(binaryMat)

binaryMat2<-matrix(sample(c(0,1),60,replace=TRUE),nrow=3,ncol=20)
corall_bc<-cor(sBFs[2,],sBFs[3,])
responsive<-binaryMat2[1,]|binaryMat2[2,]
corpair_ab<-cor(sBFs[1,responsive],sBFs[2,responsive])
responsive<-binaryMat2[1,]|binaryMat2[3,]
corpair_ac<-cor(sBFs[1,responsive],sBFs[3,responsive])
responsive<-binaryMat2[2,]|binaryMat2[3,]
corpair_bc<-cor(sBFs[2,responsive],sBFs[3,responsive])
dimnames(binaryMat2)<-dimnames(binaryMat)



cmp<-NULL
CancerType<-"PANCAN"
removeCore<-FALSE
ncores<-1

testCM<-matrix(1,3,3)

testCM[1,2]<-corpair_ab
testCM[2,1]<-corpair_ab
testCM[1,3]<-corpair_ac
testCM[3,1]<-corpair_ac
testCM[2,3]<-corpair_bc
testCM[3,2]<-corpair_bc
dimnames(testCM)<-list(paste0("gene",letters[1:3]),paste0("gene",letters[1:3]))

testCM2<-cor(t(sBFs))
subless<-abs(testCM2)<abs(testCM)
testCM2[subless]<-testCM[subless]

PaircorT<-0.1
testCM3<-testCM2
diag(testCM3)<-NA
testCM3[upper.tri(testCM3)]<-NA
Cmat<-testCM2
SCheck<-abs(Cmat)>=PaircorT+0
SCheck<-as.matrix(SCheck)
SCheck[upper.tri(SCheck)] = NA
diag(SCheck)<-NA

pairVals<-reshape2::melt(SCheck, na.rm=T)
pairCors<-reshape2::melt(testCM3,na.rm=T)
pairCors<-cbind(pairCors,pairVals[,3])
pairCors<-pairCors[pairCors[,4],]
rownames(pairCors)<-paste0("Group",1:nrow(pairCors))
name<-paste0("Group",1:nrow(pairCors))
pairList<-pairCors[,1:2]

#extra_sBF example for extra profiles needed to get pair l1/l2 priority scores:
extra_sBF<-apply(sBFs[1:2,],2,mean)
#for pos pick min FC and get the binary mat as either:
extra_bDep<-(binaryMat[1,]&binaryMat[2,])+0

#extra_sBF example for extra profiles needed to get pair l1/l2 priority scores:
maxvals<-apply(sBFs[1:2,],2,which.max)
extra_sBFN<-sapply(1:ncol(sBFs),function(x) sBFs[maxvals[x],x])
#for pos pick min FC and get the binary mat as either:
extra_bDepN<-(binaryMat[1,]|binaryMat[2,])+0
