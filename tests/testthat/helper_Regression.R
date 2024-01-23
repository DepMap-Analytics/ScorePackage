
BEMDep2<-BEMDep
BEMDep2[1,7]<-NA
BEMDep2[2,8]<-NA
BEMDepNA<-rbind(BEMDep2,c(1,1,rep(NA,18)))
rownames(BEMDepNA)<-paste0("gene",letters[1:3],"_mut")


lmout<-list()
for(i in 1:nrow(BEMDepNA)){

  output<-limmaRegress(FCs[,!is.na(BEMDepNA[i,])],BEMDepNA[i,])
  output$FEATURE<-paste0(rownames(BEMDepNA)[i],"_Expr")
  output$`Depleted Gene`<-rownames(output)
  lmout[[i]]<-output
}
limmaout<-do.call(rbind,lmout)
rownames(limmaout)<-NULL
limmaout[is.na(limmaout$P.Value),"P.Value"]<-1
limmaout$FDR<-qvalue(limmaout$P.Value)$qvalues*100
limmaout<-limmaout[order(limmaout$FDR,decreasing=FALSE),]

testData1<-rnorm(10)
testData2<-rnorm(15,-1,1.5)

testD1<-rbind(testData1,testData1)
testD2<-rbind(testData2,testData2)

g1<-ANOVA_glass_Ds(testData1,testData2)$g1
