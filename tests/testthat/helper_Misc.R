testFCs<-matrix(rnorm(20),2,10)
rownames(testFCs)<-paste0("gene",letters[1:2])
colnames(testFCs)<-c(paste0("SIDM",1:4),paste0("ACH",12:17))
testFCs2<-matrix(rnorm(20),2,10)
rownames(testFCs2)<-paste0("gene",letters[1:2])
colnames(testFCs2)<-c(paste0("SIDM",1:10))

annot<-data.frame(model_id=paste0("SIDM",1:10),BROAD_ID=paste0("ACH",10:19),stringsAsFactors = FALSE)


inputbinary<-matrix(sample(c(0,1),size=20,replace=TRUE),nrow=2,ncol=20)
outputbinary<-inputbinary
rownames(outputbinary)<-c("mut1",c("mut2, mut3"))
inputbinary<-rbind(inputbinary,inputbinary[2,])
rownames(inputbinary)<-paste0("mut",1:3)
