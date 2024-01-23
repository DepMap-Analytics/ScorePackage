set.seed(123)
geneAnnot<-data.frame(cbind(name=paste0("gene",letters[1:20]),gene_family=paste0("group",1:20),location_sortable=NA,
                            entrez_id=NA,
                            ensembl_gene_id=NA,
                            pubmed_id=NA,string_id=paste0("SID",letters[1:20])),stringsAsFactors = FALSE)

PPInet<-cbind(symbol=paste0("gene",letters[1:3]),STRING_id=paste0("SID",letters[1:3]))
rownames(PPInet)<-PPInet[,"symbol"]



#general scaled BF matrix
sBFs<-matrix(rnorm(60,0,2),nrow=3,ncol=20)
colnames(sBFs)<-paste0("cl",1:20)
rownames(sBFs)<-paste0("gene",letters[1:3])


#general test FC data
FCs<-matrix(rnorm(60,0,2),nrow=3,ncol=20)
rownames(FCs)<-paste0("gene",letters[1:3])
colnames(FCs)<-paste0("cl",1:20)

#general fully connected PPIigraph
testcon<-matrix(1,3,3)

dimnames(testcon)<-list(paste0("SID",letters[1:3]),paste0("SID",letters[1:3]))
PPIgraphT<-graph_from_adjacency_matrix(testcon)


#need a test CMP/MASTER_LIST

#need a set of bDepletions
bDepletions<-matrix(sample(c(0,1),200,replace=TRUE),nrow=20,ncol=10)
dimnames(bDepletions)<-list(rownames(Qnorm),colnames(Qnorm))

MASTER_LIST<-data.frame(model_id=paste0("cl",1:20),BROAD_ID=paste0("ACH",21:40),tissue=c(rep("Breast",5),rep("Large Intestine",10),rep("Ovary",5)),msi_status=c(rep("MSI",10),rep("",5),rep("MSS",5)))


