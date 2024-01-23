tissue<-"Breast"
refSets<-vector("list",9)
names(refSets)<-c("curated_BAGEL_essential","BAGEL_nonEssential","BAGEL_essential",
                  "EssGenes.ribosomalProteins","EssGenes.DNA_REPLICATION_cons",
                  "EssGenes.PROTEASOME_cons","EssGenes.SPLICEOSOME_cons",
                  "histones","EssGenes.KEGG_rna_polymerase")
cancerDrivers<-list()
cancerDrivers[[1]]<-c("KRAS","BRAF","WRN","CCND1")


refSets[[1]]<-c("ABCD","EFGH")
refSets[[4]]<-c("USP7")

PRIORIEss<-data.frame(NGUIDES=rep(4,5), IN_CURATED_BAGEL_REF_SET=rep(0,5), ribosomalProteins=rep(0,5),
                      dna_replication=rep(0,5), proteasome=rep(0,5), spliceosome=rep(0,5), histones=rep(0,5),
                      rna_polymerase=rep(0,5), `predicted as PANCAN CoreFitness`=rep(0,5), `predicted asBreastCoreFitness`=rep(0,5),stringsAsFactors=FALSE)

rownames(PRIORIEss)<-c("USP7","DOHH", "CBWD1", "DHX8","RPL12")
NGUIDES<-rep(4,5)
names(NGUIDES)<-c("USP7","DOHH", "CBWD1", "DHX8","RPL12")

coreFitnessGenes<-c("USP7","ABCD")
PanCancerCoreFitnessGenes<-c("RPL12","EDFG")

PRIORIEss["USP7",c("predicted.asBreastCoreFitness","ribosomalProteins")]<-1
PRIORIEss["RPL12","predicted.as.PANCAN.CoreFitness"]<-1
PRIORIEss<-PRIORIEss[,3:9]
colnames(PRIORIEss)[7]<-'predictedAsPanCancerCoreFitness'

