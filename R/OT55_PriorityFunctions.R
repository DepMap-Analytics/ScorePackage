Get_PriorityThreshold<-function(set1,set2,logratio=2){
  #set1 the  positive controls
  #set2 the negative controls
  kess<-density(set1, kernel = "gaussian")
  knon<-density(set2, kernel = "gaussian")

  x <- seq(0,100,0.01)
  nonfitx <- approx(knon$x,knon$y,x)$y

  logratio_sample <- log2( approx(kess$x,kess$y,x)$y / approx(knon$x,knon$y,x)$y )

  priority_threshold<-x[min(which(logratio_sample>=logratio & x>=mean(set2)))]
  return(priority_threshold)
}

Get_L1filename<-function(subdir,mutationNumber,ctype,omictype,dir.Results){
  if(subdir==""){
    L1filename<-paste(dir.Results,'/33_L1/',ctype,'_L0_L1',subdir,mutationNumber,'.csv',sep='')
  }else{
    if(subdir=="PairBD"){
      L1filename<-paste(dir.Results,'/33_L1/',subdir,"/",omictype,"/",ctype,"/",ctype,'_L0_L1',omictype,mutationNumber,'.csv',sep='')

    }else{
      L1filename<-paste(dir.Results,'/33_L1/',subdir,"/",ctype,"/",ctype,'_L0_L1',subdir,mutationNumber,'.csv',sep='')

    }
  }
  return(L1filename)

}

Get_L1score<-function(L1,componentNames=c("MutPrimTum","marker.Ttest.p...0.05..Class.D.","ANOVA.p...0.001..Class.C.",
                                       "ANOVA.FDR...30...Class.B.","ANOVA.FDR...25...and.p...0.001..or.MSI.FDR...5...and.p...0.001..and.Glass.Deltas...1..Class.A."),
                      componentWeights=rep(20,5)){
  L1score<-0
  if(length(componentNames)!=length(componentWeights)){
    stop('Must have number of component Weights equal number of components')
  }
  if(sum(componentWeights)!=100){
    stop('Component Weights must equal 100')
  }
  for(i in 1:length(componentNames)){
    temp<-L1[,componentNames[i]]*componentWeights[i]
    L1score<-L1score+temp
  }

  return(L1score)


}

Get_L2score<-function(L2,site){
  ns<-ncol(L2)
  nc<-ns/15

  L2score<-NULL


  for (j in 1:nc){
    #1.7.19 added in option for Combined data set that removes the mageck contribution to the L2 score
    #This is the format of L2:
    #currentBench<-
    # cbind(FILTER_MATRIX[,((i-1)*3+1):(i*3)],
    #      Mageck_10[,cellLines[i]],Mageck_5[,cellLines[i]],
    #      FOLD_1_sBF[,cellLines[i]],FOLD_2_sBF[,cellLines[i]],FOLD_3_sBF[,cellLines[i]],
    #      FOLD_minus2_FC[,cellLines[i]],FOLD_minus3_FC[,cellLines[i]],FOLD_minus5_FC[,cellLines[i]],
    #     HIGHLY_EXPRESSED[,as.character(cids[i])],
    #     IS_MUTATED[,cellLines[i]],
    #     IS_COSMIC_MUTATED[,cellLines[i]],
    #    is.element(rownames(FOLD_1_sBF),inDepletedENRpath)+0)

    currentBlock<-L2[,(15*(j-1)+1):(15*j)]
    #filt 1 is it depleted according to bagel
    #filt 2 is expressed
    #filt 3 is not homozygously deleted
    #add min sBF>1 filter 17.9.22
    currentL2filter<-(currentBlock[,1]>0 & currentBlock[,2]!=1 & currentBlock[,3]!=1&currentBlock[,5]>0)
    #changed to set to zero cell lines where it's not depleted:
    #currentL2filter<-(currentBlock[,1]>0 & currentBlock[,2]!=1 & currentBlock[,3]!=1&currentBlock[,4]==1)
    if(site=="Combined"){
      #currentL2score<-currentBlock[,c(6:8,12,13,15)]*100
      #column 4 is the mageck 10 fdr column
      #columns 6-8 are the scaled BF columns
      #column 12 is highly expressed
      #column 13 is mutated patient population
      #column 15 is in depleted path (slap enrich)
      currentL2score<-( currentBlock[,4]*22+currentBlock[,6]*16+currentBlock[,7]*16+currentBlock[,8]*16+
                  currentBlock[,12]*10+currentBlock[,13]*10+currentBlock[,15]*10)
      #currentL2score<-rowSums(currentL2score,na.rm=T)/(rowSums(!is.na(currentL2score)))
      #temp 7.12.20 just to see effect of weights on new version:
      #currentL2score<-( currentBlock[,6]*(3/8)*100+currentBlock[,7]*(1/8)*100+currentBlock[,8]*(1/8)*100+
      #                    currentBlock[,12]*(1/8)*100+currentBlock[,13]*(1/8)*100+currentBlock[,15]*(1/8)*100)

    }else{
      currentL2score<-(currentBlock[,4]*12.5+currentBlock[,5]*12.5+
                         currentBlock[,6]*12.5+currentBlock[,7]*12.5+currentBlock[,8]*12.5+
                         currentBlock[,12]*12.5+currentBlock[,13]*12.5+currentBlock[,15]*12.5)
    }

    currentL2score<-currentL2score*(currentL2filter+0)

    L2score<-cbind(L2score,currentL2score)
  }
  return(L2score)
}

Get_TRACTinfo<-function(subdir,FINAL_priority,tractability_both){
  if(subdir=="PairBD"){
    genepairs<-rownames(FINAL_priority)
    glist<-sapply(genepairs,function(x) strsplit(x,"[&|]"))
    temp<-lapply(glist,function(x) x[1])
    TRACTinfos<-lapply(glist,function(x) apply(tractability_both[x,],2,function(y)paste(y,collapse="//")))
    TRACTinfos<-do.call(rbind,TRACTinfos)
    glist<-unlist(temp)
  }else{
    glist<-rownames(FINAL_priority)
    TRACTinfos<-tractability_both[rownames(FINAL_priority),]
  }
  return(list(TRACTinfos=TRACTinfos,glist=glist))
}

Get_CompoundIndications<-function(TRACTinfos,glist,tract_sm,tract_ab){
  compounds<-rep(NA,nrow(TRACTinfos))
  indications<-rep(NA,nrow(TRACTinfos))

  if(!c("type")%in%colnames(TRACTinfos)){
    stop('Type column not in TRACTinfos')
  }
  idsm<-grep('small-molecule',TRACTinfos[,"type"])
  idab<-grep('antibody',TRACTinfos[,"type"])
  idboth<-grep("small-molecule/antibody",TRACTinfos[,"type"])

  if(!c("drug_name")%in%colnames(tract_sm)){
    stop('drug_name column not in tract_sm')
  }
  if(!c("drug_name")%in%colnames(tract_ab)){
    stop('drug_name column not in tract_ab')
  }
  if(!c("indication_efo_term")%in%colnames(tract_sm)){
    stop('indication_efo_term column not in tract_sm')
  }
  if(!c("indication_efo_term")%in%colnames(tract_ab)){
    stop('indication_efo_term column not in tract_ab')
  }
  compounds[idsm]<-tract_sm[glist[idsm],"drug_name"]
  compounds[idab]<-tract_ab[glist[idab],"drug_name"]
  compounds[idboth]<-paste(
    tract_sm[glist[idboth],"drug_name"],
    tract_ab[glist[idboth],"drug_name"],sep='|')



  indications[idsm]<-tract_sm[glist[idsm],"indication_efo_term"]
  indications[idab]<-tract_ab[glist[idab],"indication_efo_term"]
  indications[idboth]<-paste(
    tract_sm[glist[idboth],"indication_efo_term"],
    tract_ab[glist[idboth],"indication_efo_term"],sep='|')

  #add null results for targets not in any list:

  return(list(compounds=compounds,indications=indications))
}
Get_TractableBuckets<-function(subdir,tractability_both,glist,L1_NNL2,L1_NNL2L3){
  TRACTABLE<-vector("list",10)

  TRACTABLEL3<-vector("list",10)

  uu<-sort(unique(as.numeric(tractability_both[,"min_bucket"])))
  if(length(uu)!=10){
    warning('Number of tractability buckets not equal to 10')
  }
  if(subdir=="PairBD"){
    allGenesCheck<-setdiff(unlist(glist),rownames(tractability_both))

  }else{
    allGenesCheck<-setdiff(names(L1_NNL2),rownames(tractability_both))
  }
  if(length(allGenesCheck)>0){
    warning('Some priority targets not in tractability information')
    print(paste('missing targets:',allGenesCheck))
  }
  for (j in 1:10){
    if(subdir=="PairBD"){


      gcheck<-lapply(glist,function(x) length(intersect(x,rownames(tractability_both)[which(tractability_both$min_bucket==j)]))>0)

      if(sum(unlist(gcheck))>0){

        TRACTABLE[[j]]<-sort(L1_NNL2[names(glist)[unlist(gcheck)]],decreasing=TRUE)

        TRACTABLEL3[[j]]<-sort(L1_NNL2L3[names(glist)[unlist(gcheck)]],decreasing=TRUE)
      }else{
        warning(paste('No priority genes found for tractability bucket',j))

      }
    }else{
      genes<-intersect(rownames(tractability_both)[which(tractability_both$min_bucket==j)],
                       names(L1_NNL2))
      if(length(genes)>0){
        TRACTABLE[[j]]<-sort(L1_NNL2[genes],decreasing=TRUE)
        TRACTABLEL3[[j]]<-sort(L1_NNL2L3[genes],decreasing=TRUE)
      }else{
        warning(paste('No priority genes found for tractability bucket',j))

      }
    }

  }
  return(list(TRACTABLE=TRACTABLE,TRACTABLEL3=TRACTABLEL3))
}


Get_EFOIinfo<-function(EFOI_list,CTYPE,tractability_both,tract_ab,tract_sm){
  names(EFOI_list)<-tolower(make.names(names(EFOI_list)))
  antiCancerSpecificIndication<-NULL
  if(CTYPE!='PANCAN'){

    Ctest<-tolower(make.names(CTYPE))
    antiCancerSpecificIndication<-EFOI_list[[Ctest]]
  }

  OtherAntiCancerIndication<-EFOI_list$anti.cancer.generic

  OtherAntiCancerTargets<-sort(union(rownames(tract_ab)[grep(paste(OtherAntiCancerIndication,collapse = '|'),tract_ab[,"indication_efo_term"])],
                                     rownames(tract_sm)[grep(paste(OtherAntiCancerIndication,collapse = '|'),tract_sm[,"indication_efo_term"])]))

  if(CTYPE!='PANCAN' & length(antiCancerSpecificIndication)>0){
    AntiCancerSpecificTargets<-sort(union(rownames(tract_ab)[grep(paste(antiCancerSpecificIndication,collapse = '|'),tract_ab[,"indication_efo_term"])],
                                          rownames(tract_sm)[grep(paste(antiCancerSpecificIndication,collapse = '|'),tract_sm[,"indication_efo_term"])]))
  } else{
    AntiCancerSpecificTargets<-NULL
  }

  OtherAntiCancerTargets<-setdiff(OtherAntiCancerTargets,AntiCancerSpecificTargets)

  approved<-
    rownames(tractability_both)[which(as.numeric(tractability_both[,"min_bucket"])<4)]

  approvedAntiCancerSpecific<-intersect(approved,AntiCancerSpecificTargets)
  approvedAntiCancer<-intersect(approved,OtherAntiCancerTargets)
  approvedOthers<-setdiff(approved,union(approvedAntiCancer,approvedAntiCancerSpecific))

  nonApproved<-
    rownames(tractability_both)[which(as.numeric(tractability_both[,"min_bucket"])>=4)]
  return(list(approvedAntiCancerSpecific=approvedAntiCancerSpecific,approvedAntiCancer=approvedAntiCancer,
              approvedOthers=approvedOthers,nonApproved=nonApproved))
}

charPriority<-function(ctype,subdir,L1weight=20,L2weight=60,site,mutationNumber,dir.Results,omictype=NULL,L3column="MaxTargetScoreFDR",L3weight=20){

  CTYPE<-ctype
  L1filename<-Get_L1filename(subdir,mutationNumber,ctype,omictype,dir.Results)

  L1<-read.csv(L1filename, header=TRUE,
               stringsAsFactors = FALSE,row.names = 1)


  if(subdir=="PairBD"){
    load(paste(dir.Results,'/35_L2/PairBD/',ctype,'_L2.Rdata',sep=''))
    if(file.exists(file=paste(dir.Results,'/33_L3/PairBD',CTYPE,'_L0_L3_all.RData',sep=''))){
      load(file=paste(dir.Results,'/33_L3/PairBD',CTYPE,'_L0_L3_all.RData',sep=''))}else{
        L3scores<-NULL
      }

  }else{
    load(paste(dir.Results,'/35_L2/',ctype,'_L2.Rdata',sep=''))
    if(file.exists(file=paste(dir.Results,'/33_L3/',CTYPE,'_L0_L3_all.RData',sep=''))){
      load(file=paste(dir.Results,'/33_L3/',CTYPE,'_L0_L3_all.RData',sep=''))}else{
        L3scores<-NULL
      }

  }

  #changed 10.4.20 to intersect of L1 and L2 genes, there can be a discrepancy, guide library (e.g ANOVA used to define genes in L1, ~18K where L2 only has number in overlap fcs data)
  #problem with the guide library - might not have same hugo gene symbol as the L2 fcs. Need to fix this.
  #update 4.11.21 checked new guide_librarycombined contains all genes in the CorrectedData FC set

  rownames(L1)<-gsub("[()]","",rownames(L1))
  rownames(L2)<-gsub("[()]","",rownames(L2))

  usegenes<-intersect(rownames(L1),rownames(L2))
  if(length(usegenes)>0){
    #L1<-L1[rownames(L2),]
    L1<-L1[usegenes,]
    #temporary fix:
    L2<-L2[usegenes,]



    if(!dir.exists(paste0(dir.Results,'/35_PrioritizedHits/',subdir))){
      dir.create(paste0(dir.Results,'/35_PrioritizedHits/',subdir))
    }
    if(subdir=="PairBD"){
      if(!dir.exists(paste0(dir.Results,'/35_PrioritizedHits/',subdir,"/",omictype,"/"))){
        dir.create(paste0(dir.Results,'/35_PrioritizedHits/',subdir,"/",omictype,"/"))
      }
      if(!dir.exists(paste0(dir.Results,'/35_PrioritizedHits/',subdir,"/",omictype,"/",ctype,"/"))){
        dir.create(paste0(dir.Results,'/35_PrioritizedHits/',subdir,"/",omictype,"/",ctype,"/"))
      }

    }else{

      if(!dir.exists(paste0(dir.Results,'/35_PrioritizedHits/',subdir,"/",ctype,"/"))){
        dir.create(paste0(dir.Results,'/35_PrioritizedHits/',subdir,"/",ctype,"/"))
      }

    }

    L2score<-Get_L2score(L2,site)
    if(subdir=="PairBD"){
      save(L2score,file=paste(dir.Results,'/35_PrioritizedHits/',subdir,"/",omictype,"/",ctype,"/",ctype,'_L2Score.RData',sep=''))
    }else{
      save(L2score,file=paste(dir.Results,'/35_PrioritizedHits/',subdir,"/",ctype,"/",ctype,'_L2Score.RData',sep=''))
    }
    NonNullL2score<-NULL
    NonNullL2score<-rowSums(L2score>0)

    AvgL2Score<-rowMeans(L2score)
    AvgNN2Score<-rowSums(L2score)/NonNullL2score
    AvgNN2Score[is.na(AvgNN2Score)]<-0



    if(!is.null(L3scores)){
    L1L3<-Get_L1L3score(L1=L1,L3scores=L3scores,L3col=L3column)}else{
      L1L3<-L1
    }

    if(CTYPE!='PANCAN'){
      L1score<-Get_L1score(L1=L1L3)
    }else{
      L1score<-Get_L1score(L1=L1L3,componentNames = c("MutPrimTum","marker.Ttest.p...0.05...Glass.Deltas...1..Class.D.",
                                                      "ANOVA.p...0.001...at.least.one.Glass.Delta...1..Class.C.",
                                                      "ANOVA.FDR...30....at.least.one.Glass.Delta...1..Class.B.",
                                                      "ANOVA.FDR...25...and.p...0.001..or.MSI.FDR...5...and.p...0.001..and.Glass.Deltas...1..Class.A."))
    }

    names(L1score)<-rownames(L1)
    L1filter<-(rowSums(L1L3[,2:10]) == 0 & L1L3[,1]>1)


    if(CTYPE!='PANCAN'){
      L1scoreO<-Get_L1score(L1=L1)
    }else{
      L1scoreO<-Get_L1score(L1=L1,componentNames = c("MutPrimTum","marker.Ttest.p...0.05...Glass.Deltas...1..Class.D.",
                                                      "ANOVA.p...0.001...at.least.one.Glass.Delta...1..Class.C.",
                                                      "ANOVA.FDR...30....at.least.one.Glass.Delta...1..Class.B.",
                                                      "ANOVA.FDR...25...and.p...0.001..or.MSI.FDR...5...and.p...0.001..and.Glass.Deltas...1..Class.A."))
    }

    names(L1scoreO)<-rownames(L1)
    #changed 6.7.19 to allow for flexible weights
    #L1_avgL2<-(L1filter+0)*(NonNullL2score>2+0)*(L1score*30+AvgL2Score*70)/100
    #L1_NNL2<-(L1filter+0)*(NonNullL2score>2+0)*(L1score*30+AvgNN2Score*70)/100
    #if statement added 19.3.20 to adjust if the input weights are percentage e.g. 0.3 instead of 30
    if(L1weight<1&L2weight<1){
      L1weight<-L1weight*100
      L2weight<-L2weight*100
      L3weight<-L3weight*100
    }
    L1_avgL2<-(L1filter+0)*(NonNullL2score>2+0)*(L1score*L1weight+AvgL2Score*L2weight)/100
    L1_NNL2<-(L1filter+0)*(NonNullL2score>2+0)*(L1scoreO*L1weight+AvgNN2Score*L2weight)/100
    if(!is.null(L3scores)){
      #L1_NNL2L3<-(L1_NNL2*70+unlist(L3scores[names(L1_NNL2),"MaxNS"])*30)/100
      L1_NNL2L3<-(L1filter+0)*(NonNullL2score>2+0)*(L1score*L1weight+AvgNN2Score*L2weight+
                                                      as.numeric(L3scores[names(L1_NNL2),L3column])*L3weight)/100
    }else{
      L1_NNL2L3<-(L1filter+0)*(NonNullL2score>2+0)*(L1score*20+AvgNN2Score*60)/100
      L3scores<-matrix(0,nrow=length(L1_NNL2),ncol=1,dimnames=list(names(L1_NNL2),L3column))
    }

    L1_avgL2<-L1_avgL2[L1_avgL2>0]
    L1_NNL2<-L1_NNL2[L1_NNL2>0]
    L1_NNL2L3<-L1_NNL2L3[L1_NNL2L3>0]
    ##error handling added 25.9.19 as L2[names(L1_NNL2),] giving error ##
    genesinc<-intersect(rownames(L1),rownames(L2))
    if(length(genesinc)!=nrow(L1)){
      warnings('Some genes in L1 not in L2 and will be discarded')
    }
    if(length(genesinc)!=nrow(L2)){
      warnings('Some genes in L2 not in L1 and will be discarded')
    }
    L1_NNL2<-L1_NNL2[intersect(names(L1_NNL2),genesinc)]
    L1_NNL2L3<-L1_NNL2L3[intersect(names(L1_NNL2L3),genesinc)]
    ##end error handling
    if(subdir=="PairBD"){
      write.table(cbind(L1L3,L2),quote=FALSE,sep='\t',
                  file=paste(dir.Results,'/35_PrioritizedHits/',subdir,"/",omictype,"/",ctype,"/",ctype,mutationNumber,'_L1_L2.tsv',sep=''))


    }else{
      write.table(cbind(L1L3,L2),quote=FALSE,sep='\t',
                  file=paste(dir.Results,'/35_PrioritizedHits/',subdir,"/",ctype,"/",ctype,mutationNumber,'_L1_L2.tsv',sep=''))

    }
    FINAL_priority<-cbind(L1L3[names(L1_NNL2L3),],L2[names(L1_NNL2L3),],unlist(L3scores[names(L1_NNL2L3),L3column]))

    if(CTYPE!='PANCAN'){
      FINAL_priority<-FINAL_priority[,-11]

    }


    Output<-Get_TRACTinfo(subdir,FINAL_priority,tractability_both)
    TRACTinfos<-Output$TRACTinfos
    glist<-Output$glist

    GENEinfos<-OT15_retrievegeneInfo(glist)

    CIoutput<-Get_CompoundIndications(TRACTinfos,glist,tract_sm,tract_ab)
    compounds<-CIoutput$compounds
    indications<-CIoutput$indications

    FINAL_priority<-cbind(GENEinfos,TRACTinfos,compounds,indications,FINAL_priority)

    if(subdir=="PairBD"){
      write.table(FINAL_priority,quote=FALSE,sep='\t',
                  file=paste(dir.Results,'/35_PrioritizedHits/',subdir,"/",omictype,"/",ctype,"/",ctype,mutationNumber,'_FINAL_priority.tsv',sep=''))

      save(L1_avgL2,file=paste(dir.Results,'/35_PrioritizedHits/',subdir,"/",omictype,"/",ctype,"/",ctype,mutationNumber,'_L1_avgL2.RData',sep=''))
      save(L1_NNL2,file=paste(dir.Results,'/35_PrioritizedHits/',subdir,"/",omictype,"/",ctype,"/",ctype,mutationNumber,'_L1_avgNNL2.RData',sep=''))
      save(L1_NNL2L3,file=paste0(dir.Results,'/35_PrioritizedHits/',subdir,"/",omictype,"/",ctype,"/",ctype,mutationNumber,'_L1_avgNNL2L3.RData'))

    }else{
      write.table(FINAL_priority,quote=FALSE,sep='\t',
                  file=paste(dir.Results,'/35_PrioritizedHits/',subdir,"/",ctype,"/",ctype,mutationNumber,'_FINAL_priority.tsv',sep=''))

      save(L1_avgL2,file=paste(dir.Results,'/35_PrioritizedHits/',subdir,"/",ctype,"/",ctype,mutationNumber,'_L1_avgL2.RData',sep=''))
      save(L1_NNL2,file=paste(dir.Results,'/35_PrioritizedHits/',subdir,"/",ctype,"/",ctype,mutationNumber,'_L1_avgNNL2.RData',sep=''))
      save(L1_NNL2L3,file=paste0(dir.Results,'/35_PrioritizedHits/',subdir,"/",ctype,"/",ctype,mutationNumber,'_L1_avgNNL2L3.RData'))
    }

    Tout<-Get_TractableBuckets(subdir,tractability_both,glist,L1_NNL2,L1_NNL2L3)
    TRACTABLE<-Tout$TRACTABLE
    TRACTABLEL3<-Tout$TRACTABLEL3

    EFOIout<-Get_EFOIinfo(EFOI_list,CTYPE,tractability_both,tract_ab,tract_sm)
    approvedAntiCancerSpecific<-EFOIout$approvedAntiCancerSpecific
    approvedAntiCancer<-EFOIout$approvedAntiCancer
    approvedOthers<-EFOIout$approvedOthers
    nonApproved<-EFOIout$nonApproved
    if(subdir=="PairBD"){
      # gcheck<-unlist(lapply(glist,function(x) length(intersect(x,nonApproved))>0))
      missing<-setdiff(unlist(glist),rownames(tractability_both))
      if(length(missing)>0){
        warning('Some targets missing from tractability information, set to non Approved status')
        nonApproved<-c(nonApproved,missing)
      }


      return(list(AACS=L1_NNL2[unlist(lapply(glist,function(x) length(intersect(x,approvedAntiCancerSpecific))>0))],
                  AACO=L1_NNL2[unlist(lapply(glist,function(x) length(intersect(x,approvedAntiCancer))>0))],
                  AOD=L1_NNL2[unlist(lapply(glist,function(x) length(intersect(x,approvedOthers))>0))],
                  NAP=L1_NNL2[unlist(lapply(glist,function(x) length(intersect(x,nonApproved))>0))],
                  AACSL3=L1_NNL2L3[unlist(lapply(glist,function(x) length(intersect(x,approvedAntiCancerSpecific))>0))],
                  AACOL3=L1_NNL2L3[unlist(lapply(glist,function(x) length(intersect(x,approvedAntiCancer))>0))],
                  AODL3=L1_NNL2L3[unlist(lapply(glist,function(x) length(intersect(x,approvedOthers))>0))],
                  NAPL3=L1_NNL2L3[unlist(lapply(glist,function(x) length(intersect(x,nonApproved))>0))],
                  PRIORITIES=TRACTABLE, PRIORITIESL3=TRACTABLEL3))
    }else{
      missing<-setdiff(names(L1_NNL2),rownames(tractability_both))
      if(length(missing)>0){
        warning('Some targets missing from tractability information, set to non Approved status')
        nonApproved<-c(nonApproved,missing)
      }

      return(list(AACS=L1_NNL2[intersect(approvedAntiCancerSpecific,names(L1_NNL2))],
                  AACO=L1_NNL2[intersect(approvedAntiCancer,names(L1_NNL2))],
                  AOD=L1_NNL2[intersect(approvedOthers,names(L1_NNL2))],
                  NAP=L1_NNL2[intersect(nonApproved,names(L1_NNL2))],
                  AACSL3=L1_NNL2L3[intersect(approvedAntiCancerSpecific,names(L1_NNL2L3))],
                  AACOL3=L1_NNL2L3[intersect(approvedAntiCancer,names(L1_NNL2L3))],
                  AODL3=L1_NNL2L3[intersect(approvedOthers,names(L1_NNL2L3))],
                  NAPL3=L1_NNL2L3[intersect(nonApproved,names(L1_NNL2L3))],
                  PRIORITIES=TRACTABLE,PRIORITIESL3=TRACTABLEL3))
    }
  }else{
    return(NULL)
  }
}
Get_L1L3score<-function(L1,L3scores,componentNames=c("MutPrimTum","marker.Ttest.p...0.05..Class.D.","ANOVA.p...0.001..Class.C.",
                                                     "ANOVA.FDR...30...Class.B.","ANOVA.FDR...25...and.p...0.001..or.MSI.FDR...5...and.p...0.001..and.Glass.Deltas...1..Class.A.",
                                                     componentWeights=rep(20,5)),L3col="MaxTargetScoreFDR",BMclassCol="BMclass"
){
  if(!is.null(L3scores)){
    #reset L1 scores, where a network score avail:
    withL3score<-rownames(L3scores)[L3scores[,L3col]>0]
    L1[withL3score,componentNames[2:5]]<-0
    L3score<-rep(0,nrow(L1))
    names(L3score)<-rownames(L1)
    L3score[withL3score]<-L3scores[withL3score,L3col]
    #L3scores<-L3scores[withL3score,]
    #get biomarker classes from L3scores:

    classA<-rownames(L3scores)[L3scores[,BMclassCol]=="A"]
    classB<-rownames(L3scores)[L3scores[,BMclassCol]=="B"]
    classC<-rownames(L3scores)[L3scores[,BMclassCol]=="C"]
    classD<-rownames(L3scores)[L3scores[,BMclassCol]=="D"]
    L1[classA,componentNames[2:5]]<-1
    L1[classB,componentNames[2:4]]<-1
    L1[classC,componentNames[2:3]]<-1
    L1[classD,componentNames[2]]<-1
    L1<-cbind(L1,L3score)

  }else{
    L3score<-rep(0,nrow(L1))
    L1<-cbind(L1,L3score)
  }

  return(L1)
}

Get_BMPPIScores<-function(L3list,GENES,subdirlist,L3thresh=c(75,50,25,0),NetworkScore="FDRScore",L3names){
  nullres<-L3list[[1]][1,]
  nullres<-rep(0,length(nullres))
  names(nullres)<-names(L3list[[1]][1,])

  NetworkClasses<-matrix(0,nrow=length(GENES),ncol=length(subdirlist))
  sdlistcheck<-intersect(subdirlist,L3names)
  if(length(sdlistcheck)!=length(L3names)){

    stop("L3 priority list not in the subdir list, check inputs")
  }
  colnames(NetworkClasses)<-subdirlist
  rownames(NetworkClasses)<-GENES
  for(i in L3names){
    input<-L3list[[i]]
    if(!is.null(L3thresh)){
      NetworkClasses[,i]<-apply(input,1,function(y) sum(y[NetworkScore]>L3thresh)*25)
    }else{
      NetworkClasses[,i]<-input[,c(NetworkScore,"MaxTargetScoreFDR")]
    }

  }

  MaxNS<-apply(NetworkClasses,1,max)
  NetworkClasses<-cbind(NetworkClasses,MaxNS)
  OutMat<-list()
  for(i in 1:nrow(NetworkClasses)){
    mns<-NetworkClasses[i,"MaxNS"]
    if(mns!=0){
      maxset<-subdirlist[which(NetworkClasses[i,subdirlist]==mns)]
      target<-rownames(NetworkClasses)[i]
      outputs<-sapply(maxset,function(x) L3list[[x]][target,])

      if(length(maxset)>1){
        outputs<-apply(outputs,1,function(x) paste(x,collapse="||"))
      }
      OutMat[[i]]<-t(outputs)
    }else{
      OutMat[[i]]<-nullres
    }
  }
  Output<-do.call(rbind,OutMat)
  rownames(Output)<-rownames(NetworkClasses)

  L3scores<-cbind(NetworkClasses,Output,NetworkClasses[,"MaxTargetScoreFDR"])
  colnames(L3scores)[1]<-"MaxNS"
  colnames(L3scores)[ncol(L3scores)]<-"MaxTargetScoreFDR"
  return(L3scores)
}
Get_BMPPIScoresCombined<-function(CombinedData,GENES,subdirlist,L3thresh=c(75,50,25,0),NetworkScore="FDRScore",L3names){
  NetworkClasses<-matrix(0,nrow=length(GENES),ncol=4)

  colnames(NetworkClasses)<-c("FDRScore","MaxTargetScoreFDR","BMclass","RWRscore")
  rownames(NetworkClasses)<-GENES

  if(!is.null(L3thresh)){


  NetworkClasses[CombinedData$Target,1]<-apply(CombinedData,1,function(y) sum(y[NetworkScore]>L3thresh)*25)

  }else{
    NetworkClasses[CombinedData$Target,1]<-CombinedData[,"FDRScore"]
    NetworkClasses[CombinedData$Target,2]<-CombinedData[,"MaxTargetScoreFDR"]
    NetworkClasses[CombinedData$Target,3]<-CombinedData[,"BMclass"]
    NetworkClasses[CombinedData$Target,4]<-CombinedData[,"RWRscore"]
  }


  return(NetworkClasses)
}

charPriorityO<-function(ctype,dir.Results)
{

  CTYPE<-ctype

  L1<-read.csv(paste(dir.Results,'/33_L1/',ctype,'_L0_L1.csv',sep=''), header=TRUE,
               stringsAsFactors = FALSE,row.names = 1)
  load(paste(dir.Results,'/35_L2/',ctype,'_L2.Rdata',sep=''))
  usegenes<-intersect(rownames(L1),rownames(L2))
  L1<-L1[usegenes,]
  L2<-L2[usegenes,]
  #L1<-L1[rownames(L2),]
  print(head(L1))
  if(CTYPE!='PANCAN'){
    L1score<-(L1$MutPrimTum*20+
                L1$marker.Ttest.p...0.05..Class.D.*20+
                L1$ANOVA.p...0.001..Class.C.*20+
                L1$ANOVA.FDR...30...Class.B.*20+
                L1$ANOVA.FDR...25...and.p...0.001..or.MSI.FDR...5...and.p...0.001..and.Glass.Deltas...1..Class.A.*20)


  }else{
    L1score<-(L1$MutPrimTum*20+
                L1$marker.Ttest.p...0.05...Glass.Deltas...1..Class.D.*20+
                L1$ANOVA.p...0.001...at.least.one.Glass.Delta...1..Class.C.*20+
                L1$ANOVA.FDR...30....at.least.one.Glass.Delta...1..Class.B.*20+
                L1$ANOVA.FDR...25...and.p...0.001..or.MSI.FDR...5...and.p...0.001..and.Glass.Deltas...1..Class.A.*20)
  }

  names(L1score)<-rownames(L1)

  L1filter<-(rowSums(L1[,2:10]) == 0 & L1[,1]>1)
  #if they exist get the BMPPI scores as well:
  if(sum(c("Score","PvalScore","RWRscore","NumberBMs")%in%colnames(L1))==4){
    L3<-L1$Score
    names(L3)<-rownames(L1)
    L3<-L3[usegenes]
    #have BMPPI data:
    ScoreBoundaries<-seq(from=20,to=80,by=20)
    L3score<-sapply(L3,function(x) sum(x>=ScoreBoundaries)*25)
  }else{
    L3score<-NULL
  }
  ns<-ncol(L2)
  nc<-ns/15

  L2score<-NULL
  NonNullL2score<-NULL

  for (j in 1:nc){

    currentBlock<-L2[,(15*(j-1)+1):(15*j)]

    currentL2filter<-(currentBlock[,1]>0 & currentBlock[,2]!=1 & currentBlock[,3]!=1)
    if(site=="Combined"){

      currentL2score<-currentBlock[,c(6:8,12,13,15)]*100
      #currentL2score<-c( currentBlock[,6]*(1/6)*100,currentBlock[,7]*(1/6)*100,currentBlock[,8]*(1/6)*100,
      #            currentBlock[,12]*(1/6)*100,currentBlock[,13]*(1/6)*100,currentBlock[,15]*(1/6)*100)
      currentL2score<-rowSums(currentL2score,na.rm=T)/(rowSums(!is.na(currentL2score)))
      #temp 7.12.20 just to see effect of weights on new version:
      #currentL2score<-( currentBlock[,6]*(3/8)*100+currentBlock[,7]*(1/8)*100+currentBlock[,8]*(1/8)*100+
      #                    currentBlock[,12]*(1/8)*100+currentBlock[,13]*(1/8)*100+currentBlock[,15]*(1/8)*100)

    }else{
      currentL2score<-(currentBlock[,4]*12.5+currentBlock[,5]*12.5+
                         currentBlock[,6]*12.5+currentBlock[,7]*12.5+currentBlock[,8]*12.5+
                         currentBlock[,12]*12.5+currentBlock[,13]*12.5+currentBlock[,15]*12.5)
    }


    currentL2score<-currentL2score*(currentL2filter+0)

    L2score<-cbind(L2score,currentL2score)
  }

  NonNullL2score<-rowSums(L2score>0)

  AvgL2Score<-rowMeans(L2score)
  AvgNN2Score<-rowSums(L2score)/NonNullL2score
  AvgNN2Score[is.na(AvgNN2Score)]<-0

  L1_avgL2<-(L1filter+0)*(NonNullL2score>2+0)*(L1score*30+AvgL2Score*70)/100
  L1_NNL2<-(L1filter+0)*(NonNullL2score>2+0)*(L1score*30+AvgNN2Score*70)/100
  L1_NNL2L3<-(L1_NNL2*70+L3score[names(L1_NNL2)]*30)/100
  L1_avgL2<-L1_avgL2[L1_avgL2>0]
  L1_NNL2<-L1_NNL2[L1_NNL2>0]
  L1_NNL2L3<-L1_NNL2L3[L1_NNL2L3>0]

  #added error handling 25.9.19 because CERES throwing error on L2[names(L1_NNL2),]) below:
  genesinc<-intersect(rownames(L1),rownames(L2))
  L1_NNL2<-L1_NNL2[intersect(names(L1_NNL2),genesinc)]
  #end of error handling addition
  write.table(cbind(L1,L2),quote=FALSE,sep='\t',
              file=paste(dir.Results,'/35_PrioritizedHits/',ctype,'_L1_L2.tsv',sep=''))

  FINAL_priority<-cbind(L1[names(L1_NNL2),],L2[names(L1_NNL2),],L3score[names(L1_NNL2)])

  if(CTYPE!='PANCAN'){
    FINAL_priority<-FINAL_priority[,-11]
  }
  glist<-rownames(FINAL_priority)

  TRACTinfos<-tractability_both[rownames(FINAL_priority),]
  GENEinfos<-OT15_retrievegeneInfo(rownames(FINAL_priority))

  compounds<-rep(NA,nrow(TRACTinfos))
  indications<-rep(NA,nrow(TRACTinfos))

  idsm<-grep('small-molecule',TRACTinfos[,"type"])
  idab<-grep('antibody',TRACTinfos[,"type"])
  idboth<-grep("small-molecule/antibody",TRACTinfos[,"type"])


  compounds[idsm]<-tract_sm[glist[idsm],"drug_name"]
  compounds[idab]<-tract_ab[glist[idab],"drug_name"]
  compounds[idboth]<-paste(
    tract_sm[glist[idboth],"drug_name"],
    tract_ab[glist[idboth],"drug_name"],sep='|')



  indications[idsm]<-tract_sm[glist[idsm],"indication_efo_term"]
  indications[idab]<-tract_ab[glist[idab],"indication_efo_term"]
  indications[idboth]<-paste(
    tract_sm[glist[idboth],"indication_efo_term"],
    tract_ab[glist[idboth],"indication_efo_term"],sep='|')

  FINAL_priority<-cbind(GENEinfos,TRACTinfos,compounds,indications,FINAL_priority)


  write.table(FINAL_priority,quote=FALSE,sep='\t',
              file=paste(dir.Results,'/35_PrioritizedHits/',ctype,'_FINAL_priority.tsv',sep=''))

  save(L1_avgL2,file=paste(dir.Results,'/35_PrioritizedHits/',ctype,'_L1_avgL2.RData',sep=''))
  save(L1_NNL2,file=paste(dir.Results,'/35_PrioritizedHits/',ctype,'_L1_avgNNL2.RData',sep=''))

  uu<-sort(unique(as.numeric(tractability_both[,"min_bucket"])))

  TRACTABLE<-list()
  TRACTABLEL3<-list()

  for (j in 1:length(uu)){

    genes<-intersect(rownames(tractability_both)[which(tractability_both$min_bucket==uu[j])],
                     names(L1_NNL2))
    TRACTABLE[[j]]<-sort(L1_NNL2[genes],decreasing=TRUE)
    TRACTABLEL3[[j]]<-sort(L1_NNL2L3[genes],decreasing=TRUE)
  }
  names(EFOI_list)<-tolower(make.names(names(EFOI_list)))
  antiCancerSpecificIndication<-NULL
  if(CTYPE!='PANCAN'){

    Ctest<-tolower(make.names(CTYPE))
    antiCancerSpecificIndication<-EFOI_list[[Ctest]]
  }

  OtherAntiCancerIndication<-EFOI_list$anti.cancer.generic

  OtherAntiCancerTargets<-sort(union(rownames(tract_ab)[grep(paste(OtherAntiCancerIndication,collapse = '|'),tract_ab[,"indication_efo_term"])],
                                     rownames(tract_sm)[grep(paste(OtherAntiCancerIndication,collapse = '|'),tract_sm[,"indication_efo_term"])]))

  if(CTYPE!='PANCAN'& length(antiCancerSpecificIndication)>0){
    AntiCancerSpecificTargets<-sort(union(rownames(tract_ab)[grep(paste(antiCancerSpecificIndication,collapse = '|'),tract_ab[,"indication_efo_term"])],
                                          rownames(tract_sm)[grep(paste(antiCancerSpecificIndication,collapse = '|'),tract_sm[,"indication_efo_term"])]))
  } else{
    AntiCancerSpecificTargets<-NULL
  }

  OtherAntiCancerTargets<-setdiff(OtherAntiCancerTargets,AntiCancerSpecificTargets)

  approved<-
    rownames(tractability_both)[which(as.numeric(tractability_both[,"min_bucket"])<4)]

  approvedAntiCancerSpecific<-intersect(approved,AntiCancerSpecificTargets)
  approvedAntiCancer<-intersect(approved,OtherAntiCancerTargets)
  approvedOthers<-setdiff(approved,union(approvedAntiCancer,approvedAntiCancerSpecific))

  nonApproved<-
    rownames(tractability_both)[which(as.numeric(tractability_both[,"min_bucket"])>=4)]

  return(list(AACS=L1_NNL2[intersect(approvedAntiCancerSpecific,names(L1_NNL2))],
              AACO=L1_NNL2[intersect(approvedAntiCancer,names(L1_NNL2))],
              AOD=L1_NNL2[intersect(approvedOthers,names(L1_NNL2))],
              NAP=L1_NNL2[intersect(nonApproved,names(L1_NNL2))],
              AACSL3=L1_NNL2L3[intersect(approvedAntiCancerSpecific,names(L1_NNL2L3))],
              AACOL3=L1_NNL2L3[intersect(approvedAntiCancer,names(L1_NNL2L3))],
              AODL3=L1_NNL2L3[intersect(approvedOthers,names(L1_NNL2L3))],
              NAPL3=L1_NNL2L3[intersect(nonApproved,names(L1_NNL2L3))],
              PRIORITIES=TRACTABLE,PRIORITIESL3=TRACTABLEL3))

}

comboBoxes<-function(targets,CFEs,tissues,Mclass){

  ScatterData<-list()
  allCol<-NULL

  TissueTypeColors<-c(TissueTypeColors,'darkgray')
  names(TissueTypeColors)[length(TissueTypeColors)]<-'PANCAN'
  for (i in 1:length(CFEs)){
    currentScatterData<-defScatter(CFE = CFEs[i],
                                   target = targets[i],
                                   tissues= tissues[[i]])
    ScatterData<-c(ScatterData,currentScatterData,list(0))

    allCol<-c(allCol,
              c(rbind(TissueTypeColors[tissues[[i]]],TissueTypeColors[tissues[[i]]])),
              'black')

  }

  ScatterData<-rev(ScatterData)
  allCol<-rev(allCol)

  beeswarm(ScatterData,horizontal = TRUE,corral = 'wrap',xlim=c(range(unlist(ScatterData))),
           bg=makeTransparent(allCol),col=allCol,cex=3,pch=21,ylim=c(0,35))
  par(new=TRUE)
  boxplot(ScatterData,
          horizontal = TRUE,frame.plot=FALSE,ylim=c(range(unlist(ScatterData))),outline = FALSE,
          border =allCol,lwd=2,xlim=c(0,35),
          col=NA)

  locMclass<-rev(unlist(lapply(Mclass,function(x){c(c(t(cbind(x,rep('',length(x))))),'')})))

  text(min(unlist(ScatterData)),1:length(ScatterData),locMclass,cex=4)

  text(min(unlist(ScatterData)),which(allCol=='black'),paste(rev(targets),rev(CFEs),sep=' / '),cex=2,pos = 4)

}

defScatter<-function(CFE,target,tissues){


  RES<-list()

  for (i in 1:length(tissues)){
    if(tissues[i]=='PANCAN'){
      CLs<-rownames(ESSprofiles)
    }else{
      tissueSpec_CLs<-intersect(unique(as.character(manifest$COSMIC_ID[which(is.element(manifest$Tissue,tissues[i]) |
                                                                               is.element(manifest$Cancer.Type,tissues[i]))])),
                                rownames(ESSprofiles))
      CLs<-tissueSpec_CLs
    }

    sensPat<-ESSprofiles[CLs,target]

    if(CFE=='MSI_Status'){
      mutPat<-InputFeatures$MSI_VARIABLE[CLs]
    }else{
      mutPat<-InputFeatures$BEM[CFE,CLs]
    }

    RES<-c(RES,list(negPop=sensPat[mutPat==0],posPop=sensPat[mutPat==1]))
  }



  return(RES)
}

get_PriorityRes<-function(MARKERclass,ctype,mutationNumber,TRACTABLE,PCHSYM,indications,PRIORITY_vectors,PRIORITY_vectorsL3,th,sigOnly){
  XX<-NULL
  YY<-NULL
  NN<-NULL
  YYL3<-NULL
  ALLMARK<-NULL
  ALLPCH<-NULL
  ALLBUCK<-NULL
  currentMarkers<-MARKERclass[which(MARKERclass$ANALYSIS==paste0(ctype,mutationNumber)),]


  #TRACTABLE are priority scores for dependencies across different tractability buckets
  uu<-1:10


  for (j in 1:length(uu)){

    if(sigOnly){
      genes<-names(which(TRACTABLE[[j]]>=th))
    }else{
      genes<-names(TRACTABLE[[j]])
    }

    ngenes<-length(genes)
    if(ngenes>0){
      if(ctype=="PANCAN"){

        symbols<-PCHSYM[colSums(do.call(rbind,lapply(indications[2:4],function(x){is.element(genes,names(x))}))*2:4)]
      }else{
        symbols<-PCHSYM[colSums(do.call(rbind,lapply(indications,function(x){is.element(genes,names(x))}))*1:4)]}

      if(ngenes==1){
        xc<-(uu[j])
      }else{
        xc<-seq(uu[j]-0.4,uu[j]+0.4,0.8/(ngenes-1))
      }



      if(ctype!='PANCAN'){
        PP<-PRIORITY_vectors[[paste0(ctype,mutationNumber)]][genes]
        inL3<-intersect(names(PRIORITY_vectorsL3),genes)

        noL3<-setdiff(genes,names(PRIORITY_vectorsL3))
        if(length(noL3)>0){
          add<-rep(0,length(noL3))
          names(add)<-noL3
          PPL3<-c(PRIORITY_vectorsL3[inL3],add)
        }else{
          PPL3<-PRIORITY_vectorsL3[genes]}

      }else{
        PP<-PRIORITY_vectors[genes]

        inL3<-intersect(names(PRIORITY_vectorsL3),genes)

        noL3<-setdiff(genes,names(PRIORITY_vectorsL3))
        if(length(noL3)>0){
          add<-rep(0,length(noL3))
          names(add)<-noL3
          PPL3<-c(PRIORITY_vectorsL3[inL3],add)
        }else{
          PPL3<-PRIORITY_vectorsL3[genes]}

      }

      ad<-names(sort(PP,decreasing=TRUE))

      indi<-do.call(rbind,lapply(indications[1:3],function(x){is.element(ad,names(x))}))+0
      colnames(indi)<-genes

      genes<-names(sort(PP,decreasing=TRUE))

      MARKERclass<-NULL
      for (k in 1:length(genes)){
        #this alphabetically sorts the marker class for all biomarkers of selected gene dependency
        currentC<-sort(unlist(currentMarkers[which(currentMarkers$Depleted.Gene==genes[k]),'CLASS']))[1]
        if(length(currentC)==0){currentC<-'NA'}
        MARKERclass<-c(MARKERclass,currentC)
      }
      #removed a final NA from the newSymbols vector 22.8.20
      newSymbols<-c(8,3,4,2)
      names(newSymbols)<-c('A','B','C','D')
      MARKERclass<-newSymbols[MARKERclass]
      names(MARKERclass)<-genes

      ALLMARK<-c(ALLMARK,MARKERclass)
      XX<-c(XX,xc)
      YY<-c(YY,sort(PP,decreasing=TRUE))
      YYL3<-c(YYL3,PPL3[order(PP,decreasing=TRUE)])
      NN<-c(NN,names(sort(PP,decreasing=TRUE)))
      ALLPCH<-c(ALLPCH,symbols[order(PP,decreasing=TRUE)])
      ALLBUCK<-c(ALLBUCK,rep(j,length(xc)))
    }
  }



  priorSet<-unlist(lapply(TRACTABLE,function(x){names(which(x>=th))[1:3]}))
  priorSet<-priorSet[!is.na(priorSet)]

  if(length(priorSet)){

    ii<-match(priorSet,NN)

    # text(XX[ii]+0.4,YY[ii],NN[ii])
  }


  tmp<-lapply(TRACTABLE,function(x){
    ht<-names(which(x>=th))
    scores<-x[which(x>=th)]
  })
  rescheck<-sum(unlist(lapply(tmp,length)))
  if(rescheck>0){
    res<-do.call('rbind',lapply(1:10,function(i){cbind(rep(i,length(tmp[[i]])),tmp[[i]])}))

    colnames(res)<-c('bucket','priority')


    RES<-data.frame(ctype=paste0(ctype,mutationNumber),TARGET=NN,BUCKET=ALLBUCK,
                    PRIORITY=as.numeric(YY),
                    MARKERsymbol=ALLMARK,
                    indiSymbol=ALLPCH,
                    PRIORITYL3=as.numeric(YYL3),
                    stringsAsFactors = FALSE)
  }else{
    RES<-NULL
  }
  return(RES)
}
priorityPlot<-function(ctype,PAB,th,sigOnly=TRUE,IDENTIFY=FALSE,priorSet=NULL,MARKERclass,TissueTypeColors,subdir="",mutationNumber="",PAI=NULL,Pvectors=NULL,PvectorsL3=NULL){


  TissueColors<-TissueTypeColors


    if(!is.null(PAI)){
      PRIORITISED_across_indications<-PAI
    }else{
      if(ctype=="PANCAN"){
        load(paste0(dir.Results,'/35_PrioritizedHits/',subdir,'/_PRIORITISED_across_indications_PANCAN.RData'))
      }else{
        load(paste0(dir.Results,'/35_PrioritizedHits/',subdir,'/_PRIORITISED_across_indications.RData'))
      }
    }
    if(!is.null(Pvectors)){
      PRIORITY_vectors<-Pvectors
    }else{
      if(ctype=="PANCAN"){
        load(paste0(dir.Results,'/35_PrioritizedHits/',subdir,'/_PRIORITY_vectors_PANCAN.RData'))
      }else{
        load(paste0(dir.Results,'/35_PrioritizedHits/',subdir,'/_PRIORITY_vectors.RData'))
      }
    }

    if(!is.null(PvectorsL3)){
      if(ctype=="PANCAN"){
        PRIORITY_vectorsL3<-PvectorsL3
      }else{
        PRIORITY_vectorsL3<-PvectorsL3[[paste0(ctype,mutationNumber)]]
      }

    }else{
      if(ctype=="PANCAN"){
        load(paste0(dir.Results,'/35_PrioritizedHits/',subdir,'/_PRIORITY_vectorsL3_PANCAN.RData'))
      }else{
        if(file.exists(paste0(dir.Results,'/35_PrioritizedHits/',subdir,'/_PRIORITY_vectorsL3.RData'))){
          load(paste0(dir.Results,'/35_PrioritizedHits/',subdir,'/_PRIORITY_vectorsL3.RData'))

          PRIORITY_vectorsL3<-PRIORITY_vectorsL3[[paste0(ctype,mutationNumber)]]
        }else{
          PRIORITY_vectorsL3<-rep(0,length(PRIORITY_vectors[[paste0(ctype,mutationNumber)]]))
          names(PRIORITY_vectorsL3)<-names(PRIORITY_vectors[[paste0(ctype,mutationNumber)]])
        }
      }

    }




  TissueColors_tr<-makeTransparent(TissueColors,100)

  if(paste0(ctype,mutationNumber)%in%names(PRIORITISED_across_indications)){
    #changes 15/4/20
    indications<-PRIORITISED_across_indications[[paste0(ctype,mutationNumber)]]
  }else{
    indications<-PRIORITISED_across_indications
  }

  PCHSYM<-c(23,24,22,21)
  names(PCHSYM)<-names(indications)



  if(paste0(ctype,mutationNumber)%in%names(PRIORITISED_across_indications)){
    TRACTABLE<-PAB[[paste0(ctype,mutationNumber)]]
  }else{
    TRACTABLE<-PAB
  }



  #if(sigOnly){
   # YL<-c(th-1,80)
  #}else{
   # YL<-c(1,80)
  #}



  RES<-get_PriorityRes(MARKERclass,ctype,mutationNumber,TRACTABLE,PCHSYM,indications,PRIORITY_vectors,PRIORITY_vectorsL3,th,sigOnly=sigOnly)


  return(RES)
}

AnnotInfo<-function(TOTRES,tractability_both){
  GROUP<-TOTRES$BUCKET

  GROUP[which(TOTRES$BUCKET>7)]<-3
  GROUP[which(TOTRES$BUCKET>3 & TOTRES$BUCKET<=7)]<-2
  GROUP[which(TOTRES$BUCKET<=3)]<-1

  MARKERCLASS<-rep('N/A',length(TOTRES$MARKERsymbol))


  MARKERCLASS[which(TOTRES$MARKERsymbol==8)]<-'A'
  MARKERCLASS[which(TOTRES$MARKERsymbol==3)]<-'B'
  MARKERCLASS[which(TOTRES$MARKERsymbol==4)]<-'C'
  MARKERCLASS[which(TOTRES$MARKERsymbol==2)]<-'D'

  INDICATION<-rep('N/A',length(TOTRES$indiSymbol))
  INDICATION[TOTRES$indiSymbol==23]<-'Anti-Cancer Specific'
  INDICATION[TOTRES$indiSymbol==24]<-'other Anti-Cancer'
  INDICATION[TOTRES$indiSymbol==22]<-'other Disease'
  INDICATION[TOTRES$indiSymbol==21]<-'no compound available'

  MOLTYPE<-rep('N/A',length(TOTRES$BUCKET))
  MOLTYPE<-tractability_both[TOTRES$TARGET,"type"]

  return(list(GROUP=GROUP,MARKERCLASS=MARKERCLASS,INDICATION=INDICATION,MOLTYPE=MOLTYPE))
}

chunk <- function(x,n) split(x, factor(sort(rank(x)%%n)))
Annotate_PPIdistance<-function(allMarkers,subdir,ncores,PPIigraph,PPInet){
  if(subdir%in%c("","Cont","CExpr","CMet","CCN","CProt","all","Binary")){
    if("BMstringID"%in%colnames(allMarkers)&"string_id"%in%colnames(allMarkers)){
      if(nrow(allMarkers)==1){
      anovaMarkers<-matrix(unlist(allMarkers[,c("FEATURE","BMstringID","string_id")]),ncol=3,byrow=TRUE)
      colnames(anovaMarkers)<-c("FEATURE","BMstringID","string_id")
      }else{
        anovaMarkers<-allMarkers[,c("FEATURE","BMstringID","string_id")]
      }
      if(nrow(anovaMarkers)>5*ncores){
        nMarkers<-nrow(anovaMarkers)

        idxs<-chunk(1:nMarkers,ncores)

        test<-foreach(i=1:ncores,.combine=cbind)%dopar%{

          apply(anovaMarkers[idxs[[i]],],1,function(x) GetSPFromIDs(x,PPIigraph,subdir=subdir))
        }
      }else{
        test<-apply(anovaMarkers,1,function(x) GetSPFromIDs(x,PPIigraph,subdir=subdir))
      }
    }else{
      anovaMarkers<-matrix(allMarkers[,c("Depleted.Gene","FEATURE")],ncol=2,byrow=TRUE)
      colnames(anovaMarkers)<-c("Depleted.Gene","FEATURE")
      nMarkers<-nrow(anovaMarkers)

      if(nMarkers>5*ncores){
        idxs<-chunk(1:nMarkers,ncores)

        test<-foreach(i=1:ncores,.combine=cbind)%dopar%{
          inputdata<-anovaMarkers[idxs[[i]],]
          apply(inputdata,1,function(x) GetSP(PPInet,x,PPIigraph,subdir=subdir))
        }
      }else{
        test<-apply(anovaMarkers,1,function(x) GetSP(PPInet,x,PPIigraph,subdir=subdir))
      }
    }
    ppis<-test["ppi",]
    markers<-test["markers",]
  }


  if(subdir=="Interaction"){
    anovaMarkers<-allMarkers[,c("Depleted.Gene","FEATURE_expression")]
    nMarkers<-nrow(allMarkers)

    idxs<-chunk(1:nMarkers,ncores)
    test<-foreach(i=1:ncores,.combine=cbind)%dopar%{
      inputdata<-anovaMarkers[idxs[[i]],]
      apply(inputdata,1,function(x) GetSP(PPInet,x,PPIigraph,subdir=subdir))
    }
    ppis<-test["ppi",]
    markers<-test["markers",]


  }

  if(subdir%in%c("PairBD","Compound","CompoundME","ME")){
    return(list(ppi_target=rep(NA,nrow(allMarkers)),markers=rep(NA,nrow(allMarkers))))

  }else{
    return(list(ppi_target=ppis,markers=markers))


  }

}

collapseMarkers<-function(TOTRES,allMarkers,PPIigraph,PPInet,gmtfile,subdir="",ifPPI=TRUE,BMPPIdir=NULL,CommonEssentials,ncores,AInfo){
  MARKER<-rep('N/A',length(TOTRES$MARKERCLASS))
  ASSOCIATION_EFFECT<-rep('N/A',length(TOTRES$MARKERCLASS))
  ANOVA_table_entry<-rep('N/A',length(TOTRES$MARKERCLASS))
  PPI_distance<-rep('N/A',length(TOTRES$MARKERCLASS))
  PPI_info<-rep('N/A',length(TOTRES$MARKERCLASS))
  PPI_enrich<-rep('N/A',length(TOTRES$MARKERCLASS))
  PPI_enrichFn<-rep('N/A',length(TOTRES$MARKERCLASS))
  MARKER_TYPE<-rep('N/A',length(TOTRES$MARKERCLASS))
  PPI_min<-rep('N/A',length(TOTRES$MARKERCLASS))
  PPI_max<-rep('N/A',length(TOTRES$MARKERCLASS))
  CE<-rep('N/A',length(TOTRES$MARKERCLASS))
  RWRscore<-rep('N/A',length(TOTRES$MARKERCLASS))
  NetworkScore<-rep('N/A',length(TOTRES$MARKERCLASS))
  IQRscore<-rep('N/A',length(TOTRES$MARKERCLASS))
  #targetsNF<-allMarkersNF[,"Depleted.Gene"]
  CE<-rep("FALSE",length(TOTRES$MARKERCLASS))
  names(CE)<-TOTRES$TARGET
  CE[which(names(CE)%in%CommonEssentials)]<-"TRUE"

  if(!is.null(BMPPIdir)){
    RWRcombined<-read.table(paste0(BMPPIdir,"/RWR_ScoreTargetsCombined.txt"),header=T,stringsAsFactors = FALSE,sep="\t")
    RWRcombined<-RWRcombined[RWRcombined[,"ScoreType"]=="Avg",]
    #rownames(RWRcombined)<-RWRcombined$IDS
  }else{
    RWRcombined<-vector(mode='list',length=1)
    names(RWRcombined)<-"IDS"
  }
  for (i in 1:nrow(TOTRES)){
    idR<-paste0(TOTRES[i,"TARGET"],TOTRES[i,"ctype"])
    if(idR%in%RWRcombined$IDS){
      sel<-which(RWRcombined$IDS==idR)
      RWRbiomarkers<-RWRcombined[sel,"BMSymbol"]
      RWRtype<-RWRcombined[sel,"BMType"]
      RWRclass<-RWRcombined[sel,"BMclass"]
      NetworkScore[i]<-RWRcombined[sel,"MaxTargetScoreFDR"]
      RWRscore[i]<-RWRcombined[sel,"RWRscore"]
    }else{
      RWRbiomarkers<-""
      RWRtype<-""
      RWRclass<-""
      RWRscore[i]<-0
      NetworkScore[i]<-0
    }
    if(TOTRES$MARKERCLASS[i]!='N/A'){
      #this is where only markers of the maximum class (TOTRES$MARKERCLASS) for depleted gene meeting tractability thresholds are selected

        id<-which(allMarkers$ANALYSIS==TOTRES$ctype[i] &
                    allMarkers$Depleted.Gene==TOTRES$TARGET[i] &
                    allMarkers$CLASS==TOTRES$MARKERCLASS[i])
        if(RWRclass!=""){
          id2<-which(allMarkers$ANALYSIS==TOTRES$ctype[i] &
                       allMarkers$Depleted.Gene==TOTRES$TARGET[i] &
                       allMarkers$CLASS==RWRclass)
        }else{id2<-NULL}
        if(subdir=="Interaction"){
          currentMar<-allMarkers$FEATURE_expression[id]
          currentMar2<-allMarkers$FEATURE_expression[id2]
        }
        if(subdir%in%c("Compound","Cell","CompoundME")){
          currentMar<-allMarkers$FEATURE[id]
          currentMar2<-allMarkers$FEATURE[id2]
        }

        if(!subdir%in%c("Interaction","Compound","Cell","CompoundME")){
          currentMar<-allMarkers$FEATURE[id]
          currentMar2<-allMarkers$FEATURE[id2]

        }

        if(RWRbiomarkers==""){
          MARKER[i]<-paste(currentMar,collapse=' // ')
        }else{
          output<-Get_MarkerIDs(currentMar,currentMar2,id,id2,RWRbiomarkers)
          currentMar<-output$currentMar
          id<-output$id
          MARKER[i]<-paste(currentMar,collapse="//")
          TOTRES$MARKERCLASS[i]<-RWRclass
        }
      #get marker types:
      markersummary<-getMarkerSummary(currentMar)
      markertype<-getMarkerType(markersummary)



      #get expression markers for enrichment
      expressionMarkers<-unique(getSingleExpression(currentMar,subdir))
      expressionMarkers<-PPImarkerCheck(expressionMarkers,PPInet)

      if(length(expressionMarkers)>2){
        expressionEnrich<-tryCatch(ppi_enrichment(expressionMarkers,PPIigraph)$enrichment,error=function(e){return('N/A')})
        }else{
        expressionEnrich<-'N/A'
        }

      if(exists("expressionEnrich")){
        if(length(expressionEnrich)>0){

          if(expressionEnrich!='N/A'){
            hyperres<-lapply(gmtfile,function(x) testGmt(x,expressionMarkers,5000))

            adjhyper<-p.adjust(unlist(hyperres))
            tophyperres<-adjhyper[adjhyper<0.05]
            if(length(tophyperres)>5){
              tophyperres<-sort(tophyperres)[1:5]
            }
            if(length(tophyperres)==1){
              tophyperres<-'NoSig'
            }
          }else{
            tophyperres<-'NoEnrich'
          }
         }else{expressionEnrich<-'N/A'}
      }else{tophyperres<-"NoEnrichT"}

      if(ifPPI){
        ##added in the PPI information 6.8.19

        if(length(id)>1){
          inputdata<-allMarkers[id,]
        }else{

          inputdata<-matrix(unlist(allMarkers[id,]),nrow=1,ncol=ncol(allMarkers))
          colnames(inputdata)<-colnames(allMarkers)
        }

        currentPPIDist<-Annotate_PPIdistance(inputdata,subdir,ncores,PPIigraph,PPInet)$ppi_target
        #currentPPIDist<-allMarkers$ppi_target[id]
        numericPPI<-getPPIsummary(currentPPIDist)
        PPI_min[i]<-numericPPI$minM
        PPI_max[i]<-numericPPI$maxM
        #currentPPIinfo<-allMarkers$ppi_target_detailed[id]
        PPI_distance[i]<-paste(unlist(currentPPIDist),collapse='//')
      }
      #PPI_info[i]<-paste(currentPPIinfo,collapse=" // ")
      #changed from FEATURE_pos to FEATURE_delta on 4.8.19
      ASSOCIATION_EFFECT[i]<-paste(allMarkers[id,"depAssoc"],collapse='//')
      ANOVA_table_entry[i]<-paste(allMarkers$assoc_id[id],collapse='//')
      IQRscore[i]<-paste(allMarkers[id,"IQR"],collapse = "//")
      #CE[i]<-allMarkers$CE[id[1]]
      if(length(expressionEnrich)>0){
        PPI_enrich[i]<-expressionEnrich
        PPI_enrichFn[i]<-paste(tophyperres,collapse="//")
      }
      MARKER_TYPE[i]<-markertype
    }else{
      MARKER[i]<-'N/A'
      ASSOCIATION_EFFECT[i]<-'N/A'
      ANOVA_table_entry[i]<-'N/A'
      PPI_distance[i]<-'N/A'
      PPI_info[i]<-'N/A'
      PPI_enrich[i]<-'N/A'
      PPI_enrichFn[i]<-'N/A'
      MARKER_TYPE[i]<-'N/A'
      PPI_min[i]<-'N/A'
      PPI_max[i]<-'N/A'
      #CE[i]<-'N/A'
      RWRscore[i]<-'N/A'
      IQRscore[i]<-'N/A'
      NetworkScore[i]<-'N/A'
    }

  }
  rownames(TOTRES)<-NULL
  names(CE)<-NULL
  TRACTABILITY<-TOTRES[,3]

  #GLOBAL<-cbind(TOTRES[,c(1,2,4)],GROUP,TRACTABILITY,MARKERCLASS,ASSOCIATION_EFFECT,MARKER,ANOVA_table_entry,INDICATION,PPI_info,PPI_distance)
  GLOBAL<-cbind(TOTRES[,c("ctype","TARGET","PRIORITY","PRIORITYL3")],GROUP=AInfo$GROUP,TRACTABILITY,MOLTYPE=AInfo$MOLTYPE,MARKERCLASS=TOTRES$MARKERCLASS, ASSOCIATION_EFFECT,MARKER,
                ANOVA_table_entry,INDICATION=AInfo$INDICATION,PPI_distance,PPI_enrich,PPI_enrichFn,MARKER_TYPE=MARKER_TYPE,PPI_min=PPI_min,CE=CE,RWRscore=RWRscore,NetworkScore=NetworkScore,IQR=IQRscore)

  return(GLOBAL)
}

Get_ScoreComponents<-function(collapseMarkersInfo,dir.Results,subdir,ctype,mutationNumber,omictype="Binary"){
  L1filename<-Get_L1filename(subdir,mutationNumber,ctype,omictype,dir.Results)

  L1<-read.csv(L1filename, header=TRUE,
               stringsAsFactors = FALSE,row.names = 1)
  collapseMarkersInfo<-collapseMarkersInfo[collapseMarkersInfo$ctype==ctype,]
  rownames(collapseMarkersInfo)<-collapseMarkersInfo$TARGET
  L1L3<-Get_L1L3score(L1,collapseMarkersInfo,L3col="NetworkScore",BMclassCol = "MARKERCLASS")

  write.csv(L1L3,file=paste0(dir.Results,"/",subdir,"/",ctype,"_L1L3",subdir,mutationNumber,".csv"),quote=FALSE,row.names=TRUE)
}

Get_PriorityScore<-function(dir.Results,subdir,ctype,mutationNumber,L1weight,L2weight,L3weight,site="Combined",L3column="L3score"){
  L1L3<-read.csv(file=paste0(dir.Results,"/",subdir,"/",ctype,"_L1L3",subdir,mutationNumber,".csv"),stringsAsFactors = FALSE,row.names=1)
  L2<-read.csv(paste(dir.Results,'/35_L2/',ctype,'_L2.csv',sep=''),row.names =1,stringsAsFactors = FALSE)
  usegenes<-intersect(rownames(L1L3),rownames(L2))
  L1L3<-L1L3[usegenes,]
  L2<-L2[usegenes,]

  L1score<-Get_L1score(L1L3)
  names(L1score)<-rownames(L1L3)
  L1filter<-(rowSums(L1L3[,2:10]) == 0 & L1L3[,1]>1)

  L2score<-Get_L2score(L2,site)
  NonNullL2score<-NULL
  NonNullL2score<-rowSums(L2score>0)

  AvgL2Score<-rowMeans(L2score)
  AvgNN2Score<-rowSums(L2score)/NonNullL2score
  AvgNN2Score[is.na(AvgNN2Score)]<-0



  L1_NNL2L3<-(L1filter+0)*(NonNullL2score>2+0)*(L1score*L1weight+AvgNN2Score*L2weight+
                                                  as.numeric(L1L3[,L3column])*L3weight)/100

  return(L1_NNL2L3)
}

Get_MarkerIDs<-function(currentMar,currentMar2,id,id2,RWRbiomarkers){
  Mmet<-currentMar[grep("_Met",currentMar)]
  CMS<-currentMar[grep("CMS",currentMar)]
  Breastsub<-currentMar[grep("ScMod",currentMar)]
  Breastsub2<-currentMar[grep('PAM50',currentMar)]
  Mcomposite<-currentMar[grep("_Composite",currentMar)]
  Mcomposite<-Mcomposite[grep('[&|]',Mcomposite,invert=T)]

  Mmet2<-currentMar2[grep("_Met",currentMar2)]
  CMS2<-currentMar2[grep("CMS",currentMar2)]
  Breastsub2<-currentMar2[grep("ScMod",currentMar2)]
  Breastsub22<-currentMar2[grep('PAM50',currentMar2)]
  Mcomposite2<-currentMar2[grep("_Composite",currentMar2)]
  Mcomposite2<-Mcomposite2[grep('[&|]',Mcomposite2,invert=T)]
  RWRbiomarkers<-RWRbiomarkers
  id2.1<-NULL
  m2<-NULL
  m3<-NULL
  for(i in 1:length(currentMar2)){
    x<-currentMar2[i]

    out<-unlist(sapply(unlist(strsplit(x,", ",fixed=TRUE)),function(z)
              z[toupper(z)%in%toupper(unlist(strsplit(RWRbiomarkers,"//")))]))
    if(length(out)>0){
      id2.1<-c(id2.1,id2[i])
      if(length(out)>1){
        m2<-c(m2,paste(out,collapse=", "))
      }else{
        m2<-c(m2,out)
      }
    }
    out2<-unlist(sapply(unlist(strsplit(x,", ",fixed=TRUE)),function(z)
      z[toupper(paste0(unlist(sapply(unlist(strsplit(z,"[&|]")),
                             function(p) strsplit(p,'_',fixed=TRUE)[[1]][1]))[2:1],collapse=":"))%in%toupper(gsub("[&|]",":",unlist(strsplit(RWRbiomarkers,"//"))))]))

    if(length(out2)>0){
      id2.1<-c(id2.1,id2[i])
      if(length(out2)>1){
        m3<-c(m3,paste(out2,collapse=", "))
      }else{
        m3<-c(m3,out2)
      }
    }
    out3<-unlist(sapply(unlist(strsplit(x,", ",fixed=TRUE)),function(z)
      z[toupper(strsplit(z,"_",fixed=TRUE)[[1]][1])%in%toupper(unlist(strsplit(RWRbiomarkers,"//")))]))
    if(length(out3)>0){
      id2.1<-c(id2.1,id2[i])
      if(length(out3)>1){
          m3<-c(m3,paste(out3,collapse=", "))
        }else{
          m3<-c(m3,out3)
        }
    }
  }


  m2<-unique(c(m2,m3))
  id2.1<-unique(id2.1)


  m1<-currentMar[currentMar%in%c(Mcomposite,Mmet,CMS,Breastsub,Breastsub2)]
  id1<-id[currentMar%in%c(Mcomposite,Mmet,CMS,Breastsub,Breastsub2)]
  m1.2<-currentMar2[currentMar2%in%c(Mcomposite2,Mmet2,CMS2,Breastsub2,Breastsub22)]
  id1.2<-id2[currentMar2%in%c(Mcomposite2,Mmet2,CMS2,Breastsub2,Breastsub22)]

  id1<-unique(c(id1,id2.1))
  id<-unique(c(id1,id1.2))
  m1<-unique(c(m1,m2))
  currentMar<-unique(c(m1,m1.2))
  id<-id[!is.na(id)]
  currentMar<-currentMar[!is.na(currentMar)]
  return(list(currentMar=currentMar,id=id))
}
combineRes<-function(r1,r2){
  allvals<-unique(c(names(r1),names(r2)))
  dmax<-matrix(NA,nrow=length(allvals),ncol=2)
  rownames(dmax)<-allvals
  dmax[names(r1),1]<-r1
  dmax[names(r2),2]<-r2
  return(apply(dmax,1,function(x) max(x,na.rm=T)))
}

PPImarkerCheck<-function(markers,PPInet){
  stringM<-PPInet[markers,"STRING_id"]
  stringM[!is.na(stringM)]
  return(stringM)
}
getSingleExpression<-function(currentMar,subdir=""){
  singlemarkers<-sapply(currentMar,function(x) strsplit(x,",",fixed=TRUE))
  if(subdir=="EMDE"){
    genes<-singlemarkers
  }else{
    expressionmarkers<-unlist(lapply(singlemarkers,function(x) x[grep("_Expr",x)]))
    genes<-sapply(expressionmarkers,function(x) strsplit(x,"_",fixed=TRUE)[[1]][1])
  }
  return(unlist(genes))
}

getMarkerSummary<-function(markers){
  allmarkers<-NULL
  for(i in 1:length(markers)){
    splitmarkers<-unlist(strsplit(markers[i],",",fixed=TRUE))
    expr<-grep("_Expr",splitmarkers)
    mut<-grep("_mut",splitmarkers)
    cna<-grep("gain",splitmarkers)
    cna2<-grep("loss",splitmarkers)
    cna3<-grep("_CN",splitmarkers)
    met<-grep("_Met",splitmarkers)
    prot<-grep("_Prot",splitmarkers)
    hypmet<-grep("_HypMET",splitmarkers)
    variant<-grep("_var",splitmarkers)
    composite<-grep("_Composite",splitmarkers)
    dual<-grep("[&|]",splitmarkers)
    allmarkers<-rbind(allmarkers,c(expr=length(expr),mut=length(mut),cna=length(cna)+length(cna2)+length(cna3),
                                   hypMet=length(hypmet),met=length(met),prot=length(prot),variant=length(variant),
                                   composite=length(composite),dual=length(dual)))
  }
  rownames(allmarkers)<-markers
  return(allmarkers)
}

getMarkerType<-function(markersummary){
  markercount<-colSums(markersummary)
  names(markercount)<-c("expr","mut","cna","hypMet","met","prot","variant","composite","dual")
  markerexist<-markercount!=0
  markertype<-paste0(names(markercount)[markercount!=0],collapse = ",")
  return(markertype)
}

getPPIsummary<-function(markers){
  aM<-NULL
  for(i in 1:length(markers)){
    if(is.character(markers[i])){
    splitmarkers<-unlist(strsplit(markers[i],"[,-]"))
    if(is.list(splitmarkers)){
    numericMarkers<-suppressWarnings(as.numeric(splitmarkers[[1]]))}else{
      numericMarkers<-suppressWarnings(as.numeric(splitmarkers))
    }
    aM<-c(aM,numericMarkers)}else{
      #for trying to debug:
      print(markers[i])
    }

  }

  return(list(minM=min(aM,na.rm=T),maxM=max(aM,na.rm=T)))
}

collapseTargets<-function(TOTRES,allMarkers,PPIigraph,PPInet,AnnotInfo=NULL,gmtfile=NULL,subdir="",splitmarkers=TRUE){
  #TOTRES has all the priority targets and ctypes that pass the priority score threshold

  print(colnames(allMarkers))
  markersuse<-NULL
  for (i in 1:nrow(TOTRES)){

    if(TOTRES$MARKERCLASS[i]!='N/A'){
      #this is where only markers of the maximum class (TOTRES$MARKERCLASS) for depleted gene meeting tractability thresholds are selected

      id<-which(allMarkers$ANALYSIS==TOTRES$ctype[i] &
                  allMarkers$Depleted.Gene==TOTRES$TARGET[i] &
                  allMarkers$CLASS==TOTRES$MARKERCLASS[i])

      markersuse<-rbind(markersuse,allMarkers[id,])
    }

  }
  if(splitmarkers){
    #now need to split any compound/interaction or markers:
    if(subdir=="Interaction"){
      features<-unlist(sapply(markersuse$FEATURE_expression,function(x) changeCNAsep(x)))

      }
    if(subdir%in%c("Compound","Cell","CompoundME")){
      features<-unlist(sapply(markersuse$FEATURE_child,function(x) changeCNAsep(x)))

    }
    if(!subdir%in%c("Interaction","Compound","Cell","CompoundME")){
      features<-unlist(sapply(markersuse$FEATURE,function(x) changeCNAsep(x)))

    }
    compoundmarkers<-grep(",",features,fixed=TRUE)
    alldf<-NULL
    for(i in 1:length(compoundmarkers)){
      sel<-compoundmarkers[i]
      morig<-markersuse[sel,]

      individm<-unlist(strsplit(features[sel],",",fixed=TRUE))
      individppi<-unlist(strsplit(morig$ppi_target,",",fixed=TRUE))
      nmarkers<-length(individm)
      newdf<-morig
      for(j in 2:nmarkers){

        newdf<-rbind(newdf,morig)
      }
      if(subdir=="Interaction"){
        newdf[,"FEATURE_expression"]<-unlist(individm)
      }
      if(subdir%in%c("Compound","Cell","CompoundME")){
        newdf[,"FEATURE_child"]<-unlist(individm)
      }
      if(!subdir%in%c("Interaction","Compound","Cell","CompoundME")){
        newdf[,"FEATURE"]<-unlist(individm)
      }
      if(length(unlist(individppi))==nrow(newdf)){
      newdf[,"ppi_target"]<-unlist(individppi)}else{
        newdf[,"ppi_target"]<-unlist(individppi)[[1]]
      }
      if(is.null(alldf)){
        alldf<-newdf
      }else{
      alldf<-rbind(alldf,newdf)}
    }
    markersuse<-markersuse[-compoundmarkers,]
    allm<-rbind(markersuse,alldf)
  }else{
    allm<-markersuse
  }
  if(subdir=="Interaction"){
  biomarkers<-unique(allm[,c("FEATURE_expression","ANALYSIS")])
  }
  if(subdir%in%c("Compound","Cell","CompoundME")){
    biomarkers<-unique(allm[,c("FEATURE_child","ANALYSIS")])
  }

  if(!subdir%in%c("Interaction","Compound","Cell","CompoundME")){
    biomarkers<-unique(allm[,c("FEATURE","ANALYSIS")])
  }
  nbm<-nrow(biomarkers)
  CTYPE<-rep('N/A',nbm)
  TARGET<-rep('N/A',nbm)
  SCORE<-rep('N/A',nbm)
  MARKERCLASS<-rep('N/A',nbm)
  Biomarker<-rep('N/A',nbm)
  for (i in 1:nbm){

      if(subdir=="Interaction"){
        id<-which(allm$ANALYSIS==biomarkers$ANALYSIS[i] &
                    allm$FEATURE_expression==biomarkers$FEATURE_expression[i]

        )
        bmarker<-biomarkers$FEATURE_expression[i]
      }
      if(subdir%in%c("Compound","Cell","CompoundME")){
        id<-which(allm$ANALYSIS==biomarkers$ANALYSIS[i] &
                  allm$FEATURE_child==biomarkers$FEATURE_child[i]

        )
      bmarker<-biomarkers$FEATURE_child[i]
      }

      if(!subdir%in%c("Compound","Interaction","Cell","CompoundME")){
        id<-which(allm$ANALYSIS==biomarkers$ANALYSIS[i] &
                    allm$FEATURE==biomarkers$FEATURE[i]


        )
        bmarker<-biomarkers$FEATURE[i]
      }
      currentT<-allm$Depleted.Gene[id]
      #can do enrichment on targets here:
      Biomarker[i]<-bmarker
      CTYPE[i]<-as.character(biomarkers$ANALYSIS[i])
      TARGET[i]<-paste(currentT,collapse=' // ')

      SCORE[i]<-paste(allm$PRIORITY[id],collapse=' // ')
      MARKERCLASS[i]<-paste(allm$CLASS[id],collapse=' // ')



  }

  #GLOBAL<-cbind(TOTRES[,c(1,2,4)],GROUP,TRACTABILITY,MARKERCLASS,ASSOCIATION_EFFECT,MARKER,ANOVA_table_entry,INDICATION,PPI_info,PPI_distance)
  GLOBAL<-data.frame(cbind(CTYPE,Biomarker,TARGET,SCORE,MARKERCLASS),stringsAsFactors = FALSE)

  return(GLOBAL)
}

FilterTargets<-function(GLOBALM){
  newm<-GLOBALM
  for(i in 1:nrow(GLOBALM)){
    singleclass<-unlist(strsplit(GLOBALM[i,"MARKERCLASS"]," // "))
    minclass<-sort(singleclass,decreasing=F)[1]
    toinc<-which(singleclass==minclass)
    singletarget<-unlist(strsplit(GLOBALM[i,"TARGET"]," // "))
    singlescore<-unlist(strsplit(GLOBALM[i,"SCORE"], " // "))
    newm[i,"TARGET"]<-paste(singletarget[toinc],collapse=" // ")
    newm[i,"MARKERCLASS"]<-minclass
    newm[i,"SCORE"]<-paste(singlescore[toinc],collapse=" // ")



  }
  return(newm)
}

collapseTargets2<-function(TOTRES,allMarkers,PPIigraph,PPInet,AnnotInfo,gmtfile,subdir="",splitmarkers=FALSE){
  #TOTRES has all the priority targets and ctypes that pass the priority score threshold

  print(colnames(allMarkers))
  allm<-allMarkers
  markersuse<-NULL
  if(subdir=="Interaction"){
    biomarkers<-unique(allm[,c("FEATURE_expression","ANALYSIS")])
  }
  if(subdir%in%c("Compound","Cell","CompoundME")){
    biomarkers<-unique(allm[,c("FEATURE_child","ANALYSIS")])
  }

  if(!subdir%in%c("Interaction","Compound","Cell","CompoundME")){
    biomarkers<-unique(allm[,c("FEATURE","ANALYSIS")])
  }

  nbm<-nrow(biomarkers)
  CTYPE<-rep('N/A',nbm)
  TARGET<-rep('N/A',nbm)
  SCORE<-rep('N/A',nbm)
  MARKERCLASS<-rep('N/A',nbm)
  Biomarker<-rep('N/A',nbm)
  for (i in 1:nbm){


    if(subdir=="Interaction"){
      id<-which(allm$ANALYSIS==biomarkers$ANALYSIS[i] &
                  allm$FEATURE_expression==biomarkers$FEATURE_expression[i]

      )
      bmarker<-biomarkers$FEATURE_expression[i]
    }
    if(subdir%in%c("Compound","Cell","CompoundME")){
      id<-which(allm$ANALYSIS==biomarkers$ANALYSIS[i] &
                  allm$FEATURE_child==biomarkers$FEATURE_child[i]

      )
      bmarker<-biomarkers$FEATURE_child[i]
    }

    if(!subdir%in%c("Compound","Interaction","Cell","CompoundME")){
      id<-which(allm$ANALYSIS==biomarkers$ANALYSIS[i] &
                  allm$FEATURE==biomarkers$FEATURE[i]


      )
      bmarker<-biomarkers$FEATURE[i]
    }
    currentT<-allMarkers$Depleted.Gene[id]
    #can do enrichment on targets here:
    Biomarker[i]<-bmarker
    allClass<-allMarkers$CLASS[id]
    minclass<-sort(allClass,decreasing=F)[1]
    toinc<-which(allClass==minclass)


    CTYPE[i]<-as.character(biomarkers$ANALYSIS[i])
    TARGET[i]<-paste(currentT[toinc],collapse=' // ')

    SCORE[i]<-paste(allm$PRIORITY[id[toinc]],collapse=' // ')
    MARKERCLASS[i]<-paste(allm$CLASS[id[toinc]],collapse=' // ')



  }

 GLOBAL<-data.frame(cbind(CTYPE,Biomarker,TARGET,SCORE,MARKERCLASS),stringsAsFactors = FALSE)

  return(GLOBAL)
}


harvestBiomarkers<-function(target,ctype,inputdirectoryAnova){

  dname<-paste0(inputdirectoryAnova)
  fn<-dir(dname)

  dn<-grep(ctype,fn,value=TRUE)

  if(ctype!='PANCAN'){
    res<-read.table(paste(inputdirectoryAnova,'/31_anova/',dn,
                          '/OUTPUT/usedInPrioritiz.txt',sep=''),
                    sep='\t',header=TRUE,stringsAsFactors = FALSE)

  }else{
    availdirs<-list.dirs(inputdirectoryAnova)
    Pancandir<-grep("PANCAN",availdirs,value=T)
    pcfiles<-c()
    for(f in 1:length(Pancandir)){
      pcfiles<-c(pcfiles,list.files(Pancandir[f]))
    }
    hitfile<-grep("hits.txt",pcfiles,value=TRUE)
    res<-read.table(hitfile,
                    sep='\t',header=TRUE,stringsAsFactors = FALSE)

  }



  id<-which(res$Depleted.Gene==target)


  if(length(id)>0){

    markers<-res$FEATURE[id]
    SIGN<-sign(as.numeric(res$FEATURE_deltaMEAN_ESS[id]))
    EFFsizeNeg<-as.numeric(res$FEATUREneg_Glass_delta[id])
    EFFsizePos<-as.numeric(res$FEATUREpos_Glass_delta[id])
    EFFsizeCoh<-as.numeric(res$FEATURE_ESS_effect_size[id])
    TP<-as.numeric(res$FEATURE_ESS_T_pval[id])
    FDR<-as.numeric(res$ANOVA.FEATURE.FDR..[id])
    pval<-as.numeric(res$FEATURE_ANOVA_pval[id])
    markerClass<-rep(1,length(SIGN))

    markerClass<-markerClass+(as.numeric(res$FEATURE_ESS_T_pval[id])<0.05)
    markerClass<-markerClass+(as.numeric(res$FEATURE_ANOVA_pval[id])<0.05)
    markerClass<-markerClass+(as.numeric(res$FEATURE_ANOVA_pval[id])<0.001)
    markerClass<-markerClass+(as.numeric(res$ANOVA.FEATURE.FDR..[id])<30)

    if(ctype=='Colorectal Carcinoma' | ctype=='Ovarian Carcinoma'){
      load(paste('../../20181009_PaperRevision_results/31_anova/',dn,'/MSI_HITS.rdata',sep=''))

      idmsi<-which(hits[,'Depleted Gene']==target)

      if(length(idmsi)>0){


        NmarkerClass<-2
        NmarkerClass<-NmarkerClass+(as.numeric(hits[idmsi,'FEATURE_ANOVA_pval'])<0.05)
        NmarkerClass<-NmarkerClass+(as.numeric(hits[idmsi,'FEATURE_ANOVA_pval'])<0.001)
        NmarkerClass<-NmarkerClass+(as.numeric(hits[idmsi,'ANOVA FEATURE FDR %'])<5)


        if(as.numeric(hits[idmsi,'ANOVA FEATURE FDR %'])<5){
          markers<-c('MSI',markers)
          markerClass<-c(NmarkerClass,markerClass)
          SIGN<-c(sign(as.numeric(hits[idmsi,"FEATURE_deltaMEAN_ESS"])),SIGN)
          EFFsizeNeg<-c(as.numeric(hits[idmsi,"FEATUREneg_Glass_delta"]),EFFsizeNeg)
          EFFsizePos<-c(as.numeric(hits[idmsi,"FEATUREpos_Glass_delta"]),EFFsizePos)
          EFFsizeCoh<-c(as.numeric(hits[idmsi,"FEATURE_ESS_effect_size"]),EFFsizeCoh)
          FDR<-c(as.numeric(hits[idmsi,"ANOVA FEATURE FDR %"]),FDR)
          pval<-c(as.numeric(hits[idmsi,'FEATURE_ANOVA_pval']),pval)
          TP<-c(NA,TP)
        }
      }
    }



    markerClass<-c('N/A','D','C','B','A')[markerClass]
    #markerType
    TARGET<-rep(target,length(SIGN))
    CTYPE<-rep(ctype,length(SIGN))

    res<-data.frame(ctype=CTYPE,
                    target=TARGET,
                    marker=markers,
                    sign=SIGN,
                    effectSizeC=EFFsizeCoh,
                    effectSizePos=EFFsizePos,
                    effectSizeNeg=EFFsizeNeg,
                    FDR=FDR,
                    pval=pval,
                    TP=TP,
                    markerClass=markerClass)
  }else{
    res<-NULL
  }

  return(res)
}
Get_DependencyAssociation<-function(TOTRES,subdir){

  if(subdir==""){
    minset<-apply(TOTRES[,c("FEATUREpos_ESS_MEAN","FEATUREneg_ESS_MEAN")],1,which.min)
    depAssoc<-sapply(minset,function(x) ifelse(x==1,"Increased Dep.","Decreased Dep."))
    TOTRES<-cbind(TOTRES,depAssoc)
  }
  if(subdir%in%c("Cont","CExpr","CCN","CMet","CProt","Binary")){
    depAssoc<-sapply(TOTRES[,"logFC"],function(x) ifelse(x<0,"Increased Dep.","Decreased Dep."))
    TOTRES<-cbind(TOTRES,depAssoc)
  }
  if(subdir%in%c("Interaction","Compound","CompoundME","ME")){
    minset<-apply(TOTRES[,c("FEATUREwt_ESS_MEAN","FEATUREcompound_ESS_MEAN")],1,which.min)
    depAssoc<-sapply(minset,function(x) ifelse(x==1,"Decreased Dep.","Increased Dep."))
    TOTRES<-cbind(TOTRES,depAssoc)
  }
  return(TOTRES)
}
LoadBiomarkerRes<-function(CTYPE,subdir="",mutationNumber="",fname=NULL,sortCol="ANOVA FEATURE FDR %",BMdir=NULL,loadF=TRUE,ANOVA.results.dir,pairBD=FALSE,MSI=FALSE){
  #set load to FALSE for function unit testing only

  if(!is.null(BMdir)){
    ANOVA.results.dir<-BMdir
  }
  ddPC<-""
  if(!is.null(fname)){
    print(paste("fname is",fname))
    LoadFile=fname
    if(loadF){
    load(fname)}

  }else{
    if(subdir%in%c("Expr","")){
      if(CTYPE=='PANCAN'){
        #dname<-grep(paste(CTYPE,'_2019',sep=''),dir(paste0(ANOVA.results.dir,'/27_anova/')),value=TRUE)
        ddPC<-(paste0(ANOVA.results.dir,"/27_anova/"))
        if(file.exists(paste0(ddPC,"/PANCAN/OUTPUT/ANOVA_results",subdir,mutationNumber,'.rdata'))){
          LoadFile<-paste(ddPC,'/PANCAN/OUTPUT/ANOVA_results',subdir,mutationNumber,'.rdata',sep='')
        }else{
          ddPC<-""
          if(file.exists(paste0(ANOVA.results.dir,"/PANCAN/OUTPUT/ANOVA_results",subdir,mutationNumber,'.rdata'))){
            LoadFile<-paste(ANOVA.results.dir,'/PANCAN/OUTPUT/ANOVA_results',subdir,mutationNumber,'.rdata',sep='')
          }
        }

      }else{
       LoadFile<-paste(ANOVA.results.dir,"/",CTYPE,'/OUTPUT/ANOVA_results',subdir,mutationNumber,'.rdata',sep='')
      }

    }

    if(subdir%in%c("Interaction","Compound","CompoundME")){

      LoadFile<-paste(ANOVA.results.dir,"/",CTYPE,'/OUTPUT/ANOVA_results',subdir,mutationNumber,'.rdata',sep='')

    }
    if(subdir%in%c("Cont","CExpr","CCN","CMet","CProt","Binary")){

      if(pairBD){
        LoadFile<-paste(ANOVA.results.dir,"/",CTYPE,'/OUTPUT/ANOVA_AggresultsPairBD',subdir,mutationNumber,'.rdata',sep='')
      }else{
        if(MSI){
          LoadFile<-paste(ANOVA.results.dir,"/",CTYPE,'/OUTPUT/LM_resultsMSI',subdir,mutationNumber,'.rdata',sep='')
        }else{
          LoadFile<-paste(ANOVA.results.dir,"/",CTYPE,'/OUTPUT/LM_results',subdir,mutationNumber,'.rdata',sep='')
        }
      }

    }
  }
  print("File loaded:")
  print(LoadFile)

  if(loadF){
    load(file=LoadFile)
    if(MSI){TOTRES<-TOTRES_MSI}
    checkCols<-if(!sortCol%in%colnames(TOTRES)){
      stop(paste("column",sortCol,"not available in output (TOTRES) to sort by"))
    }
    TOTRES<-Get_DependencyAssociation(TOTRES,subdir)
    TOTRES<-TOTRES[order(unlist(TOTRES[,sortCol])),]}else{TOTRES<-NULL}
  return(list(TOTRES=TOTRES,loadedFile=LoadFile))
}
LoadPriorityVectors<-function(CTYPE,Combine,dir.Results,subdir,pairBD,mutationNumber){
  if (CTYPE!='PANCAN'){
    if(Combine){
      load(paste0(dir.Results,'/PvectorsL3.Rdata'))
      PRIORITY_vectors<-PvectorsL3
      if(pairBD){
        load(paste0(dir.Results,'/36_targetDistributions/PairBD/priority_threshold_all.RData'))
      }else{
        load(paste0(dir.Results,'/36_targetDistributions/priority_threshold_all.RData'))
      }
      PRIORITY_vectors<-PRIORITY_vectors[[paste0(CTYPE,mutationNumber)]]
    }else{
      if(pairBD){
        load(paste0(dir.Results,'/35_PrioritizedHits/PairBD/',subdir,'/_PRIORITY_vectorsL3.RData'))
      }else{
        load(paste0(dir.Results,'/35_PrioritizedHits/',subdir,'/_PRIORITY_vectorsL3.RData'))
      }
      if(pairBD){
        outdir<-"/36_targetDistributions/PairBD/"
      }else{
        outdir<-"/36_targetDistributions/"
      }
      if(subdir!=""){
        load(paste0(dir.Results,outdir,'priority_threshold_',subdir,'.RData'))
      }else{
        load(paste0(dir.Results,outdir,'priority_threshold.RData'))
      }
      PRIORITY_vectors<-PRIORITY_vectorsL3[[paste0(CTYPE,mutationNumber)]]
    }


  }else{
    #PANCAN
    if(Combine){
      load(paste0(dir.Results,'/36_targetDistributions/priority_threshold_all_PANCAN.RData'))
      load(paste0(dir.Results,'/PvectorsL3_PANCAN.Rdata'))
      PRIORITY_vectors<-PvectorsL3
    }else{
      if(subdir==""){
        load(paste0(dir.Results,'/35_PrioritizedHits/_PRIORITY_vectorsL3_PANCAN.RData'))
        load(paste0(dir.Results,'/36_targetDistributions/priority_threshold_PANCAN.RData'))
      }
      if(subdir=="PairBD"){
        load(paste0(dir.Results,'/35_PrioritizedHits/PairBD/_PRIORITY_vectorsL3_PANCAN.RData'))
        load(paste0(dir.Results,'/36_targetDistributions/PairBD/priority_threshold_PANCAN.RData'))
      }
      if(subdir%in%c("Cont","CExpr","CCN","CMet","CProt","Binary")){
        load(paste0(dir.Results,'/35_PrioritizedHits/',subdir,'/_PRIORITY_vectorsL3_PANCAN.RData'))
        load(paste0(dir.Results,'/36_targetDistributions/',subdir,'/priority_threshold_PANCAN.RData'))
      }
      PRIORITY_vectors<-PRIORITY_vectorsL3
    }
  }
  return(list(PRIORITY_vectors=PRIORITY_vectors,priority_threshold=priority_threshold))
}
GetClassIDs<-function(subdir,TOTRES,fdra=5,fdrb=5,fdrc=10,fdrd=15,t_val=1.5,f2_val=1.5){
  if(subdir%in%c("","PairBD")){
    # class D
    class_D_id<-which(unlist(TOTRES[,"ANOVA FEATURE FDR %"])<15)
    #class C
    class_C_id<-which(unlist(TOTRES[,"ANOVA FEATURE FDR %"])<10)
    #class B
    class_B_id<-which(unlist(TOTRES[,"ANOVA FEATURE FDR %"])<5)
    #class A
    class_A_id<-which(unlist(TOTRES[,"ANOVA FEATURE FDR %"])<5 &
                        unlist(TOTRES[,"FEATURE_ANOVA_pval"])<0.001 &
                        unlist(TOTRES[,"FEATUREpos_Glass_delta"])>1 &
                        unlist(TOTRES[,"FEATUREneg_Glass_delta"])>1)

  }
  if(subdir%in%c("Cont","CExpr","CCN","CMet","CProt","Binary")){


    colcheck<-which(colnames(TOTRES)=="Depleted Gene")
    if(length(colcheck)==0){
      colcheck<-which(colnames(TOTRES)=="Depleted.Gene")
      colnames(TOTRES)[colcheck]<-"Depleted Gene"
    }

    # class D
    class_D_id<-which(unlist(TOTRES[,'FDR'])<fdrd)

    #class C
    class_C_id<-which(unlist(TOTRES[,"FDR"])<fdrc)

    #class B
    class_B_id<-which(unlist(TOTRES[,"FDR"])<fdrb)

    #class A
    if(subdir=="CExpr"){
      class_A_id<-which(unlist(TOTRES[,"FDR"])<fdra & abs(unlist(TOTRES[,"logFC"]))>t_val &
                          unlist(TOTRES[,"f2"])>1.5)
    }
    if(subdir=="CMet"){
      class_A_id<-which(unlist(TOTRES[,"FDR"])<fdra & abs(unlist(TOTRES[,"logFC"]))>t_val &
                          unlist(TOTRES[,"f2"])>1.5)
    }
    if(subdir=="CProt"){
      class_A_id<-which(unlist(TOTRES[,"FDR"])<fdra & abs(unlist(TOTRES[,"logFC"]))>t_val &
                          unlist(TOTRES[,"f2"])>1.5)
    }
    if(subdir=="CCN"){
      class_A_id<-which(unlist(TOTRES[,"FDR"])<fdra & abs(round(unlist(TOTRES[,"logFC"])))>=t_val &
                          unlist(TOTRES[,"f2"])>0.9)
    }

    if(subdir=="Binary"){
      #class_A_id<-which(unlist(TOTRES[,"FDR"])<5 & unlist(TOTRES[,"P.Value"])<0.001 &
      #                     unlist(TOTRES[,"IQR"])>1&unlist(TOTRES[,"Range"]>1))
      class_A_id<-which(unlist(TOTRES[,"FDR"])<fdra & abs(unlist(TOTRES[,"logFC"]))>t_val &
                          unlist(TOTRES[,"f2"])>0.9)
    }

  }
  if(subdir=="Interaction"){
    # class D
    class_D_id<-which(unlist(TOTRES[,'FEATURE_INT_pval'])<0.05)

    #class C
    class_C_id<-which(unlist(TOTRES[,"FEATURE_INT_pval"])<0.001)

    #class B
    class_B_id<-which(unlist(TOTRES[,"ANOVA FEATURE FDR %"])<30|unlist(TOTRES[,"FEATURE_INT_pval"]<0.00001))

    #class A
    class_A_id<-which(unlist(TOTRES[,"ANOVA FEATURE FDR %"])<25 &
                        unlist(TOTRES[,"FEATURE_INT_pval"])<0.00001 &
                        ((unlist(TOTRES[,"FEATUREmutP_Epos_Glass_delta"])>1 &
                            unlist(TOTRES[,"FEATUREmutP_Eneg_Glass_delta"])>1)|(unlist(TOTRES[,"FEATUREmutN_Epos_Glass_delta"])>1 &
                                                                                  unlist(TOTRES[,"FEATUREmutN_Eneg_Glass_delta"])>1)))
  }
  if(subdir%in%c("Compound","CompoundME","ME")){

    # class D
    class_D_id<-which(unlist(TOTRES[,'LRT FEATURE FDR %'])<15&unlist(TOTRES[,"ANOVA FEATURE FDR %"])<15)

    #class C
    class_C_id<-which(unlist(TOTRES[,"LRT FEATURE FDR %"])<10&unlist(TOTRES[,"ANOVA FEATURE FDR %"])<10)

    #class B
    class_B_id<-which(unlist(TOTRES[,"ANOVA FEATURE FDR %"])<5&unlist(TOTRES[,"LRT FEATURE FDR %"])<5)

    #class A
    class_A_id<-which(unlist(TOTRES[,"ANOVA FEATURE FDR %"])<5 &
                        unlist(TOTRES[,"LRT FEATURE FDR %"])<5&
                        unlist(TOTRES[,"FEATUREcompound_pos_Glass_delta"])>1 &
                        unlist(TOTRES[,"FEATUREcompound_neg_Glass_delta"])>1
                 )




  }
  if(subdir=="Cell"){
    # class D
    class_D_id<-which(unlist(TOTRES[,'FEATURE_compound_pval'])<0.05&unlist(TOTRES[,"LRT FEATURE FDR %"])<30)

    #class C
    class_C_id<-which(unlist(TOTRES[,"FEATURE_compound_pval"])<0.01&unlist(TOTRES[,"LRT FEATURE FDR %"])<20)

    #class B
    class_B_id<-which(unlist(TOTRES[,"ANOVA FEATURE FDR %"])<15&unlist(TOTRES[,"LRT FEATURE FDR %"])<15)

    #class A
    class_A_id<-which(unlist(TOTRES[,"ANOVA FEATURE FDR %"])<10&unlist(TOTRES[,"FEATURE_compound_pva"])<0.001&unlist(TOTRES[,"LRT FEATURE FDR %"])<10
    )
  }

  if(subdir=="Expr"){
    # class D
    class_D_id<-which(unlist(TOTRES[,'FEATURE_ESS_T_pval'])<0.05)

    #class C
    class_C_id<-which(unlist(TOTRES[,"FEATURE_ANOVA_pval"])<0.001)

    #class B
    class_B_id<-which(unlist(TOTRES[,"ANOVA FEATURE FDR %"])<30)

    #class A
    class_A_id<-which(unlist(TOTRES[,"ANOVA FEATURE FDR %"])<25 &
                        unlist(TOTRES[,"FEATURE_ANOVA_pval"])<0.001 &
                        unlist(TOTRES[,"FEATUREpos_Glass_delta"])>1 &
                        unlist(TOTRES[,"FEATUREneg_Glass_delta"])>1)
  }
  return(list(class_A_id=class_A_id,class_B_id=class_B_id,class_C_id=class_C_id,class_D_id=class_D_id))
}
GetClassIDsPANCAN<-function(subdir,TOTRES,fdra=1,fdrb=1,fdrc=1,fdrd=1,t_val=1.25,f2_val=0.15){
  if(subdir==""){

    # class D
    class_D_id<-which(unlist(TOTRES[,'ANOVA FEATURE FDR %'])<10 )

    #class C
    class_C_id<-which(unlist(TOTRES[,"ANOVA FEATURE FDR %"])<5 )

    #class B
    class_B_id<-which(unlist(TOTRES[,"ANOVA FEATURE FDR %"])<1)

    #class A
    class_A_id<-which(unlist(TOTRES[,"ANOVA FEATURE FDR %"])<1 &
                        unlist(TOTRES[,"FEATUREpos_Glass_delta"])>1 &
                        unlist(TOTRES[,"FEATUREneg_Glass_delta"])>1)
  }
  if(subdir=="PairBD"){

    # class D
    class_D_id<-which(unlist(TOTRES[,'ANOVA FEATURE FDR %'])<15)
    #class C
    class_C_id<-which(unlist(TOTRES[,"ANOVA FEATURE FDR %"])<10)
    #class B
    class_B_id<-which(unlist(TOTRES[,"ANOVA FEATURE FDR %"])<5)
    #class A
    class_A_id<-which(unlist(TOTRES[,"ANOVA FEATURE FDR %"])<5 &
                        unlist(TOTRES[,"FEATUREpos_Glass_delta"])>1 &
                        unlist(TOTRES[,"FEATUREneg_Glass_delta"])>1)

  }
  if(subdir%in%c("Cont","CExpr","CCN","CMet","CProt","Binary")){


    class_D_id<-which(unlist(TOTRES[,'FDR'])<fdrd)

    #class C
    class_C_id<-which(unlist(TOTRES[,"FDR"])<fdrc &
                        unlist(TOTRES[,"f2"])>0.07)

    #class B
    class_B_id<-which(unlist(TOTRES[,"FDR"])<fdrb& abs(unlist(TOTRES[,"logFC"]))>1 &
                        unlist(TOTRES[,"f2"])>0.08)

    #class A
    if(subdir=="CExpr"){
      class_A_id<-which(unlist(TOTRES[,"FDR"])<fdra & abs(unlist(TOTRES[,"logFC"]))>t_val &
                          unlist(TOTRES[,"f2"])>f2_val)
    }
    if(subdir=="CMet"){
      class_A_id<-which(unlist(TOTRES[,"FDR"])<fdra & abs(unlist(TOTRES[,"logFC"]))>t_val &
                          unlist(TOTRES[,"f2"])>f2_val)
    }
    if(subdir=="CProt"){
      class_A_id<-which(unlist(TOTRES[,"FDR"])<fdra &abs(unlist(TOTRES[,"logFC"]))>t_val &
                          unlist(TOTRES[,"f2"])>f2_val)
    }
    if(subdir=="CCN"){
      class_A_id<-which(unlist(TOTRES[,"FDR"])<fdra & abs(unlist(TOTRES[,"logFC"]))>t_val &
                          unlist(TOTRES[,"f2"])>f2_val)
    }
    if(subdir=="Binary"){
      #class_A_id<-which(unlist(TOTRES[,"FDR"])<1 & unlist(TOTRES[,"P.Value"])<0.001 &
       #                   unlist(TOTRES[,"IQR"])>1&unlist(TOTRES[,"Range"]>1))
      class_A_id<-which(unlist(TOTRES[,"FDR"])<fdra & abs(unlist(TOTRES[,"logFC"]))>t_val &
                          unlist(TOTRES[,"f2"])>f2_val)
    }

  }

  if(subdir%in%c("Compound","CompoundME","ME")){
    class_D_id<-which(unlist(TOTRES[,"ANOVA FEATURE FDR %"])<10 &
                                                unlist(TOTRES[,"LRT FEATURE FDR %"])<10 )

    #class C
    class_C_id<-which(unlist(TOTRES[,"ANOVA FEATURE FDR %"])<5 &
                                                unlist(TOTRES[,"LRT FEATURE FDR %"])<5 )

    #class B
    class_B_id<-which(unlist(TOTRES[,"ANOVA FEATURE FDR %"])<1 &
                                                 unlist(TOTRES[,"LRT FEATURE FDR %"])<1)

    #class A
    class_A_id<-which(unlist(TOTRES[,"ANOVA FEATURE FDR %"])<1 &
                                                             unlist(TOTRES[,"LRT FEATURE FDR %"])<1 &
                                                             unlist(TOTRES[,"FEATUREcompound_pos_Glass_delta"])>1 &
                                                             unlist(TOTRES[,"FEATUREcompound_neg_Glass_delta"])>1)

  }
  return(list(class_A_id=class_A_id,class_B_id=class_B_id,class_C_id=class_C_id,class_D_id=class_D_id))

}
saveFilesForPrioritization<-function(ANOVA.results.dir,CTYPE,subdir="",mutationNumber="",fname=NULL,sortCol="ANOVA FEATURE FDR %",
                                     keepCols=c("ANALYSIS","CLASS","assoc_id","FEATURE","Depleted Gene","Name","FEATURE_deltaMEAN_ESS","ANOVA FEATURE FDR %","BMstringID","string_id","depAssoc"),
                                     DepCol="FEATURE_deltaMEAN_ESS",BMdir=NULL,Combine=FALSE,pairBD=FALSE,dir.Results,MSI=FALSE,
                                     fa=5,fb=5,fc=10,fd=15,t_val=1.5,f2_val=0.75){
  TOTRES<-LoadBiomarkerRes(CTYPE,subdir,mutationNumber,fname,sortCol,BMdir,ANOVA.results.dir=ANOVA.results.dir,pairBD=pairBD,MSI=MSI)$TOTRES

  PriorityRes<-LoadPriorityVectors(CTYPE,Combine,dir.Results,subdir=subdir,pairBD,mutationNumber=mutationNumber)
  PRIORITY_vectors<-PriorityRes$PRIORITY_vectors
  priority_threshold<-PriorityRes$priority_threshold
  if(!Combine){
    #try to keep more marker information in case passes on combined threshold if not on individual threshold
    priority_threshold<-min(20,priority_threshold)
  }

    if(subdir=="EMDE"){
      fn<-paste(ANOVA.results.dir,"/",CTYPE,'/OUTPUT/EMfilt_results.txt',sep='')
      print(fn)
      if(file.exists(fn)){

       EMRESfilt<-read.delim(fn,header=TRUE,stringsAsFactors = F,sep="\t")
       maxGSEA<-sapply(EMRESfilt[,"gsea"],function(x) max(as.numeric(unlist(strsplit(as.character(x),"//",fixed=TRUE)))))
       hyperPval<-sapply(EMRESfilt[,"hyper"],function(x) min(as.numeric(unlist(strsplit(as.character(x),"//",fixed=TRUE)))))
        # class D
        class_D_id<-which(unlist(EMRESfilt[,'nnmd'])>4&!unlist(EMRESfilt[,'gsea'])%in%c("NoDE","NoBimodal"))

        #class C
        class_C_id<-which(unlist(EMRESfilt[,"nnmd"])>5&!unlist(EMRESfilt[,'gsea'])%in%c("NoDE","NoBimodal"))

        #class B
        class_B_id<-which(unlist(EMRESfilt[,"nnmd"])>5&maxGSEA>0.3)

        #class A
        class_A_id<-which(unlist(EMRESfilt[,"nnmd"])>5 &maxGSEA>0.3&
                          hyperPval<0.001 )
        TOTRES<-EMRESfilt
        colnames(TOTRES)[colnames(TOTRES)=="dep"]<-"Depleted Gene"
        colnames(TOTRES)[colnames(TOTRES)=="DEgenes"]<-"FEATURE"
        temp<-gsub("//",",",TOTRES[,"FEATURE"],fixed=TRUE)
        TOTRES[,"FEATURE"]<-temp
        colnames(TOTRES)[colnames(TOTRES)=="em_id"]<-"assoc_id"
      }else{
        TOTRES<-NULL
        RESALL<-NULL

        return(list(TOTRES=TOTRES,RESALL=RESALL))
      }

    }else{
      if(CTYPE!="PANCAN"){
        PriorityClassIDs<-GetClassIDs(subdir,TOTRES,fdra=fa,fdrb=fb,fdrc=fc,fdrd=fd,t_val=t_val,f2_val=f2_val)
      }else{
        if(MSI){
          PriorityClassIDs<-GetClassIDsPANCAN(subdir,TOTRES,fdra=1,fdrb=1,fdrc=5,fdrd=10,t_val=1,f2_val=0.1)
        }else{
          PriorityClassIDs<-GetClassIDsPANCAN(subdir,TOTRES)
        }
      }
      class_A_id<-PriorityClassIDs$class_A_id
      class_B_id<-PriorityClassIDs$class_B_id
      class_C_id<-PriorityClassIDs$class_C_id
      class_D_id<-PriorityClassIDs$class_D_id
    }
  #get class definitions to add to output:
  CLASS<-rep('NA',nrow(TOTRES))

  CLASS[class_D_id]<-'D'
  CLASS[class_C_id]<-'C'
  CLASS[class_B_id]<-'B'
  CLASS[class_A_id]<-'A'
  PRIORITY<-rep('NA',nrow(TOTRES))

  ANALYSIS<-rep(paste0(CTYPE,mutationNumber),nrow(TOTRES))
  TOTRES<-cbind(ANALYSIS,CLASS,TOTRES)
  id<-unique(c(class_A_id,class_B_id,class_C_id,class_D_id))
  PRIORITY<-PRIORITY_vectors[unlist(TOTRES[,"Depleted Gene"])]

  print(length(id))
  TOTRES<-cbind(TOTRES,PRIORITY)


  print(colnames(TOTRES))
  print(dim(TOTRES))

  if(subdir!="EMDE"){

    print('starting check hpT')
    print(CTYPE)

    hpT<-names(which(PRIORITY_vectors>=priority_threshold))
    id<-which(is.element(unlist(TOTRES[,"Depleted Gene"]),hpT))
    print(paste("priority genes",head(hpT)))
    if(length(id)>0){
      if(length(id)>1){
        TOTRES<-TOTRES[id,]
        }else{
          cnames<-colnames(TOTRES)
          TOTRES<-TOTRES[id,]
          TOTRES<-as.data.frame(TOTRES,nrow=1)
          colnames(TOTRES)<-cnames
        }
        TOTRES<-TOTRES[,intersect(keepCols,as.character(colnames(TOTRES)))]

      for(i in 1:ncol(TOTRES)){
        TOTRES[,i]<-unlist(TOTRES[,i])
      }

      print('finished getting main results')
      print(colnames(TOTRES))
      if(subdir==""){
        #not doing separate MSI analysis
        ddPC<-NULL
        #changed temporarily because MSI analysis not coded yet
        if(CTYPE=='PANCAN' |CTYPE=='Colorectal.Carcinoma' |CTYPE=='Gastric.Carcinoma' |CTYPE=='Ovarian.Carcinoma'){
            if(CTYPE=="PANCAN"){
              fn<-paste(ANOVA.results.dir,'/27_anova/',ddPC,'_MSI_allTests.rdata',sep='')
            }else{
              fn<-paste(ANOVA.results.dir,"/",CTYPE,'_MSI_allTests.rdata',sep='')
            }
          if(file.exists(fn)){
            load(fn)
            print('getting MSI results')
            id<-which(as.numeric(unlist(hits[,"ANOVA FEATURE FDR %"]))<5 &
                    as.numeric(unlist(hits[,"FEATURE_ANOVA_pval"]))<0.001 &
                    as.numeric(unlist(hits[,"FEATUREpos_Glass_delta"]))>1 &
                    as.numeric(unlist(hits[,"FEATUREneg_Glass_delta"])>1))
            hits<-matrix(hits[id,],length(id),ncol(hits),dimnames = list(NULL,colnames(hits)))
            hpT<-names(which(PRIORITY_vectors>=priority_threshold))
            id<-which(is.element(unlist(hits[,"Depleted Gene"]),hpT))
            hits<-matrix(hits[id,c(1:4,14)],length(id),5,dimnames = list(NULL,colnames(hits)[c(1:4,14)]))
            hits[,1]<-paste('msi_',hits[,1],sep='')
            #hits[which(as.numeric(unlist(hits[,5]))<0),5]<-'Increased dep.'
            #hits[which(unlist(hits[,5])!='Increased dep.'),5]<-'Decreased dep.'
            #hits<-cbind(rep(CTYPE,nrow(hits)),rep('A',nrow(hits)),hits)
            PRIORITY<-PRIORITY_vectors[unlist(hits[,"Depleted Gene"])]
            hits<-cbind(hits,PRIORITY)
            print('finished getting MSI results')
            colnames(hits)[1:2]<-c("ANALYSIS","CLASS")
            print(colnames(hits))
            print(colnames(TOTRES))
            TOTRES<-rbind(TOTRES,hits)
            for(i in 1:ncol(TOTRES)){
              TOTRES[,i]<-unlist(TOTRES[,i])
            }
          }
        }
      }
    }else{
      #dont have any high priority targets hpT=0
      TOTRES<-NULL

    }
  }else{
    #looking at EMDE output so different file format etc:
    hpT<-names(which(PRIORITY_vectors>=priority_threshold))
    id<-which(is.element(unlist(TOTRES[,"Depleted Gene"]),hpT))
    if(length(id)>0){

      TOTRES<-TOTRES[id,]
      for(i in 1:ncol(TOTRES)){
        TOTRES[,i]<-unlist(TOTRES[,i])
      }
      PRIORITY<-PRIORITY_vectors[unlist(TOTRES[,"Depleted Gene"])]

      TOTRES<-cbind(TOTRES,PRIORITY)

    }else{TOTRES<-NULL}


  }


  #changed to allow use of foreach
  print(head(TOTRES))
  return(TOTRES)
}
getStringID<-function(mkr,geneAnnotations){
  splitC<-changeCNAsep(mkr)
  splitm<-sapply(unlist(splitC),function(x) strsplit(x,",",fixed=TRUE))
  splitm<-unlist(splitm)
  markers<-getMarkers(splitm)
  stringM<-geneAnnotations[markers,"string_id"]
  ids<-paste(stringM,collapse=",")
  return(ids)
}
changeCNAsep<-function(mkr){
  splitC<-unlist(sapply(unlist(mkr),function(x) strsplit(x,"[()]")))
  if(length(splitC)>1){
    #have a cna marker with associated genes
    #find things that are not full markers:
    cnaList<-!grepl("mut|Expr",splitC)
    splitC[cnaList]<-gsub(",",";",splitC[cnaList])
    splitC<-paste0(splitC,collapse="")
  }
  return(splitC)
}
getMarkers<-function(splitm){
  splitM<-sapply(splitm,function(x) strsplit(x,"_",fixed=TRUE))
  type<-unlist(lapply(splitM,function(x) x[length(x)]))
  type[!type%in%c("mut","Expr","var")]<-"CNA"
  nmarkers<-length(type)
  markers<-c()
  j=1
  for(i in 1:nmarkers){
    if(type[i]%in%c("mut","Expr","var")){
      markers[j]<-splitM[[i]][1]

      j=j+1
    }else{
      possCNAmarker<-splitM[[i]][length(splitM[[i]])]
      if(possCNAmarker!=""){

        CNAmarker<-unlist(strsplit(possCNAmarker,";",fixed=TRUE))
        nM<-length(CNAmarker)
        markers[j:(j+nM-1)]<-CNAmarker
        j=j+nM
      }

    }
  }
  markers<-gsub(" ","",markers, fixed = T)
  return(markers)
}
GetSP<-function(genemap,anovamarkers,PPIigraph,maxdist=2,subdir="",distances=NULL){
  depgene<-unlist(anovamarkers[1])
  mkr<-anovamarkers[2]
  splitC<-changeCNAsep(mkr)

  splitm<-sapply(unlist(splitC),function(x) strsplit(x,",",fixed=TRUE))
  splitm<-unlist(splitm)

  if(subdir=="EMDE"){
    markers<-gsub(" ","",splitm,fixed=T)
  }else{
    markers<-getMarkers(splitm)

  }
  newmarkers<-markers

    stringM<-genemap[markers,"STRING_id"]
    stringD<-genemap[depgene,"STRING_id"]

    sp<-vector("numeric",length=length(markers))
    for(i in 1:length(markers)){
      if(!is.na(stringM[i])){
        if(!is.na(stringD)){
          if(is.null(distances)){
            temp<-tryCatch(distances(PPIigraph,stringD,stringM[i]),error=function(e){NA})
          }else{
            temp<-tryCatch(distances[stringD,stringM[i]])
          }
          if(!is.na(temp)){
            if(temp>maxdist){sp[i]<-Inf}
            if(length(temp)==0){sp[i]<-"No path"}else{sp[i]<-temp}
          }else{sp[i]<-"No path"}
        }else{sp[i]<-"Target not in Network"}
      }else{sp[i]<-"Marker not in Network"}
    }

    if(sum(sp=="Inf")==length(sp)){
      #all markers greater than max PPI distance:
      sp<-"No Marker"
    }

    newmarkers<-newmarkers[sp!="Inf"]
    sp<-sp[sp!="Inf"]

  markers<-paste(newmarkers,collapse=",")
  ppi<-paste(sp,collapse = ",")
  return(c(markers=markers,ppi=ppi))
}

GetSPFromIDs<-function(anovamarkers,PPIigraph,maxdist=Inf,subdir="",distances=NULL){
  stringD<-anovamarkers["string_id"]
  markers<-unlist(sapply(anovamarkers["BMstringID"],function(x) strsplit(x,",",fixed=T)))
  markers<-gsub(" ","",markers)
  sp<-vector("numeric",length=length(markers))
  for(i in 1:length(markers)){
  stringM<-markers[i]

  stringM<-unlist(strsplit(stringM,"[&|]",fixed=F))
  mp<-vector("numeric",length=length(stringM))
  for(j in 1:length(stringM)){
    if(!is.na(stringM[j])){
      if(!is.na(stringD)){
        if(is.null(distances)){

          temp<-tryCatch(distances(PPIigraph,stringD,stringM[j]),error=function(e){NA})
        }else{
          temp<-tryCatch(distances[stringD,stringM[j]],error=function(e){NA})
        }
        if(!is.na(temp)){
          if(temp>maxdist){mp[j]<-Inf}
          if(length(temp)==0){mp[j]<-"No path"}else{mp[j]<-temp}
        }else{mp[j]<-"No path"}
      }else{mp[j]<-"Target not in Network"}
    }else{mp[j]<-"Marker not in Network"}
  }

  if(sum(mp=="Inf")==length(mp)){
    #all markers greater than max PPI distance:
    sp[i]<-"No Path"
  }

  if(length(mp)>1){

      sp[i]<-paste0(mp,collapse=":")

  }else{
    sp[i]<-mp
  }


  }
  ppi<-paste(sp,collapse = ",")

  markers<-unlist(anovamarkers["FEATURE"])
  names(markers)<-NULL
  return(c(markers=markers,ppi=ppi))
}


shortestPathFn<-function(x,y,z,output){
  tryCatch(
    expr={
      temp<-shortest_paths(x,y,z,output=output)

    },
    error={
      temp<-NA

    },
    finally = {return(temp)}

  )
}

filterPPI<-function(markeroutput,PPIigraph,PPInet){

  marker<-markeroutput$MARKER

  PPImin<-markeroutput$PPI_min
  mDist<-markeroutput$PPI_distance
  if(!is.na(marker)){
  if(marker!=""){
    if(markeroutput$Mtype=="EMDE"){
      #split by comma
      marker<-unlist(strsplit(marker,",",fixed=TRUE))
      mDist<-unlist(strsplit(mDist,",",fixed=TRUE))
    }else{
      #split by //
      marker<-unlist(strsplit(marker,"//",fixed=TRUE))
      mDist<-unlist(strsplit(mDist,"//",fixed=TRUE))
    }
    #remove all those less than min PPI unless min PPI is target? Then allow 1?
    PPImin<-max(c(1,as.numeric(PPImin)))
    mDistN<-as.numeric(mDist)
    for(i in 1:length(mDistN)){
      if(is.na(mDistN[i])){
        mDistN[i]<-min(as.numeric(unlist(strsplit(mDist[i],",",fixed=TRUE))))
      }
    }
    mDistN[mDistN>PPImin]<-NA
    marker<-marker[!is.na(mDistN)]
    mDistN<-mDistN[!is.na(mDistN)]
    #do enrichment
    subdir<-markeroutput$Mtype
    expressionMarkers<-unique(getSingleExpression(marker,subdir))
    expressionMarkers<-PPImarkerCheck(expressionMarkers,PPInet)

    if(length(expressionMarkers)>2){
      expressionEnrich<-tryCatch(ppi_enrichment(expressionMarkers,PPIigraph)$enrichment,error=function(e){return(NA)})
    }else{
      expressionEnrich<-NA
    }
  }else{
    expressionEnrich<-NA
  }
  }else{
    expressionEnrich<-NA
  }
  return(list(FilterMarker=paste0(marker,collapse="//"),FilterDist=PPImin,FilterEnrich=expressionEnrich))
}
densityThresholding<-function(set1,set2,x=seq(0,100,length.out=100)){
  #x default to range of priority scores 0-100. Change for other types of data.
  #set2 should be positive set (e.g. known targets or biomarkers)
  #set1 rest of data.
  kess<-density(set1, kernel = "gaussian")
  knon<-density(set2, kernel = "gaussian")


  nonfitx <- approx(knon$x,knon$y,x)$y
  incx<-x[which(!is.na(nonfitx))]
  if(length(incx)>5){
    logratio_sample <- log2( approx(kess$x,kess$y,x)$y / approx(knon$x,knon$y,x)$y )

    priority_threshold<-x[min(which(logratio_sample>=1 & x>=mean(set2)))]
  }else{
    priority_threshold<-median(set2,na.rm=T)
  }
  return(priority_threshold)
}

getGenesFromBiomarkers<-function(mkr){

splitC<-unlist(sapply(unlist(mkr),function(x) strsplit(x,"[()]")))
if(length(splitC)>1){
  #have a cna marker with associated genes
  #find things that are not full markers:
  cnaList<-!grepl("mut|Expr",splitC)
  splitC[cnaList]<-gsub(",",";",splitC[cnaList])
  splitC<-paste0(splitC,collapse="")
}
splitm<-sapply(unlist(splitC),function(x) strsplit(x,",",fixed=TRUE))
splitm<-unlist(splitm)


  splitM<-sapply(splitm,function(x) strsplit(x,"_",fixed=TRUE))
  type<-unlist(lapply(splitM,function(x) x[length(x)]))
  type[!type%in%c("mut","Expr")]<-"CNA"
  nmarkers<-length(type)
  markers<-c()
  j=1
  for(i in 1:nmarkers){
    if(type[i]%in%c("mut","Expr")){
      markers[j]<-splitM[[i]][1]

      j=j+1
    }else{
      possCNAmarker<-splitM[[i]][length(splitM[[i]])]
      if(possCNAmarker!=""){

        CNAmarker<-unlist(strsplit(possCNAmarker,";",fixed=TRUE))
        nM<-length(CNAmarker)
        markers[j:(j+nM-1)]<-CNAmarker
        j=j+nM
      }

    }
  }
  markers<-gsub(" ","",markers, fixed = T)



return(markers)
}
