collapseMarkersSingle<-function(TOTRES,allMarkers){
  MARKER<-rep('N/A',length(TOTRES$MARKERCLASS))
  ASSOCIATION_EFFECT<-rep('N/A',length(TOTRES$MARKERCLASS))
  ANOVA_table_entry<-rep('N/A',length(TOTRES$MARKERCLASS))
  
  MARKER_TYPE<-rep('N/A',length(TOTRES$MARKERCLASS))
  
  for (i in 1:nrow(TOTRES)){
    
    if(TOTRES$MARKERCLASS[i]!='N/A'){
      #this is where only markers of the maximum class (TOTRES$MARKERCLASS) for depleted gene meeting tractability thresholds are selected
      
      id<-which(allMarkers$ANALYSIS==TOTRES$ctype[i] & 
                  allMarkers$Depleted.Gene==unlist(TOTRES$TARGET[i]) &
                  allMarkers$CLASS==TOTRES$MARKERCLASS[i])
      
      currentMar<-allMarkers$FEATURE[id]
      
      MARKER[i]<-paste(currentMar,collapse=' // ')
      
      #changed from FEATURE_pos to FEATURE_delta on 4.8.19
      ASSOCIATION_EFFECT[i]<-paste(allMarkers[id,6],collapse=' // ')
      ANOVA_table_entry[i]<-paste(allMarkers$assoc_id[id],collapse=' // ')
      
      
    }
    
  }
  TRACTABILITY<-TOTRES[,3]
  #GLOBAL<-cbind(TOTRES[,c(1,2,4)],GROUP,TRACTABILITY,MARKERCLASS,ASSOCIATION_EFFECT,MARKER,ANOVA_table_entry,INDICATION,PPI_info,PPI_distance)
  GLOBAL<-cbind(TOTRES[,c(4,5,26,27,28,29:32)],ASSOCIATION_EFFECT,MARKER,ANOVA_table_entry)
  
  return(GLOBAL)
}


processNetMod<-function(directoryInput,outputSuffix,AOVRes=NULL,FDRthresh=15,AOVcol='ANOVA FEATURE FDR %',CMethod="Single",cancerTypes,SingleAnova=NULL){
  ModPPI<-NULL
  for(i in 1:length(directoryInput)){

    temp<-paste0(directoryInput[i],"/",outputSuffix)
    if(file.exists(temp)){
      load(temp)
      temp<-as.data.frame(TOTRES)
      temp$CTYPE<-cancerTypes[i]
      temp$MARKERCLASS<-"A"
      temp$PPI<-"TRUE"
      temp$CorMethod<-CMethod
      temp$individ_effectSize<-"N/A"
      temp$individ_pvals<-"N/A"
      temp$individ_fdrs<-"N/A"
      temp<-temp[temp[,AOVcol]<FDRthresh,]
      aovOut<-paste0(SingleAnova,cancerTypes[i],"/OUTPUT/ANOVA_results.rdata")
      if(file.exists(aovOut)&nrow(temp)>0){
        load(aovOut)
        TOTRES<-as.data.frame(TOTRES)
        for(i in 1:nrow(temp)){
          markers<-temp[i,"gene family"]
          markers<-unlist(strsplit(unlist(markers),"//"))
          rownames(TOTRES)<-paste0(TOTRES$FEATURE,TOTRES$`Depleted Gene`)
          combo<-paste0(temp[i,"FEATURE"],markers)
          temp[i,"individ_pvals"]<-paste(signif(as.numeric(unlist(TOTRES[combo,"FEATURE_ANOVA_pval"])),2),collapse="//")
          temp[i,"individ_fdrs"]<-paste(signif(as.numeric(unlist(TOTRES[combo,"ANOVA FEATURE FDR %"])),2),collapse="//")
          temp[i,"individ_effectSize"]<-paste(signif(as.numeric(unlist(TOTRES[combo,"FEATURE_ESS_effect_size"])),2),collapse="//")
          
        }
        
        
      }
      ModPPI<-rbind(ModPPI,temp)}
  }
  return(ModPPI)
}

compareCTNetMod<-function(NetModResult){
  #compare same results between the different cancer types - just looking at existence of regulatory networks/overlap:
  allCT<-unique(NetModResult[,"CTYPE"])
  comparepairs<-combn(allCT,2)
  regOver<-apply(comparepairs,2,function(x) compareNM(NetModResult,x))
  names(regOver)<-apply(comparepairs,2,function(x) paste0(x[1],x[2],collapse="-"))
  return(regOver)
}
compareNM<-function(NMRes,x){
  print(x[1])
  print(x[2])
  N1<-NMRes[NMRes[,"CTYPE"]==x[1],]
  N2<-NMRes[NMRes[,"CTYPE"]==x[2],]
  ResMat<-matrix(0,nrow=nrow(N1),ncol=nrow(N2))
  overlap<-apply(N1,1,function(x) setoverlap(unlist(x["gene family"]),unique(N2[,"gene family"])))
  if(is.matrix(overlap)){
  colnames(overlap)<-unlist(N1[,"gene family"])
  rownames(overlap)<-unlist(unique(N2[,"gene family"]))}
  return(overlap)
}
jaccardSets <- function(set1, set2){
  in_length = length(intersect(set1, set2))
  un_length = length(union(set1, set2))
  jaccard = in_length/un_length
  return(jaccard)
}
setoverlap<-function(invec,compareMat){
  set1<-unlist(strsplit(as.character(invec),split='//',fixed=TRUE))
  set2<-lapply(compareMat,function(x) strsplit(as.character(x),split='//',fixed=TRUE))
  output<-unlist(lapply(set2,function(x) jaccardSets(set1,unlist(x))))

  return(output)
}
compare2NetMod<-function(NetMod1,NetMod2){
  
}

annotateNetMod<-function(NetModRes,PPIgraph,PPInet,ProteinComplex){
  features<-unlist(NetModRes[,"FEATURE"])
  processfeat<-unlist(sapply(features,function(x) splitMarkers(x)))
  targets<-unlist(NetModRes[,"gene family"])
  allfeatures<-unique(sapply(processfeat,function(x) strsplit(as.character(x),"_",fixed=TRUE)[[1]][1]))
  alltargets<-unique(unlist(sapply(targets,function(x) strsplit(as.character(x),"//",fixed=TRUE))))
  ft<-unique(c(allfeatures,alltargets))
  print(ft)
  stringids<-PPInet[ft,"STRING_id"]
  stringids<-stringids[!is.na(stringids)]
  vertices<-which(as_ids(V(PPIgraph))%in%stringids)
  subnet<-induced_subgraph(PPIgraph,vertices)
  distancesG<-distances(PPIgraph)
  PPInet<-PPInet[PPInet$STRING_id%in%as_ids(V(subnet)),]
  PPIsp<-apply(NetModRes,1,function(x) shortestNM(x,distancesG,PPInet,ComplexMemb = ProteinComplex))
  return(PPIsp)
}
splitMarkers<-function(mkr){
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
    
  
  newmarkers<-markers
  return(newmarkers)
}

shortestNM<-function(NetModRes,distances,PPInet,ComplexMemb){
  set1<-strsplit(as.character(NetModRes['gene family']),split='//',fixed=TRUE)
  markers<-splitMarkers(NetModRes['FEATURE'])
  set1id<-PPInet[set1[[1]],"STRING_id"]
  set1id<-set1id[!is.na(set1id)]
  markerid<-PPInet[markers,"STRING_id"]
  markerid<-markerid[!is.na(markerid)]
  set1id<-intersect(set1id,rownames(distances))
  markerid<-intersect(markerid,colnames(distances))
  if(length(set1id)>0&length(markerid)>0){
    distances<-sapply(set1id,function(x) min(distances[x,markerid],na.rm=T))
    mindist<-min(distances,na.rm=T)
  }else{
    mindist<-NA
  }
  if(is.list(ComplexMemb)){
    InMemb<-unlist(lapply(ComplexMemb,function(x) compareComplex(set1,x,markers)))
  }else{
    InMemb<-compareComplex(set1,ComplexMemb,markers)
  }
  InMemb<-c(InMemb,minPPI=mindist)
  return(data.frame(InMemb))
}

compareComplex<-function(set1,ComplexMemb,markers){
  set1id<-intersect(set1[[1]],rownames(ComplexMemb))
  markerid<-intersect(markers,colnames(ComplexMemb))

  if(length(set1id)>0&length(markerid)>0){
    subset<-ComplexMemb[set1id,markerid]
    if(is.matrix(subset)){values<-subset[lower.tri(subset)]}else{values<-subset}
    vals<-sum(values=="0")
    if(vals<length(values)){
      InMemb<-"True"
    }else{InMemb<-"False"}
  }else{
    InMemb<-"N/A" 
  }
  return(InMemb)
}
NetModScores<-function(NetModRes,CTYPE,L2dir){
  allL2<-getL2scores(L2dir,CTYPE)
  targetscores<-apply(NetModRes,1,function(x) targetScores(x,allL2))
  NetModRes$PriorityL2<-unlist(targetscores)
  return(NetModRes)  
}
targetScores<-function(NetModRes,L2scores){
  targets<-strsplit(as.character(NetModRes["gene family"]),"//",fixed=TRUE)[[1]]
  L2s<-L2scores[targets]
  return(paste0(L2s,collapse="//"))
}
getL2scores<-function(L2dir,CTYPE,site="Combined"){
  load(file=paste(L2dir,CTYPE,'_L2.Rdata',sep=''))
 # paste(cellLines[i],c('_Mgk_10percFDR','_Mgk_5percFDR',
      #                 '_FOLD_1_sBF','_FOLD_2_sBF','_FOLD_3_sBF',
      #                 '_FOLD_minus2_FC','_FOLD_minus3_FC','_FOLD_minus5_FC',
       #                '_HIGHLY_EXP','_MUT','_COSMIC_MUT',
        #               '_IN_DEP_PATH'),sep=''))
  
  ns<-ncol(L2)
  nc<-ns/15
  
  L2score<-NULL
  NonNullL2score<-NULL
  
  for (j in 1:nc){

    
    currentBlock<-L2[,(15*(j-1)+1):(15*j)]
    
    currentL2filter<-(currentBlock[,1]>0 & currentBlock[,2]==0 & currentBlock[,3]==0)
    if(site=="Combined"){
      currentL2score<-( currentBlock[,6]*(1/6)*100+currentBlock[,7]*(1/6)*100+currentBlock[,8]*(1/6)*100+
                          currentBlock[,12]*(1/6)*100+currentBlock[,13]*(1/6)*100+currentBlock[,15]*(1/6)*100)
      #temp 7.12.20 just to see effect of weights on new version:
      #currentL2score<-( currentBlock[,6]*(3/8)*100+currentBlock[,7]*(1/8)*100+currentBlock[,8]*(1/8)*100+
      #                    currentBlock[,12]*(1/8)*100+currentBlock[,13]*(1/8)*100+currentBlock[,15]*(1/8)*100)
      #just use BF for the moment:
      currentL2score<-( currentBlock[,6]*(1/3)*100+currentBlock[,7]*(1/3)*100+currentBlock[,8]*(1/3)*100)
      
    }else{
      currentL2score<-(currentBlock[,4]*12.5+currentBlock[,5]*12.5+
                         currentBlock[,6]*12.5+currentBlock[,7]*12.5+currentBlock[,8]*12.5+
                         currentBlock[,12]*12.5+currentBlock[,13]*12.5+currentBlock[,15]*12.5)
    }
    
    currentL2score<-currentL2score*(currentL2filter+0)
    
    L2score<-cbind(L2score,currentL2score)
  }
  NonNullL2score<-rowSums(L2score>0)

  AvgNN2Score<-rowSums(L2score)/NonNullL2score
  AvgNN2Score[is.na(AvgNN2Score)]<-0
  FinalScore<-(NonNullL2score>2+0)*(AvgNN2Score)
  FinalScore<-round(FinalScore)
  return(FinalScore)
}
NetModData<-function(DynClust,Group,FCs,corMat,corL=0.6,BinaryDep,Amat=NULL,Gexp=NULL,ScaledBF=NULL,annotation){
  temp<-CLnameMapping(ScaledBF,Gexp,"col",annotation=cmp)
  ScaledBF<-temp$inputdata
  Gexp<-temp$refdata
  temp<-CLnameMapping(ScaledBF,BinaryDep,"col",annotation=cmp)
  BinaryDep<-temp$refdata
  
  
  genes<-unlist(DynClust[[1]][[Group]])
  if(is.null(genes)){
 
    #have pairs output:
    genes<-DynClust[[1]][Group,1:2]
    genes<-as.character(unlist(genes))
  }
  if(length(intersect(genes,rownames(FCs)))==length(genes)){
  alldata<-FCs[genes,]

  cm<-corMat[genes,genes]
  diag(cm)<-0
  cm[upper.tri(cm)]<-0

  if(length(genes)>2){
    aggData<-aggData(alldata,cm,BinaryDep,corL=corL,Amat,scaledBFs=ScaledBF,expressionProfile = Gexp)
  }else{
    direction<-sign(cm[2,1])
    outputdata<-combineGeneProfiles(direction,alldata,BinaryDep,Gexp=Gexp,BF=ScaledBF)
    #return(list(subdata=subdata,subbinary=subbinary,subgexp=Gexp,subBF=BF))
    aggData<-outputdata
    
    names(aggData)<-c("newdata","binarydata","Extra_Gexp","Extra_BF","namesdata")
    aggData$OutputMat<-cm
    aggData$newdataplus=NULL
    aggData$binarydataplus=NULL
  }
#  return(list(Extra_BF=Extra_BF,Extra_Gexp=Extra_Gexp,newdata=newdata,binarydata=binarydata,OutputMat=OutputMat,newdataplus=newdataplus,binarydataplus=binarydataplus))

  return(aggData)}else{return(NULL)}
}

aggData<-function(InputData,CorMat,BinaryDep,corL,Amat,scaledBFs=NULL,expressionProfile=NULL){
  subdata<-list()
  subbinary<-list()
  j=1
  namesdata<-NULL
  newdata<-InputData
  binarydata<-BinaryDep
  start=TRUE
  Finished<-FALSE
  subAmat<-Amat[rownames(CorMat),colnames(CorMat)]
  diag(CorMat)<-0
  CorMat<-CorMat*subAmat
  
  while(!is.null(newdata)&!Finished){
    if(!start){
      CorMat<-corMatBD(newdata,binarydata)
      #set variables to return:
      if(!is.null(subAmat)){
        CorMat<-CorMat*subAmat[rownames(CorMat),colnames(CorMat)]
      }
      #check to make sure have one new correlation above threshold or return multiple entries:
      diag(CorMat)<-0
      if(!max(abs(CorMat)>corL)){
        #this will stop the while loops
        OutputMat<-CorMat
        CorMat<-NULL
        Finished<-TRUE
        removesingletons<-rownames(keepdata)
        OutputMat<-OutputMat[!rownames(OutputMat)%in%removesingletons,!colnames(OutputMat)%in%removesingletons]
        newdata<-newdata[!rownames(newdata)%in%removesingletons,]
        binarydata<-binarydata[!rownames(binarydata)%in%removesingletons,]
      }
    }
    
    j=1
    subdata<-list()
    subbinary<-list()
    Extra_Gexp<-list()
    Extra_BF<-list()
    namesdata<-NULL
    keepdata<-newdata
    keepbinary<-binarydata
    while(!is.null(CorMat)){
      maxpair<-arrayInd(which.max(abs(as.matrix(CorMat))), dim(CorMat))
      maxCor<-max(abs(CorMat))
      
      if(maxCor>corL){
        r1<-rownames(CorMat)[maxpair[1]]
        r2<-colnames(CorMat)[maxpair[2]]
        direction<-sign(CorMat[r1,r2])
        pairdata<-newdata[c(r1,r2),]
        keepdata<-keepdata[!rownames(keepdata)%in%c(r1,r2),]
        bdep<-binarydata[c(r1,r2),]
        keepbinary<-binarydata[!rownames(binarydata)%in%c(r1,r2),]
        if(!is.null(expressionProfile)&!is.null(scaledBFs)){
          pairgexp<-expressionProfile[c(r1,r2),]
          pairBF<-scaledBFs[c(r1,r2),]
          combinedData<-combineGeneProfiles(direction,pairdata,bdep,Gexp=pairgexp,BF=pairBF)
          Extra_BF[[j]]<-combinedData$subBF
          Extra_Gexp[[j]]<-combinedData$subgexp
        }else{
          combinedData<-combineGeneProfiles(direction,pairdata,bdep,Gexp=NULL,BF=NULL)
        }
        
        subdata[[j]]<-combinedData$subdata
        subbinary[[j]]<-combinedData$subbinary
        
        if(direction==(-1)){
          
          newrow<-subAmat[r1,]|subAmat[r2,]
          subAmat<-subAmat[!rownames(subAmat)%in%c(r1,r2),!colnames(subAmat)%in%c(r1,r2)]
          incnames<-colnames(subAmat)
          subAmat<-rbind(subAmat,newrow[incnames])
          subAmat<-cbind(subAmat,c(newrow[incnames],0))
          namesdata[j]<-paste(paste0("(",r1),paste0(r2,")"),sep="|")
          rownames(subAmat)[nrow(subAmat)]<-namesdata[j]
          colnames(subAmat)[ncol(subAmat)]<-namesdata[j]
        }else{
          
          newrow<-subAmat[r1,]&subAmat[r2,]
          subAmat<-subAmat[!rownames(subAmat)%in%c(r1,r2),!colnames(subAmat)%in%c(r1,r2)]
          incnames<-colnames(subAmat)
          subAmat<-rbind(subAmat,newrow[incnames])
          subAmat<-cbind(subAmat,c(newrow[incnames],0))
          namesdata[j]<-paste(paste0("(",r1),paste0(r2,")"),sep="&")
          rownames(subAmat)[nrow(subAmat)]<-namesdata[j]
          colnames(subAmat)[ncol(subAmat)]<-namesdata[j]
        }
        
        CorMat<-CorMat[!rownames(CorMat)%in%c(r1,r2),!colnames(CorMat)%in%c(r1,r2)]
        
        j=j+1
        if(is.matrix(CorMat)){
          if(nrow(CorMat)<2){CorMat<-NULL}
        }else{CorMat<-NULL}
      }else{
        CorMat<-NULL
      }
      
    }
    if(!Finished){
      names(subdata)<-namesdata
      names(subbinary)<-namesdata
      newdata<-do.call(rbind,subdata)
      binarydata<-do.call(rbind,subbinary)
      newdataplus<-rbind(newdata,keepdata)
      binarydataplus<-rbind(binarydata,keepbinary)
      Extra_Gexp<-do.call(rbind,Extra_Gexp)
      Extra_BF<-do.call(rbind,Extra_BF)
      start=FALSE
    }
    
  }
  return(list(Extra_BF=Extra_BF,Extra_Gexp=Extra_Gexp,newdata=newdata,binarydata=binarydata,OutputMat=OutputMat,newdataplus=newdataplus,binarydataplus=binarydataplus))
}

corMatBD<-function(FCs,BinaryDep,method="pearson"){
  idxs<-rownames(FCs)
  tempout<-c()
  
  inidx<-length(idxs)
  
  for(k in 1:inidx){
    temp<-subcor(FCs,idxs[k],BinaryDep,method=method)
    temp$gene1<-idxs[k]
    tempout<-rbind(tempout,temp)
    
  }
  tempCor<-cor(t(FCs),method=method)
  corMat<-acast(tempout,gene1~gene2, value.var="corr")
  corMat<-corMat[,rownames(corMat)]
  temp<-corMat
  corMat<-t(corMat)
  corMat[upper.tri(corMat)]<-temp[upper.tri(temp)]
  diag(corMat)<-1
  corMat[is.na(corMat)]<-0
  
  corMatMax<-pmax(abs(tempCor),abs(corMat),na.rm=T)
  #sign of correlations:
  corSigns<-sapply(1:nrow(tempCor),function(x) which.pmax(abs(tempCor[x,]),abs(corMat[x,])))
  corSigns<-t(corSigns)
  s1<-sign(tempCor)
  s2<-sign(corMat)
  signMat<-s1
  signMat[corSigns==2]<-s2[corSigns==2]
  corMat<-corMatMax*signMat
  return(corMat)
}

combineGeneProfiles<-function(direction,FC,Binary,Gexp=NULL,BF=NULL){
  #  load(paste0(inputdirectoryAnova,"/Extra_sBFs.Rdata"))
  # load(paste0(inputdirectoryAnova,"/Extra_correctedFCs.Rdata"))
  #  load(paste0(inputdirectoryAnova,"/Extra_GExp.Rdata"))
  
  r1<-rownames(FC)[1]
  r2<-rownames(FC)[2]
  Binary<-Binary[c(r1,r2),]
  if(!direction){
    #correlation is negative
    subdata<-apply(FC,2,min)
    #for gene expression and BF select for each cell line the one corresponding to the gene used for the FC:
    
    if(!is.null(Gexp)){
      cl<-intersect(colnames(FC),colnames(Gexp))
      FCu<-FC[c(r1,r2),cl]
        if(length(intersect(rownames(Gexp),c(r1,r2)))==2){
          Gexp<-Gexp[c(r1,r2),cl]
          whichgene<-apply(FCu,2,which.min)
          Gexp<-sapply(1:ncol(FCu),function(x) Gexp[whichgene[x],x])
        }else{
          Gexp<-rep(NA,ncol(FCu))
          names(Gexp)<-colnames(FCu)
        }
    }
    if(!is.null(BF)){
      cl<-intersect(colnames(FC),colnames(BF))
      BFu<-BF[c(r1,r2),cl]
      whichgene<-apply(BFu,2,which.max)

      BF<-sapply(1:ncol(BFu),function(x) BFu[whichgene[x],x])}
    
    subbinary<-apply(Binary,2,function(x) as.numeric(x[1]|x[2]+0))
    
    namesdata<-paste(paste0("(",r1),paste0(r2,")"),sep="|")
    
  }else{
    #correlation is positive
    subdata<-apply(FC,2,mean)
    subbinary<-apply(Binary,2,function(x) as.numeric(x[1]&x[2]+0))
    if(!is.null(Gexp)){
      cl<-intersect(colnames(FC),colnames(Gexp))
      if(length(intersect(rownames(Gexp),c(r1,r2)))==2){
        Gexp<-Gexp[c(r1,r2),cl]
        Gexp<-apply(Gexp,2, mean)}else{
          Gexp<-rep(NA,length(cl))
          names(Gexp)<-cl
        }
    }
    if(!is.null(BF)){
      cl<-intersect(colnames(FC),colnames(BF))
   
      BF<-BF[c(r1,r2),cl]
      BF<-apply(BF,2, mean)
 
    }
    
    namesdata<-paste(paste0("(",r1),paste0(r2,")"),sep="&")
    
  }

  return(list(subdata=subdata,subbinary=subbinary,subgexp=Gexp,subBF=BF,namesdata=namesdata))
}



