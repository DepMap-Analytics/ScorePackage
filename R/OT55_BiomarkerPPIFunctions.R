

Combine_SingleDualBMPPI<-function(AllBM_Single,AllBM_Dual){
  ctypes<-union(unique(AllBM_Single$CancerType),unique(AllBM_Dual$CancerType))
  AllOut<-NULL
  for(i in ctypes){
    if(i%in%AllBM_Single$CancerType&i%in%AllBM_Dual$CancerType){
      SubSingle<-AllBM_Single[AllBM_Single$CancerType==i,]
      SubDual<-AllBM_Dual[AllBM_Dual$CancerType==i,]
      bothtargets<-intersect(SubSingle$Target,SubDual$Target)
      if(length(bothtargets)>0){
        if(length(!SubSingle$Target%in%bothtargets)>0){
          AllOut<-rbind(AllOut,SubSingle[!SubSingle$Target%in%bothtargets,])}
        if(length(!SubDual$Target%in%bothtargets)>0){
          AllOut<-rbind(AllOut,SubDual[!SubDual$Target%in%bothtargets,])
        }
        for(j in bothtargets){
          classes<-c(SubSingle[SubSingle$Target==j,"MaxTargetScoreFDR"],SubDual[SubDual$Target==j,"MaxTargetScoreFDR"])
          if(classes[1]==classes[2]){
            combineRes<-SubSingle[SubSingle$Target==j,]
            combineRes$BMSymbol<-paste0(c(SubSingle[SubSingle$Target==j,"BMSymbol"],SubDual[SubDual$Target==j,"BMSymbol"]),collapse="//")
            combineRes$BMType<-paste0(c(SubSingle[SubSingle$Target==j,"BMType"],SubDual[SubDual$Target==j,"BMType"]),collapse="//")
            combineRes$BMs<-paste0(c(SubSingle[SubSingle$Target==j,"BMs"],SubDual[SubDual$Target==j,"BMs"]),collapse="//")

            AllOut<-rbind(AllOut,combineRes)
          }else{
            maxbm<-which.max(as.numeric(classes))
            ifelse(maxbm==1,temp<-SubSingle[SubSingle$Target==j,],temp<-SubDual[SubDual$Target==j,])
            AllOut<-rbind(AllOut,temp)
          }
        }
      }else{
        AllOut<-rbind(AllOut,SubSingle)
        AllOut<-rbind(AllOut,SubDual)
      }
    }else{
      AllOut<-rbind(AllOut,AllBM_Single[AllBM_Single$CancerType==i,])
      AllOut<-rbind(AllOut,AllBM_Dual[AllBM_Dual$CancerType==i,])
    }
  }
  return(AllOut)
}



splitDependencies<-function(res){
  depgenes<-unlist(strsplit(unlist(res["Depleted Gene"]),"[|&]"))
  string_id<-unlist(strsplit(unlist(res["string_id"]),"//"))
  out<-rbind(res,res)
  out[,"Depleted Gene"]<-depgenes
  out[,"string_id"]<-string_id
  pairlist<-c(depgenes,res["assoc_id"],res["Depleted Gene"])
  names(pairlist)<-c("Dep1","Dep2","assoc_id","pair")
  return(list(out=out,pairlist=pairlist))
}

splitMultiFeature<-function(res,type){

  out<-tryCatch({depgenes<-unlist(strsplit(unlist(res["FEATURE"]),"[|&]"))
    string_id<-unlist(strsplit(unlist(res["BMstringID"]),"[|&]"))
    out<-rbind(res,res)
    out[,"FEATURE"]<-depgenes
    out[,"BMstringID"]<-string_id
    if(length(string_id)==2){
      pairlist<-c(depgenes,res["assoc_id"],res["FEATURE"],string_id,type)
    }else{
      pairlist<-c(depgenes,res["assoc_id"],res["FEATURE"],rep(string_id,2),type)
    }
    names(pairlist)<-c("Feat1","Feat2","assoc_id","pair","BM1","BM2","DualBM")
  })

  if(class(out)!="try-error"){
    return(list(out=out,pairlist=pairlist))
  }else{
    return(list(out=NULL,pairlist=NULL))
  }
}

getBMPPIdata<-function(BMPPIdir,CTYPE,ScoreType="Avg",GENES,biomarkertype="mut",omictype="",TargetData=NULL){
  if(!is.null(BMPPIdir)){
    TargetData<-read.table(paste0(BMPPIdir,"/RWR_ScoreTargets.txt"),sep="\t",stringsAsFactors = FALSE,header=T)
  }
  TargetData<-TargetData[TargetData$CancerType==CTYPE,]
  TargetData<-TargetData[TargetData$ScoreType==ScoreType,]
  TargetSplit<-split.data.frame(TargetData,f=TargetData$BMType)
  numberBMs<-names(TargetSplit)
  if(biomarkertype=="CCN"){biomarkertype<-"CN"}
  if(biomarkertype=="CProt"){biomarkertype<-"Prot"}
  if(biomarkertype=="CExpr"){biomarkertype<-"Expr"}
  if(biomarkertype=="CMet"){biomarkertype<-"Met"}
    TD<-matrix(0,nrow=length(GENES),ncol=11)
    colnames(TD)<-c("Score","MaxTargetScoreFDR","FDRScore","PvalScore","RWRscore","NumberBMs","StringScore","BMs","RWRclasses","PvalClasses","FDRClasses")
    rownames(TD)<-GENES
    TD<-data.frame(TD,stringsAsFactors = FALSE)
    Inputdata<-data.frame(TargetSplit[[biomarkertype]],stringsAsFactors = FALSE)

  if(!is.null(Inputdata)){
    incg<-intersect(Inputdata$Target,GENES)
    rownames(Inputdata)<-Inputdata$Target
    if(length(incg)>0){
      TD[incg,]<-Inputdata[incg,colnames(TD)]
      TD[,"biomarkertype"]<-biomarkertype
    }



  return(TD)}else{return(NULL)}
}

getBMPPIannot<-function(BMPPIdir=NULL,CTYPE,ScoreType="Avg",GENES,biomarkertypes,selScore="FDRScore",TargetData=NULL){
  if(!is.null(BMPPIdir)){
    TargetData<-read.table(paste0(BMPPIdir,"/RWR_ScoreTargets.txt"),sep="\t",stringsAsFactors = FALSE,header=T)
  }
  TargetData<-TargetData[TargetData$CancerType==CTYPE,]
  TargetData<-TargetData[TargetData$ScoreType==ScoreType,]
  TargetSplit<-split.data.frame(TargetData,f=TargetData$BMType)
  numberBMs<-names(TargetSplit)
  output<-NULL

  for(i in 1:length(biomarkertypes)){
    biomarkertype<-biomarkertypes[i]
    if(biomarkertype=="CCN"){biomarkertype<-"CN"}
    if(biomarkertype=="CProt"){biomarkertype<-"Prot"}
    if(biomarkertype=="CExpr"){biomarkertype<-"Expr"}
    if(biomarkertype=="CMet"){biomarkertype<-"Met"}


    Inputdata<-data.frame(TargetSplit[[biomarkertype]],stringsAsFactors = FALSE)

    if(!is.null(Inputdata)){
      incg<-intersect(Inputdata$Target,GENES)
      rownames(Inputdata)<-Inputdata$Target

      colnamesTD<-c("Score","MaxTargetScoreFDR",selScore,"PvalScore","RWRscore","NumberBMs","StringScore","BMs","PvalClasses","RWRclasses")


      if(length(incg)>0){
        TD<-Inputdata[incg,colnamesTD]
        TD[,"biomarkertype"]<-biomarkertype

        output<-rbind(output,TD)
      }


    }
  }
  if(!is.null(output)){
    if(nrow(output)>1){
      maxScore<-max(output[,selScore])

      output<-output[output[,selScore]==maxScore,]

    }
  }
    return(output)
}


getWeightedScore<-function(Res,PPIigraphWeighted,PPIgraph){
  Target<-unlist(Res["TargetID"])
  Targets<-unlist(strsplit(Target,"//",fixed=TRUE)[[1]])
  BMs<-unlist(Res["BMs"])
  BMs<-unlist(strsplit(BMs,'//',fixed=TRUE)[[1]])
  BMs<-intersect(BMs,rownames(PPIigraphWeighted))
  avgweight<-list()
  avgperT<-list()
  print(Res)
  for(k in 1:length(Targets)){
    Target<-Targets[k]
    for(i in 1:length(BMs)){
      BMID<-BMs[i]
      sps<-igraph:::shortest_paths(PPIgraph,Target,BMID,output="vpath")

      vpath<-sps$vpath[[1]]
      vpath<-V(PPIgraph)$name[as.numeric(vpath)]
      weights<-list()
      if(length(vpath)>1){
        for(j in 1:(length(vpath)-1)){

          if(sum(vpath[j:(j+1)]%in%colnames(PPIigraphWeighted))==2){
            weights[[j]]<-PPIigraphWeighted[vpath[j],vpath[j+1]]
          }else{
            weights[[j]]<-NA
          }
        }
        avgweight[[i]]<-mean(unlist(weights),na.rm=T)
      }else{
        avgweight[[i]]<-1000
      }
    }
    avgperT[[k]]<-mean(unlist(avgweight))
  }
  return(mean(unlist(avgperT)))
}
normalizeW<-function(adjmat,degreelist,normBy=c("Column","Both","Row")){
  normBy <- match.arg(normBy)
  if(normBy=="Row"){
    dmat<-adjmat/degreelist
  }
  if(normBy=="Column"){
    degreemat<-matrix(degreelist,nrow=length(degreelist),ncol=length(degreelist),byrow=TRUE)
    dmat<-adjmat/degreemat
  }
  if(normBy=="Both"){
    degreemat<-matrix(degreelist,nrow=length(degreelist),ncol=length(degreelist),byrow=TRUE)
    denom<-degreelist*degreemat
    denom<-sqrt(denom)
    dmat<-adjmat/denom
  }

  return(dmat)
}
rwrset<-function(startingmatrix,adjacencymatrix,restartprob,degreelist,PPInet=NULL,normBy=c("Column","Both","Row","none"),idxs=NULL){
  normBy<-match.arg(normBy)

  if(normBy!="none"){
    normW<-normalizeW(adjacencymatrix,degreelist,normBy=normBy)
  }else{
      normW<-adjacencymatrix
  }
  ni<-length(idxs)
  if(!is.null(idxs)){
    allRankList<-matrix(NA,nrow=nrow(startingmatrix),ncol=ncol(startingmatrix))
    for(i in 1:ni){
      input<-startingmatrix[,idxs[[i]]]
      temp<-apply(input,2,function(x) rwr(startingvec = x,normW,restartprob = restartprob,PPInet = PPInet))
      allRankList[,idxs[[i]]]<-temp
    }
  }else{
    allRankList<-apply(startingmatrix,2,function(x) rwr(startingvec = x,normW,restartprob = restartprob,PPInet = PPInet))
  }
  rownames(allRankList)<-rownames(adjacencymatrix)
  colnames(allRankList)<-colnames(startingmatrix)
  return(allRankList)
}
rwr<-function(startingvec,normW,restartprob,convergetol=1e-012,maxiter=100,PPInet=NULL){
  #BIOKDD Can 2005 Protein-Protein Interaction Networks using Random Walk

  iter<-1
  startval<-restartprob*startingvec
  CurrVal<-startingvec
  diff<-1
  while(diff>convergetol&iter<maxiter){
    NewVal<-(1-restartprob)*normW%*%CurrVal+startval
    diff<-sum(abs(NewVal-CurrVal))
    CurrVal<-NewVal
    iter<-iter+1

  }

  #NewVal<-NewVal[order(NewVal[,1],decreasing=TRUE),1]

  if(!is.null(PPInet)){
    rownames(NewVal)<-PPInet[match(rownames(NewVal),PPInet$STRING_id),"symbol"]
  }
  return(NewVal)
}

rwrcost<-function(finalprob,pvaluemat){

}

heatdiff<-function(startvec,graph,time){
  Laplacian<-laplacian_matrix(graph)
  return(startvec*exp(-Laplacian*time))
}

heatstart<-function(endvec,graph,time){
  Laplacian<-laplacian_matrix(graph)
  return(endvec*inv(exp(-Laplacian*time)))
}
msplit<-function(markerrow,subdir="",PPInet){
  mkr<-markerrow["FEATURE"]
  ftype<-NULL
  splitC<-unlist(sapply(unlist(mkr),function(x) strsplit(x,"[()]")))
  if(length(splitC)>1){
    #have a cna marker with associated genes
    #find things that are not full markers:
    cnaList<-!grepl("mut|Expr|var",splitC)
    splitC[cnaList]<-gsub(",",";",splitC[cnaList])
    splitC<-paste0(splitC,collapse="")
  }
  splitm<-sapply(unlist(splitC),function(x) strsplit(x,",",fixed=TRUE))
  splitm<-unlist(splitm)

    if(subdir=="EMDE"){
      markers<-gsub(" ","",splitm,fixed=T)
    }else{
      splitM<-sapply(splitm,function(x) strsplit(x,"_",fixed=TRUE))
      type<-unlist(lapply(splitM,function(x) x[length(x)]))
      type[!type%in%c("mut","Expr","var")]<-"CNA"
      nmarkers<-length(type)
      markers<-c()
      j=1
      for(i in 1:nmarkers){
        if(type[i]%in%c("mut","Expr","var")){
          markers[j]<-splitM[[i]][1]
          ftype[j]<-type[i]
          j=j+1
        }else{
          possCNAmarker<-splitM[[i]][length(splitM[[i]])]
          if(possCNAmarker!=""){

            CNAmarker<-unlist(strsplit(possCNAmarker,";",fixed=TRUE))
            nM<-length(CNAmarker)
            markers[j:(j+nM-1)]<-paste0(splitM[[i]][1],"(",CNAmarker,")")
            ftype[j:(j+nM-1)]<-type[i]
            j=j+nM
          }

        }
      }
      markers<-gsub(" ","",markers, fixed = T)

    }
  if(length(markers)>1){
    mgene<-markers
    ppi<-unlist(sapply(markerrow["ppi_target"],function(x) strsplit(x,",",fixed=TRUE)))
    newmarkers<-splitm
    newRes<-matrix(markerrow,nrow=1)
    colnames(newRes)<-colnames(markerrow)
    for(i in 2:length(newmarkers)){
      newRes<-rbind(newRes,markerrow)
    }
    newRes[,"FEATURE"]<-newmarkers
    newRes[,"markers"]<-mgene
    newRes[,"ppi_target"]<-ppi
    newRes[,"TstringID"]<-unlist(PPInet[unlist(newRes[,"Depleted.Gene"]),"STRING_id"])
    newRes[,"StringID"]<-unlist(PPInet[unlist(mgene),"STRING_id"])
    newRes[,"type"]<-ftype

  }else{
    #don't have a multi "," biomarker
    markerrow["markers"]<-markers
    markerrow["TstringID"]<-PPInet[unlist(markerrow["Depleted.Gene"]),"STRING_id"]
    markerrow["StringID"]<-PPInet[unlist(markerrow["markers"]),"STRING_id"]
    type<-strsplit(as.character(markerrow["FEATURE"]),"_",fixed=TRUE)[[1]]
    type<-type[length(type)]
    if(!type%in%c("mut","Expr")){type<-"CNA"}
    markerrow["type"]<-type
    newRes<-matrix(markerrow,nrow=1)



  }
  return(newRes)
}
msplitAov<-function(markerrow,subdir="",PPInet,filename){

  mkr<-markerrow["FEATURE"]
  ftype<-NULL

  singleMarker<-c()
  markers<-c()
  splitC<-unlist(sapply(unlist(mkr),function(x) strsplit(x,"[()]")))
  if(length(splitC)>1){
    #have a cna marker with associated genes
    #find things that are not full markers:
    cnaList<-!grepl("mut|Expr|Met|Copy|Prot",splitC)
    splitC[cnaList]<-gsub(",",";",splitC[cnaList])
    splitC<-paste0(splitC,collapse="")
  }
  splitm<-sapply(unlist(splitC),function(x) strsplit(x,",",fixed=TRUE))
  splitm<-unlist(splitm)
  splitid<-sapply(unlist(markerrow["BMstringID"]),function(x) strsplit(x,",",fixed=TRUE))
  splitid<-unlist(splitid)
    #have more than one biomarker
    if(subdir=="EMDE"){
      markers<-gsub(" ","",splitm,fixed=T)
    }else{
      splitM<-sapply(splitm,function(x) strsplit(x,"_",fixed=TRUE))
      type<-unlist(lapply(splitM,function(x) x[length(x)]))
      type[!type%in%c("mut","Expr","Met","Copy","Prot")]<-"CNA"
      nmarkers<-length(type)
      markers<-c()
      j=1
      if(nmarkers==1){
        temp<-splitM
        splitM[[1]]<-unlist(temp)
      }
      for(i in 1:nmarkers){
        if(type[i]%in%c("mut","Expr","Met","CN","Prot")){
          singleMarker[j]<-splitm[i]
          ftype[j]<-type[i]

          markers[j]<-unlist(splitM[[i]])[1]

          j=j+1
        }else{
          possCNAmarker<-splitM[[i]][length(splitM[[i]])]
          if(possCNAmarker!=""){

            CNAmarker<-unlist(strsplit(possCNAmarker,";",fixed=TRUE))
            nM<-length(CNAmarker)
            singleMarker[j:(j+nM-1)]<-paste0(splitM[[i]][1],"(",CNAmarker,")")
            markers[j:(j+nM-1)]<-CNAmarker
            ftype[j:(j+nM-1)]<-type[i]
            j=j+nM
          }

        }
      }
      markers<-gsub(" ","",markers, fixed = T)
      singleMarker<-gsub(" ","",singleMarker,fixed=T)
    }

  output<-markerrow
  for(i in 2:length(markers)){
    output<-rbind(output,output)
  }
  output[,"markers"]<-markers
  output[,'FEATURE']<-splitm
  output[,"BMstringID"]<-splitid
  type<-ftype
  ppi_target<-rep("",nrow(output))
  output<-data.frame(output,ppi_target,type,stringsAsFactors = FALSE)
  return(output)
}

annotateAovPPI<-function(AnovaRes,PPInet,pvalCol="FEATURE_ANOVA_pval",getPPIid=TRUE,filename=NULL){
  if(getPPIid){
    for(i in 1:ncol(AnovaRes)){
      AnovaRes[,i]<-unlist(AnovaRes[,i])
    }
    splitM<-apply(AnovaRes,1,function(x) msplitAov(x,PPInet=PPInet,filename=filename))
    if(is.null(filename)){
      allMarkers<-do.call(rbind,splitM)
    }else{
      allMarkers<-NULL
    }
  }else{
    colnames(AnovaRes)[colnames(AnovaRes)=="string_id"]<-"TstringID"
    singlemarkers<-AnovaRes[grep(",",unlist(AnovaRes[,"BMstringID"]),fixed=T,invert=T),]
    singlemarkers<-data.frame(singlemarkers,ppi_target="",stringsAsFactors = FALSE)
    markers<-unlist(sapply(unlist(singlemarkers[,"FEATURE"]),function(x) strsplit(as.character(x),"_",fixed=T)[[1]][1]))
    gsub("(","",markers,fixed=T)
    gsub(")","",markers,fixed=T)
    singlemarkers<-data.frame(singlemarkers,markers,stringsAsFactors = FALSE)
    splitM<-sapply(unlist(singlemarkers[,"FEATURE"]),function(x) strsplit(as.character(x),"_",fixed=TRUE))
    type<-unlist(lapply(splitM,function(x) x[length(x)]))
    singlemarkers<-data.frame(singlemarkers,type,stringsAsFactors = FALSE)


    multimarkers<-grep(",",unlist(AnovaRes[,"BMstringID"]),fixed=T)
    if(length(multimarkers)>0){
      splitmarkers<-list()
      count<-1
      for(k in multimarkers){
        splitmarkers[[count]]<-msplitAov(AnovaRes[k,],PPInet=PPInet)
        count<-count+1
      }

      allmulti<-do.call(rbind,splitmarkers)

      onames<-colnames(singlemarkers)
      colnames(singlemarkers)<-make.names(colnames(singlemarkers))
      allMarkers<-rbind(singlemarkers,allmulti[,colnames(singlemarkers)])
      colnames(allMarkers)<-onames
      colnames(allMarkers)[colnames(allMarkers)=="BMstringID"]<-"StringID"

    }else{
      allMarkers<-singlemarkers
      colnames(allMarkers)[colnames(allMarkers)=="BMstringID"]<-"StringID"
    }

  }
  if(is.null(filename)){
  allMarkers<-allMarkers[!(is.na(allMarkers[,"StringID"])|is.na(allMarkers[,"TstringID"])),]
  colnames(allMarkers)[colnames(allMarkers)==pvalCol]<-"AnovaPval"
  colnames(allMarkers)<-make.names(colnames(allMarkers))
  allMarkers$FEATURE<-gsub(" ","",allMarkers$FEATURE)
  }
  return(allMarkers)
}
splitmultimarker<-function(bmres,subdir=""){
  mkr<-bmres["FEATURE"]
  ftype<-NULL

  singleMarker<-c()
  markers<-c()
  splitC<-unlist(sapply(unlist(mkr),function(x) strsplit(x,"[()]")))
  if(length(splitC)>1){
    #have a cna marker with associated genes
    #find things that are not full markers:
    cnaList<-!grepl("mut|Expr|Met|Copy|Prot",splitC)
    splitC[cnaList]<-gsub(",",";",splitC[cnaList])
    splitC<-paste0(splitC,collapse="")
  }
  splitm<-sapply(unlist(splitC),function(x) strsplit(x,",",fixed=TRUE))
  splitm<-unlist(splitm)

  #have more than one biomarker
  if(subdir=="EMDE"){
    markers<-gsub(" ","",splitm,fixed=T)
  }else{
    splitM<-sapply(splitm,function(x) strsplit(x,"_",fixed=TRUE))
    type<-unlist(lapply(splitM,function(x) x[length(x)]))
    type[!type%in%c("mut","Expr","Met","Copy","Prot")]<-"CNA"
    nmarkers<-length(type)
    markers<-c()
    j=1
    if(nmarkers==1){
      temp<-splitM
      splitM[[1]]<-unlist(temp)
    }
    for(i in 1:nmarkers){
      if(type[i]%in%c("mut","Expr","Met","CN","Prot")){
        singleMarker[j]<-splitm[i]
        ftype[j]<-type[i]

        markers[j]<-unlist(splitM[[i]])[1]

        j=j+1
      }else{
        possCNAmarker<-splitM[[i]][length(splitM[[i]])]
        if(possCNAmarker!=""){

          CNAmarker<-unlist(strsplit(possCNAmarker,";",fixed=TRUE))
          nM<-length(CNAmarker)
          singleMarker[j:(j+nM-1)]<-paste0(splitM[[i]][1],"(",CNAmarker,")")
          markers[j:(j+nM-1)]<-CNAmarker
          ftype[j:(j+nM-1)]<-type[i]
          j=j+nM
        }

      }
    }
    markers<-gsub(" ","",markers, fixed = T)
    singleMarker<-gsub(" ","",singleMarker,fixed=T)
  }

    newRes<-matrix(bmres,nrow=length(singleMarker),byrow=T)


    newRes[,"FEATURE"]<-singleMarker
    newFeat<-data.frame(markers=markers,ppi_target="",type=ftype,stringsAsFactors = FALSE)



    newRes[,"BMstringID"]<-unlist(strsplit(bmres[,"BMstringID"],",",fixed=TRUE))
    newRes<-cbind(newRes,newFeat)


  return(newRes)
}

graphMatricesRWR<-function(PPIigraph,PPInet,subdir="Genomic",single=TRUE,AnovaRes,pvalCol="FEATURE_ANOVA_pval",fdrCol="ANOVA FEATURE FDR %",
                           deltaPosCol="FEATUREpos_Glass_delta",deltaNegCol="FEATUREneg_Glass_delta",Int=FALSE,dir.Results=NULL,IntSplit=FALSE,HubCorrect=FALSE){
  allMarkers<-AnovaRes

  colnames(allMarkers)[colnames(allMarkers)=='Depleted Gene']<-"Depleted.Gene"

  #PPI prior matrix:
  allMarkers[,"Depleted.Gene"]<-as.character(allMarkers[,"Depleted.Gene"])
  allMarkers[,"FEATURE"]<-as.character(unlist(allMarkers[,"FEATURE"]))
  #allMarkers$markers<-as.character(allMarkers$markers)
  #allMarkers[,"CLASS"]<-as.character(allMarkers[,"CLASS"])
  allMarkers[,"assoc_id"]<-as.character(allMarkers[,"assoc_id"])
  type<-rep("",nrow(allMarkers))
  DualBM<-rep("",nrow(allMarkers))
  markers<-rep(NA,nrow(allMarkers))
  ppi_target<-rep(" ",nrow(allMarkers))
  allMarkers<-data.frame(allMarkers,type,DualBM,markers,ppi_target,stringsAsFactors = FALSE)
  #if(subdir%in%c("Genomic","")){

   # allMarkers$TstringID<-""
   # allMarkers$StringID<-""

  #}else{
    colnames(allMarkers)[which(colnames(allMarkers)=="BMstringID")]<-"StringID"
    colnames(allMarkers)[which(colnames(allMarkers)=="string_id")]<-"TstringID"

  #}
    #allMarkers$AnovaFdr<-AnovaRes[match(paste0(allMarkers$Depleted.Gene,allMarkers$FEATURE,sep="-"),paste0(AnovaRes[,"Depleted Gene"],AnovaRes[,"FEATURE"],sep="-")),fdrCol]
    #allMarkers$AnovaPval<-AnovaRes[match(paste0(allMarkers$Depleted.Gene,allMarkers$FEATURE,sep="-"),paste0(AnovaRes[,"Depleted Gene"],AnovaRes[,"FEATURE"],sep="-")),pvalCol]
    #allMarkers$AnovaDeltaP<-AnovaRes[match(paste0(allMarkers$Depleted.Gene,allMarkers$FEATURE,sep="-"),paste0(AnovaRes[,"Depleted Gene"],AnovaRes[,"FEATURE"],sep="-")),deltaPosCol]
    #allMarkers$AnovaDeltaN<-AnovaRes[match(paste0(allMarkers$Depleted.Gene,allMarkers$FEATURE,sep="-"),paste0(AnovaRes[,"Depleted Gene"],AnovaRes[,"FEATURE"],sep="-")),deltaNegCol]

    AnovaFdr<-unlist(AnovaRes[,fdrCol])
    AnovaPval<-unlist(AnovaRes[,pvalCol])
    AnovaDeltaP<-unlist(AnovaRes[,deltaPosCol])
    AnovaDeltaN<-unlist(AnovaRes[,deltaNegCol])

    allMarkers<-cbind(allMarkers,AnovaFdr,AnovaPval,AnovaDeltaP,AnovaDeltaN)
   if(!Int){
    #need to un collapse markers here:
    #composite markers:
    cmarkers<-grep(",",allMarkers[,"FEATURE"],fixed=T)
    smarkers<-grep(",",allMarkers[,"FEATURE"],fixed=T,invert = T)

    if(length(cmarkers)>0){
      splitM<-list()
      count<-1
      for(k in cmarkers){
        splitM[[count]]<-msplit(allMarkers[k,],PPInet=PPInet)
        count<-count+1
      }

      if(is.list(splitM)){
        allCMarkers<-do.call(rbind,splitM)
        allMarkers<-rbind(allMarkers[smarkers,],allCMarkers)
        }else{
        allMarkers<-allMarkers[smarkers,]
          for(k in 1:ncol(splitM)){
            temp<-splitM[,k]
            allMarkers<-rbind(allMarkers,unique(temp))
          }

        }
    }
   }else{
     load(file=paste0(dir.Results,"/PairList_Compound.Rdata"))
     #pairlist has STRING IDS for two biomarkers in BM1 and BM2
     pairList<-as.matrix(Allpair[,c("pair","BM1","BM2","DualBM")])

     if(is.matrix(pairList)){
       if(length(unique(pairList))==1){
         pairList<-matrix(unique(pairList),nrow=1,ncol=4,byrow = T)
         colnames(pairList)<-c("pair","BM1","BM2","DualBM")
       }else{
        pairList<-unique(pairList)
       }
     }
     print(dim(pairList))

   }

  if(Int){
    childFeat<-unlist(sapply(Allpair[,"Feat2"],
                  function(x) strsplit(x,"_",fixed=TRUE)[[1]][length(strsplit(x,"_",fixed=TRUE)[[1]])]))


    nm<-length(childFeat)
    type<-sapply(1:nm,function(x) paste0(Allpair[x,"DualBM"],childFeat[x],collapse="_"))
    allMarkers[,"type"]<-type

    allMarkers[,"DualBM"]<-unlist(Allpair[,"DualBM"])
  }else{
    allMarkers[,"type"]<-unlist(sapply(allMarkers[,"FEATURE"],
                                       function(x) strsplit(x,"_",fixed=TRUE)[[1]][length(strsplit(x,"_",fixed=TRUE)[[1]])]))
    allMarkers[,"DualBM"]<-"None"
  }

  if(Int){
    keep<-!(is.na(allMarkers[,'TstringID'])|(rowSums(Allpair[,c("BM1","BM2")]=="NA"))!=0)
    allMarkers<-allMarkers[keep,]

    Allpair<-Allpair[keep,]

  }else{
    allMarkers<-allMarkers[!(is.na(allMarkers[,"StringID"])|is.na(allMarkers[,"TstringID"])),]

  }

  NullCheck<-unlist(apply(allMarkers,1,function(x) sum(is.null(unlist(x[c("AnovaFdr","AnovaPval","AnovaDeltaP","AnovaDeltaN")])))))
  allMarkers<-allMarkers[NullCheck==0,]
  if(Int){
    Allpair<-Allpair[NullCheck==0,]
    check2<-sum(pairList[,"pair"]%in%Allpair[,"pair"])
    pairList<-pairList[pairList[,"pair"]%in%Allpair[,"pair"],]
    if(check2==1){
      pairList<-matrix(pairList,nrow=1,ncol=4,byrow = T)
      colnames(pairList)<-c("pair","BM1","BM2","DualBM")
    }
    print(dim(pairList))
  }


    UniqueBiomarker<-unique(allMarkers[,"FEATURE"])

  UniqueTargets<-unique(allMarkers[,"Depleted.Gene"])

  graphIDs<-vertex_attr(PPIigraph, "name")
  allAdjMat<-as.matrix(as_adjacency_matrix(PPIigraph))
  diag(allAdjMat)<-0
  if(Int){
      if(is.matrix(pairList)){
      BiomarkerString<-unique(as.vector(unlist(pairList[,c("BM1","BM2")])))
      }else{
        BiomarkerString<-pairList[c("BM1","BM2")]
      }
  }else{
    BiomarkerString<-unique(allMarkers[,"StringID"])
  }

  BMandTarget<-c(BiomarkerString,unique(AnovaRes[,"string_id"]))
  cts<-components(PPIigraph,"strong")
  maxct<-which.max(cts$csize)
  incvec<-names(cts$membership[cts$membership==maxct])

  adjmat<-allAdjMat[incvec,incvec]
  BiomarkerString<-intersect(BiomarkerString,incvec)

  inputID<-na.omit(unique(c(incvec,BiomarkerString)))
  names(inputID)<-PPInet[match(inputID,PPInet$STRING_id),"symbol"]
  inputID<-inputID[which(!is.na(names(inputID)))]
  #inputMat<-induced_subgraph(PPIigraph,which(as.vector(graphIDs)%in%inputID))
  WeightMat<-as.matrix(as_adjacency_matrix(PPIigraph))
  checkID<-rownames(WeightMat)%in%PPInet$STRING_id
  WeightMat<-WeightMat[checkID,]
  checkID<-colnames(WeightMat)%in%PPInet$STRING_id
  WeightMat<-WeightMat[,checkID]

  useID<-intersect(colnames(WeightMat),rownames(WeightMat))
  WeightMat<-WeightMat[useID,useID]
  degreelist<-colSums(WeightMat)

  incvec<-names(degreelist)[degreelist!=0]

  degreelist<-degreelist[incvec]

  WeightMat<-WeightMat[incvec,incvec]

  BiomarkerString<-intersect(BiomarkerString,incvec)

  featname<-"FEATURE"
  if(Int){
    keep<-unlist(apply(pairList,1,function(x) sum(unlist(x[c("BM1","BM2")])%in%rownames(WeightMat))==2))
    allMarkers<-allMarkers[keep,]
    if(sum(keep)==1){
      pairList<-matrix(pairList[keep,],nrow=1,ncol=4,byrow = T)
      colnames(pairList)<-c("pair","BM1","BM2","DualBM")
    }else{
      pairList<-pairList[keep,]
    }
    print(dim(pairList))
  }else{
  allMarkers<-allMarkers[allMarkers[,"StringID"]%in%rownames(WeightMat),]}
  allMarkers<-allMarkers[allMarkers[,"TstringID"]%in%rownames(WeightMat),]
  #allMarkers<-allMarkers[,c("Depleted.Gene","FEATURE","StringID","TstringID","CLASS","assoc_id","type","AnovaPval","AnovaFdr","AnovaDeltaP","AnovaDeltaN")]
  allMarkers<-allMarkers[,c("Depleted.Gene","FEATURE","StringID","TstringID","assoc_id","type","AnovaPval","AnovaFdr","AnovaDeltaP","AnovaDeltaN")]


  if(single){
    startingvec<-rep(0,nrow(WeightMat))
    names(startingvec)<-rownames(WeightMat)
    startingvec[names(startingvec)%in%BiomarkerString]<-1
  }else{
    #generate matrix, each column is a single startingvec for a single biomarker

      if(Int&IntSplit){
        if(is.matrix(pairList)){


          if(nrow(pairList)==1){
            pairList<-matrix(unique(pairList[,c("pair","BM1","BM2")]),nrow=1,ncol=3,byrow = T)
            colnames(pairList)<-c("pair","BM1","BM2")
          }else{
            pairList<-unique(pairList[,c("pair","BM1","BM2")])
          }
        print(dim(pairList))
        startingvec<-matrix(0,nrow=nrow(WeightMat),ncol=nrow(pairList))
        rownames(startingvec)<-rownames(WeightMat)
        colnames(startingvec)<-pairList[,'pair']
        for(i in 1:nrow(pairList)){
          startingvec[unlist(pairList[i,c("BM1","BM2")]),i]<-0.5
        }
        }else{
          startingvec<-matrix(0,nrow=nrow(WeightMat),ncol=1)
          rownames(startingvec)<-rownames(WeightMat)
          colnames(startingvec)<-pairList['pair']
          startingvec[unlist(pairList[c("BM1","BM2")]),1]<-0.5
        }
      }else{
        startingvec<-matrix(0,nrow=nrow(WeightMat),ncol=length(BiomarkerString))
        rownames(startingvec)<-rownames(WeightMat)
        colnames(startingvec)<-BiomarkerString
        for(i in 1:length(BiomarkerString)){

          startingvec[BiomarkerString[i],i]<-1
        }
      }
  }
  #also generate a weighted startingvec according to the number/strength of target p-values associated to biomarker
  #stringscore<-unlist(sapply(BiomarkerString,function(x) biomarkerPriorWeight(allMarkers,x)))

  weightstart<-rep(1,nrow(WeightMat))
  names(weightstart)<-rownames(WeightMat)
  #for(i in 1:length(stringscore)){
  #  weightstart[names(weightstart)==names(stringscore)[i]]<-stringscore[i]
  #}
  wsname<-weightstart
  names(wsname)<-PPInet[match(names(weightstart),PPInet$STRING_id),"symbol"]

  allMarkers$Depleted.Gene<-unlist(allMarkers$Depleted.Gene)
  allMarkers$assoc_id<-unlist(allMarkers$assoc_id)
  allMarkers$AnovaPval<-unlist(allMarkers$AnovaPval)
  allMarkers$AnovaFdr<-unlist(allMarkers$AnovaFdr)
  allMarkers$AnovaDeltaP<-unlist(allMarkers$AnovaDeltaP)
  allMarkers$AnovaDeltaN<-unlist(allMarkers$AnovaDeltaN)
  allMarkers$FEATURE<-gsub(" ","",allMarkers$FEATURE)
  rownames(allMarkers)<-NULL

  if(HubCorrect){

    mnames<-dimnames(WeightMat)
    degreerow<-rowSums(WeightMat)
    DegMat<-t(sapply(degreerow,function(x) x/degreelist[colnames(WeightMat)]))
    minMat<-matrix(1,nrow=nrow(DegMat),ncol=nrow(DegMat))
    DegMat<-pmin(minMat,DegMat)
    WeightMat<-t(apply(DegMat,2,function(x) x/degreerow))
    dimnames(WeightMat)<-mnames
    degreelist<-colSums(WeightMat)
  }

  return(list(wsname=wsname,weightstart=weightstart,startingvec=startingvec,AllM=allMarkers,WeightMat=WeightMat,degreelist=degreelist))
}
randRWRinput<-function(graphRWR,nbiomarkers=100,nNetworks=10,AnovaFull){
  adjmat<-graphRWR$WeightMat
  allFeat<-unlist(unique(AnovaFull[,"FEATURE"]))


  randNetworks<-list()
  BSmat<-AnovaFull[,c("FEATURE","StringID")]
  BSmat<-BSmat[!is.na(BSmat[,"StringID"]),]
  if(!is.character(BSmat)){
    BSmat[,"FEATURE"]<-as.character(BSmat[,"FEATURE"])
    BSmat[,"StringID"]<-as.character(BSmat[,"StringID"])
  }
  BSmat<-unique(BSmat)

  BSmat<-BSmat[BSmat[,"StringID"]%in%rownames(adjmat),]
  nbiomarkers<-min(nrow(BSmat),nbiomarkers)
  #randomly select nbiomarkers from the ANOVA list
  selbm<-sample(unlist(BSmat[,"FEATURE"]),nbiomarkers,replace=FALSE)
  rownames(BSmat)<-unlist(BSmat[,1])
  startingvecs<-matrix(0,nrow=nrow(adjmat),ncol=nbiomarkers)
  rownames(startingvecs)<-rownames(adjmat)
  BiomarkerString<-unlist(BSmat[selbm,"StringID"])

  colnames(startingvecs)<-selbm
  for(i in 1:length(BiomarkerString)){
    startingvecs[BiomarkerString[i],i]<-1
  }
  #create nNetworks of random networks - where the node labels from the adjacency matrix are perturbed.
  for(i in 1:nNetworks){
    temp<-adjmat
    tnames<-sample(colnames(adjmat),ncol(adjmat))
    rownames(temp)<-tnames
    colnames(temp)<-tnames
    randNetworks[[i]]<-temp
  }
  return(list(networks=randNetworks,biomarkers=selbm,svecs=startingvecs))
}
rwrRand<-function(randRWRinput,rprob,PPInet){
  rnetworks<-randRWRinput$networks
  bms<-randRWRinput$biomarkers
  svecs<-randRWRinput$svecs[,bms]
  nNetworks<-length(rnetworks)
  output<-foreach(i=1:nNetworks)%dopar%{
    #strength has the option of degree values if a weighted network.
    deglist<-strength(graph_from_adjacency_matrix(rnetworks[[i]]))
    incvec<-names(deglist)[deglist!=0]
    svecs<-svecs[incvec,]
    deglist<-deglist[incvec]
    innet<-rnetworks[[i]]
    innet<-innet[incvec,incvec]
    rwrset(svecs,innet,rprob,deglist)

  }
  return(output)
}
GetDeltaThresh<-function(netsizes,RWRres){

  deltavals<-c()
  for(i in 1:length(netsizes)){
    deltavals[i]<-min(unlist(sapply(1:ncol(RWRres),function(x) sort(RWRres[,x],decreasing=TRUE)[netsizes[i]])),na.rm=T)
  }
  return(deltavals)
}
rwrRandstats<-function(rwrRand,AllM,PPInet,hthresh=FALSE,useBMscore=FALSE,netsizes=c(5,10,15)){
  res<-list()
  nnet<-length(rwrRand)
  #for each random network{
    #sapply to delta range
    #for each delta to each anova
  #}
  res<-foreach(i=1:nnet)%dopar%{

    if(hthresh){
      RWRthresh<-sigRWR(rwrRand[[i]])
    }else{
      RWRthresh<-rwrRand[[i]]
    }
    if(useBMscore){
      #ScoreRWRbiomarker(RWRthresh,AllM,PPInet)[[1]]
      temp<-GetDeltaThresh(netsizes,RWRthresh)
    }else{
      #run for a set of deltas find number of signif ANOVA res
      deltarange<-seq(from=0.00001,to=0.001,length.out=100)
      anovaranges<-c(0.2,0.15,0.1,0.05,0.025,0.01,0.005,0.001)
       anovascore<-as.data.frame(AllM,stringsAsFactors=FALSE)
      anovascore$AnovaPval<-as.character(anovascore$AnovaPval)
      anovascore<-anovascore[anovascore[,"AnovaPval"]!="NULL",]
      if(nrow(anovascore)>0){
        genes<-intersect(rownames(RWRthresh),anovascore[,"Depleted.Gene"])
        if(length(genes)==0){
          rownames(RWRthresh)<-PPInet[match(rownames(RWRthresh),PPInet$STRING_id),"symbol"]
       }
        anovascore$AnovaPval<-unlist(anovascore$AnovaPval)
        anovascore$AnovaPval<-as.numeric(anovascore$AnovaPval)
        anovascore$Depleted.Gene<-unlist(anovascore$Depleted.Gene)
        anovascore$StringID<-unlist(anovascore$StringID)
        AnovaPval<-acast(anovascore,Depleted.Gene~StringID,value.var="AnovaPval",fill=1,fun.aggregate = function(x) min(unlist(x)))
        LogPval<-matrix(1,nrow=nrow(RWRthresh),ncol=ncol(RWRthresh))
        dimnames(LogPval)<-dimnames(RWRthresh)
        usegenes<-intersect(rownames(AnovaPval),rownames(LogPval))
        useT<-intersect(colnames(LogPval),colnames(AnovaPval))
        if(length(useT)==0){
          #LogPval has columns as biomarkers not, stringIDs, convert:
          colnames(LogPval)<-anovascore[match(colnames(LogPval),anovascore[,"FEATURE"]),"StringID"]
          colnames(RWRthresh)<-anovascore[match(colnames(RWRthresh),anovascore[,"FEATURE"]),"StringID"]
          useT<-intersect(colnames(LogPval),colnames(AnovaPval))
        }
        if(length(usegenes)==0){
          #LogPval has rows as stringIDs convert to symbols:
          rownames(LogPval)<-PPInet[match(rownames(LogPval),PPInet$STRING_id),"symbol"]
          rownames(RWRthresh)<-PPInet[match(rownames(RWRthresh),PPInet$STRING_id),"symbol"]
          usegenes<-intersect(rownames(LogPval),rownames(AnovaPval))
        }
        if(length(usegenes)==0){
          warning('No genes found in intersection between background RWR data sets. Check inputs.')
        }
        if(length(useT)==0){
          warning('No biomarkers overlapping background ANOVA and RWR. Check inputs.')
        }
        if(length(useT)<ncol(RWRthresh)){
          warning('Not all biomarkers overlapping background ANOVA and RWR. Check inputs.')

        }
        LogPval[usegenes,useT]<- AnovaPval[usegenes,useT]
        if(!all.equal(dimnames(RWRthresh),dimnames(LogPval))){
          warning("RWR and P-value backgrounds not same dim names")
          #  SignifMat<-LogPval*RWRthreshold[rownames(LogPval),colnames(LogPval)]
          print(setdiff(rownames(LogPval),rownames(RWRthresh)))
          print(setdiff(colnames(LogPval),colnames(RWRthresh)))
        }
        checkR<-!is.na(rownames(LogPval))
        LogPval<-LogPval[checkR,]
        checkC<-!is.na(colnames(LogPval))
        LogPval<-LogPval[,checkC]
        temp<-sapply(deltarange,function(x) threshRWRdelta(RWRthresh,x,anovaranges,LogPval))
        rownames(temp)<-paste("Anova p-value",anovaranges)
      }else{
        temp<-NA
      }

    }
    temp
  }
  return(res)
}
threshRWRdelta<-function(RWRthreshold,delta,anovapvals,LogPval){
  RWRthreshold[RWRthreshold<=delta]<-0
  #for each p-value threshold find number signif

  res<-list()


    for(i in 1:length(anovapvals)){
      out<-ScoreRWRrand(RWRthreshold,LogPval,anovapvals[i])
      #do total number of connections summed over all biomarkers tested as well.
      #change to output number signifcant connections as changing delta threshold changes underlying network size as well.
      res[[i]]<-out$numbersignif
    }


  return(res)
}

ScoreRWRrand<-function(RWRthreshold,LogPval,pvalthresh){


  #split so have an anova table for each type of biomarker and then there can't be any overlapping
  #have a type column for anovascore and for the all markers (anovascore)

 numbersignif<-c()
 propsignif<-c()
  outnames<-c()


    #system.time(temp<-GetBMrand(anovascore,RWRthreshold,PPInet,pvalthresh))
  temp<-GetSignifRand(LogPval,RWRthreshold,pvalthresh)
    if(!is.null(temp)){
      temp1<-temp$SignifMat
      temp1<-temp1[,!colnames(temp1)=="NULL"]
      numbersignif<-sum(temp1)/ncol(temp$RWRthreshold)

      temp2<-temp$RWRthreshold
     propsignif<-numbersignif/sum(temp2!=0)

     }


  return(list(numbersignif=numbersignif,propsignif=propsignif))
}
GetSignifRand<-function(LogPval,RWRthreshold,pvalthresh){
  LogPval[LogPval>=pvalthresh]<-2
  LogPval[LogPval!=2]<-1
  LogPval[LogPval==2]<-0
  SignifMat<-LogPval*RWRthreshold[rownames(LogPval),colnames(LogPval)]
  SignifMat<-SignifMat!=0+0
  return(list(SignifMat=SignifMat,RWRthreshold=RWRthreshold))
}
GetBMrand<-function(anovascore,RWRthreshold,PPInet,pvalthresh){
  anovascore<-as.data.frame(anovascore,stringsAsFactors=FALSE)
  anovascore$AnovaPval<-as.character(anovascore$AnovaPval)
  anovascore<-anovascore[anovascore[,"AnovaPval"]!="NULL",]

  if(nrow(anovascore)>0){
    anovascore$AnovaPval<-unlist(anovascore$AnovaPval)
    anovascore$AnovaPval<-as.numeric(anovascore$AnovaPval)
    anovascore$Depleted.Gene<-unlist(anovascore$Depleted.Gene)
    anovascore$StringID<-unlist(anovascore$StringID)
    AnovaPval<-acast(anovascore,Depleted.Gene~StringID,value.var="AnovaPval",fill=1,fun.aggregate = function(x) min(unlist(x)))
    LogPval<-matrix(0,nrow=nrow(RWRthreshold),ncol=ncol(RWRthreshold))
    SigScores<-matrix(0,nrow=nrow(RWRthreshold),ncol=ncol(RWRthreshold))
    dimnames(LogPval)<-dimnames(RWRthreshold)
    dimnames(SigScores)<-dimnames(RWRthreshold)
    usegenes<-intersect(rownames(AnovaPval),rownames(LogPval))
    useT<-intersect(colnames(LogPval),colnames(AnovaPval))
    LogPval[usegenes,useT]<- AnovaPval[usegenes,useT]
    LogPval<-LogPval[usegenes,useT]
    LogPval[LogPval>=pvalthresh]<-2
    LogPval[LogPval!=2]<-1
    LogPval[LogPval==2]<-0
    SignifMat<-LogPval*RWRthreshold[usegenes,useT]
    SignifMat<-SignifMat!=0+0

   return(list(SignifMat=SignifMat,RWRthreshold=RWRthreshold))
   }else{
      return(NULL)
    }
}
getDelta<-function(rwrOut,pval,boundaries=c(1,2,3),dranges=seq(from=0.000001,to=0.0001,length.out=100)){
  rname<-paste("Anova p-value",pval)
  deltas<-sapply(boundaries,function(x) dranges[which(unlist(rwrOut[rname,])<x)[1]])
  return(deltas)
}
biomarkerPriorWeight<-function(allMarkers,markers){
  allscore<-0
  scoretable<-c(1,2,3,4)
  names(scoretable)<-c("D","C","B","A")
  subM<-allMarkers[allMarkers$StringID==markers,]
  subM<-subM[subM$CLASS!="D",]
  tempt<-table(subM$CLASS)
  classscore<-as.vector(tempt)
  names(classscore)<-names(tempt)
  for(i in 1:length(classscore)){
    allscore<-allscore+scoretable[names(classscore)[i]]*classscore[i]
  }
  names(allscore)<-NULL
  names(classscore)<-NULL
  return(allscore)
}
GSEAfunction<-function(rankedlist,geneset){
  guse<-intersect(rankedlist,geneset)
  Inset<-(rankedlist%in%guse)*1
  sizeset<-length(guse)
  N_r<-sizeset
  PosScores<-Inset/N_r
  NegScores<-(!rankedlist%in%guse)*1
  N_Nh<-length(rankedlist)-sizeset
  NegScores<-NegScores/N_Nh
  Allscores<-PosScores-NegScores
  RunningSum<-cumsum(Allscores)
  maxdev<-RunningSum[which.max(abs(RunningSum))]
  return(list(ESscore=maxdev,RunningSum=RunningSum))

}
compareRWRanova<-function(rankRWRlist,anovaset,plotname,dir.Results,randomNull=100){
  output<-GSEAfunction(rankRWRlist,anovaset)
  #random sets:
  randomES<-sapply(1:randomNull,function(x) GSEAfunction(sample(rankRWRlist,size=length(rankRWRlist)),anovaset)$ESscore)
  pvalNull<-sum(unlist(randomES)>=output$ESscore)/randomNull
  if(pvalNull==0){
    pvalNull<-"<0.01"
  }
  pdf(paste0(dir.Results,"/",plotname,".pdf"))
  plot(output$RunningSum,ylim=c(min(output$RunningSum)-0.05,max(output$RunningSum)+0.05),ylab="Running Sum")
  text(5000,max(output$RunningSum)+0.03,paste0("Enrichment Score ",round(output$ESscore,2)," p-value ",pvalNull))

  dev.off()
}
sigRWR<-function(RWRresult,delta=NULL){
  if(is.matrix(RWRresult)){
    res<-apply(RWRresult,2,function(x) threshRWR(x,delta))
  }else{
    res<-threshRWR(RWRresult)
  }
  return(res)
}
threshRWR<-function(RWRresult,delta=NULL){


  if(is.null(delta)){
    med<-median(RWRresult)
    madr<-mad(RWRresult)
    thresh<-med+2*madr
    RWRresult[RWRresult<=thresh]<-0
  }else{
    RWRresult[RWRresult<=delta]<-0
  }

  return(RWRresult)
}
minmax<-function(x){
  minv<-min(x)
  maxv<-max(x)
  mm<-(x-minv)/(maxv-minv)
  return(mm)
}
getIndividScores<-function(RankWeight){
  Targets<-sapply(1:ncol(RankWeight),function(x) rownames(RankWeight)[RankWeight[,x]!=0])
  Targets<-unlist(lapply(Targets,function(x) paste(x,collapse="//")))
  Scores<-sapply(1:ncol(RankWeight),function(x) RankWeight[RankWeight[,x]!=0,x])
  Scores<-unlist(lapply(Scores,function(x) paste(x,collapse="//")))
  return(list(Scores=Scores,Targets=Targets))
}
processBMPPICombine<-function(ctype,inputRWR,methodName=c("_Compound","_CompoundME","_CCN","_CExpr","_Genomic","_Int","_Binary","_CProt"),outputName,type=c("n1","Avg","NP","NPAvg"),
                              classMat,biomarkerMat=NULL,TPset=NULL,RWRthreshMat,PvalthreshMat=NULL,FDRscoreMat=NULL,PPInet=NULL,
                              pair=FALSE,DualBM=FALSE,
                              RegressScoreMat=NULL,RegressEffectMat=NULL,IntSplit=FALSE){
  AllBM<-NULL
  BMPriority<-NULL
  BMScores<-NULL
  NumberBiomarkers<-NULL
  allTSdf<-NULL
  aovoutAll<-NULL
  j=1
  type<-match.arg(type)
  if(pair){
    inputRWR<-paste0(inputRWR,"/PairBD/")
  }
  if(DualBM){
    inputRWR<-paste0(inputRWR,"/Int/")
  }
  for(i in 1:length(ctype)){
    inputAll<-NULL
    rwrThreshAll<-NULL
    print(ctype[i])
    for(k in 1:length(methodName)){

        if(file.exists(paste0(inputRWR,ctype[i],"/RWRinputAll",methodName[k],".Rdata"))&file.exists(paste0(inputRWR,ctype[i],"/RWRthresh",methodName[k],".Rdata"))){
          pairList<-NULL
          if(pair){
            load(file=paste0(inputRWR,"/",ctype[i],"/PairList",methodName[k],".Rdata"))
          }
          if(DualBM){

            load(file=paste0(inputRWR,"/",ctype[i],"/PairList_Compound.Rdata"))
            pairList<-Allpair
          }

          load(file=paste0(inputRWR,"/",ctype[i],"/RWRinputAll",methodName[k],".Rdata"))
          load(file=paste0(inputRWR,"/",ctype[i],"/RWRthresh",methodName[k],".Rdata"))
          inputAll<-rbind(inputAll,RWRinputAll$AllM)
          rwrThreshAll<-cbind(rwrThreshAll,RWRthresh[,setdiff(colnames(RWRthresh),colnames(rwrThreshAll))])


      }
    }

    if(!is.null(rwrThreshAll)&sum(colnames(rwrThreshAll)!="")>0){

      BMscoresWeight<-NULL
      BMscoresWeightAvg<-NULL
      BMindTargets<-NULL
      BMindScores<-NULL
      TargetScoresWeight<-NULL
      TargetScoresWeightAvg<-NULL

      #BMS1<-ScoreRWRbiomarkerCombined(rwrThreshAll,inputAll,PPInet,pval=0.1,RWRthreshScore=RWRthreshMat,PvalthreshScore=PvalthreshMat,FDRscoreMat=FDRscoreMat,pairList=pairList,pair=pair,DualBM=DualBM,RegressScoreMat=RegressScoreMat,RegressEffectMat=RegressEffectMat)
      BMS2<-ScoreRWRbiomarkerCombined(rwrThreshAll,inputAll,PPInet,pval=0.1,RWRthreshScore=RWRthreshMat,PvalthreshScore=PvalthreshMat,
                                      FDRscoreMat=FDRscoreMat,scoring="avg",pairList=pairList,pair=pair,DualBM=DualBM,
                                      RegressScoreMat=RegressScoreMat,RegressEffectMat=RegressEffectMat,IntSplit=IntSplit,ctype=ctype[i])
      BMS1<-BMS2
      BMrw<-BMS1$RankWeight
      BMPriority<-AnnotateBiomarkers(BMS1)
      BiomarkerScores<-BMS1
      BMscores<-BMS1$BiomarkerScore
      aovout<-BMS1$anovascore

      if(!is.null(BMPriority)){
        aovout$ctype<-ctype[i]
        if(length(BMrw)==1){
          BMrw<-BMrw[[1]]

          BMSindivid<-getIndividScores(BMrw)
          BMindTargets<-BMSindivid$Targets
          BMindScores<-BMSindivid$Scores
          names(BMindTargets)<-colnames(BMrw)
          names(BMindScores)<-colnames(BMrw)
          TargetScoresWeight<-BMS1$TargetScore[[1]]
          TargetScoresWeightAvg<-BMS2$TargetScore[[1]]

          BMscoresWeight<-BMS1$BiomarkerScore
          BMscoresWeightAvg<-BMS2$BiomarkerScore
          BMscoresWeight<-BMscoresWeight[names(BMPriority)]
          BMscoresWeightAvg<-BMscoresWeightAvg[names(BMPriority)]
          BMindTargets<-BMindTargets[names(BMPriority)]
          BMindScores<-BMindScores[names(BMPriority)]

        }
      }
      AnovaPval<-BiomarkerScores$AnovaPval
      AnovaFdr<-BiomarkerScores$AnovaFdr
      AnovaDeltaP<-BiomarkerScores$AnovaDeltaP
      AnovaDeltaN<-BiomarkerScores$AnovaDeltaN
      if(length(BMindScores)==0){
        BMindScores<-rep(NA,length(BMPriority))
        BMindTargets<-rep(NA,length(BMPriority))
        BMscoresWeight<-rep(NA,length(BMPriority))
        BMscoresWeightAvg<-rep(NA,length(BMPriority))
      }

      #PriorityBiomarkers<-BMscores[round(BMscores,2)>0]
      #PriorityBiomarkers<-PriorityBiomarkers[!is.na(names(PriorityBiomarkers))]
      if(!is.null(BMscores)){
        PriorityBiomarkers<-BMscores[!is.na(names(BMscores))]
        PriorityBiomarkers<-PriorityBiomarkers[names(BMPriority)]
        aovoutAll<-rbind(aovoutAll,aovout)
        if(length(PriorityBiomarkers)>0){
          BMScores<-data.frame(BM=names(PriorityBiomarkers),BMScore=PriorityBiomarkers,BMScoreClass=BMscoresWeight,BMScoreClassAvg=BMscoresWeightAvg,BMindTargets=BMindTargets,BMindScores=BMindScores,stringsAsFactors = FALSE)
          BMScores$CTYPE<-ctype[i]
          BMtargets<-lapply(BMPriority,function(x) paste0(names(sort(x,decreasing=TRUE)),collapse="//"))
          #BMclass<-sapply(names(BMPriority),function(x) getBMclass(x,BMPriority,AnovaFdr,AnovaDeltaP,AnovaDeltaN,AnovaPval,classMat))
          BMScores$ScoringTargets<-unlist(BMtargets[BMScores$BM])
          #BMScores$ScoringClass<-unlist(BMclass[BMScores$BM])
          BMScores$NumberTargets<-unlist(sapply(BMScores$ScoringTargets,function(x) length(unlist(strsplit(x,"//",fixed=TRUE)))))
          if(!dir.exists(paste0(dir.Results,ctype[i]))){dir.create(paste0(dir.Results,ctype[i]))}

          AllBM<-rbind(AllBM,BMScores)
        }
        SigTS<-TargetScoresWeight[TargetScoresWeight$Score!=0,]
        SigTSavg<-TargetScoresWeightAvg[TargetScoresWeightAvg$Score!=0,]
        SigTS<-NULL
        tdf1<-NULL
        tdf2<-NULL
        if(!is.null(SigTS)){
          tdf1<-data.frame(Target=rownames(SigTS),MaxTargetScoreFDR=SigTS$MaxTargetScoreFDR,Score=SigTS$Score,FDRScore=SigTS$FDRScore,PvalScore=SigTS$PvalScore,RWRscore=SigTS$RWRscore,BMType=SigTS$BMtype,BMSymbol=SigTS$BMsymbol,BMclass=SigTS$BMclass,CancerType=ctype[i],ScoreType="Sum",BMs=SigTS$BMs,NumberBMs=SigTS$NumberBMs,RWRclasses=SigTS$RWRclasses,PvalClasses=SigTS$PvalClasses,FDRClasses=SigTS$FDRClasses,stringsAsFactors = FALSE)
        }
        if(!is.null(SigTSavg)){
          tdf2<-data.frame(Target=rownames(SigTSavg),MaxTargetScoreFDR=SigTSavg$MaxTargetScoreFDR,Score=SigTSavg$Score,FDRScore=SigTSavg$FDRScore,PvalScore=SigTSavg$PvalScore,RWRscore=SigTSavg$RWRscore,BMType=SigTSavg$BMtype,BMSymbol=SigTSavg$BMsymbol,BMclass=SigTSavg$BMclass,CancerType=ctype[i],ScoreType="Avg",BMs=SigTSavg$BMs,NumberBMs=SigTSavg$NumberBMs,RWRclasses=SigTSavg$RWRclasses,PvalClasses=SigTSavg$PvalClasses,FDRClasses=SigTSavg$FDRClasses,stringsAsFactors = FALSE)
        }
        allTSdf<-rbind(allTSdf,tdf1,tdf2)
        j=j+1
      }

    }

  }
  if(!is.null(AllBM)){
    AllBM<-AllBM[order(AllBM$BMScore,decreasing=TRUE),]
    BMclass<-rep("none",nrow(AllBM))
    if(!is.null(biomarkerMat)){
      biomarkerMat<-biomarkerMat[order(biomarkerMat$class,decreasing=TRUE),]
      for(i in 1:nrow(biomarkerMat)){
        check<-which(AllBM$BMScore>biomarkerMat[i,"Score"])
        BMclass[check]<-biomarkerMat[i,"class"]
      }
      AllBM$BMClass<-BMclass
    }

    #determine significant results here:
    BMscores<-AllBM$BMScore
    if(is.null(TPset)){
      PriorityBiomarkers<-which(round(BMscores,2)>0)
    }else{
      #have a set of biomarkers to use to define a tp distribution of scores:
      BMscores<-BMscores[BMscores!=0]
      BMgenes<-names(BMscores)
      BMgenes<-sapply(BMgenes,getGenesFromBiomarkers)
      tpset<-which(BMgenes%in%TPset)
      if(length(tpset)>5){
        dthresh<-densityThresholding(set1=BMscores[-tpset],set2=BMscores[tpset],x=seq(min(BMscores),max(BMscores),length.out=100))
        PriorityBiomarkers<-which(BMscores>dthresh)}else{
          PriorityBiomarkers<-which(BMscores>median(BMscores,na.rm=T))
        }
    }
    AllBM<-AllBM[PriorityBiomarkers,]
    NumberBiomarkers<-as.vector(table(AllBM$CTYPE))
    names(NumberBiomarkers)<-names(table(AllBM$CTYPE))
    write.table(AllBM,file=paste0(dir.Results,"/",outputName,".tsv"),row.names=FALSE,quote=F,sep="\t")
    if(!is.null(aovoutAll)){
      aovoutAll<-as.data.frame(aovoutAll)
      for(i in 1:ncol(aovoutAll)){
        aovoutAll[,i]<-unlist(aovoutAll[,i])
      }
      write.table(aovoutAll,file=paste0(dir.Results,"/",outputName,"_aovout.tsv"),row.names=FALSE,quote=F,sep="\t")
    }
    return(list(AllBM,NumberBiomarkers,allTSdf))
  }else{
    return(NULL)
  }
}

processBMPPI<-function(ctype,inputRWR,methodName,outputName,type=c("n1","Avg","NP","NPAvg"),classMat,biomarkerMat=NULL,TPset=NULL,RWRthreshMat,PvalthreshMat=NULL,FDRscoreMat=NULL,PPInet=NULL,pair=FALSE,DualBM=FALSE){
  AllBM<-NULL
  BMPriority<-NULL
  BMScores<-NULL
  NumberBiomarkers<-NULL
  allTSdf<-NULL
  j=1
  type<-match.arg(type)
  if(pair){
    inputRWR<-paste0(inputRWR,"/PairBD/")}
  if(DualBM){inputRWR<-paste0(inputRWR,"/Int/")}
  for(i in 1:length(ctype)){


    if(file.exists(paste0(inputRWR,ctype[i],"/RWRinputAll",methodName,".Rdata"))&file.exists(paste0(inputRWR,ctype[i],"/RWRthresh",methodName,".Rdata"))){
      pairList<-NULL
      if(pair){
      load(file=paste0(inputRWR,"/",ctype[i],"/PairList",methodName,".Rdata"))
      }
      if(DualBM){

        load(file=paste0(inputRWR,"/",ctype[i],"/PairList_Compound.Rdata"))
        pairList<-Allpair
      }


      BMscoresWeight<-NULL
      BMscoresWeightAvg<-NULL
      BMindTargets<-NULL
      BMindScores<-NULL
      TargetScoresWeight<-NULL
      TargetScoresWeightAvg<-NULL

        load(file=paste0(inputRWR,"/",ctype[i],"/RWRinputAll",methodName,".Rdata"))
        load(file=paste0(inputRWR,"/",ctype[i],"/RWRthresh",methodName,".Rdata"))

        BMS1<-ScoreRWRbiomarker(RWRthresh,RWRinputAll$AllM,PPInet,pval=0.1,RWRthreshScore=RWRthreshMat,PvalthreshScore=PvalthreshMat,FDRscoreMat=FDRscoreMat,pairList=pairList,pair=pair,DualBM=DualBM)
        BMS2<-ScoreRWRbiomarker(RWRthresh,RWRinputAll$AllM,PPInet,pval=0.1,RWRthreshScore=RWRthreshMat,PvalthreshScore=PvalthreshMat,FDRscoreMat=FDRscoreMat,scoring="avg",pairList=pairList,pair=pair,DualBM=DualBM)

        BMrw<-BMS1$RankWeight
        BMPriority<-AnnotateBiomarkers(BMS1)
        BiomarkerScores<-BMS1
        BMscores<-BMS1$BiomarkerScore
        if(!is.null(BMPriority)){
          if(length(BMrw)==1){
            BMrw<-BMrw[[1]]

            BMSindivid<-getIndividScores(BMrw)
            BMindTargets<-BMSindivid$Targets
            BMindScores<-BMSindivid$Scores
            names(BMindTargets)<-colnames(BMrw)
            names(BMindScores)<-colnames(BMrw)
            TargetScoresWeight<-BMS1$TargetScore[[1]]
            TargetScoresWeightAvg<-BMS2$TargetScore[[1]]

            BMscoresWeight<-BMS1$BiomarkerScore
            BMscoresWeightAvg<-BMS2$BiomarkerScore
            BMscoresWeight<-BMscoresWeight[names(BMPriority)]
            BMscoresWeightAvg<-BMscoresWeightAvg[names(BMPriority)]
            BMindTargets<-BMindTargets[names(BMPriority)]
            BMindScores<-BMindScores[names(BMPriority)]

          }
        }
        AnovaPval<-BiomarkerScores$AnovaPval
        AnovaFdr<-BiomarkerScores$AnovaFdr
        AnovaDeltaP<-BiomarkerScores$AnovaDeltaP
        AnovaDeltaN<-BiomarkerScores$AnovaDeltaN
        if(length(BMindScores)==0){
          BMindScores<-rep(NA,length(BMPriority))
          BMindTargets<-rep(NA,length(BMPriority))
          BMscoresWeight<-rep(NA,length(BMPriority))
          BMscoresWeightAvg<-rep(NA,length(BMPriority))
        }

        #PriorityBiomarkers<-BMscores[round(BMscores,2)>0]
        #PriorityBiomarkers<-PriorityBiomarkers[!is.na(names(PriorityBiomarkers))]
        if(!is.null(BMscores)){
          PriorityBiomarkers<-BMscores[!is.na(names(BMscores))]
          PriorityBiomarkers<-PriorityBiomarkers[names(BMPriority)]

          if(length(PriorityBiomarkers)>0){
            BMScores<-data.frame(BM=names(PriorityBiomarkers),BMScore=PriorityBiomarkers,BMScoreClass=BMscoresWeight,BMScoreClassAvg=BMscoresWeightAvg,BMindTargets=BMindTargets,BMindScores=BMindScores,stringsAsFactors = FALSE)
            BMScores$CTYPE<-ctype[i]
            BMtargets<-lapply(BMPriority,function(x) paste0(names(sort(x,decreasing=TRUE)),collapse="//"))
            #BMclass<-sapply(names(BMPriority),function(x) getBMclass(x,BMPriority,AnovaFdr,AnovaDeltaP,AnovaDeltaN,AnovaPval,classMat))
            BMScores$ScoringTargets<-unlist(BMtargets[BMScores$BM])
            #BMScores$ScoringClass<-unlist(BMclass[BMScores$BM])
            BMScores$NumberTargets<-unlist(sapply(BMScores$ScoringTargets,function(x) length(unlist(strsplit(x,"//",fixed=TRUE)))))
            if(!dir.exists(paste0(dir.Results,ctype[i]))){dir.create(paste0(dir.Results,ctype[i]))}

            AllBM<-rbind(AllBM,BMScores)
          }
          SigTS<-TargetScoresWeight[TargetScoresWeight$Score!=0,]
          SigTSavg<-TargetScoresWeightAvg[TargetScoresWeightAvg$Score!=0,]

          if(!is.null(SigTS)){
            tdf1<-data.frame(Target=rownames(SigTS),MaxTargetScoreFDR=SigTS$MaxTargetScoreFDR,Score=SigTS$Score,FDRScore=SigTS$FDRScore,PvalScore=SigTS$PvalScore,RWRscore=SigTS$RWRscore,BMType=SigTS$BMtype,CancerType=ctype[i],ScoreType="Sum",BMs=SigTS$BMs,NumberBMs=SigTS$NumberBMs,RWRclasses=SigTS$RWRclasses,PvalClasses=SigTS$PvalClasses,FDRClasses=SigTS$FDRClasses,stringsAsFactors = FALSE)
          }
          if(!is.null(SigTSavg)){
            tdf2<-data.frame(Target=rownames(SigTSavg),MaxTargetScoreFDR=SigTSavg$MaxTargetScoreFDR,Score=SigTSavg$Score,FDRScore=SigTSavg$FDRScore,PvalScore=SigTSavg$PvalScore,RWRscore=SigTSavg$RWRscore,BMType=SigTSavg$BMtype,CancerType=ctype[i],ScoreType="Avg",BMs=SigTSavg$BMs,NumberBMs=SigTSavg$NumberBMs,RWRclasses=SigTSavg$RWRclasses,PvalClasses=SigTSavg$PvalClasses,FDRClasses=SigTSavg$FDRClasses,stringsAsFactors = FALSE)
          }
          allTSdf<-rbind(allTSdf,tdf1,tdf2)
          j=j+1
        }

      }

  }
  if(!is.null(AllBM)){
    AllBM<-AllBM[order(AllBM$BMScore,decreasing=TRUE),]
    BMclass<-rep("none",nrow(AllBM))
    if(!is.null(biomarkerMat)){
      biomarkerMat<-biomarkerMat[order(biomarkerMat$class,decreasing=TRUE),]
      for(i in 1:nrow(biomarkerMat)){
        check<-which(AllBM$BMScore>biomarkerMat[i,"Score"])
        BMclass[check]<-biomarkerMat[i,"class"]
      }
      AllBM$BMClass<-BMclass
    }

    #determine significant results here:
    BMscores<-AllBM$BMScore
    if(is.null(TPset)){
      PriorityBiomarkers<-which(round(BMscores,2)>0)
    }else{
      #have a set of biomarkers to use to define a tp distribution of scores:
      BMscores<-BMscores[BMscores!=0]
      BMgenes<-names(BMscores)
      BMgenes<-sapply(BMgenes,getGenesFromBiomarkers)
      tpset<-which(BMgenes%in%TPset)
      if(length(tpset)>5){
        dthresh<-densityThresholding(set1=BMscores[-tpset],set2=BMscores[tpset],x=seq(min(BMscores),max(BMscores),length.out=100))
        PriorityBiomarkers<-which(BMscores>dthresh)}else{
        PriorityBiomarkers<-which(BMscores>median(BMscores,na.rm=T))
      }
    }
    AllBM<-AllBM[PriorityBiomarkers,]
    NumberBiomarkers<-as.vector(table(AllBM$CTYPE))
    names(NumberBiomarkers)<-names(table(AllBM$CTYPE))
    write.table(AllBM,file=paste0(dir.Results,"/",outputName,".tsv"),row.names=FALSE,quote=F,sep="\t")
    return(list(AllBM,NumberBiomarkers,allTSdf))
  }else{
    return(NULL)
  }
}
getBMclass<-function(BM_target,BMPriority,AnovaFdr,AnovaDeltaP,AnovaDeltaN,AnovaPval,classMat){
  BM<-BM_target
  BMtype<-strsplit(BM,"_",fixed=TRUE)[[1]]
  BMtype<-BMtype[length(BMtype)]
  if(BMtype=="Expr"){
    A1<-AnovaFdr$Expr[,BM]

    A3<-AnovaDeltaN$Expr[,BM]
    A4<-AnovaPval$Expr[,BM]
    inputdata<-cbind(A4,A1,A3)
    classMat<-classMat[order(classMat$class,decreasing=TRUE),c("Pval","Fdr","Range","class")]
  }
  if(BMtype=="CN"){
    A1<-AnovaFdr$CN[,BM]

    A3<-AnovaDeltaN$CN[,BM]
    A4<-AnovaPval$CN[,BM]
    inputdata<-cbind(A4,A1,A3)
    classMat<-classMat[order(classMat$class,decreasing=TRUE),c("Pval","Fdr","Range","class")]
  }
  if(BMtype=="mut"){
    A1<-AnovaFdr$mut[,BM]
    A2<-AnovaDeltaP$mut[,BM]
    A3<-AnovaDeltaN$mut[,BM]
    A4<-AnovaPval$mut[,BM]
    inputdata<-cbind(A4,A1,A2,A3)
    classMat<-classMat[order(classMat$class,decreasing=TRUE),c("Pval","Fdr","DeltaP","DeltaN","class")]
  }
  if(BMtype=="Prot"){
    A1<-AnovaFdr$Prot[,BM]

    A3<-AnovaDeltaN$Prot[,BM]
    A4<-AnovaPval$Prot[,BM]
    inputdata<-cbind(A4,A1,A3)
    classMat<-classMat[order(classMat$class,decreasing=TRUE),c("Pval","Fdr","Range","class")]
  }
  if(!BMtype%in%c("Expr","mut","CN","Prot")){
    A1<-AnovaFdr$CNA[,BM]
    A2<-AnovaDeltaP$CNA[,BM]
    A3<-AnovaDeltaN$CNA[,BM]
    A4<-AnovaPval$CNA[,BM]
    inputdata<-cbind(A4,A1,A2,A3)
    classMat<-classMat[order(classMat$class,decreasing=TRUE),c("Pval","Fdr","DeltaP","DeltaN","class")]
  }


  if(nrow(inputdata)>0){
    targetclass<-rep('none',nrow(inputdata))}else{
      targetclass<-"none"
    }
  if(BMtype=="mut"){
    for(i in 1:nrow(classMat)){

      if(nrow(inputdata)>1){

        check<-which(inputdata[,1]<classMat[i,1]&inputdata[,2]<classMat[i,2]&inputdata[,3]>classMat[i,3]&inputdata[,4]>classMat[i,4])}else{
          check<-which(inputdata[1]<classMat[i,1]&inputdata[2]<classMat[i,2]&inputdata[3]>classMat[i,3]&inputdata[4]>classMat[i,4])
        }
      targetclass[check]<-classMat[i,"class"]
    }
    names(targetclass)<-rownames(inputdata)
    output<-paste0(targetclass[names(sort(BMPriority[BM_target][[1]],decreasing=TRUE))],collapse="//")
  }
  if(BMtype%in%c("Expr","CN","Prot")){
    for(i in 1:nrow(classMat)){

      if(nrow(inputdata)>1){

        check<-which(inputdata[,1]<classMat[i,1]&inputdata[,2]<classMat[i,2]&inputdata[,3]>classMat[i,3])}else{
          check<-which(inputdata[1]<classMat[i,1]&inputdata[2]<classMat[i,2]&inputdata[3]>classMat[i,3])
        }
      targetclass[check]<-classMat[i,"class"]
    }
    names(targetclass)<-rownames(inputdata)
    output<-paste0(targetclass[names(sort(BMPriority[BM_target][[1]],decreasing=TRUE))],collapse="//")
  }
  return(output)
}
ScoreRWRbiomarker<-function(RWRthreshold,anovascore,PPInet,scoring=c("sum","avg"),normPval=c("none","minmax"),pval=1,RWRthreshScore=NULL,PvalthreshScore=NULL,FDRscoreMat=NULL,pairList=NULL,pair=FALSE,DualBM=FALSE){
  scoring<-match.arg(scoring)
  normPval<-match.arg(normPval)
  genes<-intersect(rownames(RWRthreshold),anovascore[,"Depleted.Gene"])
  if(length(genes)==0){
    rownames(RWRthreshold)<-PPInet[match(rownames(RWRthreshold),PPInet$STRING_id),"symbol"]
    genes<-intersect(rownames(RWRthreshold),anovascore[,"Depleted.Gene"])
   # RWRthreshold<-RWRthreshold[!is.na(rownames(RWRthreshold)),]
  }
  #split so have an anova table for each type of biomarker and then there can't be any overlapping
  #have a type column for anovascore and for the all markers (anovascore)
  utypes<-unique(unlist(anovascore[,'type']))
  utypes<-utypes[utypes%in%c("Expr","mut","CNA","CN","Met","Prot")]
  BiomarkerScore<-NULL
  RWRscorevec<-NULL
  SigScores<-NULL
  AnovaPval<-list()
  AnovaFdr<-list()
  AnovaDeltaP<-list()
  AnovaDeltaN<-list()
  RankWeight<-list()
  RWRprob<-list()
  TargetScore<-list()
  outnames<-c()
  for(i in 1:length(utypes)){
    subanova<-anovascore[anovascore[,"type"]==utypes[i],]
    if(is.matrix(RWRthreshold)){
    temp<-GetBMscores(subanova,RWRthreshold,PPInet,normPval,scoring,pval=pval,RWRthreshScore=RWRthreshScore,PvalthreshScore=PvalthreshScore,FDRthreshScore=FDRscoreMat,pairList=pairList,pair=pair,DualBM=DualBM)
    if(sum(temp$BiomarkerScore,na.rm=T)!=0){
      BiomarkerScore<-c(BiomarkerScore,temp$BiomarkerScore)
      SigScores<-cbind(SigScores,temp$SigScores)
      tR<-temp$RWRscorevec
      tR$bmtype<-utypes[i]
      RWRscorevec<-rbind(RWRscorevec,tR)

      AnovaPval[[i]]<-temp$AnovaPval
      AnovaFdr[[i]]<-temp$AnovaFdr
      AnovaDeltaP[[i]]<-temp$AnovaDeltaP
      AnovaDeltaN[[i]]<-temp$AnovaDeltaN
      RWRprob[[i]]<-temp$RWRprob
      RankWeight[[i]]<-temp$RankWeight
      tset<-names(temp$TargetScore)
      PvalscoreClass<-temp$PvalscoreClass
      RWRscoreClass<-temp$RWRscoreClass
      TargetBMs<-temp$TargetBMs
      NumberBMs<-temp$NumberBMs
      MaxRWR<-temp$MaxRWR
      RWRclasses<-temp$RWRclasses
      PvalClasses<-temp$Pvalclasses
      FDRClasses<-temp$FDRclasses
      Tout<-data.frame(MaxTargetScoreFDR=temp$MaxTargetScoreFDR,Score=temp$TargetScore,FDRScore=temp$TargetScoreFDR,BMtype=utypes[i],BMsymbol=temp$BMSymbol,PvalScore=PvalscoreClass[tset],RWRscore=MaxRWR,BMs=TargetBMs[tset],NumberBMs=NumberBMs[tset],RWRclasses=RWRclasses[tset],PvalClasses=PvalClasses[tset],FDRClasses=FDRClasses[tset],stringsAsFactors = FALSE)
      rownames(Tout)<-tset
      TargetScore[[i]]<-Tout
      outnames<-c(outnames,utypes[i])

    }
    }
  }
  names(AnovaPval)<-outnames
  names(AnovaFdr)<-outnames
  names(AnovaDeltaP)<-outnames
  names(AnovaDeltaN)<-outnames
  names(RWRprob)<-outnames
  return(list(RWRscorevec=RWRscorevec,BiomarkerScore=BiomarkerScore,SigScores=SigScores,RWRprob=RWRprob,AnovaPval=AnovaPval,AnovaFdr=AnovaFdr,AnovaDeltaN=AnovaDeltaN,AnovaDeltaP=AnovaDeltaP,RankWeight=RankWeight,TargetScore=TargetScore))
}
getMaxAovMatrix<-function(anovascore){
  idset<-paste(anovascore[,"Depleted.Gene"],anovascore[,"StringID"],sep="_")
  anovascore[,"idset"]<-idset
  freqv<-table(idset)
  idvals<-names(freqv)
  singleRes<-names(which(freqv==1))
  if(length(singleRes)>0){
    singleaov<-anovascore[idset%in%singleRes,]}else{
    singleaov<-NULL
  }
  multiRes<-names(which(freqv!=1))
  if(length(multiRes)>0){
    selAovs<-sapply(multiRes,function(x) getMaxBM(anovascore[idset==x,]))
    mAovs<-t(selAovs)
    colnames(mAovs)<-colnames(singleaov)
    allanova<-rbind(singleaov,mAovs)
  }else{allanova<-singleaov}
  return(allanova)
}
getMaxBM<-function(inputres,selMethod=c("FDR","Pvalue")){
  scoring<-match.arg(selMethod)
  inputres[,4]<-unlist(inputres[,4])
  inputres[,3]<-unlist(inputres[,3])
  if(scoring=="FDR"){
    sel<-which.min(inputres[,"AnovaFdr"])
    output<-as.matrix(inputres[sel,],nrow=1)
  }
  if(scoring=="Pvalue"){
    sel<-which.min(inputres[,"AnovaPval"])
    output<-as.matrix(inputres[sel,],nrow=1)
  }
  return(output)
}
ScoreRWRbiomarkerCombined<-function(RWRthreshold,anovascore,PPInet,scoring=c("sum","avg"),normPval=c("none","minmax"),pval=1,
                                    RWRthreshScore=NULL,PvalthreshScore=NULL,FDRscoreMat=NULL,pairList=NULL,pair=FALSE,DualBM=FALSE,
                                    RegressScoreMat=NULL,RegressEffectMat=NULL,IntSplit=FALSE,ctype){
  scoring<-match.arg(scoring)
  normPval<-match.arg(normPval)
  genes<-intersect(rownames(RWRthreshold),anovascore[,"Depleted.Gene"])
  if(length(genes)==0){
    rownames(RWRthreshold)<-PPInet[match(rownames(RWRthreshold),PPInet$STRING_id),"symbol"]
    genes<-intersect(rownames(RWRthreshold),anovascore[,"Depleted.Gene"])
    # RWRthreshold<-RWRthreshold[!is.na(rownames(RWRthreshold)),]
  }
  #split so have an anova table for each type of biomarker and then there can't be any overlapping
  #have a type column for anovascore and for the all markers (anovascore)
  utypes<-unique(unlist(anovascore[,'type']))
  utypes<-utypes[utypes%in%c("Expr","mut","CNA","CN","Met","Prot","var","Epi")]
  BiomarkerScore<-NULL
  SigScores<-NULL
  RWRscorevec<-NULL
  aovout<-NULL
  AnovaPval<-list()
  AnovaFdr<-list()
  AnovaDeltaP<-list()
  AnovaDeltaN<-list()
  RankWeight<-list()
  RWRprob<-list()
  TargetScore<-list()
  outnames<-c()
  #anovause<-getMaxAovMatrix(anovascore)
  anovause<-anovascore
  for(i in 1:1){

    temp<-GetBMscores(anovause,RWRthreshold,PPInet,normPval,scoring,pval=pval,RWRthreshScore=RWRthreshScore,PvalthreshScore=PvalthreshScore,
                      FDRthreshScore=FDRscoreMat,pairList=pairList,pair=pair,DualBM=DualBM,
                      RegressScoreMat=RegressScoreMat,RegressEffectMat=RegressEffectMat,IntSplit=IntSplit,ctype=ctype)
    if(sum(temp$BiomarkerScore,na.rm=T)!=0){
      BiomarkerScore<-c(BiomarkerScore,temp$BiomarkerScore)
      SigScores<-cbind(SigScores,temp$SigScores)
      aovout<-temp$anovascore
      AnovaPval[[i]]<-temp$AnovaPval
      AnovaFdr[[i]]<-temp$AnovaFdr
      AnovaDeltaP[[i]]<-temp$AnovaDeltaP
      AnovaDeltaN[[i]]<-temp$AnovaDeltaN
      RWRprob[[i]]<-temp$RWRprob
      RankWeight[[i]]<-temp$RankWeight
      tset<-names(temp$TargetScore)
      PvalscoreClass<-temp$PvalscoreClass
      RWRscoreClass<-temp$RWRscoreClass
      TargetBMs<-temp$TargetBMs
      NumberBMs<-temp$NumberBMs
      MaxRWR<-temp$MaxRWR
      RWRclasses<-temp$RWRclasses
      PvalClasses<-temp$Pvalclasses
      FDRClasses<-temp$FDRclasses
      BMSymbol<-temp$BMSymbol
      BMtype<-temp$BMTypeClasses
      BMclass<-temp$BMclass
      Tout<-data.frame(MaxTargetScoreFDR=temp$MaxTargetScoreFDR,Score=temp$TargetScore,FDRScore=temp$TargetScoreFDR,BMtype=BMtype[tset],BMsymbol=BMSymbol[tset],BMclass=BMclass[tset],PvalScore=PvalscoreClass[tset],RWRscore=MaxRWR,BMs=TargetBMs[tset],NumberBMs=NumberBMs,RWRclasses=RWRclasses[tset],PvalClasses=PvalClasses[tset],FDRClasses=FDRClasses[tset],stringsAsFactors = FALSE)
      rownames(Tout)<-tset
      TargetScore[[i]]<-Tout

      outnames<-"Combined"
      names(AnovaPval)<-outnames
      names(AnovaFdr)<-outnames
      names(AnovaDeltaP)<-outnames
      names(AnovaDeltaN)<-outnames
      names(RWRprob)<-outnames
    }
  }

  return(list(anovascore=aovout,BiomarkerScore=BiomarkerScore,SigScores=SigScores,RWRprob=RWRprob,AnovaPval=AnovaPval,AnovaFdr=AnovaFdr,AnovaDeltaN=AnovaDeltaN,AnovaDeltaP=AnovaDeltaP,RankWeight=RankWeight,TargetScore=TargetScore))
}
combineRWRprobs<-function(RWRmatrix,pairList){
  allpairs<-unique(pairList[,c("Dep1","Dep2","pair")])

  bms<-colnames(RWRmatrix)

  RWRpair<-apply(allpairs,1,function(x) apply(RWRmatrix[x[1:2],],2,max))
  RWRm<-t(RWRpair)

  rownames(RWRm)<-allpairs[,"pair"]
  colnames(RWRm)<-bms
  return(RWRm)
}
combineAnovaRes<-function(anovares,pairList){
  associds<-unique(pairList[,"assoc_id"])
  allpairs<-unique(pairList[,"pair"])
  outAR<-NULL
  for(i in associds){
    inaov<-anovares[anovares[,"assoc_id"]==i,]
    newaov<-inaov[1,]
    newaov[,"AnovaPval"]<-min(inaov[,"AnovaPval"],na.rm=T)
    newaov[,"AnovaFdr"]<-min(inaov[,"AnovaFdr"],na.rm=T)
    newaov[,"Depleted.Gene"]<-pairList[pairList[,"assoc_id"]==i,"pair"]
    outAR<-rbind(outAR,newaov)

  }

  return(outAR)
}
combineRWRprobsDual<-function(RWRmatrix,pairList){
  #pairlist<-c(depgenes,res["assoc_id"],res["Depleted Gene"],string_id,type)
  #names(pairlist)<-c("Feat1","Feat2","assoc_id","pair","BM1","BM2","DualBM")

  outOr<-tryCatch({
    allpairs<-unique(pairList[,c("BM1","BM2","pair","DualBM")])

    allpairs<-split.data.frame(allpairs,f=unlist(allpairs[,"DualBM"]))
    allpairsOr<-allpairs[["Or"]]
    allpairsAnd<-allpairs[["And"]]

    bms<-colnames(RWRmatrix)
    deps<-rownames(RWRmatrix)
    if(nrow(allpairsOr)>1){
       allpairsOr<-allpairsOr[unlist(apply(allpairsOr[,1:2],1,function(x) sum(x%in%bms)==2)),]

      RWRdualOr<-apply(allpairsOr,1,function(x) apply(RWRmatrix[,unlist(x[1:2])],1,min))

      colnames(RWRdualOr)<-unlist(apply(allpairsOr[,c("BM1","BM2")],1,function(x) paste0(x,collapse="|")))
      rownames(RWRdualOr)<-deps
      RWRdualOr
    }else{
      if(sum(allpairsOr[1:2]%in%bms)==2){
        RWRdualOr<-matrix(min(RWRmatrix[,allpairsAnd[1:2]]),nrow=1,ncol=1)
        colnames(RWRdualOr)<-paste(allpairsOr[c("BM1","BM2")],collapse="|")
        rownames(RWRdualOr)<-deps
        RWRdualOr
      }
    }


  },
  error=function(cond) {
    return(NULL)
  })
  outAnd<-tryCatch({
    allpairs<-unique(pairList[,c("BM1","BM2","pair","DualBM")])

    allpairs<-split.data.frame(allpairs,f=unlist(allpairs[,"DualBM"]))
    allpairsOr<-allpairs[["Or"]]
    allpairsAnd<-allpairs[["And"]]

    bms<-colnames(RWRmatrix)
    deps<-rownames(RWRmatrix)

    if(nrow(allpairsAnd)>1){


      allpairsAnd<-allpairsAnd[unlist(apply(allpairsAnd[,1:2],1,function(x) sum(x%in%bms)==2)),]

      RWRdualAnd<-apply(allpairsAnd,1,function(x) apply(RWRmatrix[,unlist(x[1:2])],1,function(x) min(x)))

      colnames(RWRdualAnd)<-unlist(apply(allpairsAnd[,c("BM1","BM2")],1,function(x) paste0(x,collapse="&")))
      rownames(RWRdualAnd)<-deps
      RWRdualAnd
    }else{
      if(sum(allpairsAnd[1:2]%in%bms)==2){
        RWRdualAnd<-matrix(mean(RWRmatrix[,allpairsAnd[1:2]]),nrow=1,ncol=1)
        colnames(RWRdualAnd)<-paste(allpairsAnd[c("BM1","BM2")],collapse="&")
        rownames(RWRdualAnd)<-deps
        RWRdualAnd
      }


    }


  },
  error=function(cond) {
    return(NULL)
  })
  out<-cbind(outOr,outAnd)
  return(out)
}
combineAnovaResDual<-function(anovares,pairList){
  associds<-unique(pairList[,"assoc_id"])
  allpairs<-unique(pairList[,"pair"])
  outAR<-NULL
  for(i in associds){
    inaov<-anovares[anovares[,"assoc_id"]==i,]
    newaov<-inaov[1,]
    #how to combine anova fdrs for multiple biomarkers to single dependency?
    newaov[,"AnovaPval"]<-min(inaov[,"AnovaPval"],na.rm=T)
    newaov[,"AnovaFdr"]<-min(inaov[,"AnovaFdr"],na.rm=T)
    #
    newaov[,"FEATURE"]<-pairList[pairList[,"assoc_id"]==i,"pair"]
    outAR<-rbind(outAR,newaov)

  }

  return(outAR)
}
GetBMscores<-function(anovascore,RWRthreshold,PPInet,normPval,scoring,pval,RWRthreshScore=NULL,
                      PvalthreshScore=NULL,FDRthreshScore=NULL,pairList=NULL,pair,DualBM,
                      RegressScoreMat=NULL,RegressEffectMat=NULL,IntSplit=FALSE,ctype){


  anovascore<-as.data.frame(anovascore,stringsAsFactors=FALSE)
  anovascore$AnovaPval<-as.character(anovascore$AnovaPval)
  anovascore<-anovascore[anovascore[,"AnovaPval"]!="NULL",]
  anovascore<-anovascore[anovascore[,"AnovaFdr"]!="NULL",]
  anovascore<-anovascore[anovascore[,"StringID"]!="",]

  if(nrow(anovascore)>0){

    if((!is.null(pairList))&pair){
      #want to combine RWRscores across dependencies. The pvalscore  is the same for each dependency
      #combine the RWRthreshold and anovascore matrices back into the combined dependencies
      #RWRthreshold is RWR probs rows are gene symbol (targets) columns are biomarkers as string ids
      pairList<-pairList[pairList[,"assoc_id"]%in%anovascore[,"assoc_id"],]
      if(!IntSplit){
        RWRthreshold<-combineRWRprobs(RWRthreshold,pairList)
      }


      anovascore<-combineAnovaRes(anovascore,pairList)

    }
    if((!is.null(pairList))&DualBM){
      #want to combine RWRscores across biomarkers.
      pairList<-pairList[pairList[,"assoc_id"]%in%anovascore[,"assoc_id"],]
      #RWRthreshold<-combineRWRprobsDual(RWRthreshold,pairList)

      #anovascore<-combineAnovaResDual(anovascore,pairList)

    }
if(!is.null(RWRthreshold)){
  RWRthreshold<-matrix(RWRthreshold[,colnames(RWRthreshold)!=""],ncol=sum(colnames(RWRthreshold)!=""),dimnames=list(rownames(RWRthreshold),colnames(RWRthreshold)[colnames(RWRthreshold)!=""]))

  anovascore$AnovaPval<-unlist(anovascore$AnovaPval)
  anovascore$AnovaPval<-as.numeric(anovascore$AnovaPval)
  anovascore$AnovaFdr<-unlist(anovascore$AnovaFdr)
  anovascore$AnovaFdr<-as.numeric(anovascore$AnovaFdr)
  anovascore$AnovaDeltaP<-unlist(anovascore$AnovaDeltaP)
  anovascore$AnovaDeltaP<-as.numeric(anovascore$AnovaDeltaP)
  anovascore$AnovaDeltaN<-unlist(anovascore$AnovaDeltaN)
  anovascore$AnovaDeltaN<-as.numeric(anovascore$AnovaDeltaN)
  anovascore$Depleted.Gene<-unlist(anovascore$Depleted.Gene)
  anovascore$StringID<-unlist(anovascore$StringID)
  if(DualBM){
    AnovaPval<-acast(anovascore,Depleted.Gene~FEATURE,value.var="AnovaPval",fill=1,fun.aggregate = function(x) min(unlist(x)))
    AnovaFdr<-acast(anovascore,Depleted.Gene~FEATURE,value.var="AnovaFdr",fill=100,fun.aggregate = function(x) min(unlist(x)))
    AnovaDeltaP<-acast(anovascore,Depleted.Gene~FEATURE,value.var="AnovaDeltaP",fill=0,fun.aggregate = function(x) max(unlist(x)))
    AnovaDeltaN<-acast(anovascore,Depleted.Gene~FEATURE,value.var="AnovaDeltaN",fill=0,fun.aggregate = function(x) max(unlist(x)))

  }else{
    AnovaPval<-acast(anovascore,Depleted.Gene~StringID,value.var="AnovaPval",fill=1,fun.aggregate = function(x) min(unlist(x)))
    AnovaFdr<-acast(anovascore,Depleted.Gene~StringID,value.var="AnovaFdr",fill=100,fun.aggregate = function(x) min(unlist(x)))
    AnovaDeltaP<-acast(anovascore,Depleted.Gene~StringID,value.var="AnovaDeltaP",fill=0,fun.aggregate = function(x) max(unlist(x)))
    AnovaDeltaN<-acast(anovascore,Depleted.Gene~StringID,value.var="AnovaDeltaN",fill=0,fun.aggregate = function(x) max(unlist(x)))
  }

  AnovaPval[AnovaPval>=pval]=1
  PvalScore<-matrix(0,nrow=nrow(RWRthreshold),ncol=ncol(RWRthreshold))
  FDRScore<-matrix(0,nrow=nrow(RWRthreshold),ncol=ncol(RWRthreshold))
  PvalMatrix<-matrix(1,nrow=nrow(RWRthreshold),ncol=ncol(RWRthreshold))
  FDRMatrix<-matrix(100,nrow=nrow(RWRthreshold),ncol=ncol(RWRthreshold))
  LogPval<-matrix(0,nrow=nrow(RWRthreshold),ncol=ncol(RWRthreshold))
  SigScores<-matrix(0,nrow=nrow(RWRthreshold),ncol=ncol(RWRthreshold))
  ClassScore<-matrix(0,nrow=nrow(RWRthreshold),ncol=ncol(RWRthreshold))
  dimnames(LogPval)<-dimnames(RWRthreshold)
  dimnames(SigScores)<-dimnames(RWRthreshold)
  dimnames(PvalScore)<-dimnames(RWRthreshold)
  dimnames(PvalMatrix)<-dimnames(RWRthreshold)
  dimnames(FDRMatrix)<-dimnames(RWRthreshold)
  dimnames(FDRScore)<-dimnames(RWRthreshold)
  dimnames(ClassScore)<-dimnames(RWRthreshold)
  usegenes<-intersect(rownames(AnovaPval),rownames(LogPval))
  useT<-intersect(colnames(LogPval),colnames(AnovaPval))
  LogPval[usegenes,useT]<- -1*log(AnovaPval[usegenes,useT])
  PvalMatrix[usegenes,useT]<- AnovaPval[usegenes,useT]
  FDRMatrix[usegenes,useT]<-AnovaFdr[usegenes,useT]
  if(normPval=="minmax"){
    LogPval<-apply(LogPval,2,function(x) minmax(x))
  }
  if(!is.null(RWRthreshScore)){
    #also score using threshold groups for RWR probabilities:
    RWRthreshScore<-RWRthreshScore[order(RWRthreshScore[,"Score"],decreasing = FALSE),]

    RWRscore<-matrix(0,nrow=nrow(RWRthreshold),ncol=ncol(RWRthreshold))
    dimnames(RWRscore)<-dimnames(RWRthreshold)
    for(k in 1:nrow(RWRthreshScore)){
      RWRscore[RWRthreshold>RWRthreshScore[k,"RWRthresh"]]<-RWRthreshScore[k,"Score"]
    }
    PvalthreshScore<-PvalthreshScore[order(PvalthreshScore[,"Score"],decreasing=FALSE),]
    FDRthreshScore<-FDRthreshScore[order(FDRthreshScore[,"Score"],decreasing=FALSE),]


    for(k in 1:nrow(PvalthreshScore)){
      PvalScore[PvalMatrix<PvalthreshScore[k,"Pvalthresh"]]<-PvalthreshScore[k,"Score"]
    }

    for(k in 1:nrow(FDRthreshScore)){
      FDRScore[FDRMatrix<FDRthreshScore[k,"FDRthresh"]]<-FDRthreshScore[k,"Score"]
    }
    anovascore$AnovaClass<-"None"

    if(!is.null(RegressScoreMat)){

      AnovaClass<-apply(anovascore,1,function(x) Get_BMClass(x,RegressScoreMat,RegressEffectMat,DualBM=DualBM,ctype=ctype))
      anovascore$AnovaClass<-AnovaClass

    }
    anovascore$AnovaScore<-0

    anovascore[anovascore$AnovaClass=="A","AnovaScore"]<-100
    anovascore[anovascore$AnovaClass=="B","AnovaScore"]<-75
    anovascore[anovascore$AnovaClass=="C","AnovaScore"]<-50
    anovascore[anovascore$AnovaClass=="D","AnovaScore"]<-25
    RWRscorevec<-reshape2::melt(RWRscore)
    RWRscorevec$id<-paste(RWRscorevec$Var1,RWRscorevec$Var2,sep="-")
    if(DualBM){
      anovascore$idrwr<-paste(anovascore$Depleted.Gene,anovascore$FEATURE,sep="-")
    }else{
      anovascore$idrwr<-paste(anovascore$Depleted.Gene,anovascore$StringID,sep="-")
    }
    anovascore$RWRscore<-RWRscorevec[match(anovascore$idrwr,RWRscorevec$id),"value"]
    if(DualBM){
      ClassScoreV<-acast(anovascore,Depleted.Gene~FEATURE,value.var="AnovaScore",fill=0,fun.aggregate = function(x) max(unlist(x)))

    }else{
      ClassScoreV<-acast(anovascore,Depleted.Gene~StringID,value.var="AnovaScore",fill=0,fun.aggregate = function(x) max(unlist(x)))
    }

    ClassScore[usegenes,useT]<-ClassScoreV[usegenes,useT]




    RankWeight<-(RWRscore+PvalScore[rownames(RWRscore),colnames(RWRscore)])/2
    RankWeight<-RankWeight*((RWRscore!=0)*1)
    RankWeight<-RankWeight*((PvalScore!=0)*1)

    #RankWeightFDR<-(RWRscore+FDRScore[rownames(RWRscore),colnames(RWRscore)])/2
    RankWeightFDR<-(RWRscore+ClassScore[rownames(RWRscore),colnames(RWRscore)])/2
    RankWeightFDR<-RankWeightFDR*((RWRscore!=0)*1)
    #RankWeightFDR<-RankWeightFDR*((FDRScore!=0)*1)
    RankWeightFDR<-RankWeightFDR*((ClassScore!=0)*1)
    RWRscoreClass<-rowSums(RWRscore*((FDRScore!=0)*1))
    PvalscoreClass<-rowSums(PvalScore*((RWRscore!=0)*1))
    names(RWRscoreClass)<-rownames(RankWeight)
    names(PvalscoreClass)<-rownames(RankWeight)

  }else{
    RankWeight<-RWRthreshold*LogPval[rownames(RWRthreshold),colnames(RWRthreshold)]
  }

    if(scoring=="sum"){
      #this one will favour larger networks as long as there are 'some' high scoring p-values associated
      BiomarkerScore<-colSums(RankWeight)
      TargetScore<-rowSums(RankWeight)
      TargetScoreFDR<-rowSums(RankWeightFDR)

    }else{
      #average over so it is an overall measure of the strength between ppi network and anova p-value
      BiomarkerScore<-sapply(1:ncol(RankWeight),function(x) mean(RankWeight[RankWeight[,x]!=0,x]))
      BiomarkerScore<-unlist(BiomarkerScore)
      names(BiomarkerScore)<-colnames(RankWeight)
      TargetScore<-sapply(1:nrow(RankWeight),function(x) mean(RankWeight[x,RankWeight[x,]!=0]))
      TargetScore<-unlist(TargetScore)
      TargetScore[is.na(TargetScore)]<-0
      names(TargetScore)<-rownames(RankWeight)
      TargetScoreFDR<-sapply(1:nrow(RankWeightFDR),function(x) mean(RankWeightFDR[x,RankWeightFDR[x,]!=0]))
      TargetScoreFDR<-unlist(TargetScoreFDR)
      TargetScoreFDR[is.na(TargetScoreFDR)]<-0
      names(TargetScoreFDR)<-rownames(RankWeightFDR)

    }
  #MaxTargetScoreFDR<-apply(RankWeightFDR,1,max)
  #names(MaxTargetScoreFDR)<-rownames(RankWeightFDR)
  #TargetBMs<-sapply(1:nrow(RankWeightFDR),function(x) paste0(colnames(RankWeightFDR)[RankWeightFDR[x,]!=0],collapse="//"))
  #names(TargetBMs)<-rownames(RankWeightFDR)
  #NumberBMs<-sapply(1:nrow(RankWeightFDR),function(x) sum(RankWeightFDR[x,]!=0))
  #names(NumberBMs)<-rownames(RankWeightFDR)

  #BMInfo<-sapply(1:nrow(RankWeightFDR),function(x) Get_BMInfo(RankWeightFDR,
  #                                                            anovascore,PPInet,gene=rownames(RankWeightFDR)[x],DualBM=DualBM,RWRscore))

  BMInfo<-sapply(1:nrow(RankWeightFDR),function(x) Get_BMmax(RankWeightFDR,
                                                              anovascore,PPInet,gene=rownames(RankWeightFDR)[x],DualBM=DualBM,RWRscore))
  BMTypeClasses<-unlist(BMInfo["BMtypeSub",])
  BMSymbol<-unlist(BMInfo["BMsymbolSub",])
  BMclass<-unlist(BMInfo["BMclass",])
  MaxTargetScoreFDR<-unlist(BMInfo["MaxFDRwithR",])
  TargetBMs<-unlist(BMInfo["overlapMarkers",])

  NumberBMs<-unlist(BMInfo["numberBM",])
  MaxRWR<-unlist(BMInfo["MaxRWR",])
  names(BMTypeClasses)<-rownames(RankWeightFDR)
  names(BMSymbol)<-rownames(RankWeightFDR)
  names(BMclass)<-rownames(RankWeightFDR)
  names(MaxTargetScoreFDR)<-rownames(RankWeightFDR)
  names(TargetBMs)<-rownames(RankWeightFDR)
  names(NumberBMs)<-rownames(RankWeightFDR)
  names(MaxRWR)<-rownames(RankWeightFDR)


  useSymb<-intersect(names(TargetBMs),rownames(RWRscore))
  #RWRclasses<-sapply(1:nrow(RWRscore),function(x) paste0(RWRscore[x,FDRScore[x,]!=0&RWRscore[x,]!=0],collapse="//"))
  RWRclasses<-sapply(useSymb,function(x) paste0(RWRscore[x,intersect(unlist(strsplit(TargetBMs[x],"//",fixed=T)),colnames(RWRscore))],collapse="//"))
  names(RWRclasses)<-useSymb
  Pvalclasses<-sapply(1:nrow(PvalScore),function(x) paste0(PvalScore[x,RWRscore[x,]!=0&PvalScore[x,]!=0],collapse="//"))
  names(Pvalclasses)<-rownames(PvalScore)
  #FDRclasses<-sapply(1:nrow(FDRScore),function(x) paste0(FDRScore[x,RWRscore[x,]!=0&FDRScore[x,]!=0],collapse="//"))
  FDRclasses<-sapply(useSymb,function(x) paste0(FDRScore[x,intersect(unlist(strsplit(TargetBMs[x],"//",fixed=T)),colnames(FDRScore))],collapse="//"))
  names(FDRclasses)<-useSymb

  if(!DualBM){

    names(BiomarkerScore)<-anovascore[match(names(BiomarkerScore),anovascore[,"StringID"]),"FEATURE"]
  }

  SigScores[usegenes,useT]<-RankWeight[usegenes,useT]
  RWRprob<-RWRthreshold*as.logical(RankWeight!=0)
  if(!DualBM){
  colnames(RankWeight)<-anovascore[match(colnames(RankWeight),anovascore[,"StringID"]),"FEATURE"]
  colnames(RWRprob)<-anovascore[match(colnames(RWRprob),anovascore[,"StringID"]),"FEATURE"]
  colnames(SigScores)<-anovascore[match(colnames(SigScores),anovascore[,"StringID"]),"FEATURE"]
  colnames(AnovaPval)<-anovascore[match(colnames(AnovaPval),anovascore[,"StringID"]),"FEATURE"]
  colnames(AnovaFdr)<-anovascore[match(colnames(AnovaFdr),anovascore[,"StringID"]),"FEATURE"]
  colnames(AnovaDeltaP)<-anovascore[match(colnames(AnovaDeltaP),anovascore[,"StringID"]),"FEATURE"]
  colnames(AnovaDeltaN)<-anovascore[match(colnames(AnovaDeltaN),anovascore[,"StringID"]),"FEATURE"]
  }
  #remove any missing targets i.e not in network
  keepcols<-which(!is.na(colnames(RankWeight)))
  anovaout<-anovascore[!is.na(anovascore$RWRscore),]
  anovaout<-anovaout[anovaout$RWRscore!="0",]
  anovaout<-anovaout[anovaout$RWRscore!=0,]
  return(list(RWRclasses=RWRclasses,Pvalclasses=Pvalclasses,BiomarkerScore=BiomarkerScore,FDRclasses=FDRclasses,
              MaxTargetScoreFDR=MaxTargetScoreFDR,MaxRWR=MaxRWR,TargetScore=TargetScore,TargetScoreFDR=TargetScoreFDR,TargetBMs=TargetBMs,
              NumberBMs=NumberBMs,RWRscoreClass=RWRscoreClass,PvalscoreClass=PvalscoreClass,
              SigScores=SigScores,RWRprob=RWRprob,RankWeight=RankWeight,AnovaPval=AnovaPval,anovascore=anovaout,
              AnovaFdr=AnovaFdr,AnovaDeltaP=AnovaDeltaP,AnovaDeltaN=AnovaDeltaN,BMTypeClasses=BMTypeClasses,BMSymbol=BMSymbol,BMclass=BMclass))}else{
    return(list(BiomarkerScore=0,FDRclasses=0))
              }
  }else{
    return(list(BiomarkerScore=0,FDRclasses=0))
  }
}



Get_BMClass<-function(anovares,ScoreMat,EffectMat,DualBM=FALSE,BinaryMut=TRUE,ctype){
  anovares<-unlist(anovares)

  rownames(EffectMat)<-EffectMat[,"Omic"]
  if(DualBM){bmtype="Int"}else{
    bmtype<-anovares["type"]
  }
  if(!BinaryMut){
    if(bmtype%in%c("mut","var")){
      #make it use standard anova glass deltas otherwise use the regression t, f2 values
      bmtype<-"Genomic"
    }
  }else{
    if(bmtype%in%c("mut","var")){
      bmtype<-"Binary"
    }
  }
  if(bmtype%in%c("Int","Genomic")){
    #The AnovaFDR for Int is the maximum of Anova FDR and LRT fdr so both are covered.
    #could using the AnovaPval column if necessary.
    FDRcheck<-as.numeric(anovares["AnovaFdr"])<ScoreMat[,"FDRthresh"]
    GlassDelta<-as.numeric(anovares["AnovaDeltaP"])>1&as.numeric(anovares["AnovaDeltaN"])>1
    #if(bmtype=="Int"){
    #  IntCheck<-as.numeric(anovares["LRTFdr"]<ScoreMat[,"LRTthresh"])
    #  FDRcheck<-FDRcheck&IntCheck
    #}
    FDRcheck[1]<-FDRcheck[1]&GlassDelta




  }else{
    FDRcheck<-as.numeric(anovares["AnovaFdr"])<ScoreMat[,"FDRthresh"]
    tstat<-abs(as.numeric(anovares["AnovaDeltaN"]))>EffectMat[bmtype,"EffectSize"]
    Rsq<-as.numeric(anovares["AnovaDeltaP"])>EffectMat[bmtype,"GoF"]
    if(ctype=="PANCAN"){

    RsqB<-as.numeric(anovares["AnovaDeltaP"])>0.08
    tstatB<-abs(as.numeric(anovares["AnovaDeltaN"]))>1
    RsqC<-as.numeric(anovares["AnovaDeltaP"])>0.07
   # tstatD<-round(abs(as.numeric(anovares["AnovaDeltaN"])))>=5
    #RsqD<-as.numeric(anovares["AnovaDeltaP"])>0.06
    FDRcheck[1]<-FDRcheck[1]&tstat&Rsq
    FDRcheck[2]<-FDRcheck[2]&tstatB&RsqB
    FDRcheck[3]<-FDRcheck[3]&RsqC
   # FDRcheck[4]<-FDRcheck[4]&tstatD&RsqD
    }else{
      FDRcheck[1]<-FDRcheck[1]&tstat&Rsq
    }

  }
  classRes<-FDRcheck
  names(classRes)<-ScoreMat[,"ScoringClass"]
  classID<-ifelse(sum(classRes,na.rm=T)>0,sort(names(classRes)[classRes])[1],"None")
  return(classID)
}
Get_BMmax<-function(RankWeightFDR,anovares,PPInet,gene,DualBM,RWRscore){
  anovares<-anovares[anovares[,"Depleted.Gene"]==gene,]
  inputBMs<-colnames(RankWeightFDR)[RankWeightFDR[gene,]!=0]
  nres<-nrow(anovares)
  nna<-sum(is.na(anovares[,"AnovaClass"]))
  #the anovares is just related to the specific depleted gene.
  if(length(inputBMs)>0&(nres!=nna)){
    print(gene)

    overlapMarkers<-inputBMs

    MaxFDRwithR<-max(RankWeightFDR[gene,overlapMarkers])
    MaxRWR<-RWRscore[gene,overlapMarkers[which.max(RankWeightFDR[gene,overlapMarkers])]]
      #need to extract the best biomarker class and then select the RWR class to go with:
      #will only be one combo that satisfies this as we restrict to a single biomarker class
      overlapMarkers<-overlapMarkers[which(RankWeightFDR[gene,overlapMarkers]==MaxFDRwithR)]
      if(DualBM){
        BiomarkerClasses<-sapply(overlapMarkers,function(x) anovares[anovares[,"FEATURE"]==x,"AnovaClass"])
      }else{
        BiomarkerClasses<-sapply(overlapMarkers,function(x) anovares[anovares[,"StringID"]==x,"AnovaClass"])
      }
      if(length(overlapMarkers)>1){
        BestClassAll<-min(unlist(BiomarkerClasses))
        overlapMarkers<-overlapMarkers[which(unlist(lapply(BiomarkerClasses,function(x) min(x)))==BestClassAll)]
      }else{
        BestClassAll<-min(BiomarkerClasses)
      }

      if(DualBM){
        #bminfoall<-sapply(inputBMs,function(x) paste0(PPInet[which(PPInet[,"STRING_id"]%in%strsplit(x,"&|\\|")[[1]]),"symbol"],collapse=":"))
        bminfoall<-overlapMarkers
        bmtypeall<-sapply(overlapMarkers,function(x) unlist(anovares[which(anovares[,"FEATURE"]==x),"type"]))
      }else{
        bminfoall<-sapply(overlapMarkers,function(x) PPInet[which(PPInet[,"STRING_id"]==x),"symbol"])
        bmtypeall<-sapply(overlapMarkers,function(x) unlist(anovares[which(anovares[,"StringID"]==x),"type"]))
      }
      bminfo<-paste0(bminfoall,collapse="//")
      bmtype<-paste0(bmtypeall,collapse="//")
      BMclass<-BestClassAll
      overlapIDs<-paste0(overlapMarkers,collapse="//")
      ntypePerBM<-length(unique(unlist(bmtypeall)))

        bminfosub<-paste0(bminfoall,collapse="//")
        bminfobest<-paste0(bminfoall,collapse="//")

      bmtypesub<-paste0(unlist(bmtypeall),collapse="//")
      bmtypebest<-paste0(unlist(bmtypeall),collapse="//")
      numberBMs<-length(unlist(bmtypeall))

    if(MaxFDRwithR>0){

      return(data.frame(BMtype=bmtypesub,BMsymbol=bminfosub,BMtypeSub=bmtypebest,BMsymbolSub=bminfobest,BMclass=BMclass,MaxFDRwithR=MaxFDRwithR,overlapMarkers=overlapIDs,numberBM=numberBMs,MaxRWR=MaxRWR,stringsAsFactors = FALSE))
    }else{
      return(data.frame(BMtype="",BMsymbol="",BMtypeSub="",BMsymbolSub="",BMclass=NA,MaxFDRwithR=0,overlapMarkers="",numberBM=0,MaxRWR=0,stringsAsFactors = FALSE))

    }

  }else{
    return(data.frame(BMtype="",BMsymbol="",BMtypeSub="",BMsymbolSub="",BMclass=NA,MaxFDRwithR=0,overlapMarkers="",numberBM=0,MaxRWR=0,stringsAsFactors = FALSE))
  }
}
Get_BMInfo<-function(RankWeightFDR,anovares,PPInet,gene,DualBM,RWRscore){
  #MaxTargetScoreFDR<-apply(RankWeightFDR,1,max)


  anovares<-anovares[anovares[,"Depleted.Gene"]==gene,]
  inputBMs<-colnames(RankWeightFDR)[RankWeightFDR[gene,]!=0]
  nres<-nrow(anovares)
  nna<-sum(is.na(anovares[,"AnovaClass"]))
  #the anovares is just related to the specific depleted gene.
  if(length(inputBMs)>0&(nres!=nna)){
    print(gene)


    if(DualBM){
      #bminfoall<-sapply(inputBMs,function(x) paste0(PPInet[which(PPInet[,"STRING_id"]%in%strsplit(x,"&|\\|")[[1]]),"symbol"],collapse=":"))
      bminfoall<-inputBMs
      bmtypeall<-sapply(inputBMs,function(x) unlist(anovares[which(anovares[,"FEATURE"]==x),"type"]))
    }else{
      bminfoall<-sapply(inputBMs,function(x) PPInet[which(PPInet[,"STRING_id"]==x),"symbol"])
      bmtypeall<-sapply(inputBMs,function(x) unlist(anovares[which(anovares[,"StringID"]==x),"type"]))
    }
    bminfo<-paste0(bminfoall,collapse="//")
    bmtype<-paste0(bmtypeall,collapse="//")

    #can use AnovaClass to find the best marker class for a target
    BestClassAll<-sort(unique(anovares$AnovaClass))[1]
    #look for intersection first between BestClassAll markers and PPI markers:
    if(DualBM){
      overlapMarkers<-intersect(inputBMs,anovares[which(anovares[,"AnovaClass"]==BestClassAll),"FEATURE"])
    }else{
      overlapMarkers<-intersect(inputBMs,anovares[which(anovares[,"AnovaClass"]==BestClassAll),"StringID"])
    }
    if(length(overlapMarkers)>0){
      #have overlap between regression markers with best marker class and scoring PPI markers:
      #calculate max RankWeightFDR based on overlapping markers:

      overlapIDs<-paste0(overlapMarkers,collapse="//")
      MaxFDRwithR<-max(RankWeightFDR[gene,overlapMarkers])
      MaxRWR<-RWRscore[gene,overlapMarkers[which.max(RankWeightFDR[gene,overlapMarkers])]]
      BMclass<-BestClassAll
      ntypePerBM<-unlist(lapply(bmtypeall[inputBMs%in%overlapMarkers],length))
      if(max(ntypePerBM)==1){
        bminfosub<-paste0(bminfoall[which(inputBMs%in%overlapMarkers)],collapse="//")
        bminfobest<-paste0(bminfoall[which(inputBMs%in%overlapMarkers&RankWeightFDR[gene,inputBMs]==MaxFDRwithR)],collapse="//")
      }else{
        idx<-which(inputBMs%in%overlapMarkers)
        idx2<-which(inputBMs%in%overlapMarkers&RankWeightFDR[gene,inputBMs]==MaxFDRwithR)
        bminfosub<-paste0(unlist(sapply(1:length(idx),function(x) rep(bminfoall[idx[x]],ntypePerBM[x]))),collapse="//")
        bminfobest<-paste0(unlist(sapply(1:length(idx2),function(x) rep(bminfoall[idx2[x]],ntypePerBM[x]))),collapse="//")
      }
      bmtypesub<-paste0(unlist(bmtypeall[inputBMs%in%overlapMarkers]),collapse="//")
      bmtypebest<-paste0(unlist(bmtypeall[inputBMs%in%overlapMarkers&RankWeightFDR[gene,inputBMs]==MaxFDRwithR]),collapse="//")
      numberBMs<-length(unlist(bmtypeall[inputBMs%in%overlapMarkers]))
    }else{
      #can then calculate the AnovaClass for markers with PPI results
      if(DualBM){
        bmclasses<-sapply(inputBMs,function(x) unlist(anovares[which(anovares[,"FEATURE"]==x),"AnovaClass"]),simplify=FALSE)

      }else{
        bmclasses<-sapply(inputBMs,function(x) unlist(anovares[which(anovares[,"StringID"]==x),"AnovaClass"]),simplify=FALSE)

      }
       ntypePerBM<-unlist(lapply(bmclasses,length))
      BMclass<-sort(unlist(bmclasses))[1]
      if(max(ntypePerBM)!=1){
        inputBMs<-sapply(1:length(bmclasses),function(x) rep(inputBMs[x],ntypePerBM[x]),simplify=FALSE)
        inputBMs<-unlist(inputBMs)
        bmclasses<-unlist(bmclasses)
      }

      bmmax<-inputBMs[which(bmclasses==BMclass)]

      MaxFDRwithR<-max(RankWeightFDR[gene,bmmax],na.rm = T)
      MaxRWR<-RWRscore[gene,bmmax[which.max(RankWeightFDR[gene,bmmax])]]
      bmmax2<-bmmax[which(RankWeightFDR[gene,bmmax]==MaxFDRwithR)]
      #if(length(MaxFDRwithR)>0){
        overlapIDs<-paste0(bmmax,collapse="//")
        if(DualBM){
          #bminfosub<-paste0(unlist(sapply(bmmax,function(x) paste0(PPInet[which(PPInet[,"STRING_id"]%in%strsplit(x,"&|\\|")[[1]]),"symbol"],collapse=":"))),collapse="//")
          #bminfobest<-paste0(unlist(sapply(bmmax2,function(x) paste0(PPInet[which(PPInet[,"STRING_id"]%in%strsplit(x,"&|\\|")[[1]]),"symbol"],collapse=":"))),collapse="//")
          bminfosub<-paste0(bmmax,collapse="//")
          bminfobest<-paste0(bmmax2,collapse="//")

        }else{
          bminfosub<-paste0(PPInet[match(bmmax,PPInet$STRING_id),"symbol"],collapse="//")
          bminfobest<-paste0(PPInet[match(bmmax2,PPInet$STRING_id),"symbol"],collapse="//")
        }
        bmtypesub<-paste0(bmtypeall[which(bmclasses==BMclass)],collapse="//")
        bmtypebest<-paste0(bmtypeall[which(bmclasses==BMclass&RankWeightFDR[gene,inputBMs]==MaxFDRwithR)],collapse="//")
        numberBMs<-length(which(bmclasses==BMclass&RankWeightFDR[gene,inputBMs]==MaxFDRwithR))
    }
      if(MaxFDRwithR>0){

        return(data.frame(BMtype=bmtypesub,BMsymbol=bminfosub,BMtypeSub=bmtypebest,BMsymbolSub=bminfobest,BMclass=BMclass,MaxFDRwithR=MaxFDRwithR,overlapMarkers=overlapIDs,numberBM=numberBMs,MaxRWR=MaxRWR,stringsAsFactors = FALSE))
      }else{
        return(data.frame(BMtype="",BMsymbol="",BMtypeSub="",BMsymbolSub="",BMclass=NA,MaxFDRwithR=0,overlapMarkers="",numberBM=0,MaxRWR=0,stringsAsFactors = FALSE))

      }

  }else{
    return(data.frame(BMtype="",BMsymbol="",BMtypeSub="",BMsymbolSub="",BMclass=NA,MaxFDRwithR=0,overlapMarkers="",numberBM=0,MaxRWR=0,stringsAsFactors = FALSE))
  }
}
AnnotateBiomarkers<-function(BiomarkerScores,TPset=NULL){
  BMscores<-BiomarkerScores$BiomarkerScore
  if(!is.null(BMscores)){
  #AnovaPval<-BiomarkerScores$AnovaPval
  SigScores<-BiomarkerScores$SigScores
  if(is.null(TPset)){
  PriorityBiomarkers<-BMscores[round(BMscores,2)>0]}else{
    #have a set of biomarkers to use to define a tp distribution of scores:
    BMscores<-BMscores[BMscores!=0]
    BMgenes<-names(BMscores)
    BMgenes<-sapply(BMgenes,getGenesFromBiomarkers)
    tpset<-which(BMgenes%in%TPset)
    if(length(tpset)>5){
    dthresh<-densityThresholding(set1=BMscores[-tpset],set2=BMscores[tpset],x=seq(min(BMscores),max(BMscores),length.out=100))
    PriorityBiomarkers<-BMscores[BMscores>dthresh]}else{
      PriorityBiomarkers<-BMscores[BMscores>median(BMscores,na.rm=T)]
    }
  }
  PriorityBiomarkers<-PriorityBiomarkers[!is.na(names(PriorityBiomarkers))]
  if(length(PriorityBiomarkers)>1){
    SigPriority<-SigScores[,names(PriorityBiomarkers)]
    Biomarkerlists<-list()
    for(i in 1:ncol(SigPriority)){
      temp<-sort(SigPriority[SigPriority[,i]!=0,i],decreasing=TRUE)
      names(temp)<-names(sort(SigPriority[,i],decreasing = T))[sort(SigPriority[,i],decreasing = T)!=0]
      Biomarkerlists[[i]]<-temp
    }
    names(Biomarkerlists)<-colnames(SigPriority)
    return(Biomarkerlists)
  }else{
    if(length(PriorityBiomarkers)==1){
      Biomarkerlists<-list()
      SigPriority<-SigScores[,names(PriorityBiomarkers)]
      temp<-sort(SigPriority[SigPriority!=0],decreasing=TRUE)
      names(temp)<-names(sort(SigPriority,decreasing = T))[sort(SigPriority,decreasing = T)!=0]
      Biomarkerlists[[1]]<-temp
      names(Biomarkerlists)<-names(PriorityBiomarkers)
      return(Biomarkerlists)
    }else{
    return(NULL)}
  }
  }else{
    return(NULL)
  }
}

plotBiomarkerNetwork<-function(PriorityList,PPIgraph,PPInet,graphDist,vertexCex=10,tractability,savePDF=FALSE,dir.Results="."){
  Biomarkers<-names(PriorityList)
  for(i in 1:length(Biomarkers)){
    Genes<-unique(c(Biomarkers[i],names(PriorityList[[i]])))
    stringids<-PPInet[Genes,"STRING_id"]
    names(Genes)<-stringids
    if(length(Genes)>2){
      #need to have a reference set of shortest paths between all nodes
      #then we create a subgraph based on those shortest paths, with edges weighted according to path distances
      subgraph<-graphDist[stringids,stringids]
      subgraph<-1/subgraph
      subgraph[subgraph=="Inf"]<-0
      subgraph<-graph_from_adjacency_matrix(subgraph,weighted=TRUE)
      E(subgraph)$weight<-E(subgraph)$weight*3
      E(subgraph)$width<-E(subgraph)$weight
      id_name<-PPInet[names(PriorityList[[i]]),"STRING_id"]
      StringScore<-PriorityList[[i]]
      names(StringScore)<-id_name

      Vsizes<-StringScore[V(subgraph)$name]
      Vsizes<-(Vsizes-min(Vsizes,na.rm=T))/(max(Vsizes,na.rm=T)-min(Vsizes,na.rm=T))+1

      Vsizes[is.na(Vsizes)]<-min(Vsizes,na.rm=T)+1
      Vsizes<-Vsizes*vertexCex
      V(subgraph)$size<-Vsizes
      V(subgraph)$name<-Genes[V(subgraph)$name]
      #need to weight the edges and change the biomarker vertex shape
      vShape<-rep("circle",length(Vsizes))
      names(vShape)<-Genes
      vShape[Biomarkers[i]]<-"csquare"
      V(subgraph)$shape=vShape[V(subgraph)$name]
      tractabilityscores<-tractability[V(subgraph)$name,"min_bucket"]
      names(tractabilityscores)<-V(subgraph)$name
      tractabilityColor<-brewer.pal(n = 11, name = "Set3")
      names(tractabilityColor)<-as.character(1:11)
      vcolors<-tractabilityColor[as.character(tractabilityscores)]
      names(vcolors)<-names(tractabilityscores)
      if(!Biomarkers[i]%in%names(PriorityList[[i]])){
        vcolors[Biomarkers[i]]<-"white"
      }

      V(subgraph)$color<-vcolors[V(subgraph)$name]
      if(savePDF){
        try({
          pdf(paste0(dir.Results,"/",Biomarkers[i],"_BiomarkerNetwork.pdf"))
          {  plot.igraph(subgraph,main=Biomarkers[i],directed=F,edge.arrow.size=0,vertex.label.dist=3)

            # Add a legend
            legend("bottomleft", legend=levels(as.factor(tractabilityscores) ) , col = tractabilityColor[levels(as.factor(tractabilityscores) )] , bty = "n", pch=20 , pt.cex = 1, cex = 1,  horiz = FALSE, inset = c(0.1, 0.1))
          }
          dev.off()
        })
      }else{
        plot.igraph(subgraph,main=Biomarkers[i],directed=F,edge.arrow.size=0,vertex.label.dist=3)

        # Add a legend
        legend("bottomleft", legend=levels(as.factor(tractabilityscores) ) , col = tractabilityColor[levels(as.factor(tractabilityscores) )] , bty = "n", pch=20 , pt.cex = 1, cex = 1,  horiz = FALSE, inset = c(0.1, 0.1))
      }
    }
  }
}

getStationaryMatrix<-function(PPIigraph,vertexWeights=NULL,beta=0.6,PPInet=NULL,degreelist){
  #vertexWeights are those included in the depedencies so for which we want the Smat.

  #the weighted vertex matrix:
  #Fmat=diag(vertexWeights)

  #rownames(Fmat)<-stringNames
  #colnames(Fmat)<-stringNames

  #the weighted adjacency matrix:
  Amat<-as_adjacency_matrix(PPIigraph)+0
  #transition matrix, W:
  W<-diag(0,nrow=length(degreelist))
  W<-normalizeW(Amat,degreelist,normBy="Column")


  #Fmat<-Fmat[colnames(W),colnames(W)]
  #Pmatrix:
  Pmat<-beta*solve(diag(1,nrow=length(degreelist))-(1-beta)*W)
  #the stationary similarity matrix:
  #Smat<-Pmat%*%Fmat
  Smat<-Pmat
  return(Smat)
}

getRestartProb<-function(graph,useMatrix=FALSE){
  set.seed(1234)
  cts<-components(graph,"strong")
  maxct<-which.max(cts$csize)
  incvec<-names(cts$membership[cts$membership==maxct])
  adjacencymatrix<-as_adjacency_matrix(graph)
  adjmat<-adjacencymatrix[incvec,incvec]
  #graph<-induced_subgraph(graph,incvec)
  graph<-graph_from_adjacency_matrix(adjmat)
  betweenCentr<-betweenness(graph,directed=FALSE)
  betweenCentr<-betweenCentr[betweenCentr>0]
  quantileBC<-round(quantile(betweenCentr,probs=c(0.25,0.5,0.75)))
  quantileBC<-c(min(betweenCentr),quantileBC,max(betweenCentr))
  quantileBC<-unlist(quantileBC)
  selrandv<-sapply(quantileBC,function(x) selRand(x,betweenCentr))
  possBeta<-seq(from=0.05,to=1,length.out=20)
  testsets<-list()
  for(i in 1:length(selrandv)){
    testsets[[i]]<-getTestVertices(selrandv[i],graph)
  }

  #try using the probabilities from the random walk output as the 'influence'

  degreelist<-colSums(as.matrix(adjmat))
  rwrV<-list()
  if(useMatrix){

    for(j in 1:length(possBeta)){
      rwrV[[j]]<-getStationaryMatrix(graph,beta=possBeta[j],degreelist=degreelist)
    }
  }else{
    normBy<-"Column"
    if(normBy!="none"){
     normW<-normalizeW(adjmat,degreelist,normBy=normBy)}else{
      normW<-adjmat
      }
    #select biggest connected component:

    startvec<-rep(0,length(incvec))
    names(startvec)<-incvec

    for(i in 1:length(selrandv)){
      rwrV[[i]]<-list()
      svec<-startvec
      svec[selrandv[i]]<-1
      for(j in 1:length(possBeta)){
        rwrV[[i]][[j]]<-rwr(svec,normW,possBeta[j],PPInet=NULL)
      }
      names(rwrV[[i]])<-paste0("Beta: ",possBeta)
    }
    names(rwrV)<-selrandv
  }
 return(list(rwrres=rwrV,testsets=testsets))
}
selRand<-function(selval,inputvals){
  setposs<-which(inputvals==selval)
  if(length(setposs)<2){
    rangev<-c(0.9*selval,1.1*selval)
    setposs<-which(inputvals>rangev[1]&inputvals<rangev[2])
  }
  if(length(setposs)>1){
    getv<-sample(length(setposs),1)
    setposs<-setposs[getv]

  }
  out<-unlist(names(setposs))
  return(out)
}
getTestVertices<-function(startv,graph){
  neighbours<-neighbors(graph,startv)
  nb<-V(graph)$name[as.vector(neighbours)]
  nb<-c(nb,startv)
  spset<-distances(graph,startv,V(graph))
  dist2<-colnames(spset)[spset[1,]==2]
  rest<-V(graph)$name[!V(graph)$name%in%c(nb,dist2)]
  return(list(neighbour=nb,dist2=dist2,remain=rest))


}
summaryRestartProbs<-function(rwrRestart){
  rwrlist<-rwrRestart[[1]]
  testsets<-rwrRestart[[2]]
  possBeta<-seq(from=0.05,to=1,length.out=20)
  perGeneRes<-list()
  perGeneInflect<-list()
  rVals<-list()
  for(i in 1:length(rwrlist)){
    rwrbyBeta<-rwrlist[[i]]
    perGeneRes[[i]]<-list()
    perGeneInflect[[i]]<-list()
   testgene<-testsets[[i]]
    for(j in 1:length(possBeta)){
      outvals<-rwrbyBeta[[j]]
      outvec<-as.vector(outvals)
      names(outvec)<-rownames(outvals)
      valsbysets<-lapply(testgene,function(x) outvec[x])
      rangev<-seq(from=0,to=max(outvals),length.out=50)
      res<-list()

      for(k in 1:length(valsbysets)){
        res[[k]]<-unlist(sapply(rangev,function(x) sum(valsbysets[[k]]>x)))

      }
      perGeneInflect[[i]][[j]]<-rangev[which(res[[1]]<max(res[[1]]))[1]]
      perGeneRes[[i]][[j]]<-res
    }
   input<-unlist(perGeneInflect[[i]])
   overBeta<-input[2:length(input)]-input[1:(length(input)-1)]
   selBeta<-which(overBeta<0)[1]+1
   rVals[[i]]<-possBeta[selBeta]
   names(perGeneInflect[[i]])<-paste0("Restart Prob:",possBeta)

  }
  Beta<-as.numeric(names(which.max(table(unlist(rVals)))))
  return(list(InflectionData=perGeneInflect,RestartVals=rVals,Beta=Beta))
}
