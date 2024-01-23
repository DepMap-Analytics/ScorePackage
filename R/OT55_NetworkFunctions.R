
annotateBM<-function(geneset,gmtfile,ExprList){
  hyperres<-lapply(gmtfile,function(x) testGmt(x,geneset,ExprList))

  adjhyper<-p.adjust(unlist(hyperres),method="fdr")
  names(adjhyper)<-names(hyperres)
  smallestpval<-names(sort(unlist(hyperres),decreasing = FALSE))[1:5]
    tophyperres<-adjhyper[smallestpval]

  return(tophyperres)
}
getBiomarkerActivity<-function(ContMat,Module){
  useFC<-ContMat[unlist(Module),]
  temp<-prcomp(t(useFC),scale. = TRUE,center = TRUE)
  ModuleData<-temp$x[,1]
  return(ModuleData)
}
prepareBMInputs<-function(ContMat=NULL,PPIigraph=NULL,genemap,cmp,cancer_type="PANCAN",corMethod="spearman",topVar=1000,BEMDep=NULL){
  cosmicIDs<-intersect(colnames(ContMat),cmp$COSMIC_ID)
  if(length(cosmicIDs)==0){
    colnames(ContMat)<-cmp[match(colnames(ContMat),cmp$model_id),"COSMIC_ID"]
  }
  if(cancer_type!="PANCAN"){
    incCL<-intersect(colnames(ContMat),cmp[cmp$cancer_type==cancer_type,"COSMIC_ID"])
    ContMat<-ContMat[,incCL]

  }else{
    incCL<-colnames(ContMat)
  }
  if(!is.null(topVar)){
    bmVars<-apply(ContMat,1,function(x) var(x))
    genesinc<-names(sort(bmVars,decreasing=TRUE))[1:topVar]
    ContMat<-ContMat[genesinc,]
  }
  if(!is.null(BEMDep)){
    BEMDep<-BEMDep[,intersect(colnames(BEMDep),incCL)]
    #also need to check that number depletions >2 and less than n-2:
    BEMDep<-BEMDep[rowSums(BEMDep)>2&rowSums(BEMDep)<(ncol(BEMDep-2)),]
    gexpmarkers<-rownames(BEMDep)[grep('_Expr',rownames(BEMDep))]
    gexpmarkers<-unlist(sapply(gexpmarkers,function(x) strsplit(x,"_",fixed=TRUE)[[1]][1]))
    gexpmarkers<-unique(gexpmarkers)
    ContMat<-ContMat[rownames(ContMat)%in%gexpmarkers,]

  }
  corMat<-cor(t(ContMat),method=corMethod,use="complete.obs")
  return(list(corMat=corMat,ContMat=ContMat))

}

chunk <- function(x,n) split(x, factor(sort(rank(x)%%n)))
getBDcorr<-function(ScaledBFs,BinaryDep,set){
  subBinary<-BinaryDep[set,]
  incCL<-colSums(subBinary)!=0
  temp<-0
  if(sum(incCL)>9){

    temp<-cor(ScaledBFs[set[1],incCL],ScaledBFs[set[2],incCL],method="spearman",use="complete.obs")}
  return(temp)

}
f3 <- function(mat,ref) {mat*matrix(rep(ref,ncol=ncol(mat)),nrow=nrow(mat),ncol=ncol(mat),byrow=FALSE)}
subcor<-function(BFmatrix,rowIdx,binarymat,minCL=10,method="pearson"){
  binarymat<-binarymat[rownames(BFmatrix),intersect(colnames(BFmatrix),colnames(binarymat))]
  nullmatrix<-!binarymat+0
  idxvec<-nullmatrix[rowIdx,]
  useinput<-f3(t(nullmatrix),idxvec)
  #need to exclude the one doing from the useinput vec
  useinput[,rowIdx]<-0
  BFmatrix[t(useinput)==1]<-NA
  #also want to add in a minimum cell line filter later:
  BFmatrix[rowSums(!is.na(BFmatrix))<minCL,]<-NA
  corout<-cor(matrix(BFmatrix[rowIdx,],ncol = 1),t(BFmatrix),use="pairwise.complete.obs",method=method)
  outdf<-data.frame(gene1=rep(rownames(BFmatrix)[rowIdx],nrow(BFmatrix)),gene2=rownames(BFmatrix),corr=as.vector(corout))
  return(outdf)
}
prepareInputs<-function(ScaledBFs,BEMInput,PPIigraph,BinaryDep,genemap,cmp,cancer_type="PANCAN",maxCL=2,scaledBFMax=50,removeCore=NULL,FCs,method="pearson",allObs=TRUE,cellLineID="model_id",ncores,minCL=10){
  #need to convert inputs to make sure that we have the cell line matches.

  if(cancer_type!="PANCAN"){
    incCL<-intersect(colnames(ScaledBFs),cmp[cmp$cancer_type==cancer_type,cellLineID])
    ScaledBFs<-ScaledBFs[,incCL]

  }else{incCL<-colnames(ScaledBFs)}
  #useCL<-intersect(colnames(ScaledBFs),colnames(BEMInput))
  #changed 12.1.21 to get correlations using all possible cell lines not just those with entries in BEM
  useCL<-incCL
  BinaryDep<-BinaryDep[,useCL]
  BinaryDep<-BinaryDep[intersect(rownames(BinaryDep),rownames(ScaledBFs)),]
  #if(allObs){
    usegenes<-rownames(BinaryDep)[rowSums(BinaryDep)>2&rowSums(BinaryDep)<(ncol(BinaryDep)-maxCL)]
  #}else{
  #  usegenes<-rownames(ScaledBFs)
  #}
  if(!is.null(removeCore)){
    usegenes<-setdiff(usegenes,removeCore)
    print(head(usegenes))
  }
  print(head(usegenes))
  mutBem<-BEMInput[grep("_mut",rownames(BEMInput)),]
  Mgenes<-unlist(sapply(rownames(mutBem),function(x) strsplit(x,"_mut")[[1]][1]))
  stringM<-genemap[Mgenes,"STRING_id"]
  stringD<-genemap[usegenes,"STRING_id"]
  vertices<-which(as_ids(V(PPIigraph))%in%c(stringM,stringD))
  graphg<-c(Mgenes,usegenes)
  graphgenes<-graphg[which(c(stringM,stringD)%in%as_ids(V(PPIigraph)))]
  PPIsubgraph<-induced_subgraph(PPIigraph,vertices)
  ScaledBFs<-ScaledBFs[usegenes,useCL]
  ScaledBFs[ScaledBFs==Inf]<-scaledBFMax
  FCs<-FCs[usegenes,useCL]
  tempCor<-cor(t(ScaledBFs),method=method)
  if(allObs){

    corMat<-tempCor
    subCM<-corMat
  }else{
    #genepairs<-combn(rownames(ScaledBFs),2)
    corMat<-matrix(0,nrow=nrow(ScaledBFs),ncol=nrow(ScaledBFs))
    rownames(corMat)<-rownames(ScaledBFs)
    colnames(corMat)<-rownames(ScaledBFs)
    ngenes<-nrow(ScaledBFs)

    #nMarkers<-ncol(genepairs)

    idxs<-chunk(1:ngenes,ncores)


    res<-foreach(j=1:length(idxs),.combine=rbind)%dopar%{
      #for(j in 1:length(idxs)){
      #inputdata<-genepairs[,idxs[[j]]]
      #corrvals<-apply(inputdata,2,function(x) getBDcorr(ScaledBFs,BinaryDep,x))
      inidx<-length(idxs[[j]])
      tempout<-NULL
      pb <- txtProgressBar(min=1,max=inidx,style=3)
      for(k in 1:inidx){
        temp<-subcor(ScaledBFs,idxs[[j]][k],BinaryDep,method=method,minCL=minCL)
        tempout<-rbind(tempout,temp)
        setTxtProgressBar(pb, k)
      }


      #tempout<-data.frame(gene1=inputdata[1,],gene2=inputdata[2,],corr=unlist(corrvals))

      tempout
    }

    #diagres<-data.frame(gene1=rownames(ScaledBFs),gene2=rownames(ScaledBFs),corr=1)
    #res<-rbind(res,diagres)
    corMat<-acast(res,gene1~gene2, value.var="corr")
    corMat<-corMat[,rownames(corMat)]
    temp<-corMat
    corMat<-t(corMat)
    corMat[upper.tri(corMat)]<-temp[upper.tri(temp)]
    diag(corMat)<-1
    corMat[is.na(corMat)]<-0
    subCM<-corMat
    #tempCor
    corMatMax<-pmax(abs(tempCor),abs(corMat),na.rm=T)
    #sign of correlations:
    corSigns<-sapply(1:nrow(tempCor),function(x) which.pmax(abs(tempCor[x,]),abs(corMat[x,])))
    corSigns<-t(corSigns)
    s1<-sign(tempCor)
    s2<-sign(corMat)
    signMat<-s1
    signMat[corSigns==2]<-s2[corSigns==2]
    corMat<-corMatMax*signMat

  }
  RefGraphDegree<-degree(PPIsubgraph)
  return(list(ScaledBFs=ScaledBFs,PPIigraph=PPIsubgraph,BEMmatrix=mutBem,graphgenes=graphgenes,DepMatrix=BinaryDep[usegenes,useCL],corMat=corMat,RefGraphDegree=RefGraphDegree,FCs=FCs,subCM=subCM))
}



DynamicClusters<-function(Smat,simThreshold,minClust,PPIigraph=NULL,PPInet,lout=10,PPIasPrior=FALSE,lambda=NULL,Aggregate=FALSE,splitClust=10){
  SCClist<-list()
  CorLevel<-list()
  AnnotationDF<-NULL
  if(!is.null(PPIigraph)){
    stringNames<-PPInet[rownames(Smat),"STRING_id"]
    names(stringNames)<-rownames(Smat)
    stringNames<-stringNames[!is.na(stringNames)]
    Smat<-Smat[names(stringNames),names(stringNames)]
    BioMarkerOnly<-setdiff(V(PPIigraph)$name,stringNames)
    PPIigraph<-delete.vertices(PPIigraph,which((V(PPIigraph)$name)%in%BioMarkerOnly))
    Amat<-as_adjacency_matrix(PPIigraph)+0

    incString<-intersect(stringNames,colnames(Amat))
    IncGenes<-PPInet[PPInet[,"STRING_id"]%in%incString,"symbol"]
    IncGenes<-intersect(IncGenes,rownames(Smat))
    Smat<-Smat[IncGenes,IncGenes]
    stringNames<-stringNames[IncGenes]
    rownames(Smat)<-stringNames
    colnames(Smat)<-stringNames

    if(PPIasPrior){
      #invert the Amat, for prior network, then modfiy the Smat correlations accordingly
      #have a new function for choosing level of regularization and selecting prior corrected thresholds
      penaltymatrix<-!Amat+0
      Smat<-PriorRegGraph(Smat,penaltymatrix,lambda)
    }else{

      #then  want to filter correlation matrix using prior graph:

      Smat<-Smat*Amat[stringNames,stringNames]
      rownames(Smat)<-IncGenes
      colnames(Smat)<-IncGenes
      print(head(IncGenes))
      Smat<-as.matrix(Smat)
    }
  }
  simVector<-seq(from=simThreshold[1],to=simThreshold[2],length.out=lout)
  print(simVector)
  j=1
  for(i in 1:lout){
    if(is.matrix(Smat)){
      if(nrow(Smat)>2){
        SCheck<-abs(Smat)>=simVector[i]+0
        Sgraph<-graph.adjacency(SCheck)
        SCCs<-components(Sgraph,mode="weak")
        csizes<-SCCs$csize
        incClusters<-which(csizes>=minClust)

        if(length(incClusters)>0){
          for(k in 1:length(incClusters)){
            temp<-names(SCCs$membership[SCCs$membership==incClusters[k]])
            print(temp)
            if(length(temp)>splitClust){
              #split large cluster into subclusters:
              Smat[abs(Smat)<simVector[i]]<-0
              subclust<-SubClustGraph(Smat[temp,temp])
              j=j-1
              for(a in 1:length(subclust)){
                SCClist[[j+a]]<-subclust[[a]]
                CorLevel<-c(CorLevel,simVector[i])
                name<-NULL
                gene_family<-NULL
                string_id<-NULL
                name<-c(name,paste0("Group",j+a))
                gene_family<-c(gene_family,paste0(subclust[[a]],collapse="//"))
                string_id<-c(string_id,paste0(PPInet[subclust[[a]],"STRING_id"],collapse="//"))
                tADF<-data.frame(cbind(name=name,gene_family=gene_family,location_sortable=NA,
                                       entrez_id=NA,
                                       ensembl_gene_id=NA,
                                       pubmed_id=NA,string_id=string_id),stringsAsFactors = FALSE)
                if(is.null(AnnotationDF)){
                  AnnotationDF<-tADF
                }else{
                  AnnotationDF<-rbind(AnnotationDF,tADF)
                }

              }
              j=j+length(subclust)+1
            }else{
              SCClist[[j]]<-temp
              CorLevel<-c(CorLevel,simVector[i])
              name<-NULL
              gene_family<-NULL
              string_id<-NULL
              name<-c(name,paste0("Group",j))
              gene_family<-c(gene_family,paste0(temp,collapse="//"))
              string_id<-c(string_id,paste0(PPInet[temp,"STRING_id"],collapse="//"))
              tADF<-data.frame(cbind(name=name,gene_family=gene_family,location_sortable=NA,
                                   entrez_id=NA,
                                   ensembl_gene_id=NA,
                                   pubmed_id=NA,string_id=string_id),stringsAsFactors = FALSE)
              if(is.null(AnnotationDF)){
                AnnotationDF<-tADF
              }else{
                AnnotationDF<-rbind(AnnotationDF,tADF)
              }
              j=j+1
            }
          }

        }


          if(!Aggregate){
            Smat<-Smat[!rownames(Smat)%in%unlist(SCClist),!colnames(Smat)%in%unlist(SCClist)]
          }
          print(length(SCClist))
        }
      }
    }

  if(length(SCClist)>0){
    names(SCClist)<-paste0("Group",1:length(SCClist))
    names(CorLevel)<-paste0("Group",1:length(CorLevel))
    AnnotationDF$name<-paste0("Group",1:nrow(AnnotationDF))
    rownames(AnnotationDF)<-AnnotationDF$name

  }
  if(Aggregate){
    useGroups<-names(SCClist[!duplicated(SCClist)])
    SCClist<-SCClist[useGroups]
    CorLevel<-CorLevel[useGroups]
    AnnotationDF<-AnnotationDF[useGroups,]
  }

    return(list(CorGroups=SCClist,CorLevels=CorLevel,ADF=AnnotationDF))

}
SubClustGraph<-function(InputData){
  cmat<-1-abs(InputData)
  res<-NbClust(diss=as.dist(cmat), distance = NULL, min.nc=2, max.nc=nrow(cmat)-2,
               method = "ward.D2", index = "silhouette")
  nclust<-1:res$Best.nc[1]
  values<-sapply(nclust,function(x) names(res$Best.partition)[res$Best.partition==x])
  setlengths<-lapply(values,length)
  values<-values[setlengths>1]
  isnull<-lapply(values,function(x) !is.null(as.character(x)))
  values<-values[unlist(isnull)]
  return(values)
}

corPairs<-function(Smat,simThreshold,PPIigraph=NULL,PPInet){
  SCClist<-list()
  CorLevel<-list()
  AnnotationDF<-NULL
  if(!is.null(PPIigraph)){
    stringNames<-PPInet[rownames(Smat),"STRING_id"]
    names(stringNames)<-rownames(Smat)
    stringNames<-stringNames[!is.na(stringNames)]
    Smat<-Smat[names(stringNames),names(stringNames)]
    BioMarkerOnly<-setdiff(V(PPIigraph)$name,stringNames)
    PPIigraph<-delete.vertices(PPIigraph,which((V(PPIigraph)$name)%in%BioMarkerOnly))
    Amat<-as_adjacency_matrix(PPIigraph)+0

    incString<-intersect(stringNames,colnames(Amat))
    IncGenes<-PPInet[PPInet[,"STRING_id"]%in%incString,"symbol"]
    IncGenes<-intersect(IncGenes,rownames(Smat))
    Smat<-Smat[IncGenes,IncGenes]
    stringNames<-stringNames[IncGenes]
    rownames(Smat)<-stringNames[IncGenes]
    colnames(Smat)<-stringNames[IncGenes]

    #then  want to filter correlation matrix using prior graph:

   Smat<-Smat*Amat[stringNames,stringNames]
   rownames(Smat)<-IncGenes
   colnames(Smat)<-IncGenes
   print(head(IncGenes))
   Smat<-as.matrix(Smat)

  }



    if(is.matrix(Smat)){
      if(nrow(Smat)>2){
        Cmat<-Smat
        SCheck<-abs(Smat)>=simThreshold+0
        SCheck<-as.matrix(SCheck)
        SCheck[upper.tri(SCheck)] = NA
        diag(SCheck)<-NA

        pairVals<-reshape2::melt(SCheck, na.rm=T)

        Cmat[upper.tri(Cmat)]=NA
        diag(Cmat)=NA
        pairCors<-reshape2::melt(Cmat,na.rm=T)
        pairCors<-cbind(pairCors,pairVals[,3])
        pairCors<-pairCors[pairCors[,4],]


        }
      }

rownames(pairCors)<-paste0("Group",1:nrow(pairCors))
name<-paste0("Group",1:nrow(pairCors))
gene_family=unlist(apply(pairCors,1,function(x) paste(x[1],x[2],sep="//")))
string_id<-unlist(apply(pairCors,1,function(x) paste(PPInet[x[1],"STRING_id"],PPInet[x[2],"STRING_id"],sep="//")))
ADF<-data.frame(cbind(name=name,gene_family=gene_family,location_sortable=NA,
                      entrez_id=NA,
                      ensembl_gene_id=NA,
                      pubmed_id=NA,string_id=string_id),stringsAsFactors = FALSE)


  return(list(CorGroups=pairCors,ADF=ADF))

}



getAovData<-function(ANOVAdirInput,DepList,OUTPUT_DIR,ModuleFCs,InputName="InputFeatures",OutputName="",cmp,BEMid="COSMIC_ID",ESSprofiles=NULL,biomarkertype="mut"){

  #  AovData<-getAovData(ANOVAdirInput,DepsByCoEssPPIAgg,OUTPUT_DIR,ModuleFCs,InputName="InputFeatures",OutputName="AggBD",cmp=cmp)

  geneAnnotations<-DepList$ADF
  rownames(geneAnnotations)<-geneAnnotations$name

  save(geneAnnotations,file=paste0(OUTPUT_DIR,"Annotation.Rdata"))
  #load the input BEM features for the ANOVA want to run:

  load(paste0(ANOVAdirInput,InputName,".rdata"))
  if(tolower(biomarkertype)%in%c("mut","genomic")){
  BEM<-InputFeatures$BEM}else{
    InputFeatures<-list()
    InputFeatures$BEM<-BEM
  }



  #get module data into correct form and check overlaps the BEM features:

  #ESSprofiles<-matrix(unlist(ModuleFCs),ncol=length(ModuleFCs[[1]]),byrow = TRUE)
  if(is.null(ESSprofiles)){
  if(is.list(ModuleFCs)){
    ESSprofiles<-do.call(rbind,ModuleFCs)
    colnames(ESSprofiles)<-names(ModuleFCs[[1]])
  }else{
    ESSprofiles<-t(ModuleFCs)
  }
  }
  if(BEMid!="model_id"){
    temp<-cmp[match(colnames(BEM),cmp[,BEMid]),"model_id"]
    matchMIDs<-cmp[cmp[,BEMid]%in%colnames(BEM),c("model_id",BEMid)]
    matchMIDs<-matchMIDs[matchMIDs$model_id%in%colnames(ESSprofiles),]
    colnames(BEM)<-matchMIDs[match(colnames(BEM),matchMIDs[,2]),1]

    colnames(BEM)<-make.names(colnames(BEM))
  }else{

    temp<-CLnameMapping(ESSprofiles,BEM,"col",annotation=cmp)
    ESSprofiles<-temp$inputdata
    BEM<-temp$refdata
  }
  colnames(ESSprofiles)<-make.names(colnames(ESSprofiles))
  useCL<-intersect(colnames(ESSprofiles),colnames(BEM))

  InputFeatures$BEM<-BEM[,useCL]
  ESSprofiles<-ESSprofiles[,useCL]
  if(is.null(dim(DepList[[1]]))){
    inputnames<-names(DepList[[1]])
    inputnames<-inputnames[unlist(lapply(ModuleFCs,function(x) !is.null(x)))]
    rownames(ESSprofiles)<-inputnames
  }else{
    genepairn<-apply(DepList[[1]],1,function(x) ifelse(as.numeric(x[3])<0,paste(x[1],x[2],sep = "|"),paste(x[1],x[2],sep="&")))
    rownames(ESSprofiles)<-unlist(genepairn)
    refnames<-apply(DepList[[1]],1,function(x) paste(x[1],x[2],sep = "//"))
    geneAnnotations<-geneAnnotations[match(refnames,geneAnnotations$gene_family),]
    rownames(geneAnnotations)<-genepairn
  }
  ESSprofiles<-t(ESSprofiles)
  bDepletions<-t(bDepletions)
  save(ESSprofiles,file=paste0(OUTPUT_DIR,"/ESSProfiles",OutputName,".Rdata"))
  save(InputFeatures,file=paste0(OUTPUT_DIR,"/InputFeatures",OutputName,".Rdata"))
  return(list(ESSprofiles=ESSprofiles,InputFeatures=InputFeatures,geneAnnotations=geneAnnotations))
}

getPairCombFC<-function(FCs,Pair,corV){
  corType<-ifelse(as.numeric(corV)<0,"neg","pos")
  ModuleFC<-NULL
  if(corType=="neg"){
    ModuleFC<-apply(FCs[Pair,],2,min)
  }
  if(corType=="pos"){
    ModuleFC<-apply(FCs[Pair,],2,mean)
  }
  return(ModuleFC)
}
getModuleFC<-function(FCs,Module){
  inputs<-unlist(Module)
  if(length(inputs)>0){
  useFC<-FCs[unlist(Module),]
  temp<-prcomp(t(useFC),scale. = TRUE,center = TRUE)
  ModuleData<-temp$x[,1]
  return(ModuleData)}else{
    return(NULL)
  }
}
dataPC1<-function(data){

  estpca<-prcomp(t(data),scale.=FALSE,center=TRUE)
  pca<-prcomp(t(data),scale.=FALSE,center=FALSE)
  firstPC<-estpca$rotation[,1]
  #firstPC<-pca$rotation[,1]
  unscaledVal<-t(data)%*%firstPC
  singledata<-scale(t(data),scale=FALSE)%*%firstPC
  #unscale data:
  centreAvg<-estpca$center%*%firstPC
  #singledata<-singledata+centreAvg
  #return(singledata)
  #return(unscaledVal)
  return(estpca$x[,1])
}
getEigenGenes<-function(FCs,Module){
  X<-FCs[unlist(Module),]
  datModule=t(scale(t(X)))
  eigengenes <- svd(datModule)$v
  scaledExpr = scale(t(datModule))
  averExpr= rowMeans(scaledExpr, na.rm = TRUE)
  corAve = cor(averExpr, eigengenes, use = "p")
  return(list(eigenG=eigengenes[,1],avgExpr=averExpr))

}
RemovePC<-function(data,droppcanumber=1,perfCheck=TRUE){
  estpca<-prcomp(t(data),scale.=TRUE)
  npcas<-1:ncol(data)
  pcause<-npcas[!npcas%in%droppcanumber]
  df.denoised <- estpca$x[,pcause] %*% t(estpca$rotation[,pcause])
  df.denoised<-t(df.denoised)
  correctedData<-df.denoised*estpca$scale+estpca$center
  #test data
  if(perfCheck){
    res<-classPerfCP(correctedData)
    return(list(correctedData=correctedData,Res=res))
  }else{
    return(correctedData)
  }
}

getSimilarityMatrix<-function(PPIigraph,vertexWeights,beta=0.6,PPInet,corThreshold=0.5,corMat,degreelist,useDegree=FALSE){
  #vertexWeights are those included in the depedencies so for which we want the Smat.
  #degreelist<-degree(PPIigraph)
  Isolated = which(degreelist==0)
  PPIigraph<-delete.vertices(PPIigraph,Isolated)
  #remove from graph those not in depedencies, i.e. potential biomarkers only:
  stringNames<-PPInet[names(vertexWeights),"STRING_id"]
  BioMarkerOnly<-setdiff(V(PPIigraph)$name,stringNames)
  PPIigraph<-delete.vertices(PPIigraph,which((V(PPIigraph)$name)%in%BioMarkerOnly))
  PPInames<-V(PPIigraph)$name
  if(!useDegree){
    degreelist<-degreelist[degreelist!=0]
    degreelist<-degreelist[PPInames]
    degreelist<-range01(degreelist)
    Dmat<-diag(degreelist,nrow=length(degreelist))
  }else{
    #use the above if want to penalise central "hubs" genes, alternative don't use the network topology to determine "subnetwork":
    Dmat<-diag(1,nrow=length(degreelist))
  }
  #vertexWeights<-range01(vertexWeights)
  colnames(Dmat)<-names(degreelist)
  rownames(Dmat)<-names(degreelist)
  #the weighted vertex matrix:
  Fmat=diag(vertexWeights)

  rownames(Fmat)<-stringNames
  colnames(Fmat)<-stringNames
  #the weighted adjacency matrix:
  Amat<-as_adjacency_matrix(PPIigraph)+0
  #transition matrix, W:
  W<-diag(0,nrow=length(degreelist))
  EdgeWeights<-corMat[names(vertexWeights),names(vertexWeights)]
  rownames(EdgeWeights)<-stringNames
  colnames(EdgeWeights)<-stringNames
  EdgeWeights[abs(EdgeWeights)<corThreshold]<-0
  #Amat<-Amat*EdgeWeights[PPInames,PPInames]
  W=Amat/degreelist

  Fmat<-Fmat[colnames(W),colnames(W)]
  #Pmatrix:
  Pmat<-beta*solve(diag(1,nrow=length(degreelist))-(1-beta)*W)
  #the stationary similarity matrix:
  Smat<-Pmat%*%Fmat
  return(Smat)
}


SCCbyCLlist<-function(SmatList){
  #create a list of SCCs per cell line - excluding singletons
  SCC_list<-list()
  for(i in 1:length(SmatList)){
    SCCs<-components(SmatList[[i]])
    k=1
    for(j in 1:length(SCCs)){
      if(length(SCCs[[j]])>1){
        SCCsublist[[k]]<-SCCs[[j]]
        k=k+1
      }
    }
    SCC_list[[i]]<-SCCsublist
  }
  #names should be cell line names
  names(SCC_list)<-names(SmatList)
}
getCLperModel<-function(BMgraph,BEM){
  #for a given graph identify cell lines that should be included for each biomarker

}
scoreSCC<-function(SCCList,SCC,CLuse,TPscore,FNscore){
  #SCC can be updated to combined values based on hierarchy of Biomarkers.
  perCLscore<-list()
  for(i in 1:length(CLuse)){
    SCCcl<-SCCList[[CLuse[i]]]
    #so SCCcl is a list of SCC for the CL excluding singletons.
    #for the moment include all SCC irrespective of potential different clusterings
    allNodes<-unlist(SCCcl)
    TPs<-length(intersect(allNodes,SCC))*TPscore
    FNs<-length(setdiff(SCC,allNodes))*FNscore
    perCLscore[[i]]<-c(TPs,FNs)
  }
  return(perCLscore)
}
range01 <- function(x){(x-min(x)+0.00000001)/(max(x)-min(x))}


#given a particular cluster of the score "significance" of SCC's and use these to score the models

#instead of the delta leaf structure, need to use the possible BEM entries as perturbations, hierarchy
#need a method that assigns pertubation to network subcluster - assume it should be closest to perturbation node

#start by defining all SCC's in the full similarity matrix
#match SCCs to the closest possible perturbation(BEMselected)
#use NEM framework to define hierarchy between BEMselected
#score resulting networks using existing SCC's across replicates (cell lines) and their strength? - False positives and negative rates, and use cell line replicates
#ala NEM, for scoring instead of having to create the random dendrogram clusters in the original Hierarchical HotNet
#dependencies can be further scored by connectivity in SCC and maybe distance to biomarker
matchSCCtoBEM<-function(Smat,BEMmatrix,referenceNetwork,removeNeg=TRUE,minClust=3,PPInet,simThreshold=0.02){
  #use threshold of known essential genes and their proximity =1 on PPI reference network to determine threshold simply.
  #assume all entries in BEMmatrix are possible biomarkers
  possibleMarkers<-names(BEMmatrix)
  possibleMarkers<-unlist(sapply(possibleMarkers,function(x) strsplit(x,"_")[[1]][1]))
  possibleMarkers<-PPInet[possibleMarkers,"STRING_id"]
  if(length(simThreshold)==1){
    if(removeNeg){
      Smat[abs(Smat)<simThreshold]<-0
      Smat[Smat!=0]<-1
    }
    Sgraph<-graph.adjacency(Smat)
    SCCs<-components(Sgraph,mode="strong")
    csizes<-SCCs$csize
    incClusters<-which(csizes>=minClust)
    SCClist<-sapply(incClusters,function(x) names(SCCs$membership[SCCs$membership==x]))
  }else{

    SCClist<-DynamicClusters(Smat,simThreshold,minClust)
  }
  fullMarkerName<-names(BEMmatrix)
  if(is.list(SCClist)){
    minPerSCC<-lapply(SCClist,function(x) minBEM(x,possibleMarkers,referenceNetwork,fullMarkerName))
    genelist<-lapply(SCClist,function(x) rownames(PPInet)[PPInet$STRING_id%in%x])

  }else{
    minPerSCC<-minBEM(SCClist,possibleMarkers,referenceNetwork,fullMarkerName)
    genelist<-rownames(PPInet)[PPInet$STRING_id%in%SCClist]
  }
  return(list(minMarkers=minPerSCC,genelist=genelist,StringIDlist=SCClist))
}

SCCtoBEMmap<-function(SmatList,BEMmatrix,referenceNetwork,removeNeg=TRUE,PPInet,minClust=3,simThreshold){
  #get all possible SCCs from each cell line and their associated closest biomarkers
  #restrict SCC and biomarker map to Biomarkers present for each cell line
  SCCmap<-list()
  celllines<-names(SmatList)
  for(i in 1:length(SmatList)){
    BEMmat<-as.matrix(BEMmatrix[,celllines[i]])
    BEMx<-BEMmat[BEMmat[,1]==1,]
    SCCmap[[i]]<-matchSCCtoBEM(SmatList[[i]],BEMx,referenceNetwork,removeNeg,PPInet=PPInet,minClust=minClust,simThreshold=simThreshold)
  }
  names(SCCmap)<-celllines
  return(SCCmap)
}
minBEM<-function(SCC,markerList,referenceNetwork,fullMarkerName,overDep=TRUE){
  #given a marker list find minimal marker for a group of genes in an SCC
  if(is.matrix(SCC)){
    vertices<-SCC[,1]
  }else{
    vertices<-SCC
  }
  markerList<-markerList[!is.na(markerList)]
  DepMarkerDist<-shortest.paths(referenceNetwork,v=vertices,markerList)
  minDist<-min(DepMarkerDist)
  colnames(DepMarkerDist)<-markerList
  minCheck<-DepMarkerDist==minDist
  possMarkers<-fullMarkerName[which(colSums(minCheck+0)>0)]
  if(length(possMarkers)>1&overDep){
    #work out overall distance btween potential biomarker and all dependencies in SCC
    #are we being circular?
  }
  return(possMarkers)
}

RefHclust<-function(DepData){
  distM<-dist(t(DepData))
  hclustRef<-hclust(distM)
  d2<-cophenetic(hclustRef)
  #hclustRef is dist-> dendrogram of dependencies versus cell lines. Distance between cell lines according to depdencies:
  return(d2)
}
compareToDendrogram<-function(RefCoph,InputData,penalty=NULL){
  #InputData is BEMvs Celllines
  d1<-dist(InputData,method="binary")
  #to get dist between cell lines according to values of biomarkers
  #get RefCoph from RefHclust function
  #alternative compare dendrograms not distance to dendrogram:
  #score<-cor(d1,RefCoph)
  bemHclust<-hclust(d1)
  d2<-cophenetic(bemHclust)
  score<-cor(d2,RefCoph)
  if(!is.null(penalty)){
    score<-score-penalty*ncol(InputData)
  }
  return(list(score=score,distFeature=d1))
}
plotFeat<-function(InputData){
  rdend<-hclust(dist(InputData,method="binary"))
  Heatmap(InputData,cluster_rows=rdend)
}
plotDepWithFeat<-function(InputData,DepData){
  #assume InputData has features selected by the greedy algorithm etc
  #basically finding a way to automate and identify novel features for plotting which is usually done manually by experimentalists
  #hence chance they will like/use this *fingers crossed*
  ha<-HeatmapAnnotation(df=as.data.frame(InputData),show_legend = FALSE,gap=unit(0.5,"points"),height=unit(0.1,"cm"))
  coldend<-hclust(dist(t(DepData)))
  Heatmap(as.matrix(DepData),cluster_columns=coldend,top_annotation = ha,show_heatmap_legend=FALSE)
}
splitDepByPPI<-function(Dependencies,PPIigraph,PPInet,minClust=10){

    #then also want to filter correlation matrix using prior graph:
    stringNames<-PPInet[Dependencies,"STRING_id"]
    names(stringNames)<-Dependencies
    stringNames<-stringNames[!is.na(stringNames)]

    PPIigraph<-delete.vertices(PPIigraph,which(!(V(PPIigraph)$name)%in%stringNames))
    #Amat<-as_adjacency_matrix(PPIigraph)+0

    SCCs<-components(PPIigraph,mode="strong")
    csizes<-SCCs$csize
    incClusters<-which(csizes>=minClust)
    temp<-sapply(incClusters,function(x) names(SCCs$membership[SCCs$membership==x]))
    if(length(temp)>0){
      DepNames<-lapply(temp,function(x) names(stringNames)[stringNames%in%x])
    }


  return(DepNames)

}

PriorRegGraph<-function(Smat,penaltymatrix,lambda=NULL){
  #need to set a regularization, or cross validate lambda parmeter
  if(is.null(lambda)){

  }else{
    Smat<-Smat-lambda*penaltymatrix
  }
  return(Smat)
}


