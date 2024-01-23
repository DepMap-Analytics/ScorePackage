Regression_createInputFeatures<-function(EM,
                                    oneFeatureOnly=NULL,
                                    selectedCellLines=NULL,MASTER_LIST,featFactorPopulationTh=3){



  rownames(EM)<-str_replace_all(rownames(EM),pattern = ':',replacement = '_')
  rownames(EM)<-str_replace_all(rownames(EM),pattern = ' ',replacement = '_')


  if(length(oneFeatureOnly)>0){

    for (kk in 1:length(oneFeatureOnly)){
      if (kk == 1){
        idxs<-which(str_detect(rownames(EM),oneFeatureOnly[kk]))
      }else{
        idxs<-union(idxs,which(str_detect(rownames(EM),oneFeatureOnly[kk])))
      }

    }

    if (length(idxs)==1){
      EM<-matrix(EM,1,ncol(EM),dimnames = list(rownames(EM)[idxs],colnames(EM)))
    }else{
      EM<-EM[idxs,]
    }
  }





  ##24.7.20 added in NA handling and second population size check
  #idxs<-which(rowSums(abs(EM))>ANOVA.settings.featFactorPopulationTh)
  idxs<-which(rowSums(!is.na(EM))>featFactorPopulationTh&rowSums(is.na(EM))<ncol(EM)-featFactorPopulationTh)
  if (length(idxs)==1){
    EM<-matrix(EM,1,ncol(EM),dimnames = list(rownames(EM)[idxs],colnames(EM)))
  }else{
    EM<-EM[idxs,]
  }

  if (length(selectedCellLines)){
    EM<-EM[,selectedCellLines]
  }



  return(list(EM=EM))
}
rfun<-function(values){
  rangev<-range(values,na.rm=T)
  return(rangev[2]-rangev[1])
}
AllRegression<-function(outputname,InputFC,EM,inparallel=FALSE,geneAnnotations,biomarkertype,Tissue=NULL,MSI=NULL,markerIsFactor=FALSE,pancan=FALSE,scalingFeat=FALSE){
  outputAll<-NULL
  nmarkers<-nrow(EM)
  rangeFC<-apply(InputFC,1,rfun)
  outputAll<-foreach(i=1:nmarkers,.combine=rbind)%dopar%{
    output<-NULL
    print(i)
    tempEM<-EM[i,colnames(InputFC)]
    inc<-colnames(EM)[!is.na(EM[i,])]
    tempEM<-tempEM[inc]
    if(!is.null(Tissue)|!is.null(MSI)){


        markers<-tempEM
        TissueFactor<-Tissue[inc]
        MSIfactor<-MSI[inc]
        InputData<-InputFC[,inc]
          ut<-unique(TissueFactor)
          for(k in 1:length(ut)){
            sc<-which(TissueFactor==ut[k])
            if(length(sc)<5){
              TissueFactor[sc]<-NA
            }
          }


          rmvna<-is.na(TissueFactor)|is.na(MSIfactor)
          if(sum(rmvna)>0&length(rmvna)>20){
            tempEM<-markers[!rmvna]
            TissueFactor<-TissueFactor[!rmvna]
            MSIfactor<-MSIfactor[!rmvna]
            InputData<-InputData[,!rmvna]
          }
          if(biomarkertype=="Binary"){
            check<-length(tempEM)>20&sum(tempEM)>4&sum(tempEM==0)>4
          }else{
            check<-length(tempEM)>20
            if(scalingFeat){
              tempEM<-ScorePackage::minmax(tempEM)
            }
          }
          if(check){
            if(sum(MSIfactor)>4&sum(MSIfactor==0)>4&length(unique(TissueFactor))>1){
              output<-limmaRegressPancan(InputData,tempEM,TissueFactor,MSIfactor,markerIsFactor)
            }else{

                if(length(unique(TissueFactor))>1){
                  output<-limmaRegressPancan(InputData,tempEM,TissueFactor,NULL,markerIsFactor)
                }
            }
          }
    }else{
      if(scalingFeat){
        tempEM<-ScorePackage::minmax(tempEM)
      }
        output<-limmaRegress(InputFC[,names(tempEM)],tempEM,markerIsFactor)
    }
    if(!is.null(output)){
      if(biomarkertype%in%c("Cont","CExpr")){
        output$FEATURE<-paste0(rownames(EM)[i],"_Expr")
      }
      if(biomarkertype=="CMet"){
        output$FEATURE<-paste0(rownames(EM)[i],"_Met")
      }
      if(biomarkertype=="CCN"){
        output$FEATURE<-paste0(rownames(EM)[i],"_CN")
      }
      if(biomarkertype=="CProt"){
        output$FEATURE<-paste0(rownames(EM)[i],"_Prot")
      }
      if(biomarkertype=="Binary"){
        output$FEATURE<-rownames(EM)[i]
      }
      output$`Depleted Gene`<-rownames(output)

      gA<-geneAnnotations[output$`Depleted Gene`,c('name',
                                                'gene_family',
                                                'location_sortable',
                                                'entrez_id',
                                                'ensembl_gene_id',
                                                'pubmed_id','string_id')]
      colnames(gA)<-c('Name',
                    'gene family',
                    'location',
                    'entrez id',
                    'ensembl gene id',
                    'pubmed id','string_id')
      if(markerIsFactor){
        output$BMstringID<-geneAnnotations[strsplit(rownames(EM)[i],"_",fixed=T)[[1]][1],"string_id"]
        GLASS_d<-Glass_Deltas(InputFC[,names(tempEM)[tempEM==1]],InputFC[,names(tempEM)[tempEM==0]])


        output$Range<-GLASS_d$g1[rownames(output)]
        output$IQR<-GLASS_d$g2[rownames(output)]
        output$FRange<-NA

      }else{
        output$BMstringID<-geneAnnotations[rownames(EM)[i],"string_id"]
        rangeBM<-range(EM[i,colnames(InputFC)],na.rm=T)
        rangeBM<-rangeBM[2]-rangeBM[1]
        IQRbm<-IQR(EM[i,colnames(InputFC)],na.rm=T)
        output$Range<-rangeBM
        output$IQR<-IQRbm
        if(pancan){
          inputvals<-na.omit(EM[i,colnames(InputFC)])

          q1<-quantile(inputvals,0.95)
          q2<-quantile(inputvals,0.05)
          M1<-median(inputvals[inputvals>=q1],na.rm=T)
          M2<-median(inputvals[inputvals<=q2],na.rm=T)
          output$FRange<-M1-M2

        }else{
          output$FRange<-NA

        }

      }
      output<-cbind(output,gA)
      rownames(output)<-NULL
      output
      #outputAll<-rbind(outputAll,output)
    }
  }
  assoc_id<-paste0("a:",1:nrow(outputAll))
  outputAll<-cbind(assoc_id,outputAll)
  #replace missing failing tests with p-value 1.
  if(sum(is.na(outputAll$P.Value))!=0){
    cat('warning: not all tests produced valid p-values. Set to 1 for multiple hypothesis correction')
  }
  outputAll[is.na(outputAll$P.Value),"P.Value"]<-1
  outputAll$FDR<-qvalue(outputAll$P.Value)$qvalues*100
  outputAll<-outputAll[order(outputAll$FDR,decreasing=FALSE),]
  return(outputAll)
}
Glass_Deltas<-function(x,y){
  if(is.matrix(x)){
    md<-abs(rowMeans(x)-rowMeans(y))
    sd1<-sqrt(rowVars(x))
    sd2<-sqrt(rowVars(y))
    g1<-md/sd1
    g2<-md/sd2
  }else{
    md <- abs(mean(x) - mean(y))
    g1 <- md/sd(x)
    g2 <- md/sd(y)
  }
  return(list(g1=g1,g2=g2))
}
limmaRegress<-function(Ymat,EM,markerIsFactor=FALSE){
  #add in scaling of the covariates here:

  if(markerIsFactor){
    designm<-model.matrix(~as.character(EM))
  }else{
    designm<-model.matrix(~as.vector(EM))
  }

  fit<-limma::lmFit(Ymat,designm)
  fit <- limma::eBayes(fit)
  sst<-rowSums((Ymat-rowMeans(Ymat))^2)
  ssr<-fit$df.residual*(fit$sigma^2)
  Rsquare<-1-(ssr/sst)
  tt<-topTable(fit,coef=2, number = Inf)
  output<-tt
  output$sigma<-fit$sigma
  output$sigmaB<-sqrt(fit$s2.post)
  output$Rsquare<-Rsquare[rownames(tt)]
  output$f2<-Rsquare[rownames(tt)]/(1-Rsquare[rownames(tt)])
  ordinary.t <- fit$coef[,2] / fit$stdev.unscaled[,2] / fit$sigma
  output$ord.t<-ordinary.t[rownames(output)]
  #colnames(output)[colnames(output)=="P.Value"]<-"PvalEB"
  ord.Pval<-2*pt( abs(ordinary.t), df=fit$df.residual, lower.tail=FALSE)
  output$OrdPValue<-ord.Pval[rownames(output)]
  output<-output[order(output$P.Value),]
  return(output)
}

limmaRegressPancan<-function(Ymat,EM,Tissue,MSI,markerIsFactor=FALSE){
  #want to add in scaling of the covariates and if applicable extract the MSI feature

if(!is.null(MSI)&!is.null(Tissue)){
  if(markerIsFactor){
    inputdata<-data.frame(omic=as.character(EM),Tissue=as.character(Tissue),MSI=as.character(MSI))
  }else{

    inputdata<-data.frame(omic=as.vector(EM),Tissue=as.character(Tissue),MSI=as.character(MSI))
  }

  designm<-model.matrix(~omic+Tissue+MSI,data=inputdata)
  designsub<-model.matrix(~Tissue+MSI,data=inputdata)
}else{
  if(markerIsFactor){
    inputdata<-data.frame(omic=as.character(EM),Tissue=as.character(Tissue))
  }else{

    inputdata<-data.frame(omic=as.vector(EM),Tissue=as.character(Tissue))
  }

  designm<-model.matrix(~omic+Tissue,data=inputdata)
  designsub<-model.matrix(~Tissue,data=inputdata)
}
  Ymat<-Ymat[,rowSums(is.na(inputdata))==0]
  fit<-limma::lmFit(Ymat,designm)
  fit <- limma::eBayes(fit)
  sst<-rowSums((Ymat-rowMeans(Ymat))^2)
  ssr<-fit$df.residual*(fit$sigma^2)
  Rsquare<-1-(ssr/sst)
  tt<-topTable(fit,coef=2, number = Inf)
  output<-tt
  output$Rsquare<-Rsquare[rownames(tt)]

  ordinary.t <- fit$coef[,2] / fit$stdev.unscaled[,2] / fit$sigma
  output$ord.t<-ordinary.t[rownames(output)]
  #colnames(output)[colnames(output)=="P.Value"]<-"PvalEB"
  ord.Pval<-2*pt( abs(ordinary.t), df=fit$df.residual, lower.tail=FALSE)
  output$OrdPValue<-ord.Pval[rownames(output)]

  fit2<-limma::lmFit(Ymat,designsub)
  fit2 <- limma::eBayes(fit2)
  ssr2<-fit2$df.residual*(fit2$sigma^2)
  Rsquare2<-1-(ssr2/sst)

  f2<-(Rsquare[rownames(tt)]-Rsquare2[rownames(tt)])/(1-Rsquare[rownames(tt)])
  output$f2<-f2


  output<-output[order(output$P.Value),]
  return(output)
}


DepDiffExpr<-function(BEM,geneexpression,gmtfile,OUTPUT_DIR,geneAnnotations,ncores=4){
  nmarkers<-nrow(BEM)
  output<-NULL
  geneIDs = geneAnnotations[match(rownames(geneexpression),geneAnnotations$symbol),"ensembl_gene_id"]
  tab<-vector("list",length=nmarkers)
  for(i in 1:nmarkers){
    print(i)
    alg<-rep(0,length(geneIDs))
    names(alg)<-geneIDs
      dep<-rownames(BEM)[i]
      clgroup<-as.factor(BEM[i,])
      names(clgroup)<-colnames(BEM)
      designmatrix<-model.matrix(~clgroup)
      fit<-limma::lmFit(geneexpression,designmatrix)
      fit <- limma::eBayes(fit)

      outputAll<-topTable(fit,coef=2,number=Inf,sort.by="P")

      outputAll$FDR<-qvalue(outputAll$P.Value)$qvalues*100
      outputAll<-outputAll[order(outputAll$FDR,decreasing=FALSE),]
      ttSig<-outputAll[outputAll$FDR<5,]

      if(nrow(ttSig)>20){
        logFC<-outputAll$logFC
        names(logFC)<-rownames(outputAll)
        genes<-tolower(names(sort(logFC,decreasing=TRUE)))
        gseares<-lapply(gmtfile,function(x) GSEAfunction(genes,tolower(x))$ESscore)

        hyperres<-lapply(gmtfile,function(x) testGmt(x,rownames(ttSig),rownames(geneexpression)))
        gseares<-sort(abs(unlist(gseares)),decreasing=TRUE)
        topgsea<-gseares[1:10]
        adjhyper<-qvalue(unlist(hyperres))$qvalues
        tophyperres<-adjhyper[adjhyper<0.1]
        if(length(tophyperres)>9){
          tophyperres<-sort(tophyperres)[1:10]
        }
        tophyperres<-tophyperres*100
        sigIDs<-geneAnnotations[match(rownames(ttSig),geneAnnotations$symbol),"ensembl_gene_id"]
        alg[sigIDs]<-1
        alg <- factor( alg )

        tgd <- new( "topGOdata", ontology="BP", allGenes = alg, nodeSize=10,
                    annot=annFUN.org, mapping="org.Hs.eg.db", ID = "ensembl" )

        ## run tests
        resultTopGO.elim <- runTest(tgd, algorithm = "elim", statistic = "Fisher" )
        resultTopGO.classic <- runTest(tgd, algorithm = "classic", statistic = "Fisher" )

        ## look at results
        temp<- GenTable( tgd, Fisher.elim = resultTopGO.elim,
                         Fisher.classic = resultTopGO.classic,
                         orderBy = "Fisher.classic" , topNodes = 50)
        temp$DepGene<-dep

        tab[[i]]<-temp

        temp<-data.frame(gsea=paste0(topgsea,collapse="//"),gname=paste0(names(topgsea),collapse="//"),hyper=paste0(tophyperres,collapse="//"),hypername=paste0(names(tophyperres),collapse = "//"),dep=dep,DEgenes=paste0(rownames(ttSig),collapse="//"),stringsAsFactors = FALSE)

      }else{
        temp<-data.frame(gsea="NoDE",gname="NoDE",hyper="NoDE",hypername="NoDE",dep=dep,DEgenes="NoDE",stringsAsFactors = FALSE)
      }

      output<-rbind(output,temp)
    }
  topGOResults <- rbind.fill(tab)
  return(list(fullres=output,topGOres=topGOResults))
}





