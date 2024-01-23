
ANOVA_cohens_d <- function(x, y) {
    lx <- length(x)- 1
    ly <- length(y)- 1

    md  <- abs(mean(x) - mean(y))        ## mean difference (numerator)
    csd <- lx * var(x) + ly * var(y)
    csd <- csd/(lx + ly)
    csd <- sqrt(csd)                     ## common sd computation

    cd  <- md/csd                        ## cohen's d

    return(cd)
}
ANOVA_glass_Ds<-function(x,y){

    md<-abs(mean(x)-mean(y))
    g1<-md/sd(x)
    g2<-md/sd(y)

    return(list(g1=g1,g2=g2))
}

ANOVA_individualANOVA<-function(DEP_GENE,
                                       FEATURE,
                                       display=TRUE,
                                       printTOfig=FALSE,
                                       PATH='',
                                       FN='',
                                       FDR=NA,
                                       OUTPUT_PATH='',
                                       ID='',
                                       col1=GRAY,
                                       col2=BLUE,InputFeatures,ESSprofiles,geneAnnotations,
                                ANOVA.settings.analysisType='CS',ANOVA.fakeCS=FALSE,
                                ANOVA.settings.ESSENTIALITY_domain="logFCs",ANOVA.settings.featFactorPopulationTh=ANOVA.settings.featFactorPopulationTh,
                                ANOVA.settings.MSIfactorPopulationTh=3,ANOVA.settings.includeMSI_Factor=FALSE){


  if (printTOfig){
    #png(paste(PATH,FN,'_00_Set.png',sep=''),width=793.92,height=1122.24)
    png(paste(PATH,FN,'_00_scatterSet.png',sep=''),width=400,height=450)
      #pdf(paste(PATH,FN,'_00_scatterSet.png',sep=''),width=400,height=450)
  }
  ##set up error handling to remove NAs and then do tests. So this will allow for missing data in MoBEM
  ##or when using patient subtypes e.g. intClust to compare two groups where annotation is NA for all cell lines not in either of the two groups
  commonC<-intersect(rownames(ESSprofiles),colnames(InputFeatures$BEM))
  commonC<-commonC[!is.na(ESSprofiles[commonC,DEP_GENE])]

  TISSUE_FACTOR<-InputFeatures$TISSUES[commonC]
  MSI_FACTOR<-InputFeatures$MSI_VARIABLE[commonC]

  Epattern<-ESSprofiles[commonC,DEP_GENE]

  if(ANOVA.settings.ESSENTIALITY_domain=='Mageck'){
      Epattern<- -log10(Epattern)
  }

  FEATpattern<-InputFeatures$BEM[FEATURE,commonC]

  FEATpattern<-FEATpattern[!is.na(Epattern)]
  Epattern<-Epattern[!is.na(Epattern)]
  #add in checks for NA's in features (incomplete BEMS) 24.7.20:
  Epattern<-Epattern[!is.na(FEATpattern)]
  FEATpattern<-FEATpattern[!is.na(FEATpattern)]

  FEATpattern[which(FEATpattern==0)]<-'neg'
  FEATpattern[which(FEATpattern=='1')]<-'pos'
  FEATpattern<-as.factor(FEATpattern)

  TISSUEpattern<-as.factor(TISSUE_FACTOR[names(Epattern)])
  MSIpattern<-as.factor(MSI_FACTOR[names(Epattern)])

  Y <- Epattern

  featFactorPopulationTh<-ANOVA.settings.featFactorPopulationTh
  MSIfactorPopulationTh<-ANOVA.settings.MSIfactorPopulationTh

  A<-(ANOVA.settings.includeMSI_Factor & length(which(FEATpattern=='pos'))>=featFactorPopulationTh &
        length(which(FEATpattern=='neg'))>=featFactorPopulationTh &
        length(which(MSIpattern==0))>=MSIfactorPopulationTh &
        length(which(MSIpattern==1))>=MSIfactorPopulationTh)

  B<-(!ANOVA.settings.includeMSI_Factor & length(which(FEATpattern=='pos'))>=featFactorPopulationTh &
        length(which(FEATpattern=='neg'))>=featFactorPopulationTh)


  if (A | B){

    if(ANOVA.settings.analysisType=='PANCAN' & !ANOVA.fakeCS){
            fit <- aov(Y ~ TISSUEpattern+MSIpattern+FEATpattern)
    }else{
      if(ANOVA.settings.includeMSI_Factor){
        fit <- aov(Y ~ MSIpattern+FEATpattern)
      }else{
        fit <- aov(Y ~ FEATpattern)
      }
    }
    tfit<-t.test(Epattern~FEATpattern)
    tfitEV<-t.test(Epattern~FEATpattern,var.equal=TRUE)
    Y<-anova(fit)


    if(ANOVA.settings.analysisType=='PANCAN' & !ANOVA.fakeCS){
      if(ANOVA.settings.includeMEDIA_factor){
          FEATURE_PVAL<-Y[4,5]
          MSI_PVAL<-Y[3,5]
          MEDIA_PVAL<-Y[2,5]
          tissue_PVAL<-Y[1,5]
      }else{
          FEATURE_PVAL<-Y[3,5]
          MSI_PVAL<-Y[2,5]
          MEDIA_PVAL<-NA
          tissue_PVAL<-Y[1,5]
      }
    }else{
      if(ANOVA.settings.includeMSI_Factor){
        FEATURE_PVAL<-Y[2,5]
        MSI_PVAL<-Y[1,5]
      }else{
        FEATURE_PVAL<-Y[1,5]
        MSI_PVAL<-NA
      }
      tissue_PVAL<-NA
      MEDIA_PVAL<-NA
    }

     FEATURE_ESS_WTT_pvalue<-tfit$p.value
     FEATURE_pvalT_EV<-tfitEV$p.value
     pos_ESS_MEAN<-mean(Epattern[FEATpattern=='pos'])
     neg_ESS_MEAN<-mean(Epattern[FEATpattern=='neg'])

     deltaMEAN_ESS<-pos_ESS_MEAN-neg_ESS_MEAN

     pos_ESS_sd<-sd(Epattern[FEATpattern=='pos'])
     neg_ESS_sd<-sd(Epattern[FEATpattern=='neg'])

     Npos<-length(which(FEATpattern=='pos'))
     Nneg<-length(which(FEATpattern=='neg'))

     EFFECTSIZE_ESS<-ANOVA_cohens_d(Epattern[FEATpattern=='pos'],Epattern[FEATpattern=='neg'])
     GLASS_d<-ANOVA_glass_Ds(Epattern[FEATpattern=='pos'],Epattern[FEATpattern=='neg'])

  #   inHumanKinome<-is.element(DEP_GENE,humKinome)+0
  #   uniProtKinase<-is.element(DEP_GENE,kinase_genes_and_syn)+0
  #   GSKepiTarget<-is.element(DEP_GENE,GSK_epiTargets)+0

     if (display){
         ANOVA_scatterSets(Epattern = Epattern,FEATpattern = FEATpattern,adaptive_shape = !ANOVA.settings.SCALEtoLIMIT,
                                  DEP_GENE = DEP_GENE,FEATURE = FEATURE,p = FEATURE_PVAL,ID=ID)
     }

     RES<-matrix(c(FEATURE,
                   DEP_GENE,
                   geneAnnotations[DEP_GENE,c('name',
                                      'gene_family',
                                      'location_sortable',
                                      'entrez_id',
                                      'ensembl_gene_id',
                                      'pubmed_id','string_id')],
                  Npos,
                  Nneg,
                  pos_ESS_MEAN,
                  neg_ESS_MEAN,
                  deltaMEAN_ESS,
                  pos_ESS_sd,
                  neg_ESS_sd,
                  EFFECTSIZE_ESS,
                  GLASS_d$g1,
                  GLASS_d$g2,
                  FEATURE_PVAL,
                  tissue_PVAL,
                  MEDIA_PVAL,
                  MSI_PVAL,
                  FEATURE_ESS_WTT_pvalue,FEATURE_pvalT_EV),1,25,
                dimnames=list('1',c('FEATURE',
                                    'Depleted Gene',
                                    'Name',
                                    'gene family',
                                    'location',
                                    'entrez id',
                                    'ensembl gene id',
                                    'pumbed id',
                                    'string_id',
                                    'N FEATURE pos',
                                    'N FEATURE_neg',
                                    'FEATUREpos_ESS_MEAN',
                                    'FEATUREneg_ESS_MEAN',
                                    'FEATURE_deltaMEAN_ESS',
                                    'FEATUREpos_ESS_sd',
                                    'FEATUREneg_ESS_sd',
                                    'FEATURE_ESS_effect_size',
                                    'FEATUREpos_Glass_delta',
                                    'FEATUREneg_Glass_delta',
                                    'FEATURE_ANOVA_pval',
                                    'Tissue_ANOVA_pval',
                                    'MEDIA_ANOVA_pval',
                                    'MSI_ANOVA_pval',
                                    'FEATURE_ESS_T_pval','FEATURE_T_pvalEV')))


     if (printTOfig){dev.off()}
  }else{
      RES<-NULL
  }

   return(RES)
}

ANOVA_singleDepGeneANOVA<-function(DEP_GENE,sorted=FALSE,verbose=TRUE,InputFeatures,ESSprofiles,geneAnnotations,
                                   ANOVA.settings.featFactorPopulationTh,ANOVA.settings.analysisType,ANOVA.fakeCS,
                                   ANOVA.settings.includeMSI_Factor){

    FEATURES<-names(which(rowSums(InputFeatures$BEM,na.rm=T)>=3))
    nFEATURES<-length(FEATURES)

    if(verbose){
        print(paste('Running Single Depleted Gene ANOVA: ',
                    DEP_GENE,sep=''))
        pb <- txtProgressBar(min=1,max=nFEATURES,style=3)
    }
   if(nFEATURES>0){
    for (i in 1:nFEATURES){

        if(verbose){
            setTxtProgressBar(pb, i)
        }

        RES<-ANOVA_individualANOVA(DEP_GENE=DEP_GENE,FEATURE=FEATURES[i],display = FALSE,InputFeatures=InputFeatures,
                                   ESSprofiles = ESSprofiles,geneAnnotations = geneAnnotations,ANOVA.settings.featFactorPopulationTh=ANOVA.settings.featFactorPopulationTh,
                                   ANOVA.settings.analysisType=ANOVA.settings.analysisType,ANOVA.fakeCS=ANOVA.fakeCS,
                                   ANOVA.settings.includeMSI_Factor=ANOVA.settings.includeMSI_Factor)

        if(i ==1){
            if(!is.null(RES)){
                TOTRES<-RES
            }
        }else{
            if(exists('TOTRES')){
                if(!is.null(RES)){
                    TOTRES<-rbind(TOTRES,RES)
                }
            }else{
                if(!is.null(RES)){
                    TOTRES<-RES
                }
            }
        }
    }
    if(exists('TOTRES')){
    idxs<-which(!is.na(as.numeric(TOTRES[,"FEATURE_ANOVA_pval"])))

    TOTRES<-TOTRES[idxs,]

      if (length(idxs)==1){
        TOTRES<-matrix(TOTRES,1,length(TOTRES),dimnames=list(1,names(TOTRES)))
      }

      if(verbose){
        Sys.sleep(1)
        close(pb)
      }

      if(sorted){
        TOTRES<-TOTRES[order(as.numeric(TOTRES[,"FEATURE_ANOVA_pval"])),]
      }



    }else{TOTRES=NULL}
   }else{TOTRES=NULL}
    if(!is.null(TOTRES)){
      GeneFeatures<-unlist(TOTRES[,"FEATURE"])

      BMstringID<-unlist(sapply(GeneFeatures,function(x) getStringID(x,geneAnnotations)))
      TOTRES<-cbind(TOTRES,BMstringID)
    }
  return(TOTRES)
}

ANOVA_singleDepGeneANOVACellector<-function(DEP_GENE,sorted=FALSE,verbose=FALSE,InputFeatures,ESSprofiles){

  FEATURES<-rownames(InputFeatures$cBEM)
  nFEATURES<-length(FEATURES)

  if(verbose){
    print(paste('Running Single Depleted Gene ANOVA: ',
                DEP_GENE,sep=''))
    pb <- txtProgressBar(min=1,max=nFEATURES,style=3)
  }

  if(nFEATURES>0){
    for (i in 1:nFEATURES){

      if(verbose){
        setTxtProgressBar(pb, i)
      }

      RES<-ANOVA_individualCellector(DEP_GENE=DEP_GENE,FEATURE=FEATURES[i],display = FALSE,InputFeatures=InputFeatures,ESSprofiles = ESSprofiles,featurenumber=i)
      print(RES)
      if(i ==1){
        if(!is.null(RES)){
          TOTRES<-RES
        }
      }else{
        if(exists('TOTRES')){
          if(!is.null(RES)){
            TOTRES<-rbind(TOTRES,RES)
          }
        }else{
          if(!is.null(RES)){
            TOTRES<-RES
          }
        }
      }
    }
    if(exists('TOTRES')){
      idxs<-which(!is.na(as.numeric(TOTRES[,"LRT_pval"])))

      TOTRES<-TOTRES[idxs,]

      if (length(idxs)==1){
        TOTRES<-matrix(TOTRES,1,length(TOTRES),dimnames=list(1,names(TOTRES)))
      }

      if(verbose){
        Sys.sleep(1)
        close(pb)
      }

      if(sorted){
        TOTRES<-TOTRES[order(as.numeric(TOTRES[,"LRT_pval"])),]
      }
    }else{TOTRES=NULL}
  }else{TOTRES=NULL}
  return(TOTRES)
}
ANOVA_individualCellector<-function(DEP_GENE,
                                FEATURE,
                                display=TRUE,
                                printTOfig=FALSE,
                                PATH='',
                                FN='',
                                FDR=NA,
                                OUTPUT_PATH='',
                                ID='',
                                col1=GRAY,
                                col2=BLUE,InputFeatures,ESSprofiles,featurenumber){


  commonC<-intersect(rownames(ESSprofiles),colnames(InputFeatures$cBEM))

  commonC<-commonC[!is.na(ESSprofiles[commonC,DEP_GENE])]
  ###

  #get combined single and cellector BEM entries to pass to lm as data option
  FitTerms<-InputFeatures$FitTerms[[featurenumber]]

CTerm<-FitTerms$fitC
FEATURE<-CTerm

cEntry<-InputFeatures$cBEM[CTerm,]
sTerms<-InputFeatures$FitTerms[[featurenumber]]$fitT
sTerms<-unlist(strsplit(sTerms,"+",fixed=TRUE))
sEntries<-InputFeatures$sBEM[sTerms,]
allBEM<-rbind(cEntry,sEntries)
allBEM<-allBEM[,commonC]
allBEM<-t(allBEM)

colnames(allBEM)<-c(CTerm,sTerms)
allBEMm<-allBEM
allBEM[allBEM=="0"]<-'neg'
allBEM[allBEM=="1"]<-'pos'
for(i in 1:ncol(allBEM)){

  allBEM[,i]<-as.factor(allBEM[,i])
}
sBEM<-t(sEntries)
sBEM[sBEM=="0"]<-'neg'
sBEM[sBEM=="1"]<-'pos'
for(i in 1:ncol(sBEM)){

  sBEM[,i]<-as.factor(sBEM[,i])
}

Aterms<-InputFeatures$FitTerms[[featurenumber]]$fitAll
Tterms<-InputFeatures$FitTerms[[featurenumber]]$fitT

design2<-model.matrix(~.,as.data.frame(allBEM))
design1<-model.matrix(~.,as.data.frame(sBEM))

  Epattern<-ESSprofiles[commonC,DEP_GENE]

  fit2<-lm(Epattern~design2+0,na.action=na.omit)
  fit1<-lm(Epattern~design1+0,na.action=na.omit)
  #require LRT test between fit2 and fit1 to see if the AND compound term significantly increases the fit to the dependency profile
  lmtest<-data.frame(lrtest(fit2,fit1))
  Y1<-summary(fit1)
  Y2<-summary(fit2)
  ANDpvalue<-lmtest[2,5]
  print(Y1$coefficients)
  print(Y2$coefficients)


  FEATURE_Name_PVAL<-paste(rownames(Y1$coefficients[2:nrow(Y1$coefficients),]),collapse = "\\")
  FEATURE_Single_PVAL<-paste(Y1$coefficients[2:nrow(Y1$coefficients),4],collapse="\\")
  FEATURE_Cellector_PVAL<-Y2$coefficients[nrow(Y2$coefficients),4]

  #pos_F1_ESS_MEAN is any cell line positive in any of the single terms
  pos_F1_ESS_MEAN<-mean(Epattern[colSums(sEntries)!=0])
  #wt across all terms:
  WT_ESS_MEAN<-mean(Epattern[rowSums(allBEMm)==0])
  #mean of all cell lines positive for compound pattern:
  Compound_pos_F2_ESS_MEAN<-mean(Epattern[cEntry==1])


  pos_ESS_sd<-sd(Epattern[colSums(sEntries)!=0])
  neg_ESS_sd<-sd(Epattern[rowSums(allBEMm)==0])

  Npos<-length(which(colSums(sEntries)!=0))
  Nneg<-length(which(rowSums(allBEMm)==0))

  EFFECTSIZE_ESS<-ANOVA_cohens_d(Epattern[colSums(sEntries)!=0],Epattern[rowSums(allBEMm)==0])
  EFFECTSIZE_ESSC<-ANOVA_cohens_d(Epattern[colSums(sEntries)!=0&cEntry==0],Epattern[cEntry==1])
  GLASS_dmp<-ANOVA_glass_Ds(Epattern[colSums(sEntries)!=0&rowSums(allBEMm)==0],Epattern[cEntry==1])
  GLASS_dmn<-ANOVA_glass_Ds(Epattern[cEntry!=0],Epattern[cEntry==0])

   FEATURE_child<-NA
  RES<-matrix(c(FEATURE,FEATURE_child,
                DEP_GENE,
                geneAnnotations[DEP_GENE,c('name',
                                           'gene_family',
                                           'location_sortable',
                                           'entrez_id',
                                           'ensembl_gene_id',
                                           'pubmed_id','string_id')],
                Npos,
                Nneg,
                pos_F1_ESS_MEAN,
                WT_ESS_MEAN,
                Compound_pos_F2_ESS_MEAN,
                pos_ESS_sd,
                neg_ESS_sd,
                EFFECTSIZE_ESS,
                EFFECTSIZE_ESSC,
                GLASS_dmp$g1,
                GLASS_dmp$g2,
                GLASS_dmn$g1,
                GLASS_dmn$g2,
                FEATURE_Name_PVAL,
                FEATURE_Single_PVAL,
                FEATURE_Cellector_PVAL,
                ANDpvalue

  ),1,27,
  dimnames=list('1',c('FEATURE_parent','FEATURE_child',
                      'Depleted Gene',
                      'Name',
                      'gene family',
                      'location',
                      'entrez id',
                      'ensembl gene id',
                      'pumbed id',
                      'string_id',
                      'N FEATURE pos',
                      'N FEATURE_neg',
                      'FEATUREparent_ESS_MEAN',
                      'FEATUREwt_ESS_MEAN',
                      'FEATUREcompound_ESS_MEAN',
                      'FEATUREpos_ESS_sd',
                      'FEATUREneg_ESS_sd',
                      'FEATURE_Parent_ESS_effect_size',
                      'FEATURE_Compound_effect_size',
                      'FEATUREparent_pos_Glass_delta',
                      'FEATUREparent_neg_Glass_delta',
                      'FEATUREcompound_pos_Glass_delta',
                      'FEATUREcompound_neg_Glass_delta',
                      'FEATURE_parent_pval',
                      'FEATURE_child_pval',
                      'FEATURE_compound_pval',
                      'LRT_pval'
                      )))




  return(RES)
}
ANOVA_totalANOVA<-function(fn,ESSprofiles,InputFeatures,inparallel,geneAnnotations,ANOVA.settings.featFactorPopulationTh=3,
                           gdscANOVA.settings.pval_correction_method="qvalue",ANOVA.settings.analysisType,ANOVA.fakeCS,
                           ANOVA.settings.includeMSI_Factor){

    #check names present for all BEM features:
    nameCheck<-sum(rownames(InputFeatures$BEM)=="")
    nameCheck2<-sum(is.na(rownames(InputFeatures$BEM)))
    if(nameCheck>0|nameCheck2){
      stop("Missing or NA names in Input Feature BEM, check inputs")
    }
    DEP_GENES<-colnames(ESSprofiles)

    nDEP_GENES<-length(DEP_GENES)

    print('+ Running ANOVA')

    pb <- txtProgressBar(min=1,max=nDEP_GENES,style=3)

    flag<-1
    if(inparallel){
      TOTRES<-foreach(i=1:nDEP_GENES,.combine=rbind)%dopar%{
        print(c(DEP_GENES[i],i,nDEP_GENES))
        setTxtProgressBar(pb, i)
        currentRES<-ANOVA_singleDepGeneANOVA(DEP_GENE = DEP_GENES[i],verbose=FALSE,InputFeatures=InputFeatures,ESSprofiles = ESSprofiles,
                                             geneAnnotation=geneAnnotations,ANOVA.settings.featFactorPopulationTh=ANOVA.settings.featFactorPopulationTh,
                                             ANOVA.settings.analysisType=ANOVA.settings.analysisType,ANOVA.fakeCS=ANOVA.fakeCS,
                                             ANOVA.settings.includeMSI_Factor=ANOVA.settings.includeMSI_Factor)

        #if (flag == 1 & ncol(currentRES)>0){
        ##if(flag==1&!is.null(currentRES)){

        ##    TOTRES<-currentRES
        ##    flag<-0
        ##} else{

        ##    if (!is.null(currentRES)){

        ##        TOTRES<-rbind(TOTRES,currentRES)
        ##    }
        ##}
        #if (i%%50==0){
        #    write.table(TOTRES,quote=FALSE,row.names=F,sep='\t',file=paste(fn,'.txt',sep=''))
        #}
      }
    }else{
      for(i in 1:nDEP_GENES){
        print(c(DEP_GENES[i],i,nDEP_GENES))
        setTxtProgressBar(pb, i)
        currentRES<-ANOVA_singleDepGeneANOVA(DEP_GENE = DEP_GENES[i],verbose=FALSE,InputFeatures=InputFeatures,
                                             ESSprofiles = ESSprofiles,geneAnnotation=geneAnnotations,
                                             ANOVA.settings.featFactorPopulationTh=ANOVA.settings.featFactorPopulationTh,
                                             ANOVA.settings.analysisType=ANOVA.settings.analysisType,ANOVA.fakeCS=ANOVA.fakeCS,
                                             ANOVA.settings.includeMSI_Factor=ANOVA.settings.includeMSI_Factor)

        #if (flag == 1 & ncol(currentRES)>0){
        if(flag==1&!is.null(currentRES)){

            TOTRES<-currentRES
            flag<-0
        } else{

            if (!is.null(currentRES)){

                TOTRES<-rbind(TOTRES,currentRES)
            }
        }
        #if (i%%50==0){
        #    write.table(TOTRES,quote=FALSE,row.names=F,sep='\t',file=paste(fn,'.txt',sep=''))
        #}
      }
    }
    print('tested all genes in ANOVA_totalANOVA')

    if(exists('TOTRES')){
      print(colnames(TOTRES))
      print(dim(TOTRES))
      TOTRES<-TOTRES[order(as.numeric(TOTRES[,"FEATURE_ANOVA_pval"])),]

      if(gdscANOVA.settings.pval_correction_method!='qvalue'){
        FDR<-p.adjust(as.numeric(TOTRES[,"FEATURE_ANOVA_pval"]),method='fdr')
        FDRt<-p.adjust(as.numeric(TOTRES[,"FEATURE_ESS_T_pval"]),method='fdr')
        FDRtEV<-p.adjust(as.numeric(TOTRES[,"FEATURE_T_pvalEV"]),method='fdr')
      }else{
        pvals<-as.numeric(TOTRES[,"FEATURE_ANOVA_pval"])
        print(paste("number p-vals:",length(pvals)))
        print(paste("max p-val:",max(pvals)))
        Q<-qvalue(as.numeric(TOTRES[,"FEATURE_ANOVA_pval"]))
        FDR<-Q$qvalue
        #cadded to the output 28.1.21
        Q<-qvalue(as.numeric(TOTRES[,"FEATURE_ESS_T_pval"]))
        FDRt<-Q$qvalue
        Q<-qvalue(as.numeric(TOTRES[,"FEATURE_T_pvalEV"]))
        FDRtEV<-Q$qvalue
      }

      TOTRES<-cbind(TOTRES,FDR*100)
      colnames(TOTRES)[ncol(TOTRES)]<-'ANOVA FEATURE FDR %'
      TOTRES<-cbind(TOTRES,FDRt*100)
      colnames(TOTRES)[ncol(TOTRES)]<-'t-test (unequal var) FDR %'
      #TOTRES<-cbind(TOTRES,FDRtEV*100)
      #colnames(TOTRES)[ncol(TOTRES)]<-'t-test FDR Equal Var %'
      assoc_id<-paste('a',1:nrow(TOTRES),sep='')

      TOTRES<-cbind(assoc_id,TOTRES)

      colnames(TOTRES)[1]<-'assoc_id'

      write.table(TOTRES,quote=FALSE,row.names=F,sep='\t',file=paste(fn,'.txt',sep=''))

      save(TOTRES,file=paste(fn,'.rdata',sep=''))

      if(length(range)>1){
        Sys.sleep(1)
        close(pb)
      }
    }else{
      TOTRES<-NULL
    }
    return(TOTRES)
}

ANOVA_totalANOVACellector<-function(fn,ESSprofiles,InputFeatures,inparallel){

  DEP_GENES<-colnames(ESSprofiles)

  nDEP_GENES<-length(DEP_GENES)

  print('+ Running ANOVA')

  pb <- txtProgressBar(min=1,max=nDEP_GENES,style=3)

  flag<-1
  if(inparallel){
    TOTRES<-foreach(i=1:nDEP_GENES,.combine=rbind)%dopar%{
      print(c(DEP_GENES[i],i,nDEP_GENES))
      setTxtProgressBar(pb, i)
      currentRES<-ANOVA_singleDepGeneANOVACellector(DEP_GENE = DEP_GENES[i],verbose=FALSE,InputFeatures=InputFeatures,ESSprofiles = ESSprofiles)

      currentRES

    }
  }else{
    for(i in 1:nDEP_GENES){
      print(c(DEP_GENES[i],i,nDEP_GENES))
      setTxtProgressBar(pb, i)
      currentRES<-ANOVA_singleDepGeneANOVACellector(DEP_GENE = DEP_GENES[i],verbose=FALSE,InputFeatures=InputFeatures,ESSprofiles = ESSprofiles)


      if(flag==1&!is.null(currentRES)){

        TOTRES<-currentRES
        flag<-0
      } else{

        if (!is.null(currentRES)){

          TOTRES<-rbind(TOTRES,currentRES)
        }
      }

    }
  }
  print('tested all genes in ANOVA_totalANOVA')


  if(exists('TOTRES')){
    print(colnames(TOTRES))
    print(dim(TOTRES))
    TOTRES<-TOTRES[order(as.numeric(TOTRES[,"LRT_pval"])),]

    if(gdscANOVA.settings.pval_correction_method!='qvalue'){
      FDR<-p.adjust(as.numeric(TOTRES[,"LRT_pval"]),method='fdr')
      FDRt<-p.adjust(as.numeric(TOTRES[,"FEATURE_cellector_pval"]),method='fdr')
    }else{
      pvals<-as.numeric(TOTRES[,"LRT_pval"])
      print(paste("number p-vals:",length(pvals)))
      print(paste("max p-val:",max(pvals)))
      Q<-qvalue(as.numeric(TOTRES[,"FEATURE_cellector_pval"]))
      FDR<-Q$qvalue
      Q2<-qvalue(as.numeric(TOTRES[,"LRT_pval"]))
      FDR2<-Q2$qvalue
      #commented out 15.6.19 as not included in output below
      #Q<-qvalue(as.numeric(TOTRES[,"FEATURE_ESS_T_pval"]))
      #FDRt<-Q$qvalue
    }

    TOTRES<-cbind(TOTRES,FDR*100,FDR2*100)
    colnames(TOTRES)[(ncol(TOTRES)-1):ncol(TOTRES)]<-c('ANOVA FEATURE FDR %','LRT FEATURE FDR %')

    assoc_id<-paste('a',1:nrow(TOTRES),sep='')

    TOTRES<-cbind(assoc_id,TOTRES)

    colnames(TOTRES)[1]<-'assoc_id'

    write.table(TOTRES,quote=FALSE,row.names=F,sep='\t',file=paste(fn,'.txt',sep=''))

    save(TOTRES,file=paste(fn,'.rdata',sep=''))

    if(length(range)>1){
      Sys.sleep(1)
      close(pb)
    }
  }else{
    TOTRES<-NULL
  }
  return(TOTRES)
}

DepGeneChunk<-function(GeneList,InputFeatures,ESSprofiles,OUTPUT_DIR,geneAnnotation){
  nDEP_GENES<-length(GeneList)
  pb <- txtProgressBar(min=1,max=nDEP_GENES,style=3)
  flag<-1
  for(i in 1:nDEP_GENES){
    print(c(GeneList[i],i,nDEP_GENES))
    setTxtProgressBar(pb, i)
    currentRES<-ANOVA_singleDepGeneInteraction(DEP_GENE = GeneList[i],verbose=FALSE,InputFeatures=InputFeatures,ESSprofiles = ESSprofiles,OUTPUT_DIR=OUTPUT_DIR,geneAnnotation=geneAnnotation)

    if(flag==1&!is.null(currentRES)){

      TOTRES<-currentRES
      flag<-0
    } else{

      if (!is.null(currentRES)){

        TOTRES<-rbind(TOTRES,currentRES)
      }
    }

  }
  return(TOTRES)
}

DepGeneChunkCompound<-function(GeneList,InputFeatures,ESSprofiles,OUTPUT_DIR,MutualExclusive,geneAnnotation,
                               ANOVA.settings.analysisType,ANOVA.fakeCS,
                               ANOVA.settings.includeMSI_Factor){
  nDEP_GENES<-length(GeneList)
  pb <- txtProgressBar(min=1,max=nDEP_GENES,style=3)
  flag<-1
  for(i in 1:nDEP_GENES){
    print(c(GeneList[i],i,nDEP_GENES))
    setTxtProgressBar(pb, i)
    currentRES<-ANOVA_singleDepGeneCompound(DEP_GENE = GeneList[i],verbose=FALSE,InputFeatures=InputFeatures,ESSprofiles = ESSprofiles,OUTPUT_DIR=OUTPUT_DIR,MutualExclusive=MutualExclusive,
                                            geneAnnotation=geneAnnotation,ANOVA.settings.analysisType=ANOVA.settings.analysisType,ANOVA.fakeCS=ANOVA.fakeCS,
                                            ANOVA.settings.includeMSI_Factor=ANOVA.settings.includeMSI_Factor)

    if(flag==1&!is.null(currentRES)){

      TOTRES<-currentRES
      flag<-0
    } else{

      if (!is.null(currentRES)){

        TOTRES<-rbind(TOTRES,currentRES)
      }
    }

  }
  return(TOTRES)
}

chunk <- function(x,n) split(x, factor(sort(rank(x)%%n)))
ANOVA_totalANOVAInteraction<-function(fn,ESSprofiles,InputFeatures,inparallel=TRUE,ncores=24,OUTPUT_DIR,geneAnnotation){

  DEP_GENES<-colnames(ESSprofiles)

  nDEP_GENES<-length(DEP_GENES)

  print('+ Running ANOVA')

  #pb <- txtProgressBar(min=1,max=nDEP_GENES,style=3)

  flag<-1
  idxs<-chunk(1:nDEP_GENES,ncores)

  geneChunks<-lapply(idxs,function(x) DEP_GENES[x])
  print(head(geneChunks))
  print(inparallel)
  if(inparallel){
    TOTRES<-foreach(i=1:ncores,.combine=rbind)%dopar%{
      #print(c(DEP_GENES[i],i,nDEP_GENES))

      #setTxtProgressBar(pb, i)
      #currentRES<-ANOVA_singleDepGeneInteraction(DEP_GENE = DEP_GENES[i],verbose=FALSE,InputFeatures=InputFeatures,ESSprofiles = ESSprofiles)

      DepGeneChunk(GeneList=geneChunks[[i]],InputFeatures,ESSprofiles,OUTPUT_DIR=OUTPUT_DIR)
    }
  }else{
    for(i in 1:nDEP_GENES){
      print(c(DEP_GENES[i],i,nDEP_GENES))
      setTxtProgressBar(pb, i)
      currentRES<-ANOVA_singleDepGeneInteraction(DEP_GENE = DEP_GENES[i],verbose=FALSE,InputFeatures=InputFeatures,ESSprofiles = ESSprofiles,OUTPUT_DIR=OUTPUT_DIR)

      if(flag==1&!is.null(currentRES)){

        TOTRES<-currentRES
        flag<-0
      } else{

        if (!is.null(currentRES)){

          TOTRES<-rbind(TOTRES,currentRES)
        }
      }

    }
  }
  print('tested all genes in ANOVA_totalANOVA')

  if(exists('TOTRES')){
    print(colnames(TOTRES))
    print(dim(TOTRES))
    TOTRES<-TOTRES[order(as.numeric(TOTRES[,"FEATURE_INT_pval"])),]

    if(gdscANOVA.settings.pval_correction_method!='qvalue'){
      FDR<-p.adjust(as.numeric(TOTRES[,"FEATURE_INT_pval"]),method='fdr')
      #FDRt<-p.adjust(as.numeric(TOTRES[,"FEATURE_ESS_T_pval"]),method='fdr')
    }else{
      pvals<-as.numeric(TOTRES[,"FEATURE_INT_pval"])
      print(paste("number p-vals:",length(pvals)))
      print(paste("max p-val:",max(pvals)))
      Q<-qvalue(as.numeric(TOTRES[,"FEATURE_INT_pval"]))
      FDR<-Q$qvalue
      #commented out 15.6.19 as not included in output below
      #Q<-qvalue(as.numeric(TOTRES[,"FEATURE_ESS_T_pval"]))
      #FDRt<-Q$qvalue
    }

    TOTRES<-cbind(TOTRES,FDR*100)
    colnames(TOTRES)[ncol(TOTRES)]<-'ANOVA FEATURE FDR %'

    assoc_id<-paste('a',1:nrow(TOTRES),sep='')

    TOTRES<-cbind(assoc_id,TOTRES)

    colnames(TOTRES)[1]<-'assoc_id'

    #write.table(TOTRES,quote=FALSE,row.names=F,sep='\t',file=paste(fn,'.txt',sep=''))

    #save(TOTRES,file=paste(fn,'.rdata',sep=''))

    if(length(range)>1){
      Sys.sleep(1)
      close(pb)
    }
  }else{
    TOTRES<-NULL
  }
  return(TOTRES)
}



ANOVA_singleDepGeneInteraction<-function(DEP_GENE,sorted=FALSE,verbose=TRUE,InputFeatures,ESSprofiles,OUTPUT_DIR){

nFEATURES<-nrow(InputFeatures$IntTerms)
FEATURES<-InputFeatures$IntTerms
  if(nFEATURES>0){
    for (i in 1:nFEATURES){

      if(verbose){
        setTxtProgressBar(pb, i)
      }
      RES<-ANOVA_individualANOVAInt(DEP_GENE=DEP_GENE,FEATUREmut=FEATURES[i,"mutation"],FEATUREexpr=FEATURES[i,"expression"],display = TRUE,InputFeatures=InputFeatures,ESSprofiles = ESSprofiles,OUTPUT_PATH = OUTPUT_DIR)
      #RES2<-ANOVA_individualANOVAIntM(DEP_GENE=DEP_GENE,FEATUREmut=FEATURES[i,"mutation1"],FEATUREmut2=FEATURES[i,"mutation2"],display = FALSE,InputFeatures=InputFeatures,ESSprofiles = ESSprofiles)

      if (i ==1){
        if(!is.null(RES)){
          TOTRES<-RES
        }
        #if(!is.null(RES2)){
       #   TOTRES2<-RES2
       # }
      } else{
        if(exists('TOTRES')){
          if(!is.null(RES)){
            TOTRES<-rbind(TOTRES,RES)
          }
        }else{
          if(!is.null(RES)){
            TOTRES<-RES
          }
        }
        if(exists('TOTRES2')){
          if(!is.null(RES2)){
            TOTRES2<-rbind(TOTRES2,RES2)
          }
        }else{
          #if(!is.null(RES2)){
          #  TOTRES2<-RES2
          #}
        }

      }
    }
    if(exists('TOTRES')){
      idxs<-which(!is.na(as.numeric(TOTRES[,"FEATURE_INT_pval"])))

      TOTRES<-TOTRES[idxs,]

      if (length(idxs)==1){
        TOTRES<-matrix(TOTRES,1,length(TOTRES),dimnames=list(1,names(TOTRES)))
      }

      if(verbose){
        Sys.sleep(1)
        close(pb)
      }

      if(sorted){
        TOTRES<-TOTRES[order(as.numeric(TOTRES[,"FEATURE_INT_pval"])),]
      }
    }else{TOTRES=NULL}
    if(exists('TOTRES2')){
      idxs<-which(!is.na(as.numeric(TOTRES2[,"FEATURE_INT_pval"])))

      TOTRES2<-TOTRES2[idxs,]

      if (length(idxs)==1){
        TOTRES2<-matrix(TOTRES2,1,length(TOTRES2),dimnames=list(1,names(TOTRES2)))
      }

      if(verbose){
        Sys.sleep(1)
        close(pb)
      }

      if(sorted){
        TOTRES2<-TOTRES2[order(as.numeric(TOTRES2[,"FEATURE_INT_pval"])),]
      }
    }else{TOTRES2=NULL}
  }else{TOTRES=NULL
  TOTRES2=NULL}
  #return(list(TOTRES=TOTRES,TOTRES2=TOTRES2))
return(TOTRES)
}

ANOVA_individualANOVAInt<-function(DEP_GENE,
                                FEATUREmut,FEATUREexpr,
                                display=TRUE,
                                printTOfig=FALSE,
                                PATH='',
                                FN='',
                                FDR=NA,
                                OUTPUT_PATH='',
                                ID='',
                                col1=GRAY,
                                col2=BLUE,InputFeatures,ESSprofiles){



commonC<-intersect(rownames(ESSprofiles),colnames(InputFeatures$Mut_BEM))
  commonC<-commonC[!is.na(ESSprofiles[commonC,DEP_GENE])]


  Epattern<-ESSprofiles[commonC,DEP_GENE]


  FEATpatternMut<-InputFeatures$Mut_BEM[FEATUREmut,commonC]

  FEATpatternMut<-FEATpatternMut[!is.na(Epattern)]
  FEATpatternExpr<-InputFeatures$Expr_BEM[as.character(FEATUREexpr),commonC]

  FEATpatternExpr<-FEATpatternExpr[!is.na(Epattern)]


  Epattern<-Epattern[!is.na(Epattern)]

  FEATpatternMut[which(FEATpatternMut==0)]<-'neg'
  FEATpatternMut[which(FEATpatternMut=='1')]<-'pos'
  FEATpatternMut<-as.factor(FEATpatternMut)

  FEATpatternExpr[which(FEATpatternExpr==0)]<-'neg'
  FEATpatternExpr[which(FEATpatternExpr=='1')]<-'pos'
  FEATpatternExpr<-as.factor(FEATpatternExpr)


  Y <- Epattern

      #fit <- aov(Y ~ FEATpatternMut*FEATpatternExpr)

    #Y<-anova(fit)
  fit<-lm(Y~FEATpatternMut*FEATpatternExpr)
   Y<-summary(fit)


        FEATURE_F1_PVAL<-Y$coefficients[2,4]
        FEATURE_F2_PVAL<-Y$coefficients[3,4]
        FEATURE_INT_PVAL<-Y$coefficients[4,4]
        #MSI_PVAL<-Y[2,5]
        MEDIA_PVAL<-NA
       # tissue_PVAL<-Y[1,5]
if(display){
  if(FEATURE_INT_PVAL<0.001){
    save(fit,file=paste0(OUTPUT_PATH,"/ExampleIntANOVAoutput.Rdata"))
    pdf(paste0(OUTPUT_PATH,"/IntPlotExprMut_",DEP_GENE,"_",FEATUREexpr,"_",FEATUREmut,".pdf"))
      print(plot_model(fit,type="int"))
      Sys.sleep(20)
    dev.off()
  }
}

    #FEATURE_ESS_WTT_pvalue<-fit$p.value

    pos_F1_ESS_MEAN<-mean(Epattern[FEATpatternMut=='pos'&FEATpatternExpr=='pos'])
    neg_F1_ESS_MEAN<-mean(Epattern[FEATpatternMut=='pos'&FEATpatternExpr=='neg'])
    pos_F2_ESS_MEAN<-mean(Epattern[FEATpatternMut=='neg'&FEATpatternExpr=='pos'])
    neg_F2_ESS_MEAN<-mean(Epattern[FEATpatternMut=='neg'&FEATpatternExpr=='neg'])

    deltaMEAN_ESS_F1<-pos_F1_ESS_MEAN-neg_F1_ESS_MEAN

    pos_ESS_sd<-sd(Epattern[FEATpatternMut=='pos'])
    neg_ESS_sd<-sd(Epattern[FEATpatternMut=='neg'])

    Npos<-length(which(FEATpatternMut=='pos'))
    Nneg<-length(which(FEATpatternMut=='neg'))

    EFFECTSIZE_ESS<-ANOVA_cohens_d(Epattern[FEATpatternMut=='pos'],Epattern[FEATpatternMut=='neg'])
    GLASS_dmp<-ANOVA_glass_Ds(Epattern[FEATpatternMut=='pos'&FEATpatternExpr=='pos'],Epattern[FEATpatternMut=='pos'&FEATpatternExpr=='neg'])
    GLASS_dmn<-ANOVA_glass_Ds(Epattern[FEATpatternMut=='neg'&FEATpatternExpr=='pos'],Epattern[FEATpatternMut=='neg'&FEATpatternExpr=='neg'])


    RES<-matrix(c(FEATUREmut,FEATUREexpr,
                  DEP_GENE,
                  geneAnnotations[DEP_GENE,c('name',
                                             'gene_family',
                                             'location_sortable',
                                             'entrez_id',
                                             'ensembl_gene_id',
                                             'pubmed_id','string_id')],
                  Npos,
                  Nneg,
                  pos_F1_ESS_MEAN,
                  neg_F1_ESS_MEAN,
                  pos_F2_ESS_MEAN,
                  neg_F2_ESS_MEAN,
                  deltaMEAN_ESS_F1,
                  pos_ESS_sd,
                  neg_ESS_sd,
                  EFFECTSIZE_ESS,
                  GLASS_dmp$g1,
                  GLASS_dmp$g2,
                  GLASS_dmn$g1,
                  GLASS_dmn$g2,
                  FEATURE_F1_PVAL,
                  FEATURE_F2_PVAL,
                  FEATURE_INT_PVAL
                  #tissue_PVAL,
                  #MEDIA_PVAL,
                  #MSI_PVAL,
                  #FEATURE_ESS_WTT_pvalue)
                  ),1,27,
                dimnames=list('1',c('FEATURE_mutation',"FEATURE_expression",
                                    'Depleted Gene',
                                    'Name',
                                    'gene family',
                                    'location',
                                    'entrez id',
                                    'ensembl gene id',
                                    'pumbed id',
                                    'string_id',
                                    'N FEATURE pos',
                                    'N FEATURE_neg',
                                    'FEATUREmutPexpP_ESS_MEAN',
                                    'FEATUREmutPexpN_ESS_MEAN',
                                    'FEATUREmutNexpP_ESS_MEAN',
                                    'FEATUREmutNexpN_ESS_MEAN',
                                    'FEATURE_deltaMEAN_ESS',
                                    'FEATUREpos_ESS_sd',
                                    'FEATUREneg_ESS_sd',
                                    'FEATURE_ESS_effect_size',
                                    'FEATUREmutP_Epos_Glass_delta',
                                    'FEATUREmutP_Eneg_Glass_delta',
                                    'FEATUREmutN_Epos_Glass_delta',
                                    'FEATUREmutN_Eneg_Glass_delta',
                                    'FEATURE_MUT_pval',
                                    'FEATURE_EXPR_pval',
                                    'FEATURE_INT_pval')))

  return(RES)
}

ANOVA_individualANOVAIntM<-function(DEP_GENE,
                                   FEATUREmut,FEATUREmut2,
                                   display=TRUE,
                                   printTOfig=FALSE,
                                   PATH='',
                                   FN='',
                                   FDR=NA,
                                   OUTPUT_PATH='',
                                   ID='',
                                   col1=GRAY,
                                   col2=BLUE,InputFeatures,ESSprofiles){



  commonC<-intersect(rownames(ESSprofiles),colnames(InputFeatures$Mut_BEM))
  commonC<-commonC[!is.na(ESSprofiles[commonC,DEP_GENE])]


  Epattern<-ESSprofiles[commonC,DEP_GENE]


  FEATpatternMut<-InputFeatures$Mut_BEM[FEATUREmut,commonC]

  FEATpatternMut<-FEATpatternMut[!is.na(Epattern)]
  FEATpatternMut2<-InputFeatures$Mut_BEM[as.character(FEATUREmut2),commonC]

  FEATpatternMut2<-FEATpatternMut2[!is.na(Epattern)]


  Epattern<-Epattern[!is.na(Epattern)]

  FEATpatternMut[which(FEATpatternMut==0)]<-'neg'
  FEATpatternMut[which(FEATpatternMut=='1')]<-'pos'
  FEATpatternMut<-as.factor(FEATpatternMut)

  FEATpatternMut2[which(FEATpatternMut2==0)]<-'neg'
  FEATpatternMut2[which(FEATpatternMut2=='1')]<-'pos'
  FEATpatternMut2<-as.factor(FEATpatternMut2)


  Y <- Epattern

  #fit <- aov(Y ~ FEATpatternMut*FEATpatternExpr)

  #Y<-anova(fit)
  fit<-lm(Y~FEATpatternMut*FEATpatternMut2)
  Y<-summary(fit)
  save(fit,file=paste0(dir.Results,"/ExampleIntANOVAoutput.Rdata"))

  FEATURE_F1_PVAL<-Y$coefficients[2,4]
  FEATURE_F2_PVAL<-Y$coefficients[3,4]
  FEATURE_INT_PVAL<-Y$coefficients[4,4]
  #MSI_PVAL<-Y[2,5]
  MEDIA_PVAL<-NA
  # tissue_PVAL<-Y[1,5]


  #FEATURE_ESS_WTT_pvalue<-fit$p.value

  pos_F1_ESS_MEAN<-mean(Epattern[FEATpatternMut=='pos'&FEATpatternMut2=='pos'])
  neg_F1_ESS_MEAN<-mean(Epattern[FEATpatternMut=='pos'&FEATpatternMut2=='neg'])
  pos_F2_ESS_MEAN<-mean(Epattern[FEATpatternMut=='neg'&FEATpatternMut2=='pos'])
  neg_F2_ESS_MEAN<-mean(Epattern[FEATpatternMut=='neg'&FEATpatternMut2=='neg'])

  deltaMEAN_ESS_F1<-pos_F1_ESS_MEAN-neg_F1_ESS_MEAN

  pos_ESS_sd<-sd(Epattern[FEATpatternMut=='pos'])
  neg_ESS_sd<-sd(Epattern[FEATpatternMut=='neg'])

  Npos<-length(which(FEATpatternMut=='pos'))
  Nneg<-length(which(FEATpatternMut=='neg'))

  EFFECTSIZE_ESS<-ANOVA_cohens_d(Epattern[FEATpatternMut=='pos'],Epattern[FEATpatternMut=='neg'])
  GLASS_dmp<-ANOVA_glass_Ds(Epattern[FEATpatternMut=='pos'&FEATpatternMut2=='pos'],Epattern[FEATpatternMut=='pos'&FEATpatternMut2=='neg'])
  GLASS_dmn<-ANOVA_glass_Ds(Epattern[FEATpatternMut=='neg'&FEATpatternMut2=='pos'],Epattern[FEATpatternMut=='neg'&FEATpatternMut2=='neg'])


  RES<-matrix(c(FEATUREmut,FEATUREmut2,
                DEP_GENE,
                geneAnnotations[DEP_GENE,c('name',
                                           'gene_family',
                                           'location_sortable',
                                           'entrez_id',
                                           'ensembl_gene_id',
                                           'pubmed_id','string_id')],
                Npos,
                Nneg,
                pos_F1_ESS_MEAN,
                neg_F1_ESS_MEAN,
                pos_F2_ESS_MEAN,
                neg_F2_ESS_MEAN,
                deltaMEAN_ESS_F1,
                pos_ESS_sd,
                neg_ESS_sd,
                EFFECTSIZE_ESS,
                GLASS_dmp$g1,
                GLASS_dmp$g2,
                GLASS_dmn$g1,
                GLASS_dmn$g2,
                FEATURE_F1_PVAL,
                FEATURE_F2_PVAL,
                FEATURE_INT_PVAL
                #tissue_PVAL,
                #MEDIA_PVAL,
                #MSI_PVAL,
                #FEATURE_ESS_WTT_pvalue)
  ),1,27,
  dimnames=list('1',c('FEATURE_mutation1',"FEATURE_mutation2",
                      'Depleted Gene',
                      'Name',
                      'gene family',
                      'location',
                      'entrez id',
                      'ensembl gene id',
                      'pumbed id',
                      'string_id',
                      'N FEATURE pos',
                      'N FEATURE_neg',
                      'FEATUREmut1Pmut2P_ESS_MEAN',
                      'FEATUREmut1Pmut2N_ESS_MEAN',
                      'FEATUREmut1Nmut2P_ESS_MEAN',
                      'FEATUREmut1Nmut2N_ESS_MEAN',
                      'FEATURE_deltaMEAN_ESS',
                      'FEATUREpos_ESS_sd',
                      'FEATUREneg_ESS_sd',
                      'FEATURE_ESS_effect_size',
                      'FEATUREmutP_M2pos_Glass_delta',
                      'FEATUREmutP_M2neg_Glass_delta',
                      'FEATUREmutN_M2pos_Glass_delta',
                      'FEATUREmutN_M2neg_Glass_delta',
                      'FEATURE_MUT_pval',
                      'FEATURE_MUT2_pval',
                      'FEATURE_INT_pval')))

  return(RES)
}


ANOVA_totalANOVACompound<-function(fn,ESSprofiles,InputFeatures,inparallel=TRUE,ncores=24,OUTPUT_DIR,MutualExclusive,
                                   geneAnnotation,ANOVA.settings.analysisType,ANOVA.fakeCS,
                                   ANOVA.settings.includeMSI_Factor){

  DEP_GENES<-colnames(ESSprofiles)

  nDEP_GENES<-length(DEP_GENES)

  print('+ Running ANOVA')


  flag<-1
  idxs<-chunk(1:nDEP_GENES,ncores)

  geneChunks<-lapply(idxs,function(x) DEP_GENES[x])

  if(inparallel){
    TOTRES<-foreach(i=1:ncores,.combine=rbind)%dopar%{

      temp<-DepGeneChunkCompound(GeneList=geneChunks[[i]],InputFeatures,ESSprofiles,OUTPUT_DIR=OUTPUT_DIR,MutualExclusive=MutualExclusive,
                                 geneAnnotation=geneAnnotation,ANOVA.settings.analysisType=ANOVA.settings.analysisType,ANOVA.fakeCS=ANOVA.fakeCS,
                                 ANOVA.settings.includeMSI_Factor=ANOVA.settings.includeMSI_Factor)

      temp
    }
  }else{
    for(i in 1:nDEP_GENES){
      print(c(DEP_GENES[i],i,nDEP_GENES))
      #setTxtProgressBar(pb, i)
      currentRES<-ANOVA_singleDepGeneCompound(DEP_GENE = DEP_GENES[i],verbose=FALSE,InputFeatures=InputFeatures,ESSprofiles = ESSprofiles,OUTPUT_DIR=OUTPUT_DIR,MutualExclusive=MutualExclusive,
                                              geneAnnotation=geneAnnotation,ANOVA.settings.analysisType=ANOVA.settings.analysisType,ANOVA.fakeCS=ANOVA.fakeCS,
                                              ANOVA.settings.includeMSI_Factor=ANOVA.settings.includeMSI_Factor)

      if(flag==1&!is.null(currentRES)){

        TOTRES<-currentRES
        flag<-0
      } else{

        if (!is.null(currentRES)){

          TOTRES<-rbind(TOTRES,currentRES)
        }
      }

    }
  }
  print('tested all genes in ANOVA_totalANOVA')

  if(exists('TOTRES')){

    TOTRES<-TOTRES[order(as.numeric(TOTRES[,"LRT_pval"])),]

    if(gdscANOVA.settings.pval_correction_method!='qvalue'){
      FDR<-p.adjust(as.numeric(TOTRES[,"LRT_pval"]),method='fdr')
      #FDRt<-p.adjust(as.numeric(TOTRES[,"FEATURE_ESS_T_pval"]),method='fdr')
    }else{
      pvals<-as.numeric(TOTRES[,"FEATURE_compound_pval"])
      print(paste("number p-vals:",length(pvals)))
      print(paste("max p-val:",max(pvals)))
      Q<-qvalue(as.numeric(TOTRES[,"FEATURE_compound_pval"]))
      FDR<-Q$qvalue

      Q2<-qvalue(as.numeric(TOTRES[,"LRT_pval"]))
      FDR2<-Q2$qvalue
      #commented out 15.6.19 as not included in output below
      #Q<-qvalue(as.numeric(TOTRES[,"FEATURE_ESS_T_pval"]))
      #FDRt<-Q$qvalue
    }

    TOTRES<-cbind(TOTRES,FDR*100,FDR2*100)
    colnames(TOTRES)[(ncol(TOTRES)-1):ncol(TOTRES)]<-c('ANOVA FEATURE FDR %','LRT FEATURE FDR %')

    assoc_id<-paste('a',1:nrow(TOTRES),sep='')

    TOTRES<-cbind(assoc_id,TOTRES)

    colnames(TOTRES)[1]<-'assoc_id'

    if(length(range)>1){
      Sys.sleep(1)
      close(pb)
    }
  }else{
    TOTRES<-NULL
  }
  return(TOTRES)
}

ANOVA_singleDepGeneCompound<-function(DEP_GENE,sorted=FALSE,verbose=TRUE,InputFeatures,ESSprofiles,OUTPUT_DIR,MutualExclusive,geneAnnotation
                                      ,ANOVA.settings.analysisType,ANOVA.fakeCS,
                                      ANOVA.settings.includeMSI_Factor){

  nFEATURES<-nrow(InputFeatures$IntTerms)
  FEATURES<-InputFeatures$IntTerms
  if(nFEATURES>0){
    for (i in 1:nFEATURES){

      if(verbose){
        setTxtProgressBar(pb, i)
      }
      RES<-ANOVA_individualLRTCompound(DEP_GENE=DEP_GENE,FEATUREparent=FEATURES[i,"Parent"],FEATUREchild=FEATURES[i,"Child"],display = TRUE,InputFeatures=InputFeatures,
                                       ESSprofiles = ESSprofiles,OUTPUT_PATH = OUTPUT_DIR,MutualExclusive=MutualExclusive,geneAnnotations=geneAnnotation
                                       ,ANOVA.settings.analysisType=ANOVA.settings.analysisType,ANOVA.fakeCS=ANOVA.fakeCS,
                                       ANOVA.settings.includeMSI_Factor=ANOVA.settings.includeMSI_Factor)

      if (i ==1){
        if(!is.null(RES)){
          TOTRES<-RES
        }

      } else{
        if(exists('TOTRES')){
          if(!is.null(RES)){
            TOTRES<-rbind(TOTRES,RES)
          }
        }else{
          if(!is.null(RES)){
            TOTRES<-RES
          }
        }

      }
    }
    if(exists('TOTRES')){

      idxs<-which(!is.na(as.numeric(TOTRES[,"LRT_pval"])))

      TOTRES<-TOTRES[idxs,]

      if (length(idxs)==1){
        TOTRES<-matrix(TOTRES,1,length(TOTRES),dimnames=list(1,names(TOTRES)))
      }

      if(verbose){
        Sys.sleep(1)
        close(pb)
      }

      if(sorted){
        TOTRES<-TOTRES[order(as.numeric(TOTRES[,"LRT_pval"])),]
      }
    }else{TOTRES=NULL}

  }else{TOTRES=NULL}
  #return(list(TOTRES=TOTRES,TOTRES2=TOTRES2))
  return(TOTRES)
}

ANOVA_individualLRTCompound<-function(DEP_GENE,
                                   FEATUREparent,FEATUREchild,
                                   display=TRUE,
                                   printTOfig=FALSE,
                                   PATH='',
                                   FN='',
                                   FDR=NA,
                                   OUTPUT_PATH='',
                                   ID='',
                                   col1=GRAY,
                                   col2=BLUE,InputFeatures,ESSprofiles,MutualExclusive,geneAnnotations
                                   ,ANOVA.settings.analysisType,ANOVA.fakeCS,
                                   ANOVA.settings.includeMSI_Factor){



  commonC<-intersect(rownames(ESSprofiles),colnames(InputFeatures$ParentBEM))
  commonC<-commonC[!is.na(ESSprofiles[commonC,DEP_GENE])]

  Epattern<-ESSprofiles[commonC,DEP_GENE]
  if(ANOVA.settings.analysisType=='PANCAN' & !ANOVA.fakeCS){
    print("Loading MSI and TISSUE variables ")
    TISSUE_FACTOR<-InputFeatures$TISSUE_VARIABLE[commonC]
    MSI_FACTOR<-InputFeatures$MSI_VARIABLE[commonC]

  }
  FEATpatternParent<-InputFeatures$ParentBEM[FEATUREparent,commonC]

  FEATpatternParent<-FEATpatternParent[!is.na(Epattern)]
  FEATpatternChild<-InputFeatures$ChildBEM[as.character(FEATUREchild),commonC]

  FEATpatternChild<-FEATpatternChild[!is.na(Epattern)]
  Epattern<-Epattern[!is.na(Epattern)]
  ##added NA removal for incomplete BEM 24.7.20:
  if(ANOVA.settings.analysisType=='PANCAN' & !ANOVA.fakeCS){
    rmvna<-(is.na(FEATpatternChild)|is.na(FEATpatternParent)|is.na(TISSUE_FACTOR)|is.na(MSI_FACTOR))

    MSI_FACTOR<-MSI_FACTOR[!rmvna]
    TISSUE_FACTOR<-TISSUE_FACTOR[!rmvna]
    MSI_FACTOR[which(MSI_FACTOR==0)]<-'neg'
    MSI_FACTOR[which(MSI_FACTOR=="1")]<-'pos'
    TISSUEpattern<-as.factor(TISSUE_FACTOR)
    MSIpattern<-as.factor(MSI_FACTOR)

  }else{
    rmvna<-(is.na(FEATpatternChild)|is.na(FEATpatternParent))
  }
  Epattern<-Epattern[!rmvna]
  FEATpatternParent<-FEATpatternParent[!rmvna]
  FEATpatternChild<-FEATpatternChild[!rmvna]
  ##
  FEATpatternChild[which(FEATpatternChild==0)]<-'neg'
  FEATpatternChild[which(FEATpatternChild=='1')]<-'pos'
  FEATpatternChild<-as.factor(FEATpatternChild)


  FEATpatternParent[which(FEATpatternParent==0)]<-'neg'
  FEATpatternParent[which(FEATpatternParent=='1')]<-'pos'
  FEATpatternParent<-as.factor(FEATpatternParent)

  #create OR feature here:
  FEATpatternOR<-rep('neg',length(FEATpatternChild))
  FEATpatternOR[c(which(FEATpatternParent=="pos"),which(FEATpatternChild=="pos"))]<-"pos"
  FEATpatternOR<-as.factor(FEATpatternOR)

  #create compound feature AND here:
  FEATpatternCompound<-rep('neg',length(FEATpatternChild))
  FEATpatternCompound[which(FEATpatternParent=='pos'&FEATpatternChild=='pos')]<-'pos'
  FEATpatternCompound<-as.factor(FEATpatternCompound)



  Y <- Epattern


  #fit1<-lm(Y~FEATpatternParent+FEATpatternChild+FEATpatternOR)
  #fit2<-lm(Y~FEATpatternParent+FEATpatternChild+FEATpatternOR+FEATpatternCompound)


  if(MutualExclusive){
    FEATURE<-paste0(FEATUREparent,"|",FEATUREchild)
    BMstringID<-paste0(geneAnnotations[strsplit(FEATUREparent,"_")[[1]][1],"string_id"],"&",geneAnnotations[strsplit(FEATUREchild,"_")[[1]][1],"string_id"])
    if(ANOVA.settings.analysisType=='PANCAN' & !ANOVA.fakeCS){

      fit1<-lm(Y~FEATpatternParent+MSIpattern+TISSUEpattern)
      fit2<-lm(Y~FEATpatternParent+FEATpatternOR+MSIpattern+TISSUEpattern)
      fit1.2<-lm(Y~FEATpatternChild+MSIpattern+TISSUEpattern)
      fit2.2<-lm(Y~FEATpatternChild+FEATpatternOR+MSIpattern+TISSUEpattern)

    }else{
      if(ANOVA.settings.includeMSI_Factor){
        fit1<-lm(Y~FEATpatternParent+MSIpattern)
        fit2<-lm(Y~FEATpatternParent+FEATpatternOR+MSIpattern)
        fit1.2<-lm(Y~FEATpatternChild+MSIpattern)
        fit2.2<-lm(Y~FEATpatternChild+FEATpatternOR+MSIpattern)
      }else{
        fit1<-lm(Y~FEATpatternParent)
        fit2<-lm(Y~FEATpatternParent+FEATpatternOR)
        fit1.2<-lm(Y~FEATpatternChild)
        fit2.2<-lm(Y~FEATpatternChild+FEATpatternOR)


      }
    }


    lmtest<-data.frame(lrtest(fit2,fit1))
    lmtest.2<-data.frame(lrtest(fit2.2,fit1.2))
    Y1<-summary(fit1)
    Y2<-summary(fit2)
    Y1.2<-summary(fit1.2)
    Y2.2<-summary(fit2.2)

    LRTpvalue<-max(lmtest[2,5],lmtest.2[2,5])

    FEATURE_Parent_pval<-Y2$coefficients[2,4]

    FEATURE_Child_PVAL<-Y2.2$coefficients[2,4]
    FEATURE_Compound_PVAL<-max(Y2$coefficients[3,4],Y2.2$coefficients[3,4])
    if(ANOVA.settings.analysisType=='PANCAN' & !ANOVA.fakeCS){


      FEATURE_TISSUE_pval<-min(Y2$coefficients[5:nrow(Y2$coefficients),4])
      FEATURE_MSI_pval<-Y2$coefficients[4,4]

    }else{
      if(ANOVA.settings.includeMSI_Factor){


        FEATURE_TISSUE_pval<-NA
        FEATURE_MSI_pval<-Y2$coefficients[4,4]
      }else{


        FEATURE_TISSUE_pval<-NA
        FEATURE_MSI_pval<-NA
      }
    }
    pos_F1_ESS_MEAN<-mean(Epattern[FEATpatternParent=='pos'])
    WT_ESS_MEAN<-mean(Epattern[FEATpatternParent=='neg'&FEATpatternChild=='neg'])

    Compound_pos_F2_ESS_MEAN<-mean(Epattern[FEATpatternOR=='pos'])



    deltaMEAN_ESS_F1<-Compound_pos_F2_ESS_MEAN-WT_ESS_MEAN

    pos_ESS_sd<-sd(Epattern[FEATpatternParent=='pos'])
    neg_ESS_sd<-sd(Epattern[FEATpatternParent=='neg'])

    Npos<-length(which(FEATpatternParent=='pos'))
    Nneg<-length(which(FEATpatternParent=='neg'))

    EFFECTSIZE_ESS<-ANOVA_cohens_d(Epattern[FEATpatternParent=='pos'],Epattern[FEATpatternParent=='neg'])

    EFFECTSIZE_ESSC<-ANOVA_cohens_d(Epattern[FEATpatternOR=='pos'],Epattern[FEATpatternOR=='neg'])

    GLASS_dmn<-ANOVA_glass_Ds(Epattern[FEATpatternOR=='pos'],Epattern[FEATpatternOR=='neg'])

    GLASS_dmp<-GLASS_dmn
    GLASS_dmp$g1<-NA
    GLASS_dmp$g2<-NA
  }else{
    FEATURE<-paste0(FEATUREparent,"&",FEATUREchild)
    BMstringID<-paste0(geneAnnotations[strsplit(FEATUREparent,"_")[[1]][1],"string_id"],"&",geneAnnotations[strsplit(FEATUREchild,"_")[[1]][1],"string_id"])
    if(ANOVA.settings.analysisType=='PANCAN' & !ANOVA.fakeCS){

      fit1<-lm(Y~FEATpatternParent+FEATpatternChild+MSIpattern+TISSUEpattern)
      fit2<-lm(Y~FEATpatternParent+FEATpatternChild+FEATpatternCompound+MSIpattern+TISSUEpattern)
    }else{
      if(ANOVA.settings.includeMSI_Factor){

        fit1<-lm(Y~FEATpatternParent+FEATpatternChild+MSIpattern)
        fit2<-lm(Y~FEATpatternParent+FEATpatternChild+FEATpatternCompound+MSIpattern)
      }else{
        fit1<-lm(Y~FEATpatternParent+FEATpatternChild)
        fit2<-lm(Y~FEATpatternParent+FEATpatternChild+FEATpatternCompound)
      }
    }



    #require LRT test between fit2 and fit2 to see if the AND compound term significantly increases the fit to the dependency profile
    lmtest<-data.frame(lrtest(fit2,fit1))
    Y1<-summary(fit1)
    Y2<-summary(fit2)

    LRTpvalue<-lmtest[2,5]
    FEATURE_Parent_pval<-Y2$coefficients[2,4]

    FEATURE_Child_PVAL<-Y2$coefficients[3,4]
    FEATURE_Compound_PVAL<-Y2$coefficients[4,4]

    if(ANOVA.settings.analysisType=='PANCAN' & !ANOVA.fakeCS){


      FEATURE_MSI_pval<-Y2$coefficients[5,4]
      FEATURE_TISSUE_pval<-min(Y2$coefficients[6:nrow(Y2$coefficients),4])
    }else{
      if(ANOVA.settings.includeMSI_Factor){


        FEATURE_MSI_pval<-Y2$coefficients[5,4]
        FEATURE_TISSUE_pval<-NA
      }else{

        FEATURE_MSI_pval<-NA
        FEATURE_TISSUE_pval<-NA
      }
    }

    pos_F1_ESS_MEAN<-mean(Epattern[FEATpatternParent=='pos'])
    WT_ESS_MEAN<-mean(Epattern[FEATpatternParent=='neg'&FEATpatternChild=='neg'])
    if(MutualExclusive){
      Compound_pos_F2_ESS_MEAN<-mean(Epattern[FEATpatternOR=='pos'])
    }else{
      Compound_pos_F2_ESS_MEAN<-mean(Epattern[FEATpatternCompound=='pos'])
    }


    deltaMEAN_ESS_F1<-Compound_pos_F2_ESS_MEAN-WT_ESS_MEAN

    pos_ESS_sd<-sd(Epattern[FEATpatternParent=='pos'])
    neg_ESS_sd<-sd(Epattern[FEATpatternParent=='neg'])

    Npos<-length(which(FEATpatternParent=='pos'))
    Nneg<-length(which(FEATpatternParent=='neg'))

    EFFECTSIZE_ESS<-ANOVA_cohens_d(Epattern[FEATpatternParent=='pos'],Epattern[FEATpatternParent=='neg'])

    EFFECTSIZE_ESSC<-ANOVA_cohens_d(Epattern[FEATpatternParent=='pos'&FEATpatternCompound=='neg'],Epattern[FEATpatternCompound=='pos'])
    GLASS_dmp<-ANOVA_glass_Ds(Epattern[FEATpatternParent=='pos'&FEATpatternCompound=='neg'],Epattern[FEATpatternCompound=='pos'])
    GLASS_dmn<-ANOVA_glass_Ds(Epattern[FEATpatternCompound=='pos'],Epattern[FEATpatternCompound=='neg'])

  }


  RES<-data.frame(c(FEATURE,FEATUREparent,FEATUREchild,
                DEP_GENE,
                BMstringID,
                geneAnnotations[DEP_GENE,c('name',
                                           'gene_family',
                                           'location_sortable',
                                           'entrez_id',
                                           'ensembl_gene_id',
                                           'pubmed_id','string_id')],
                Npos,
                Nneg,
                WT_ESS_MEAN,
                Compound_pos_F2_ESS_MEAN,
                deltaMEAN_ESS_F1,
                pos_ESS_sd,
                neg_ESS_sd,
                EFFECTSIZE_ESS,
                EFFECTSIZE_ESSC,
                GLASS_dmp$g1,
                GLASS_dmp$g2,
                GLASS_dmn$g1,
                GLASS_dmn$g2,
                FEATURE_Parent_pval,
                FEATURE_Child_PVAL,
                FEATURE_Compound_PVAL,
                FEATURE_TISSUE_pval,
                FEATURE_MSI_pval,
                LRTpvalue),stringsAsFactors = FALSE)

  colnames(RES)<-c('FEATURE','FEATURE_parent',"FEATURE_child",
                      'Depleted Gene',
                   "BMstringID",
                      'Name',
                      'gene family',
                      'location',
                      'entrez id',
                      'ensembl gene id',
                      'pumbed id',
                      'string_id',
                      'N FEATURE pos',
                      'N FEATURE_neg',
                      'FEATUREwt_ESS_MEAN',
                      'FEATUREcompound_ESS_MEAN',
                      'deltaMEAN_ESS_F1',
                      'FEATUREpos_ESS_sd',
                      'FEATUREneg_ESS_sd',
                      'FEATURE_Parent_ESS_effect_size',
                      'FEATURE_Compound_effect_size',
                      'FEATUREparent_pos_Glass_delta',
                      'FEATUREparent_neg_Glass_delta',
                      'FEATUREcompound_pos_Glass_delta',
                      'FEATUREcompound_neg_Glass_delta',
                      'FEATURE_parent_pval',
                      'FEATURE_child_pval',
                      'FEATURE_compound_pval',
                      'FEATURE_TISSUE_pval',
                      'FEATURE_MSI_pval',
                      'LRT_pval')

  return(RES)
}
