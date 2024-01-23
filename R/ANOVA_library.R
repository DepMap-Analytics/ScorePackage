
## data manipulation

ANOVA_createEssentialityDataInput<-function(ESSENTIALITY_DOMAIN,
                                                   CELLLINESTOREMOVE,
                                                   excludeKnownAndPredEssential=TRUE,
                                                   diagnostic=TRUE,
                                                   ExclusiveGeneList='',
                                                   scaleToLimit=FALSE,
                                            quantileNorm=FALSE,
                                            minVulnerableCellLines=2){

    fn<-paste(ANOVA.essentiality.dir,'SuperSummary.Rdata',sep='')
    load(fn)

    print('Gathering essentiality profiles...')

    if(scaleToLimit){
        ESSs<-
            SuperSummary[,seq(5,ncol(SuperSummary),5)]
    }else{
        ESSs<-
            SuperSummary[,seq(3,ncol(SuperSummary),5)]
    }


    SCALED_ESSs<-SuperSummary[,seq(5,ncol(SuperSummary),5)]


    if(quantileNorm){
        rn<-rownames(ESSs)
        cn<-colnames(ESSs)
        ESSs<-
            normalize.quantiles(ESSs)
        rownames(ESSs)<-rn
        colnames(ESSs)<-cn
    }


    colnames(ESSs)<-unlist(str_split(colnames(ESSs),' '))[seq(1,ncol(ESSs)*2,2)]
    colnames(SCALED_ESSs)<-colnames(ESSs)


    ESSs<-ESSs[,setdiff(colnames(ESSs),CELLLINESTOREMOVE)]
    SCALED_ESSs<-SCALED_ESSs[,colnames(ESSs)]

    print('')
    print('Done!')

    ESSs<-t(ESSs)
    SCALED_ESSs<-t(SCALED_ESSs)

    ndDomain<-ncol(ESSs)
    cat(paste(ndDomain,'genes in the selected domain\n'))


    if(ANOVA.settings.analysisType!='PANCAN'){
        SCALED_ESSs<-SCALED_ESSs[which(is.element(rownames(ESSs),PS_inventory$Cell_Line_Name[which(PS_inventory$GDSC_Description_2==ANOVA.settings.CELL_LINES_gdscType)])),]
        ESSs<-ESSs[which(is.element(rownames(ESSs),PS_inventory$Cell_Line_Name[which(PS_inventory$GDSC_Description_2==ANOVA.settings.CELL_LINES_gdscType)])),]
    }

    cosmicIds<-
        as.character(PS_inventory$COSMIC_ID[match(rownames(ESSs),PS_inventory$Cell_Line_Name)])

    rownames(ESSs)<-cosmicIds


    if(excludeKnownAndPredEssential){

        load(ANOVA.PANCANcoreFitnessResults.fn)
        load('../../PAPERfreeze_data/EssentialGeneSets/BAGEL/curated_BAGEL_essential.rdata')
        load('../../PAPERfreeze_data/EssentialGeneSets/BAGEL/BAGEL_v2_ESSENTIAL_GENES.rdata')

        egfn<-dir('../../PAPERfreeze_data/EssentialGeneSets/')
        egfn<-union(grep('.Rdata',egfn,value = TRUE),grep('.rdata',egfn,value = TRUE))

        priorEss<-union(curated_BAGEL_essential,BAGEL_essential)

        for (i in 1:length(egfn)){
            tmp<-ls()
            load(paste('../../PAPERfreeze_data/EssentialGeneSets/',egfn[i],sep=''))
            tmp<-setdiff(ls(),tmp)
            priorEss<-union(priorEss,eval(parse(text=tmp)))
        }

        priorEss<-union(priorEss,panCancer_cf_genes)


        if(ANOVA.settings.CELL_LINES=='BRCA'){
            load(paste(ANOVA.CScoreFitnessResults.dir,'breast_coreFitnessGenes.Rdata',sep=''))
        }
        if(ANOVA.settings.CELL_LINES=='COREAD'){
            load(paste(ANOVA.CScoreFitnessResults.dir,'large_intestine_coreFitnessGenes.Rdata',sep=''))
        }
        if(ANOVA.settings.CELL_LINES=='LUSC'){
            load(paste(ANOVA.CScoreFitnessResults.dir,'lung_NSCLC_squamous_cell_carcinoma_coreFitnessGenes.Rdata',sep=''))
        }
        if(ANOVA.settings.CELL_LINES=='OV'){
            load(paste(ANOVA.CScoreFitnessResults.dir,'ovary_coreFitnessGenes.Rdata',sep=''))
        }
        if(ANOVA.settings.CELL_LINES=='PAAD'){
            load(paste(ANOVA.CScoreFitnessResults.dir,'pancreas_coreFitnessGenes.Rdata',sep=''))
        }
        if(ANOVA.settings.CELL_LINES=='PANCAN'){
            load(ANOVA.PANCANcoreFitnessResults.fn)
            coreFitnessGenes<-panCancer_cf_genes
        }

        priorEss<-union(priorEss,coreFitnessGenes)

        load('../../PAPERfreeze_data/cancerGenes/IntoGen_cancer_drivers_13_06_2014.rdata')
        priorEss<-setdiff(priorEss,cancerDrivers$ActingDriver_Symbol)
        load('../../PAPERfreeze_data/cancerGenes/additional_9_inigo_list.rdata')
        priorEss<-setdiff(priorEss,InigoList)

        id<-setdiff(colnames(ESSs),priorEss)
        ESSs<-ESSs[,id]
        SCALED_ESSs<-SCALED_ESSs[,id]
        cat(paste(length(id),' non curated PANCAN_core-fitness, CS_core-fitness, a priori known essential genes maintained\n'))
    }

    if(length(ExclusiveGeneList)>0 & ExclusiveGeneList[1]!=''){
        id<-intersect(colnames(ESSs),ExclusiveGeneList)

        cat(paste(length(id),' exclusive genes maintained\n'))
        ESSs<-ESSs[,id]
        SCALED_ESSs<-SCALED_ESSs[,id]
    }

    id<-which(colSums(SCALED_ESSs>0)>=minVulnerableCellLines
              & colSums(SCALED_ESSs>0)<=(nrow(SCALED_ESSs)-minVulnerableCellLines))

    cat(paste(length(id),'genes singificantly depleted in at least',minVulnerableCellLines,' and at most',
              nrow(ESSs)-minVulnerableCellLines,'cell lines maintained\n'))
    ESSs<-ESSs[,id]




    return(ESSs)

}

getBiomarkerInput<-function(qnormLFC,bDepletions,minDep=2,TOTALBEM){
  bDepletions<-bDepletions[rownames(qnormLFC),]

  ESSprofiles<-qnormLFC


  ESSprofiles<-t(ESSprofiles)
  bDepletions<-t(bDepletions)


  ESSprofiles<-ESSprofiles[intersect(rownames(ESSprofiles),colnames(TOTALBEM)),]
  bDepletions<-bDepletions[rownames(ESSprofiles),colnames(ESSprofiles)]

  ESSprofiles<-ESSprofiles[,colSums(bDepletions)>minDep]
  bDepletions<-bDepletions[,colSums(bDepletions)>minDep]
  return(list(bDep=bDepletions,EssProf=ESSprofiles))
}
ANOVA_createInputFeatures<-function(additional_features=NULL,
                                           additional_features_only=FALSE,
                                           excludeHyperMetData=TRUE,
                                           oneFeatureOnly=NULL,TOTALBEM,
                                           selectedCellLines=NULL,MASTER_LIST,
                                    ANOVA.settings.includeMEDIA_factor=FALSE,
                                    ANOVA.settings.featFactorPopulationTh=2,modelID="COSMIC_ID"){
  if(additional_features_only){
    BEM<-additional_features
  }else{
    additional_features<-additional_features[,colnames(TOTALBEM)]
    BEM<-rbind(TOTALBEM,additional_features)
  }

   rownames(BEM)<-str_replace_all(rownames(BEM),pattern = ':',replacement = '_')
   rownames(BEM)<-str_replace_all(rownames(BEM),pattern = ' ',replacement = '_')

   if(excludeHyperMetData){
     idxs<-which(!str_detect(rownames(BEM),'HypMET'))
     BEM<-BEM[idxs,]
   }

  if(length(oneFeatureOnly)>0){

    for (kk in 1:length(oneFeatureOnly)){
      if (kk == 1){
        idxs<-which(str_detect(rownames(BEM),oneFeatureOnly[kk]))
      }else{
       idxs<-union(idxs,which(str_detect(rownames(BEM),oneFeatureOnly[kk])))
      }

    }

    if (length(idxs)==1){
      BEM<-matrix(BEM,1,ncol(BEM),dimnames = list(rownames(BEM)[idxs],colnames(BEM)))
    }else{
      BEM<-BEM[idxs,]
    }
  }

  idxs<-which(rowSums(abs(BEM),na.rm=T)>1)

  if (length(idxs)==1){
    BEM<-matrix(BEM,1,ncol(BEM),dimnames = list(rownames(BEM)[idxs],colnames(BEM)))
  }else{
    BEM<-BEM[idxs,]
  }

  if (length(selectedCellLines)){
      BEM<-BEM[,selectedCellLines]
  }
  ##24.7.20 added in NA handling and second population size check
  #idxs<-which(rowSums(abs(BEM))>ANOVA.settings.featFactorPopulationTh)
  idxs<-which(rowSums(abs(BEM),na.rm=T)>ANOVA.settings.featFactorPopulationTh&rowSums(1-abs(BEM),na.rm=T)>ANOVA.settings.featFactorPopulationTh)
  BEM<-BEM[idxs,]

   if(length(idxs)>1){
     cat(paste('\n\t','Features with n >',ANOVA.settings.featFactorPopulationTh,
               ' positive samples = ',nrow(BEM)))
     cat('\n\tMerging features with identical positive patterns:')
     BEM<-my.compress_identical_patterns(BEM)
   }
   BEM<-abs(BEM)
   cat(paste('\n\t','Features with n >',ANOVA.settings.featFactorPopulationTh,
' positive samples = ',nrow(BEM)))

   TISSUE_VARIABLE<-ANOVA_createTissueVariable(CID = colnames(BEM),modelID = modelID,MASTER_LIST=MASTER_LIST)

   N<-names(TISSUE_VARIABLE)
   TISSUE_VARIABLE<-as.character(TISSUE_VARIABLE)
   names(TISSUE_VARIABLE)<-N

   cat(paste('\t*d* ','Tissue variable with',length(unique(TISSUE_VARIABLE)),' factors'))

   if (length(idxs)>1){
     BEM<-BEM[,names(TISSUE_VARIABLE)]
   }else{
     BEM<-matrix(BEM,nrow = 1,ncol = length(BEM),dimnames = list(names(idxs),names(BEM)))
   }


   #MSI_VARIABLE<-as.character(MASTER_LIST[as.character(colnames(BEM)),'MMR'])
   #to be consistent with cell model passports download:
   MSI_VARIABLE<-as.character(MASTER_LIST[match(as.character(colnames(BEM)),MASTER_LIST$COSMIC_ID),'msi_status'])
   MSI_VARIABLE[which(MSI_VARIABLE!='MSI')]<-0
   MSI_VARIABLE[which(MSI_VARIABLE=='MSI')]<-1

   MSI_VARIABLE<-as.numeric(MSI_VARIABLE)
   names(MSI_VARIABLE)<-colnames(BEM)

   if(ANOVA.settings.includeMEDIA_factor){
       load(ANOVA.screeMedium.fn)
       SC_MEDIUM<-ASSAY_PROPERTIES[names(MSI_VARIABLE),"SCREEN_MEDIUM"]
       return(list(BEM=BEM,TISSUES=TISSUE_VARIABLE,SC_MEDIUM=SC_MEDIUM,MSI_VARIABLE=MSI_VARIABLE))
   }

  return(list(BEM=BEM,TISSUES=TISSUE_VARIABLE,MSI_VARIABLE=MSI_VARIABLE))
}

ANOVA_createTissueVariable<-function(CID,Limit=-Inf,modelID="COSMIC_ID",MASTER_LIST){
  if(modelID=="model_id"){
    M2<-MASTER_LIST
    M2$model_id<-M2$BROAD_ID
    MASTER_LIST<-rbind(MASTER_LIST,M2)
  }
  #TISSUE_VARIABLE<-as.character(MASTER_LIST$GDSC.description_1)
  TISSUE_VARIABLE<-as.character(MASTER_LIST$tissue)
  #names(TISSUE_VARIABLE)<-rownames(MASTER_LIST)

  names(TISSUE_VARIABLE)<-make.names(MASTER_LIST[,modelID])

  #TISSUE_VARIABLE[which(TISSUE_VARIABLE=='digestive_system')]<-MASTER_LIST$GDSC.description_2[which(TISSUE_VARIABLE=='digestive_system')]
  #TISSUE_VARIABLE[which(TISSUE_VARIABLE=='urogenital_system')]<-MASTER_LIST$GDSC.description_2[which(TISSUE_VARIABLE=='urogenital_system')]

  TISSUE_VARIABLE<-TISSUE_VARIABLE[CID]

  S<-summary(as.factor(TISSUE_VARIABLE))

  TISSUE_VARIABLE<-TISSUE_VARIABLE[is.element(TISSUE_VARIABLE,names(which(S>Limit)))]


  return(TISSUE_VARIABLE)
}
getADMfile<-function(admfiles,tissue){

  admctype<-grep(tissue,admfiles,value=TRUE,ignore.case=TRUE)

  admCF<-grep("coreFitnessGenes",admctype,value=TRUE)


  if(length(admCF)==1){

  }else{
    admCF<-admCF[grep(paste('09_ADM_',tissue,'_coreFitnessGenes.Rdata',sep=''),admCF,ignore.case=TRUE)]


  }
  return(admCF)
}

ANOVA_create_systemInfos<-function(ANOVA.settings.CELL_LINES,MAIN_DIR){

  ANALYSIS_SYSTEMS_INFOS<-list()
  DATE<-Sys.Date()
  TIME<-Sys.time()
  SYSINFO<-Sys.info()
  TIME<-str_replace_all(TIME,'[:]','-')
  TIME<-str_replace_all(TIME,'[ ]','_')
  ANALYSIS_SYSTEMS_INFOS$DATE<-DATE
  ANALYSIS_SYSTEMS_INFOS$TIME<-TIME
  ANALYSIS_SYSTEMS_INFOS$USER<-SYSINFO[7]
  ANALYSIS_SYSTEMS_INFOS$MACHINE<-SYSINFO[4]
  ANALYSIS_SYSTEMS_INFOS$ESSENTIALITY_DOMAIN<-
    ANOVA.settings.ESSENTIALITY_domain
  if(ANOVA.settings.CELL_LINES!='PANCAN'){
    analysis_type<-paste(ANOVA.settings.CELL_LINES,'specific')
  }else{
    analysis_type<-ANOVA.settings.CELL_LINES
  }

  ANALYSIS_SYSTEMS_INFOS$CELL_LINES_DOMAIN<-analysis_type
  ANALYSIS_SYSTEMS_INFOS$TISSUE_FACTOR<-ANOVA.settings.CELL_LINES=='PANCAN'
  ANALYSIS_SYSTEMS_INFOS$MSI_FACTOR<-ANOVA.settings.includeMSI_Factor

  ANALYSIS_SYSTEMS_INFOS$featFactorPopulationTh<-ANOVA.settings.featFactorPopulationTh
  if(ANOVA.settings.includeMSI_Factor){
    ANALYSIS_SYSTEMS_INFOS$MSI_FACTOR<-ANOVA.settings.MSIfactorPopulationTh
  }else{
    ANALYSIS_SYSTEMS_INFOS$MSI_FACTOR<-NA
  }

#   ANALYSIS_SYSTEMS_INFOS$N_DRUGS<-DIAGNOSTICS$N_DRUGS
#   ANALYSIS_SYSTEMS_INFOS$N_FEATURES<-DIAGNOSTICS$N_FEATURES
#   ANALYSIS_SYSTEMS_INFOS$N_COMBOS<-DIAGNOSTICS$N_COMBOS
#   ANALYSIS_SYSTEMS_INFOS$N_FEASIBLE_TESTS<-DIAGNOSTICS$N_FEASIBLE_TESTS
#   ANALYSIS_SYSTEMS_INFOS$PERC_FEAS_TESTS<-DIAGNOSTICS$PERC_FEAS_TESTS
#
  save(ANALYSIS_SYSTEMS_INFOS,file=paste(MAIN_DIR,'SYSTEM_INFOS.rdata',sep=''))
}


ANOVA_createInputFeaturesInteraction<-function(additional_features=NULL,
                                    additional_features_only=FALSE,
                                    excludeHyperMetData=TRUE,
                                    oneFeatureOnly=NULL,TOTALBEM,
                                    selectedCellLines=NULL,MASTER_LIST,ANOVA.settings.featFactorPopulationTh=5,marker2="_Expr"){

  print('getting features start')
  if(additional_features_only){
    BEM<-additional_features
  }else{
    additional_features<-additional_features[,colnames(TOTALBEM)]
    BEM<-rbind(TOTALBEM,additional_features)
  }

  rownames(BEM)<-str_replace_all(rownames(BEM),pattern = ':',replacement = '_')
  rownames(BEM)<-str_replace_all(rownames(BEM),pattern = ' ',replacement = '_')

  idxs<-which(rowSums(abs(BEM))>1)

  if (length(idxs)==1){
    BEM<-matrix(BEM,1,ncol(BEM),dimnames = list(rownames(BEM)[idxs],colnames(BEM)))
  }else{
    BEM<-BEM[idxs,]
  }

  if (length(selectedCellLines)){
    BEM<-BEM[,selectedCellLines]
  }

  Mut_BEM<-BEM[grep("_mut",rownames(BEM)),]
  Expr_BEM<-BEM[grep(marker2,rownames(BEM)),]
  #CNA_BEM<-BEM[grep(c("gain","loss"),rownames(BEM)),]
  print(paste("dim Mut_BEM",dim(Mut_BEM)))
  print(paste("dim Expr_BEM",dim(Expr_BEM)))
  idxsMut<-which(rowSums(abs(Mut_BEM))>ANOVA.settings.featFactorPopulationTh)
  idxsExpr<-which(rowSums(abs(Expr_BEM))>ANOVA.settings.featFactorPopulationTh)

  if(length(idxsMut)>2){
    #try mutation versus mutation interaction
    Mut_BEM1<-Mut_BEM[idxsMut,]
    print(Mut_BEM)
    print(Mut_BEM1)
    Mut_BEM1<-my.compress_identical_patterns(Mut_BEM1)
    print(typeof(Mut_BEM1))
    Int_Mut<-Mut_BEM1%*%t(Mut_BEM1)
    rownames(Int_Mut)<-rownames(Mut_BEM1)
    colnames(Int_Mut)<-rownames(Mut_BEM1)

    Not_Mut<-(!Mut_BEM1)*1
    Int_NM<-Not_Mut%*%t(Not_Mut)
    rownames(Int_NM)<-rownames(Mut_BEM1)
    colnames(Int_NM)<-rownames(Mut_BEM1)
    Int_N1M<-Not_Mut%*%t(Mut_BEM1)
    rownames(Int_N1M)<-rownames(Mut_BEM1)
    colnames(Int_N1M)<-rownames(Mut_BEM1)
    Int_MN1<-Mut_BEM1%*%t(Not_Mut)
    rownames(Int_MN1)<-rownames(Mut_BEM1)
    colnames(Int_MN1)<-rownames(Mut_BEM1)
    diag(Int_Mut)<-0
    diag(Int_NM)<-0
    diag(Int_N1M)<-0
    diag(Int_MN1)<-0
    check<-(Int_Mut>2)&(Int_NM>2)&(Int_N1M>2)&(Int_MN1>2)
    dimnames(check)<-dimnames(Int_Mut)

    if(sum(check)>0){
      #have interactions to test:
      #find mutations to test first:
      #number of interactions:
      output<-melt(check) %>% filter(value == 1)
      numberInt<-nrow(output)
      IntTermsMM<-data.frame(mutation1=as.character(output[,"Var2"]),mutation2=as.character(output[,"Var1"]),stringsAsFactors = FALSE)


    }else{
      IntTermsMM<-NULL
    }
  }else{IntTermsMM<-NULL}

  if(length(idxsMut)>1&length(idxsExpr)>1){
    Mut_BEM<-Mut_BEM[idxsMut,]
    Expr_BEM<-Expr_BEM[idxsExpr,]
    #my.compress_identical_patterns is in OT15.Misc

    Mut_BEM<-my.compress_identical_patterns(Mut_BEM)

    Expr_BEM<-my.compress_identical_patterns(Expr_BEM)
    #IntMat is Mutated and expressed
    IntMat<-Mut_BEM%*%t(Expr_BEM)
    rownames(IntMat)<-rownames(Mut_BEM)
    colnames(IntMat)<-rownames(Expr_BEM)
    #NullMat is not mutated and not expressed
    NullMut<-(!Mut_BEM)*1
    NullExpr<-(!Expr_BEM)*1
    NullMat<-NullMut%*%t(NullExpr)
    rownames(NullMat)<-rownames(Mut_BEM)
    colnames(NullMat)<-rownames(Expr_BEM)
    #MutNullExpr is mutated no expression
    MutNullExpr<-Mut_BEM%*%t(NullExpr)
    rownames(MutNullExpr)<-rownames(Mut_BEM)
    colnames(MutNullExpr)<-rownames(Expr_BEM)
    #NullMutExpr is not mutated and expressed
    NullMutExpr<-NullMut%*%t(Expr_BEM)
    rownames(NullMutExpr)<-rownames(Mut_BEM)
    colnames(NullMutExpr)<-rownames(Expr_BEM)
    check<-(IntMat>2)&(NullMat>2)&(MutNullExpr>2)&(NullMutExpr>2)
    dimnames(check)<-dimnames(NullMutExpr)

    if(sum(check)>0){
      #have interactions to test:
      #find mutations to test first:
      #number of interactions:
      output<-melt(check) %>% filter(value == 1)
      numberInt<-nrow(output)
      IntTerms<-data.frame(expression=as.character(output[,"Var2"]),mutation=as.character(output[,"Var1"]),stringsAsFactors = FALSE)


    }else{
      IntTerms<-NULL
    }



  }else{IntTerms<-NULL}
  return(list(Mut_BEM=Mut_BEM,Expr_BEM=Expr_BEM,IntTerms=IntTerms,IntTermsMM=IntTermsMM))
}

ANOVA_createInputFeaturesCompound<-function(additional_features=NULL,
                                               additional_features_only=FALSE,
                                               excludeHyperMetData=TRUE,
                                               oneFeatureOnly=NULL,TOTALBEM,
                                               selectedCellLines=NULL,MASTER_LIST,ANOVA.settings.featFactorPopulationTh=5,ParentBEM,MutualExclusive=FALSE){

  print('getting features start')
  if(additional_features_only){
    BEM<-additional_features
  }else{
    additional_features<-additional_features[,colnames(TOTALBEM)]
    BEM<-rbind(TOTALBEM,additional_features)
  }

  rownames(BEM)<-str_replace_all(rownames(BEM),pattern = ':',replacement = '_')
  rownames(BEM)<-str_replace_all(rownames(BEM),pattern = ' ',replacement = '_')
  rownames(ParentBEM)<-str_replace_all(rownames(ParentBEM),pattern = ':',replacement = '_')
  rownames(ParentBEM)<-str_replace_all(rownames(ParentBEM),pattern = ' ',replacement = '_')
  idxs<-which(rowSums(abs(BEM))>1)

  if (length(idxs)==1){
    BEM<-matrix(BEM,1,ncol(BEM),dimnames = list(rownames(BEM)[idxs],colnames(BEM)))
  }else{
    BEM<-BEM[idxs,]
  }

  if (length(selectedCellLines)){
    BEM<-BEM[,selectedCellLines]
    ParentBEM<-ParentBEM[,selectedCellLines]
  }


  idxsParent<-which(rowSums(abs(ParentBEM),na.rm=T)>ANOVA.settings.featFactorPopulationTh&rowSums(is.na(ParentBEM))==0)
  idxsChild<-which(rowSums(abs(BEM),na.rm=T)>ANOVA.settings.featFactorPopulationTh&rowSums(is.na(BEM))==0)
  print(length(idxsParent))
  print(length(idxsChild))

  if(length(idxsParent)>1&length(idxsChild)>1){
    BEM1<-ParentBEM[idxsParent,]
    BEM2<-BEM[idxsChild,]
    BEM1<-my.compress_identical_patterns(BEM1)
    BEM2<-my.compress_identical_patterns(BEM2)
    BEM1<-as.matrix(BEM1)
    BEM2<-as.matrix(BEM2)
    #may need to reconsider the following:
    #BEM1[is.na(BEM1)]<-0
    #BEM2[is.na(BEM2)]<-0
    Int_Mut<-BEM1%*%t(BEM2)
    rownames(Int_Mut)<-rownames(BEM1)
    colnames(Int_Mut)<-rownames(BEM2)

    NullBEM1<-(!BEM1)*1
    NullBEM2<-(!BEM2)*1

    WTmat<-NullBEM1%*%t(NullBEM2)
    rownames(WTmat)<-rownames(BEM1)
    colnames(WTmat)<-rownames(BEM2)
    Null1_2<-NullBEM1%*%t(BEM2)
    rownames(Null1_2)<-rownames(BEM1)
    colnames(Null1_2)<-rownames(BEM2)

    Null2_1<-BEM1%*%t(NullBEM2)
    rownames(Null2_1)<-rownames(BEM1)
    colnames(Null2_1)<-rownames(BEM2)


    if(MutualExclusive){
      #do OR check:
      #Mutual exclusive => not many where co-occur i.e. Int_Mut<2
      #Not all one or other otherwise cant test versus WT i.e. WTmat>2
      #Check that both individual are avail OR Check Null1_2>2 and Null2_1>2
      check<-(Int_Mut<2)&(WTmat>ANOVA.settings.featFactorPopulationTh)&(Null1_2>ANOVA.settings.featFactorPopulationTh)&(Null2_1>ANOVA.settings.featFactorPopulationTh)

      dimnames(check)<-dimnames(Int_Mut)
    }else{


      check<-(Int_Mut>ANOVA.settings.featFactorPopulationTh)&(WTmat>ANOVA.settings.featFactorPopulationTh)&(Null1_2>ANOVA.settings.featFactorPopulationTh)&(Null2_1>ANOVA.settings.featFactorPopulationTh)&(Int_Mut<ncol(BEM1)-ANOVA.settings.featFactorPopulationTh+1)&(WTmat<ncol(BEM1)-ANOVA.settings.featFactorPopulationTh+1)&(Null1_2<ncol(BEM1)-ANOVA.settings.featFactorPopulationTh+1)&(Null2_1<ncol(BEM1)-ANOVA.settings.featFactorPopulationTh+1)

      dimnames(check)<-dimnames(Int_Mut)

    }

    if(sum(check)>0){
      #have interactions to test:
      #find mutations to test first:
      #number of interactions:
      output<-reshape2::melt(check) %>% filter(value == 1)
      numberInt<-nrow(output)
      IntTerms<-data.frame(Child=as.character(output[,"Var2"]),Parent=as.character(output[,"Var1"]),stringsAsFactors = FALSE)
      check2<-IntTerms
      check2<-t(apply(check2,1,sort))
      IntTerms<-IntTerms[!duplicated(check2),]

    }else{
      IntTerms<-NULL
    }



  }else{IntTerms<-NULL}

  return(list(ParentBEM=BEM1,ChildBEM=BEM2,CompoundTerms=IntTerms))
}
decodeCNA_cp<-function(MoBEM){

  rn <- rownames(MoBEM)
  ii <- grep("cna", rownames(MoBEM))
  cnaId <- rownames(MoBEM)[ii]
  containedGenes <- unlist(lapply(str_split(cnaId, " "), function(x) {
    x[2]
  }))
  containedGenes[is.na(containedGenes)] <- ""
  segments <- unlist(lapply(str_split(unlist(lapply(str_split(cnaId,
                                                              ":"), function(x) {
                                                                x[2]
                                                              })), " "), function(x) {
                                                                x[1]
                                                              }))
  loci <- as.character(CNAdecode$locus[match(segments, CNAdecode$Identifier)])
  altType <- as.character(CNAdecode$Recurrent.Amplification..Deletion[match(segments,
                                                                            CNAdecode$Identifier)])
  altType[altType == "Amplification"] <- "G:"
  altType[altType == "Deletion"] <- "L:"
  rownames(MoBEM)[ii] <- paste(altType, loci, " ", containedGenes,
                               sep = "")
  return(MoBEM)
}


ANOVA_createInputFeaturesCellector<-function(additional_features=NULL,
                                            additional_features_only=FALSE,
                                            excludeHyperMetData=TRUE,
                                            oneFeatureOnly=NULL,TOTALBEM,
                                            selectedCellLines=NULL,MASTER_LIST,ANOVA.settings.featFactorPopulationTh=3,ParentBEM){

  print('getting features start')
  if(additional_features_only){
    BEM<-additional_features
  }else{
    additional_features<-additional_features[,colnames(TOTALBEM)]
    BEM<-rbind(TOTALBEM,additional_features)
  }
  if(length(grep("cna",rownames(BEM))>0&length(grep("cna",rownames(ParentBEM)))==0)){
    BEM<-decodeCNA_cp(BEM)
  }
  rownames(BEM)<-str_replace_all(rownames(BEM),pattern = ':',replacement = '_')
  rownames(BEM)<-str_replace_all(rownames(BEM),pattern = ' ',replacement = '')
  rownames(ParentBEM)<-str_replace_all(rownames(ParentBEM),pattern = ':',replacement = '_')
  rownames(ParentBEM)<-str_replace_all(rownames(ParentBEM),pattern = ' ',replacement = '')
  idxs<-which(rowSums(abs(BEM))>1)

  if (length(idxs)==1){
    BEM<-matrix(BEM,1,ncol(BEM),dimnames = list(rownames(BEM)[idxs],colnames(BEM)))
  }else{
    BEM<-BEM[idxs,]
  }

  if (length(selectedCellLines)){
    BEM<-BEM[,selectedCellLines]
    ParentBEM<-ParentBEM[,selectedCellLines]
  }


  idxsParent<-which(rowSums(abs(ParentBEM))>ANOVA.settings.featFactorPopulationTh)
  idxsChild<-which(rowSums(abs(BEM))>ANOVA.settings.featFactorPopulationTh)

  CellectorTerms<-list()

  sBEM<-NULL
  cBEM<-NULL
  newnamesS<-c()
  newnamesC<-c()
  if(length(idxsParent)>1&length(idxsChild)>1){
    BEM1<-ParentBEM[idxsParent,]
    BEM2<-BEM[idxsChild,]
    #BEM1<-my.compress_identical_patterns(BEM1)
    #BEM2<-my.compress_identical_patterns(BEM2)
    singleterms<-rownames(BEM2)
    for(k in 1:length(singleterms)){
      splitC<-strsplit(singleterms[k],"[()]")
      if(length(splitC[[1]])>1){
        #have a cna marker with associated genes
        #find things that are not full markers:
        cnaList<-!grepl("mut|Expr|cna",splitC[[1]])
        splitC[[1]][cnaList]<-unlist(sapply(gsub(",",";",splitC[[1]][cnaList]),function(x) paste0(paste0("(",x),")")))

        singleterms[k]<-paste0(splitC[[1]],collapse="")
      }
    }


    j=1
    for(i in 1:length(idxsParent)){
      cTerm<-rownames(BEM1)[i]

      splitC<-strsplit(cTerm,"[()]")[[1]]
      if(length(splitC)>1){
        #have a cna marker with associated genes
        #find things that are not full markers:
        cnaList<-!grepl("mut|Expr",splitC)
        splitC[cnaList]<-unlist(sapply(gsub(",",";",splitC[cnaList]),function(x) paste0(paste0("(",x),")")))
        splitC<-paste0(splitC,collapse="")
      }

      constituentparts<-sapply(unlist(splitC),function(x) gsub(" ","",strsplit(x,",",fixed=TRUE)[[1]]))

      #check to see what if any model can be run, and compared as nested model.
      parts<-unlist(sapply(unlist(constituentparts),function(x) gsub("~","",x,fixed=TRUE)))


      existsingle<-parts[parts%in%singleterms]

      orig<-unlist(constituentparts)[parts%in%singleterms]

      notterms<-grepl("\\~",orig)

      if(length(existsingle)>0){
        SingleBEM<-BEM2[existsingle,]
        if(sum(notterms)>0){
          if(length(existsingle)>1){
            SingleBEM[notterms,]<-(!SingleBEM[notterms,])+0

          }else{

            SingleBEM[notterms]<-matrix((!SingleBEM[notterms])+0,nrow=1)

          }

        }

        if(is.matrix(SingleBEM)){
        SingleBEM<-my.compress_identical_patterns(SingleBEM)
        checkmat<-rbind(SingleBEM,BEM1[i,])
        checkmat<-as.matrix(checkmat)
        rownames(checkmat)<-c(orig,cTerm)


        checkBEM<-my.compress_identical_patterns(checkmat)
          if(nrow(SingleBEM)!=nrow(checkBEM)){
            #have a nested model we can check:
            SingleBEM<-checkmat[1:(nrow(checkmat)-1),]
            if(is.null(sBEM)){

              sBEM<-SingleBEM
            }else{
              currnames<-rownames(sBEM)
              newterms<-rownames(SingleBEM)[!rownames(SingleBEM)%in%currnames]
              if(length(newterms)>0){
              newBEM<-matrix(SingleBEM[newterms,],nrow=length(newterms))
              sBEM<-rbind(sBEM,newBEM)
              newnamesS<-c(currnames,newterms)
              rownames(sBEM)<-newnamesS}
            }
            fitT<-paste(rownames(SingleBEM),collapse="+")
            fitAll<-paste(c(rownames(SingleBEM),cTerm),collapse = "+")
            if(is.null(cBEM)){
              cBEM<-matrix(BEM1[i,],nrow=1)
              colnames(cBEM)<-colnames(BEM1)
              rownames(cBEM)<-cTerm
              newnamesC<-cTerm
            }else{
              newnamesC<-rownames(cBEM)
              cBEM<-rbind(cBEM,BEM1[i,])
              newnamesC<-c(newnamesC,cTerm)
              rownames(cBEM)<-newnamesC
            }


            CellectorTerms[[j]]<-list(fitT=fitT,fitAll=fitAll,fitC=cTerm)

            j=j+1
          }
        }

      }

    }

  }else{FitTerms<-NULL}

  return(list(sBEM=sBEM,cBEM=cBEM,FitTerms=CellectorTerms))
}
getMarkerSummary<-function(markers){
  allmarkers<-NULL
  for(i in 1:length(markers)){
    splitmarkers<-unlist(strsplit(markers[i],",",fixed=TRUE))
    expr<-grep("_Expr",splitmarkers)
    mut<-grep("_mut",splitmarkers)
    cna<-grep("gain",splitmarkers)
    cna2<-grep("loss",splitmarkers)
    allmarkers<-rbind(allmarkers,c(expr=length(expr),mut=length(mut),cna=length(cna)+length(cna2)))
  }
  rownames(allmarkers)<-markers
  return(allmarkers)
}

getMarkerType<-function(markersummary){
  markercount<-colSums(markersummary)
  names(markercount)<-c("expr","mut","cna")
  markerexist<-markercount!=0
  if(markerexist["expr"]){
    if(markerexist["mut"]){
      if(markerexist["cna"]){
        markertype<-c("Expr,Mutation,CNA")
      }else{
        markertype<-c("Expr,Mutation")
      }
    }else{
      if(markerexist["cna"]){
        markertype<-c("Expr,CNA")
      }else{markertype<-c("Expr")}
    }
  }else{
    if(markerexist["mut"]){
      if(markerexist["cna"]){
        markertype<-c("Mutation,CNA")
      }else{
        markertype<-c("Mutation")
      }

    }else{
      markertype<-c("CNA")
    }
  }
  return(markertype)
}

