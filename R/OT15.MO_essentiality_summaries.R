
retrievePrioriKnownEss<-function(essData,coreFitnessGenes,panCancer_cf_genes,dir.ExternalData,nameMat=NULL,cancerDrivers=NULL,refSets){

    if(is.null(coreFitnessGenes)){
      cs_coreFitnessGenes<-NULL
    }else{

      cs_coreFitnessGenes<-list()
      cs_coreFitnessGenes[[1]]<-setdiff(coreFitnessGenes,cancerDrivers[[1]])

    }
    if(is.null(nameMat)){
      genes<-rownames(essData)
      inBagel_ess_ref<-is.element(genes,refSets[["curated_BAGEL_essential"]])
      inBagel_noness_ref<-is.element(genes,refSets[["BAGEL_nonEssential"]])


      inBagel_ess_refv2<-is.element(genes,refSets[["BAGEL_essential"]])

      ribosomalProteins<-is.element(genes,refSets[["EssGenes.ribosomalProteins"]])
      dna_replication<-is.element(genes,refSets[["EssGenes.DNA_REPLICATION_cons"]])
      proteasome<-is.element(genes,refSets[["EssGenes.PROTEASOME_cons"]])
      spliceosome<-is.element(genes,refSets[["EssGenes.SPLICEOSOME_cons"]])
      histones<-is.element(genes,refSets[["histones"]])
      rna_polymerase<-is.element(genes,refSets[["EssGenes.KEGG_rna_polymerase"]])

      predictedAsPanCancerCoreFitness<-is.element(genes,panCancer_cf_genes)

      predictedAsCS_coreFitness<-NULL
      if(!is.null(cs_coreFitnessGenes)){
        for (i in 1:length(cs_coreFitnessGenes)){
          predictedAsCS_coreFitness<-cbind(predictedAsCS_coreFitness,is.element(genes,cs_coreFitnessGenes[[i]]))
        }


      }else{
        predictedAsCS_coreFitness<-matrix(0,nrow=length(genes),ncol=1)
        colnames(predictedAsCS_coreFitness)<-"None"
      }

      res<-data.frame(cbind(inBagel_ess_ref,
               inBagel_ess_refv2,
               inBagel_noness_ref,
               ribosomalProteins,
               dna_replication,
               proteasome,
               spliceosome,
               histones,
               rna_polymerase,
               predictedAsPanCancerCoreFitness,
               predictedAsCS_coreFitness),stringsAsFactors = FALSE)


      rownames(res)<-rownames(essData)
    }else{
      inBagel_ess_ref<-unlist(lapply(nameMat,function(x) sum(is.element(x,refSets[["curated_BAGEL_essential"]]))>0))
      inBagel_noness_ref<-unlist(lapply(nameMat,function(x) sum(is.element(x,refSets[["BAGEL_nonEssential"]]))>0))


      inBagel_ess_refv2<-unlist(lapply(nameMat,function(x) sum(is.element(x,refSets[["BAGEL_essential"]]))>0))

      ribosomalProteins<-unlist(lapply(nameMat,function(x) sum(is.element(x,refSets[["EssGenes.ribosomalProteins"]]))>0))
      dna_replication<-unlist(lapply(nameMat,function(x) sum(is.element(x,refSets[["EssGenes.DNA_REPLICATION_cons"]]))>0))
      proteasome<-unlist(lapply(nameMat,function(x) sum(is.element(x,refSets[["EssGenes.PROTEASOME_cons"]]))>0))
      spliceosome<-unlist(lapply(nameMat,function(x) sum(is.element(x,refSets[["EssGenes.SPLICEOSOME_cons"]]))>0))
      histones<-unlist(lapply(nameMat,function(x) sum(is.element(x,refSets[["histones"]]))>0))
      rna_polymerase<-unlist(lapply(nameMat,function(x) sum(is.element(x,refSets[["EssGenes.KEGG_rna_polymerase"]]))>0))

      predictedAsPanCancerCoreFitness<-unlist(lapply(nameMat,function(x) sum(is.element(x,panCancer_cf_genes))>0))
      predictedAsCS_coreFitness<-NULL
      for (i in 1:length(cs_coreFitnessGenes)){
        predictedAsCS_coreFitness<-cbind(predictedAsCS_coreFitness,unlist(lapply(nameMat,function(x) sum(is.element(x,cs_coreFitnessGenes[[i]])))))
      }

      colnames(predictedAsCS_coreFitness)<-names(cs_coreFitnessGenes)
      res<-data.frame(cbind(inBagel_ess_ref,
                 inBagel_ess_refv2,
                 inBagel_noness_ref,
                 ribosomalProteins,
                 dna_replication,
                 proteasome,
                 spliceosome,
                 histones,
                 rna_polymerase,
                 predictedAsPanCancerCoreFitness,
                 predictedAsCS_coreFitness),stringsAsFactors = FALSE)


      rownames(res)<-rownames(nameMat)
    }
    return(res)
}
retrieveHCGsStatus<-function(essData,ctype,nameMat=NULL,cancerDrivers){
    CTYPE<-ctype
    cMap<-read.csv(paste0(dir.ExternalData,"/cTypeMap.csv"),stringsAsFactors=FALSE,header=F)
    m1<-cMap[match(tolower(CTYPE),tolower(cMap[,1])),3]
    m2<-cMap[match(tolower(CTYPE),tolower(cMap[,2])),3]
    ctype<-c(m1,m2)
    ctype<-ctype[!is.na(ctype)]


    cgs<-unique(cancerDrivers[[1]])

    if(!is.null(nameMat)){
      HCGs<-rep(FALSE,length(nameMat))
      HCGs<-unlist(lapply(nameMat,function(x) sum(is.element(x,cgs))>0))
      names(HCGs)<-names(nameMat)
      if(sum(is.element(ctype,cancerDrivers$Tumor_Type))>0){
        specific_HCGs<-rep(FALSE,length(nameMat))
        specific_HCGs<-unlist(lapply(nameMat,function(x) sum(is.element(x,cancerDrivers$ActingDriver_Symbol[which(is.element(cancerDrivers$Tumor_Type,ctype))]))>0))
        names(specific_HCGs)<-names(nameMat)
      }else{
        specific_HCGs<-rep(NA,nrow(essData))
        names(HCGs)<-rownames(essData)
      }
    }else{
      HCGs<-rep(FALSE,nrow(essData))
      HCGs<-is.element(rownames(essData),cgs)
      names(HCGs)<-rownames(essData)

      if(sum(is.element(ctype,cancerDrivers$Tumor_Type))>0){
        specific_HCGs<-rep(FALSE,nrow(essData))
        specific_HCGs<-is.element(rownames(essData),cancerDrivers$ActingDriver_Symbol[which(is.element(cancerDrivers$Tumor_Type,ctype))])
        names(specific_HCGs)<-rownames(essData)
      }else{
        specific_HCGs<-rep(NA,nrow(essData))
        names(HCGs)<-rownames(essData)
      }
    }
    res<-cbind(HCGs+0,specific_HCGs+0)
    colnames(res)<-c('HCG','CS_HCG')
    return(res)

}

retrievePrimaryTumourStatus<-function(essData,ctype,nameMat=NULL,dir.ExternalData){
    if(!dir.exists(paste0(dir.ExternalData,"/MultiOmic_BEMs/TUMOURS/BEMs/"))){
      stop("MultiOmic BEMs directory doesnt exist for Primary Tumour status")
    }
    CTYPE<-ctype

    cMap<-read.csv(paste0(dir.ExternalData,"/cTypeMap.csv"),stringsAsFactors=FALSE,header=F)
    m1<-cMap[match(tolower(CTYPE),make.names(tolower(cMap[,1]))),3]
    m2<-cMap[match(tolower(CTYPE),tolower(cMap[,2])),3]
    ctype<-c(m1,m2)
    ctype<-ctype[!is.na(ctype)]
    ctype<-unique(ctype)

    if (length(ctype)==1){
        load(paste(dir.ExternalData,'/MultiOmic_BEMs/TUMOURS/BEMs/',ctype,'.rdata',sep=''))
        if(!is.null(nameMat)){
          MutPrimTum<-unlist(lapply(nameMat,function(x) sum(BEM$mutFreqs[x])!=0+0))
          names(MutPrimTum)<-names(nameMat)

          MutPrimTum_Cosmic<-unlist(lapply(nameMat,function(x) sum(BEM$mutFreqs_COSMIC_filtered[x])!=0+0))
          names(MutPrimTum_Cosmic)<-names(nameMat)
        }else{
          MutPrimTum<-BEM$mutFreqs[rownames(essData)]
          names(MutPrimTum)<-rownames(essData)

          MutPrimTum_Cosmic<-BEM$mutFreqs_COSMIC_filtered[rownames(essData)]
          names(MutPrimTum_Cosmic)<-rownames(essData)
        }


    }else{
        tmpMutFreq<-NULL
        tmpMutFreqCosmic<-NULL
        MutPrimTum_Cosmic<-NULL
        for (k in 1:length(ctype)){
            print(k)
            load(paste(dir.ExternalData,'/MultiOmic_BEMs/TUMOURS/BEMs/',ctype[k],'.rdata',sep=''))
            tmpMutFreq<-cbind(tmpMutFreq,BEM$mutFreqs)
            tmpMutFreqCosmic<-cbind(tmpMutFreqCosmic,BEM$mutFreqs_COSMIC_filtered)
        }

        MutPrimTum<-apply(tmpMutFreq,MARGIN = 1,'max',na.rm=TRUE)
        names(MutPrimTum)<-rownames(essData)

        MutPrimTum_Cosmic<-apply(tmpMutFreqCosmic,MARGIN = 1,'max',na.rm=TRUE)
        names(MutPrimTum_Cosmic)<-rownames(essData)
    }

    MutPrimTum<-cbind(MutPrimTum,MutPrimTum_Cosmic)
    colnames(MutPrimTum)<-c('MutPrimTum','MutPrimTum_Cosmic')

    MutPrimTum[is.na(MutPrimTum)]<-0

    return(MutPrimTum)
}

#############################
retrieveEssData<-function(cellLine){
    mageck_pvals<-
        mgk_pvals_merged
    mageck_fdrs<-
        mgk_fdrs_merged

    tmp<-unlist(str_split(colnames(mageck_pvals),'_'))
    tmp<-tmp[seq(1,length(tmp),2)]
    colnames(mageck_pvals)<-tmp
    colnames(mageck_fdrs)<-tmp

    load(paste(bagelR_profiles.dir,cellLine,'_GeneLevelBF.rdata',sep=''))
    load(paste(bagelR_profiles.dir,cellLine,'_sgRNAs_BFs.Rdata',sep=''))

    #     nreps<-length(sgRNA_BFs_across_reps$A2058_c903R1)
    #     guides<-sort(rownames(sgRNA_BFs_across_reps[[1]]$sgRNA_BFs_inclTr))
    #
    #     for (i in 1:nreps){
    #         if (i == 1){
    #             sgBF<-sgRNA_BFs_across_reps[[1]]$sgRNA_BFs_inclTr[guides,1]
    #         }else{
    #             sgBF<-cbind(sgBF,sgRNA_BFs_across_reps[[1]]$sgRNA_BFs_inclTr[guides,1])
    #         }
    #     }
    #
    #     sgBF<-rowMeans(sgBF)

    mgk_pvals_profile<-
        mageck_pvals[,cellLine]
    mgk_fdrs_profile<-
        mageck_fdrs[,cellLine]
    bag_bfs_profile<-
        GeneLevelBFs[,1]
    bag_bfs_sd_profile<-
        GeneLevelBFs[,2]

    geneNames<-sort(unique(names(mgk_pvals_profile)))

    mgk_pvals_profile<-
        mgk_pvals_profile[geneNames]
    mgk_fdrs_profile<-100*
        mgk_fdrs_profile[geneNames]
    bag_bfs_sd_profile<-
        bag_bfs_sd_profile[geneNames]
    bag_bfs_profile<-
        bag_bfs_profile[geneNames]

    sigDeletedAccToMageck<-
        mgk_fdrs_profile<10
    sigDeletedAccToBagel<-
        bag_bfs_profile>BAGEL_BEST_PPV_th[cellLine]

    #     ug<-names(sigDeletedAccToBagel)
    #     nug<-length(ug)
    #     tmp<-lapply(ug,geneToGuideBF,sgBF)
    #     sguideBFs<-t(matrix(unlist(tmp),5,nug))
    #
    essData<-cbind(sigDeletedAccToBagel==1,
                   bag_bfs_profile,
                   bag_bfs_sd_profile,
                   sigDeletedAccToMageck==1,
                   mgk_pvals_profile,
                   mgk_fdrs_profile)

    rownames(essData)<-
        geneNames
    essData<-essData[order(essData[,'bag_bfs_profile'],decreasing=TRUE),]
    return(essData)
}

geneToGuideBF<-function(GENE,sgBF){
    sgbf<-sgBF[guideList[[GENE]]]
    sgbf<-sort(sgbf,decreasing=TRUE)
    sgbf<-sgbf[1:5]
    return(sgbf)
}



retrieveExpression<-function(essData,cellLine){

    commong<-intersect(rownames(EXPpc),rownames(essData))

    EXPpc<-EXPpc[commong,]

    cid<-
        as.character(PS_inventory$COSMIC_ID[match(cellLine,PS_inventory$Cell_Line_Name)])

    notExpressedProfile<-rep(NA,nrow(essData))
    names(notExpressedProfile)<-rownames(essData)
    higlyExpressedProfile<-rep(NA,nrow(essData))
    names(higlyExpressedProfile)<-rownames(essData)
    exp_fpkm<-rep(NA,nrow(essData))
    names(exp_fpkm)<-rownames(essData)

    if(is.element(cid,colnames(EXPpc))){
        expProf<-EXPpc[,cid]
        notExpGenes<-expProf<0.5
        notExpressedProfile[names(notExpGenes)]<-notExpGenes

        extrQuantile<-quantile(expProf,probs = 0.95)
        highExpGenes<-expProf>=extrQuantile
        higlyExpressedProfile[names(highExpGenes)]<-highExpGenes

        exp_fpkm[names(highExpGenes)]<-expProf
    }

    res<-cbind(exp_fpkm,higlyExpressedProfile+0,notExpressedProfile+0)
    return(res)
}
retrieveGenomic<-function(essData,cellLine,PS_inventory,cellLineID="Cell_Line_Name",inventoryID="COSMIC_ID",geneList=NULL,variants_catalogue, variant_types_to_exclude=c("Substitution - coding silent",
                                                                                                                                                                             "Unknown",
                                                                                                                                                                             "microRNA - Substitution",
                                                                                                                                                                             "HomDel",
                                                                                                                                                                             "microRNA - Deletion",
                                                                                                                                                                             "microRNA - Insertion")){



  checkCols<-sum(!c("DESCRIPTION","GENE_NAME","ZYGOSITY","COSMIC_Filter","AA_MUT_SYNTAX")%in%colnames(variants_catalogue))
  if(checkCols!=0){
    stop(paste("Not all required fields in inventory, missing columns:",setdiff(c("DESCRIPTION","GENE_NAME","ZYGOSITY","COSMIC_Filter","AA_MUT_SYNTAX"),colnames(variants_catalogue))))
  }
  checkCols<-sum(!cellLineID%in%colnames(PS_inventory))
  if(checkCols!=0){
    stop(paste("Cell line ID",cellLineID,"not present in inventory file"))
  }

  PS_inventory[,cellLineID]<-make.names(PS_inventory[,cellLineID])

   cid<-
        as.character(PS_inventory[match(cellLine,PS_inventory[,cellLineID]),inventoryID])



    reduced_variants<-variants_catalogue[which(variants_catalogue[,inventoryID]==cid &
                                                   !is.element(variants_catalogue$DESCRIPTION,variant_types_to_exclude)),]

    variantsProfile<-rep('wt',nrow(essData))
    names(variantsProfile)<-rownames(essData)

    variantsDescription<-rep('',nrow(essData))
    names(variantsDescription)<-rownames(essData)

    zigosity<-rep('',nrow(essData))
    names(zigosity)<-rownames(essData)

    cosmic_filter<-rep('',nrow(essData))
    names(cosmic_filter)<-rownames(essData)

    syntax<-rep('',nrow(essData))
    names(syntax)<-rownames(essData)

    genesToConsider<-intersect(reduced_variants$GENE_NAME,rownames(essData))
    ngenes<-length(genesToConsider)
    if(ngenes<500){
      message(paste("Number of genes to get Genomic information only",ngenes))
    }
    reduced_variants<-reduced_variants[is.element(reduced_variants$GENE_NAME,
                                                  genesToConsider),]
    for (i in 1:ngenes){
        id<-which(reduced_variants$GENE_NAME==genesToConsider[i])

        if (length(id)>0){
            variantsProfile[genesToConsider[i]]<-'mut'
            variantsDescription[genesToConsider[i]]<-paste(reduced_variants$DESCRIPTION[id],collapse=' \\\ ')
            zigosity[genesToConsider[i]]<-paste(reduced_variants$ZYGOSITY[id],collapse=' \\\ ')
            cosmic_filter[genesToConsider[i]]<-paste(reduced_variants$COSMIC_Filter[id]=='y',collapse=' \\\ ')
            syntax[genesToConsider[i]]<-paste(reduced_variants$AA_MUT_SYNTAX[id],collapse=' \\\ ')
        }
    }
    if(!is.null(geneList)){
      variantsProfile2<-c()
      variantsDescription2<-c()
      zigosity2<-c()
      cosmic_filter2<-c()
      syntax2<-c()
      for(i in 1:length(geneList)){
        id<-which(reduced_variants$GENE_NAME%in%geneList[[i]])

        if (length(id)>0){
          variantsProfile2[i]<-'mut'
          variantsDescription2[i]<-paste(reduced_variants$DESCRIPTION[id],collapse=' \\\ ')
          zigosity2[i]<-paste(reduced_variants$ZYGOSITY[id],collapse=' \\\ ')
          cosmic_filter2[i]<-paste(reduced_variants$COSMIC_Filter[id]=='y',collapse=' \\\ ')
          syntax2[i]<-paste(reduced_variants$AA_MUT_SYNTAX[id],collapse=' \\\ ')
        }else{
          variantsProfile2[i]<-'wt'
          variantsDescription2[i]<-''
          zigosity2[i]<-''
          cosmic_filter2[i]<-''
          syntax2[i]<-''
        }
      }
      names(variantsProfile2)<-names(geneList)
      names(variantsDescription2)<-names(geneList)
      names(zigosity2)<-names(geneList)
      names(cosmic_filter2)<-names(geneList)
      names(syntax2)<-names(geneList)
      variantsProfile<-c(variantsProfile,variantsProfile2)
      variantsDescription<-c(variantsDescription,variantsDescription2)
      zigosity<-c(zigosity,zigosity2)
      cosmic_filter<-c(cosmic_filter,cosmic_filter2)
      syntax<-c(syntax,syntax2)
    }

    res<-cbind(variantsProfile,
               variantsDescription,
               zigosity,
               cosmic_filter,
               syntax)
    return(res)

}

retrieveCNstatus<-function(essData,cellLine,PS_inventory,cellLineID="Cell_Line_Name",geneList=NULL,inventoryID="COSMIC_ID",geneLevCNA,logr=FALSE){

    commong<-intersect(rownames(geneLevCNA$CNV),rownames(essData))
    if(length(commong)<200){
      message("Number of genes with CN status annotation is less than 200")
    }
    if(sum(!c(cellLineID,inventoryID)%in%colnames(PS_inventory))){
      stop("The cell line ID and inventory ID for cell line matching are not in the inventory file")
    }
    cnas<-geneLevCNA$CNV[commong,]
    cid<-
        as.character(PS_inventory[match(cellLine,make.names(PS_inventory[,cellLineID])),inventoryID])

    cna_value<-rep(' , , , ',nrow(essData))
    names(cna_value)<-rownames(essData)
    cna_status<-rep('',nrow(essData))
    names(cna_status)<-rownames(essData)

    if(is.element(cid,colnames(cnas))){
        cna_value[commong]<-cnas[,cid]

        ntmp<-unlist(str_split(cna_value,','))

        max<-as.numeric(ntmp[seq(1,length(ntmp),4)])
        min<-as.numeric(ntmp[seq(2,length(ntmp),4)])
        HH<-ntmp[seq(3,length(ntmp),4)]
        if(logr){
          cna_status[which(min>=2.3&max>=2.3)]<-'Amplified'
          cna_status[which(min==0 | max== 0)]<-"HomDel"
        }else{
        cna_status[which(min>=8 & max>=8)]<-'Amplified'
        cna_status[which(min==0 | max==0)]<-'HomDel'}
    }

    if(!is.null(geneList)){
      cna_value2<-rep(' , , , ',length(geneList))
      names(cna_value2)<-names(geneList)
      cna_status2<-rep('',length(geneList))
      names(cna_status2)<-names(geneList)


      cnav<-lapply(geneList,function(x) paste(cna_value[intersect(names(cna_value),x)],collapse="//"))
      cnav<-unlist(cnav)
      cnas2<-lapply(geneList,function(x) combineCNstatus(x,cna_status))
      cnas2<-unlist(cnas2)
      cna_value<-c(cna_value,cna_value2)
      cna_status<-c(cna_status,cnas2)
    }

    res<-data.frame(cna_value,cna_status)
    return(res)
}
combineCNstatus<-function(x, geneCNs,priority=c("Amplified","HomDel")){
  temp<-geneCNs[intersect(names(geneCNs),x)]
  out<-""
  for(i in 1:length(priority)){
    check<-grep(priority[i],temp)
    if(length(check)>0){out<-priority[i]}
  }
  return(out)
}
updateMatrix<-function(inputmatrix,nameList,method=c("OR")){
  method<-match.arg(method)
  if(!is.null(nameList)){
  if(method=="OR"){
    newdata<-matrix(unlist(lapply(nameList,function(x) colSums(inputmatrix[x,],na.rm = T)>0)),ncol=ncol(inputmatrix),byrow=TRUE)
    rownames(newdata)<-names(nameList)
  }
  newmatrix<-rbind(inputmatrix,newdata)}else{
    newmatrix<-inputmatrix
  }
  return(newmatrix)
}
retrieveNguides<-function(essData){
    curr_lib_nGuides<-
        nGuides_old_lib[rownames(essData)]

    new_lib_nGuides<-
        nGuides_new_lib[rownames(essData)]

    res<-cbind(curr_lib_nGuides,new_lib_nGuides)
    return(res)
}
retrieveKinases_iGof_GSKtar<-function(essData){

    load(kinases.fn)
    genes<-rownames(essData)
    uniprotKinase<-is.element(genes,kinase_genes_and_syn)
    names(uniprotKinase)<-rownames(essData)

    load(humKin.fn)
    genes<-rownames(essData)
    inHumanKinome<-is.element(genes,humKinome)
    names(inHumanKinome)<-rownames(essData)


    uniprotKinase<-rep(FALSE,nrow(essData))
    names(uniprotKinase)<-rownames(essData)
    uniprotKinase<-is.element(names(uniprotKinase),kinase_genes_and_syn)
    names(uniprotKinase)<-rownames(essData)

    load(intoGen_GoF.fn)

    igGOF<-rep(FALSE,nrow(essData))
    names(igGOF)<-rownames(essData)
    igGOF<-is.element(names(igGOF),intoGen_GoF)
    names(igGOF)<-rownames(essData)

    load(GSK_epiTargets.fn)
    GSK_EPIG_TAR<-rep(FALSE,nrow(essData))
    names(GSK_EPIG_TAR)<-rownames(essData)
    GSK_EPIG_TAR<-is.element(names(GSK_EPIG_TAR),GSK_epiTargets)
    names(GSK_EPIG_TAR)<-rownames(essData)

    res<-cbind(inHumanKinome,uniprotKinase,igGOF,GSK_EPIG_TAR)
    return(res)
}
retrievegeneInfo<-function(essData){

    fc<-read.table('../../DATA/Raw/protein-coding_gene (1).txt',sep = '\t',header=TRUE,stringsAsFactors = FALSE)
    rownames(fc)<-fc$symbol

    commong<-intersect(rownames(fc),rownames(essData))

    fc<-fc[commong,]

    description<-rep('',nrow(essData))
    names(description)<-rownames(essData)

    hgnc_id<-rep('',nrow(essData))
    names(hgnc_id)<-rownames(essData)

    entrez_id<-rep('',nrow(essData))
    names(entrez_id)<-rownames(essData)

    ensemble_id<-rep('',nrow(essData))
    names(ensemble_id)<-rownames(essData)

    location<-rep('',nrow(essData))
    names(location)<-rownames(essData)

    family<-rep('',nrow(essData))
    names(family)<-rownames(essData)

    pubmed_id<-rep('',nrow(essData))
    names(pubmed_id)<-rownames(essData)

    description[commong]<-fc[,'name']
    family[commong]<-fc[,'gene_family']
    hgnc_id[commong]<-fc[,'hgnc_id']
    entrez_id[commong]<-fc[,'entrez_id']
    ensemble_id[commong]<-fc[,'ensembl_gene_id']
    location[commong]<-fc[,'location_sortable']
    pubmed_id[commong]<-fc[,'pubmed_id']

    res<-cbind(description,family,hgnc_id,entrez_id,ensemble_id,location,pubmed_id)

    return(res)
}

retrieveTractabilityData<-function(essData){
    res<-tractability[rownames(essData),]
    rownames(res)<-rownames(essData)
    return(res)
}

retrieveGDSCtractability<-function(essData){
    load('../../../../Desktop/R_MAIN_PAPER_4.0/Fi_gdsc1000_DATA/ANNOTATIONS/DRUGS/atomDT.rdata')

    gg<-rownames(essData)

    MM<-match(gg,colnames(atomDT))
    id<-which(!is.na(MM))

    MM<-MM[id]

    gdscTractability<-rep(NA,nrow(essData))

    for (i in 1:length(MM)){
        gdscTractability[id[i]]<-
            paste(sort(names(which(atomDT[,MM[i]]>0))),collapse=', ')
    }
    names(gdscTractability)<-rownames(essData)

    return(gdscTractability)
}

retrieveInEssentialPathway<-function(essData,cellLine,pathwayObject){

    currentDir<-paste(mageck_profiles.dir,cellLine,'/',sep='')

    fns<-dir(currentDir)

    fn<-paste(currentDir,grep('.gene_summary.txt_ptw.pathway_summary.txt',fns,value=TRUE),sep='')

    fc<-read.table(fn,sep='\t',header=TRUE,stringsAsFactors = FALSE)

    depletedPathways<-
        fc$id[which(fc$neg.fdr<0.10)]

    genesInDepletedPathways<-
        unique(unlist(pathwayObject[depletedPathways]))

    inPath<-is.element(rownames(essData),genesInDepletedPathways)

    names(inPath)<-rownames(essData)

    return(inPath)
}



diffExpressedGenesBasedOnMarker<-function(marker,display=TRUE){

    posSamples<-names(which(MoBEM[marker,]>0))
    negSamples<-names(which(MoBEM[marker,]==0))

    posSamples<-intersect(posSamples,PS_inventory$COSMIC_ID)
    negSamples<-intersect(negSamples,PS_inventory$COSMIC_ID)

    posSamples<-intersect(posSamples,colnames(d_basal_expression))
    negSamples<-intersect(negSamples,colnames(d_basal_expression))

    GEX<-
        d_basal_expression[,c(negSamples,posSamples)]

    GEX<-GEX[which(apply(GEX,MARGIN = 1,'sd')>1),]


    RTRES<-rowFtests(GEX,as.factor(c(rep(0,length(negSamples)),rep(1,length(posSamples)))),var.equal=FALSE)

    meanPosPop<-
        rowMeans(GEX[,posSamples])
    meanNegPop<-
        rowMeans(GEX[,negSamples])
    logFC<-log2(meanPosPop/meanNegPop)


    RTRES<-cbind(RTRES,p.adjust(RTRES$p.value,'fdr'))
    RTRES<-cbind(RTRES,qvalue(RTRES$p.value)$qvalue)
    colnames(RTRES)[ncol(RTRES)-1]<-'fdr'
    colnames(RTRES)[ncol(RTRES)]<-'qvalue'
    RTRES<-cbind(RTRES,meanPosPop,meanNegPop,logFC)

    RTRES<-RTRES[order(RTRES$p.value),]

    DEXgenes<-rownames(RTRES)[which(RTRES$qvalue<0.10 & RTRES$logFC>0)]

    if(length(DEXgenes)>0){
        if(display){
            fdrTh<-max(RTRES$p.value[which(RTRES$qvalue<0.10)])

            GRAY<-rgb(150,150,150,alpha = 180,maxColorValue = 255)
            PURPLE<-rgb(160,32,240,alpha = 180,maxColorValue = 255)
            ORANGE<-rgb(255,165,0,alpha = 180,maxColorValue = 255)

            COL<-rep(GRAY,nrow(RTRES))
            COL[which(RTRES$p.value<=fdrTh & RTRES$logFC>0.5)]<-ORANGE
            COL[which(RTRES$p.value<=fdrTh & RTRES$logFC< -0.5)]<-PURPLE

            plot(RTRES$logFC,-log10(RTRES$p.value),pch=21,bg=COL,col='white',xlab='logFC',ylab='-log10(p-value)',
                 main=paste(CTYPE,': ',marker,' vs WT',sep=''))
            abline(v=0,col='darkgray')
            #abline(v=0.5,col='black',lty=2)
            #abline(v= -0.5,col='black',lty=2)
            abline(h= -log10(fdrTh),col='black',lty=3)
        }
        return(list(DEtable=RTRES,UpRegulatedGenes=DEXgenes))
    }else{
        return(NULL)
    }
}

pathEnrichAnalysis<-function(gs){

    BG<-PATHCOM_HUMAN$backGround
    Npath<-length(PATHCOM_HUMAN$PATHWAY)


    k<-length(gs)

    flag<-1

    containedGenes<-vector()
    N_<-vector()
    n_<-vector()
    k_<-vector()
    x_<-vector()

    for (i in 1:Npath){
        N<-length(BG)
        n<-length(PATHCOM_HUMAN$HGNC_SYMBOL[[i]])
        cg<-intersect(gs,PATHCOM_HUMAN$HGNC_SYMBOL[[i]])
        x<-length(cg)

        if(x>=2){
            currentP<-my.hypTest(x,k,n,N)
            if(flag==1){
                RESTOT<-matrix(c(i,currentP),1,2)
            }else{
                RESTOT<-rbind(RESTOT,c(i,currentP))
            }
            containedGenes[flag]<-paste(cg,collapse=', ')
            N_[flag]<-N
            n_[flag]<-n
            k_[flag]<-k
            x_[flag]<-x
            flag<-flag+1
        }
    }


    if(exists('RESTOT')){
        RESTOT<-cbind(RESTOT,p.adjust(RESTOT[,2],'fdr'))
        colnames(RESTOT)<-c('pathId','p','adj.p')

        RESTOT<-data.frame(cbind(RESTOT,x_,k_,n_,N_,containedGenes,PATHCOM_HUMAN$PATHWAY[RESTOT[,1]]),stringsAsFactors = FALSE)

        ENRpath<-RESTOT$pathId[which(RESTOT$adj.p<0.05)]

        return(list(pathTAB=RESTOT,ENRpath=ENRpath))
    }else{
        return(NULL)
    }
}
familyEnrichAnalysis<-function(gs){

    BG<-GENE_FAMILIES$backGround
    Npath<-length(GENE_FAMILIES$Families)

    k<-length(gs)

    flag<-1

    containedGenes<-vector()
    N_<-vector()
    n_<-vector()
    k_<-vector()
    x_<-vector()


    for (i in 1:Npath){
        N<-length(BG)
        n<-length(GENE_FAMILIES$HGNC_SYMBOL[[i]])
        cg<-intersect(gs,GENE_FAMILIES$HGNC_SYMBOL[[i]])
        x<-length(cg)

        if(x>=2){
            currentP<-my.hypTest(x,k,n,N)
            if(flag==1){
                RESTOT<-matrix(c(i,currentP),1,2)
            }else{
                RESTOT<-rbind(RESTOT,c(i,currentP))
            }
            containedGenes[flag]<-paste(cg,collapse=', ')
            N_[flag]<-N
            n_[flag]<-n
            k_[flag]<-k
            x_[flag]<-x
            flag<-flag+1
        }
    }


    if(exists('RESTOT')){
        RESTOT<-cbind(RESTOT,p.adjust(RESTOT[,2],'fdr'))
        colnames(RESTOT)<-c('pathId','p','adj.p')

        RESTOT<-data.frame(cbind(RESTOT,x_,k_,n_,N_,containedGenes,PATHCOM_HUMAN$PATHWAY[RESTOT[,1]]),stringsAsFactors = FALSE)

        ENRpath<-RESTOT$pathId[which(RESTOT$adj.p<0.05)]

        return(list(famTAB=RESTOT,ENRfam=ENRpath))
    }else{
        return(NULL)
    }
}

ALLmarkerDiffExp<-function(ANOVAresultsDir){

    #source('Libraries/OT15.Misc.R')
    load(paste(ANOVAresultsDir,'OUTPUT/ANOVA_results.rdata',sep=''))

    id<-which((as.numeric(TOTRES[,"ANOVA FEATURE FDR %"])<30 |
                   as.numeric(TOTRES[,"FEATURE_ANOVA_pval"] < 10^-3) |
                   as.numeric(TOTRES[,"FEATURE_ANOVA_pval"] < 0.05) |
                   as.numeric(TOTRES[,"FEATURE_ESS_T_pval"] < 0.05)) &
                  (as.numeric(TOTRES[,"FEATUREpos_ESS_MEAN"]) > as.numeric(TOTRES[,"FEATUREneg_ESS_MEAN"])))

    markers<-
        unique(unlist(TOTRES[id,"FEATURE"]))

    markers<-unlist(str_split(markers,', '))

    nMarkers<-length(markers)

    flag<-0

    ugenes<-NULL
    inEnrichedPathGenes<-NULL
    incommonMarkENRpath<-NULL

    newdir<-paste(ANOVAresultsDir,'/upRegGenesPerMarker/',sep = '')
    if(!dir.exists(newdir)){
        dir.create(newdir)
    }

    for (i in 1:length(markers)){
        print(c(i,markers[i]))

        currDEPgenes<-unlist(TOTRES[which(TOTRES[,"FEATURE"]==markers[i] &
                  (as.numeric(TOTRES[,"ANOVA FEATURE FDR %"])<30 |
                       as.numeric(TOTRES[,"FEATURE_ANOVA_pval"] < 10^-3) |
                       as.numeric(TOTRES[,"FEATURE_ANOVA_pval"] < 0.05) |
                       as.numeric(TOTRES[,"FEATURE_ESS_T_pval"] < 0.05)) &
                  (as.numeric(TOTRES[,"FEATUREpos_ESS_MEAN"]) > as.numeric(TOTRES[,"FEATUREneg_ESS_MEAN"]))),"Depleted Gene"])

        if(length(currDEPgenes)>1){

                    PATHRES<-
                        pathEnrichAnalysis(gs = currDEPgenes)

                    ABpathTAB<-PATHRES$pathTAB

                    save(ABpathTAB,file=paste(newdir,str_replace(markers[i],':','_'),'_commonMarkENRpath.RData',sep=''))


                    FAMRES<-
                        familyEnrichAnalysis(gs = currDEPgenes)

                    ABfamTAB<-FAMRES$famTAB
                    save(ABfamTAB,file=paste(newdir,str_replace(markers[i],':','_'),'_commonMarkENRfam.RData',sep=''))



                    incommonMarkENRpath<-
                        union(incommonMarkENRpath,unlist(PATHCOM_HUMAN$HGNC_SYMBOL[as.numeric(PATHRES$ENRpath)]))
        }


        if(!str_detect(markers[i],'mut')){
            markers[i]<-str_replace(markers[i],pattern = '_',':')
            markers[i]<-str_replace(markers[i],pattern = '_',' ')
        }

        DEXg<-diffExpressedGenesBasedOnMarker(markers[i],display = FALSE)
        if(length(DEXg)>0){


            save(DEXg,file=paste(newdir,str_replace(markers[i],':','_'),'_diffExp.RData',sep=''))

            currentGS<-DEXg$UpRegulatedGenes
            currentGS<-intersect(currentGS,PATHCOM_HUMAN$backGround)

            if(!flag){
                ugenes<-DEXg$UpRegulatedGenes
                flag<-1
            }else{
                ugenes<-union(ugenes,DEXg$UpRegulatedGenes)
            }

            if(length(currentGS)>1){
                PATHRES<-
                    pathEnrichAnalysis(gs = currentGS)

                pathTAB<-PATHRES$pathTAB
                save(pathTAB,file=paste(newdir,str_replace(markers[i],':','_'),'_upRegPaths.RData',sep=''))

                inEnrichedPathGenes<-
                    union(inEnrichedPathGenes,unlist(PATHCOM_HUMAN$HGNC_SYMBOL[as.numeric(PATHRES$ENRpath)]))


            }
        }


    }

    UpRegulatedGenes<-matrix(0,nrow = length(markers),ncol = length(ugenes),dimnames = list(markers,ugenes))
    GenesInUpRegulatedPathways<-matrix(0,nrow = length(markers),ncol = length(inEnrichedPathGenes),dimnames = list(markers,inEnrichedPathGenes))
    GenesInENRpatwhayInCommonMark<-matrix(0,nrow = length(markers),ncol = length(incommonMarkENRpath),dimnames = list(markers,incommonMarkENRpath))

    for (i in 1:length(markers)){
        filename<-markers[i]
        filename<-str_replace(filename,pattern = ':',replacement = '_')

        filename<-paste(ANOVAresultsDir,'upRegGenesPerMarker/',filename,'_diffExp.RData',sep='')

        if(file.exists(filename)){
            load(filename)
            UpRegulatedGenes[markers[i],DEXg$UpRegulatedGenes]<-1
        }

        filename<-markers[i]
        filename<-str_replace(filename,pattern = ':',replacement = '_')
        filename<-paste(ANOVAresultsDir,'upRegGenesPerMarker/',filename,'_upRegPaths.RData',sep='')

        if(file.exists(filename)){
            load(filename)
            GenesInUpRegulatedPathways[markers[i],unique(unlist(PATHCOM_HUMAN$HGNC_SYMBOL[as.numeric(pathTAB$pathId[which(pathTAB$adj.p<0.05)])]))]<-1
        }

        filename<-markers[i]
        filename<-str_replace(filename,pattern = ':',replacement = '_')
        filename<-paste(ANOVAresultsDir,'upRegGenesPerMarker/',filename,'_commonMarkENRpath.RData',sep='')

        if(file.exists(filename)){
            load(filename)
            GenesInENRpatwhayInCommonMark[markers[i],unique(unlist(PATHCOM_HUMAN$HGNC_SYMBOL[as.numeric(ABpathTAB$pathId[which(ABpathTAB$adj.p<0.05)])]))]<-1
        }

    }


    save(UpRegulatedGenes,file=paste(ANOVAresultsDir,'UpRegulatedGenes.Rdata',sep=''))
    save(GenesInUpRegulatedPathways,file=paste(ANOVAresultsDir,'GenesInUpRegulatedPathways.Rdata',sep=''))
    save(GenesInENRpatwhayInCommonMark,file=paste(ANOVAresultsDir,'GenesInEnrichePathwaysPerCommonMarkers.Rdata',sep=''))
}

