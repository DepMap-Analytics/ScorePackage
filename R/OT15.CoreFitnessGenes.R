
bootstrappedCoreFitnessGenes<-function(CFlist,outDir){

    Sensitivity<-list()
    Precision<-list()
    Fmeasure<-list()
    PRED<-list()
    OBS<-list()

    for (i in 1:(length(CFlist)-1)){
        performances <- bootstrappedIteration(CFlist, i)

        PRED[[i]]<-performances$Predicted
        OBS[[i]]<-performances$Observed
        Sensitivity[[i]] <- performances$Sensitivity
        Precision[[i]] <- performances$Precision
        Fmeasure[[i]] <- performances$Fmeasure
    }


    pdf(paste(outDir,'02_bootstrapped_Pred_of_CFgenes_Perf.pdf',sep=''))
    par(mfrow=c(3,1))
    par(mar=c(4,6,4,1))
    boxplot(Sensitivity,names = NA,main='Sensitivity')
    boxplot(Precision,names = NA,main='Precision')
    boxplot(Fmeasure,names = c(1:length(Fmeasure)),main='Fmeasure',xlab='n. intersected cancer type cfg sets')

    perfMea<-unlist(lapply(Fmeasure,'mean'))

    id<-max(which(perfMea==max(perfMea)))

    abline(v=id,col='red',lwd=3)

    print(paste('best performing intersection cardinality =',id))
    dev.off()


    glgene<-unique(unlist(PRED[[id]]))

    appearingMat<-matrix(0,length(glgene),length(PRED[[id]]),dimnames = list(glgene,NULL))

    for (i in 1:length(PRED[[id]])){
        appearingMat[PRED[[id]][[i]],i]<-1
    }

    th<-median(rowSums(appearingMat))

    pdf(paste(outDir,'03_median_n_intersection_in_which_a_gene_occurs.pdf',sep=''))
    par(mfrow=c(1,1))
    hist(rowSums(appearingMat),ylab=paste('n. genes in at least one intersection of',id,'cancer type cfg sets'),
         xlab=paste('n. of intersections (out of ',length(PRED[[id]]),sep=''),main=paste('median = ',th))
    abline(v=th,col='red',lwd=3)
    dev.off()

    corefitness<-names(which(rowSums(appearingMat)>=th))
    return(corefitness)
}

bootstrappedIteration<-function(CFlist,k){

    nel<-length(CFlist)
    NC<-combn(names(CFlist),k)

    NN<-ncol(NC)

    PredictedPC<-list()
    ObservedPC<-list()

    TPR<-vector()
    PPV<-vector()

    for (i in 1:NN){

        Training<-NC[,i]
        Test<-setdiff(names(CFlist),Training)
        PredictedPC[[i]]<-Reduce(intersect,CFlist[Training])
        ObservedPC[[i]]<-Reduce(intersect,CFlist[Test])

        TPR[i]<-length(intersect(PredictedPC[[i]],ObservedPC[[i]]))/length(ObservedPC[[i]])
        PPV[i]<-length(intersect(PredictedPC[[i]],ObservedPC[[i]]))/length(PredictedPC[[i]])

    }

    return(list(Sensitivity=TPR,Precision=PPV,Fmeasure=(2*(PPV*TPR)/(PPV+TPR)),Predicted=PredictedPC,Observed=ObservedPC))
}

pathwayEnrichmentAnalysis<-function(geneSetTotest,PATHCOM_HUMAN,FDRth=Inf,bgGS){


   # bgGS<-intersect(bgGS,PATHCOM_HUMAN$backGround)

    npath<-length(PATHCOM_HUMAN$PATHWAY)

    k<-length(intersect(geneSetTotest,bgGS))
    N<-length(bgGS)
    p<-vector()
    x<-vector()
    n<-vector()

    includedGenes<-list()

    for (i in 1:npath){
        x[i]<-length(intersect(geneSetTotest,PATHCOM_HUMAN$HGNC_SYMBOL[[i]]))
        n[i]<-length(intersect(PATHCOM_HUMAN$HGNC_SYMBOL[[i]],bgGS))
        p[i]<-my.hypTest(x[i],k,n[i],N)
        includedGenes[[i]]<-paste(sort(intersect(geneSetTotest,PATHCOM_HUMAN$HGNC_SYMBOL[[i]])),
                                  collapse=', ')

        currentLine<-c(x[i],k,n[i],N,includedGenes[[i]])


        if (i == 1){
            RES<-currentLine
        }else{
            RES<-rbind(RES,currentLine)
        }
    }

    RES<-cbind(1:nrow(RES),RES,p,p.adjust(p,'fdr'))
    rownames(RES)<-
        PATHCOM_HUMAN$PATHWAY
    colnames(RES)<-c('id','x','k','n','N','Genes','p','adj.p')
    #RES<-RES[order(RES[,"adj.p"]),]

    id<-which(as.numeric(RES[,"adj.p"])<FDRth)

    RES<-RES[id,]

    nnames<-rownames(RES)
    nnames<-str_sub(nnames,1,60)

    rownames(RES)<-nnames
    par(mar=c(4,35,4,4))
    RES[as.numeric(RES[,'adj.p'])==0,'adj.p']<-min(as.numeric(RES[,'adj.p'])==0)

    RES<-RES[as.numeric(RES[,'adj.p'])<0.05,]

    tt<-barplot(rev(sort(-log10(as.numeric(RES[,'adj.p'])),decreasing=TRUE)),horiz = TRUE,las=2,xlab='-log10 adj.p',
            names.arg = rownames(RES)[order(-log10(as.numeric(RES[,'adj.p'])))],border = FALSE)

    text(rep(1,length(tt))-1,tt-0.06,rev(RES[,'Genes']),pos = 4,cex=0.7)

    return(RES)
}


OT15_assembleDepletionMatrix<-function(cellLines,BFths,inputFolder=NULL,BFvalues=FALSE,BFmatrix=NULL){

  if(!is.null(inputFolder)){
    print(paste(inputFolder,cellLines[1],'_GeneLevelBF.rdata',sep=''))
    if(file.exists(paste(inputFolder,cellLines[1],'_GeneLevelBF.rdata',sep=''))){
      load(paste(inputFolder,cellLines[1],'_GeneLevelBF.rdata',sep=''))
      genes<-rownames(GeneLevelBFs)

      genes<-setdiff(genes,c("",NA))

    resMat<-foreach(i=1:length(cellLines),.combine=cbind)%dopar%{


          if(file.exists(paste(inputFolder,cellLines[i],'_GeneLevelBF.rdata',sep=''))){
            load(paste(inputFolder,cellLines[i],'_GeneLevelBF.rdata',sep=''))
            if(BFvalues){
              tmp<-matrix(NA,length(genes),1,
                          dimnames = list(genes,cellLines[i]))

              inboth<-intersect(rownames(tmp),rownames(GeneLevelBFs))
              inboth<-setdiff(inboth,c(""))
              tmp[inboth,1]<-GeneLevelBFs[inboth,1]
            }else{
              tmp<-matrix(NA,length(genes),1,
                          dimnames = list(genes,cellLines[i]))
              inboth<-intersect(rownames(tmp),rownames(GeneLevelBFs))
              inboth<-setdiff(inboth,c(""))
              tmp[inboth,1]<-GeneLevelBFs[inboth,1]-BFths[cellLines[i]]


            }
          }else{
            tmp<-matrix(NA,length(genes),1,
                        dimnames = list(genes,cellLines[i]))
            tmp
          }


    }


    colnames(resMat)<-cellLines
    }
  }else{
    if(!is.null(BFmatrix)){
      resMat<-matrix(NA,nrow(BFmatrix),ncol(BFmatrix))
      dimnames(resMat)<-dimnames(BFmatrix)
      for(i in 1:ncol(resMat)){
        resMat[,i]<-BFmatrix[,i]-BFths[colnames(BFmatrix)[i]]

      }
    }
  }
    resMat<-resMat+0
    return(resMat)
}
OT15_assembleDepletionMatrixL<-function(cellLines,BFths,inputFolder=NULL,BFvalues=FALSE,BFmatrix=NULL){

  if(!is.null(inputFolder)){
    print(paste(inputFolder,cellLines[1],'_GeneLevelBFL.rdata',sep=''))
    if(file.exists(paste(inputFolder,cellLines[1],'_GeneLevelBFL.rdata',sep=''))){
      load(paste(inputFolder,cellLines[1],'_GeneLevelBFL.rdata',sep=''))
      genes<-rownames(GeneLevelBFs)
    }


    resMat<-foreach(i=1:length(cellLines),.combine=cbind)%dopar%{


      if(file.exists(paste(inputFolder,cellLines[i],'_GeneLevelBFL.rdata',sep=''))){
        load(paste(inputFolder,cellLines[i],'_GeneLevelBFL.rdata',sep=''))
        if(BFvalues){
          tmp<-matrix(NA,length(genes),1,
                      dimnames = list(genes,cellLines[i]))
          tmp[rownames(GeneLevelBFs),1]<-GeneLevelBFs[,1]
        }else{
          tmp<-matrix(NA,length(genes),1,
                      dimnames = list(genes,cellLines[i]))
          tmp[rownames(GeneLevelBFs),1]<-GeneLevelBFs[,1]-BFths[cellLines[i]]

        }
      }else{
        tmp<-matrix(NA,length(genes),1,
                    dimnames = list(genes,cellLines[i]))
        tmp
      }


    }


    colnames(resMat)<-cellLines
  }else{
    if(!is.null(BFmatrix)){
      resMat<-matrix(NA,nrow(BFmatrix),ncol(BFmatrix))
      dimnames(resMat)<-dimnames(BFmatrix)
      for(i in 1:ncol(resMat)){
        resMat[,i]<-BFmatrix[,i]-BFths[colnames(BFmatrix)[i]]

      }
    }
  }
  resMat<-resMat+0
  return(resMat)
}

OT15_assembleMageckFDRs<-function(cellLines,inputFolder,depletion=TRUE){
    mgkSummary<-read.table(paste(inputFolder,
                               cellLines[1],'.gene_summary.txt',sep=''),header=TRUE,stringsAsFactors = FALSE)

    genes<-mgkSummary$id
    genes<-setdiff(genes,c("",NA))
    mgkFDRs<-NULL
    mgkFDRs<-foreach(i=1:length(cellLines),.combine=cbind)%dopar%{
        print(i)


        mgkSummary<-read.table(paste(inputFolder,
                                     cellLines[i],'.gene_summary.txt',sep=''),header=TRUE,stringsAsFactors = FALSE)

        mgkFDRs<-rep(NA,length(genes))
        names(mgkFDRs)<-genes
        #if(i == 1){
            #genes<-sort(mgkSummary$id)
            if(depletion){
                mgkFDRs<-mgkSummary$neg.fdr[match(genes,mgkSummary$id)]
            }else{
                mgkFDRs<-mgkSummary$pos.fdr[match(genes,mgkSummary$id)]
            }

        #}else{
         #   if(depletion){
         #       mgkFDRs<-cbind(mgkFDRs,mgkSummary$neg.fdr[match(genes,mgkSummary$id)])
         #   }else{
         #       mgkFDRs<-cbind(mgkFDRs,mgkSummary$pos.fdr[match(genes,mgkSummary$id)])
        #    }
        #}

    }

    rownames(mgkFDRs)<-genes
    colnames(mgkFDRs)<-cellLines

    return(return(mgkFDRs))
}
OT15_panessprofile<-function(depMat,display=TRUE,NLIMS=1,main=''){
    depMat<-depMat[which(rowSums(depMat)>0),]
    panessprof<-rep(0,ncol(depMat))
    names(panessprof)<-as.character(1:ncol(depMat))
    paness<-summary(as.factor(rowSums(depMat)),maxsum = length(unique(as.factor(rowSums(depMat)))))
    panessprof[as.character(names(paness))]<-paness

    CUMsums<-rev(cumsum(rev(panessprof)))

    names(CUMsums)<-paste('>=',names(CUMsums),sep='')

    if(display){
        par(mfrow=c(2,1))
        par(mar=c(6,4,4,1))

        if(length(main)==0){
            main=c(nrow(depMat),'genes depleted in at least 1 cell line')
        }
        barplot(panessprof,ylab='n.genes',xlab=paste('n. cell lines with at least',NLIMS,'deleted genes'),cex.axis = 0.8,cex.names = 0.8,
                las=2,main=main)

        barplot(CUMsums,ylab='n.genes',xlab=paste('n. cell lines with at least',NLIMS,'deleted genes'),cex.axis = 0.8,cex.names = 0.6,
                las=2,main='Cumulative sums')

    }
    return(list(panessprof=panessprof,CUMsums=CUMsums))
}
OT15_randomisedepMat<-function(depMat){
    rmat<-apply(depMat,2,sample)
}
OT15_generateNullModel<-function(depMat,ntrials=100,display=TRUE){

    depMat<-depMat[which(rowSums(depMat)>0),]
    nullProf<-matrix(NA,ntrials,ncol(depMat),dimnames = list(1:ntrials,1:ncol(depMat)))
    nullCumSUM<-matrix(NA,ntrials,ncol(depMat),dimnames = list(1:ntrials,paste('≥',1:ncol(depMat),sep='')))
    print('Generating null model...')
    pb <- txtProgressBar(min=1,max=ntrials,style=3)

    for (i in 1:ntrials){
        setTxtProgressBar(pb, i)
        rMat<-
            OT15_randomisedepMat(depMat)
        Ret<-
            OT15_panessprofile(rMat,display = FALSE)
        nullProf[i,]<-Ret$panessprof
        nullCumSUM[i,]<-Ret$CUMsums
    }
    Sys.sleep(1)
    close(pb)
    print('')
    print('Done')

    if (display){
        par(mfrow=c(2,1))
        main=c(paste(ntrials,' randomised essentiality profiles of\n',nrow(depMat),' genes across ',ncol(depMat),' cell lines',
                     sep=''))
        boxplot(nullProf,las=2,xlab='n cell lines',ylab='genes depleted in n cell lines',main=main)
        colnames(nullCumSUM)<-paste(">=",1:ncol(nullCumSUM))
        boxplot(log10(nullCumSUM+1),las=2,main='Cumulative sums',xlab='n cell lines',
                ylab='log10 [number of genes + 1]',
                cex.axis=0.8)
    }

    return(list(nullProf=nullProf,nullCumSUM=nullCumSUM))
}
OT15_empiricalOdds<-function(observedCumSum,simulatedCumSum){

    nsamples<-length(observedCumSum)
    ntrials<-nrow(simulatedCumSum)

    odds<-rep(NA,1,nsamples)
    names(odds)<-paste('≥',1:nsamples,sep='')
    for (i in 1:nsamples){

        PDF<-density(simulatedCumSum[,i])


        odds[i]<- log10(observedCumSum[i]/mean(simulatedCumSum[,i]))

    }
    return(odds)
}

OT15_truePositiveRate<-function(depMat,essentialGeneSet){
    nsamples<-ncol(depMat)

    essentialGeneSet<-intersect(essentialGeneSet,rownames(depMat))

    TPR<-rep(NA,1,nsamples)
    names(TPR)<-paste('≥',1:nsamples,sep='')

    ncells<-rowSums(depMat)

    TP<-rep(NA,1,nsamples)
    names(TP)<-paste('≥',1:nsamples,sep='')

    P<-rep(NA,1,nsamples)
    names(P)<-paste('≥',1:nsamples,sep='')

    for (i in nsamples:1){
        positiveSet<-names(which(ncells>=i))
        P[i]<-length(positiveSet)
        truepositives<-intersect(positiveSet,essentialGeneSet)
        TP[i]<-length(truepositives)
        TPR[i]<-TP[i]/length(essentialGeneSet)
    }

    return(list(P=P,TP=TP,TPR=TPR))
}

OT15_tradeoffEO_TPR<-function(EO,TPR,test_set_name,point=NULL){

    if(length(point)==0){
        CCOL<-'red'
    }else{
        CCOL<-rgb(255,0,0,alpha = 100,maxColorValue = 255)
    }

    x<-EO
    x[x==Inf]<-max(x[x<Inf])
    x<-(x-min(x))/(max(x)-min(x))

    y<-TPR
    y<-(y-min(y))/(max(y)-min(y))

    orEO<-EO
    orEO[orEO==Inf]<-max(orEO[orEO<Inf])
    orTPR<-TPR

    EO<-x
    TPR<-y
    par(mar=c(4,4,4,4))
    MAIN<-c('log10 (obs/Expct) n.genes [red, left]',
            paste('% covered ',test_set_name,' [blue, right]',sep=''))
    plot(EO,type='l',xlab='genes depleted in >= # cell lines',ylab='',axes=FALSE,lwd=4,main=MAIN,col=CCOL,cex.main=0.8,
         xlim=c(0,length(EO)))
    axis(2,at = seq(0,1,0.2),format(seq(min(orEO),max(orEO),(max(orEO)-min(orEO))/5),digits=2))
    axis(1)
    par(new=TRUE)
    plot(TPR,type='l',xlab='',ylab='',axes=FALSE,lwd=4,col='blue',ylim=c(0,1),xlim=c(0,length(EO)))
    labels<-format(seq(min(orTPR),max(orTPR),(max(orTPR)-min(orTPR))/5),digits=2)
    if(length(labels!=6)){labels=rep(min(orTPR),6)}
    axis(4,at = seq(0,1,0.2),labels)

    if(length(point)==0){
        point<-min(which(!y>x))

        abline(v=point)
        abline(h=y[point],lty=2)

        points(point,y[point],pch=16,cex=2)
    }else{
        abline(v=point)
        abline(h=y[point],lty=2)
        points(point,y[point],pch=16,cex=2)
    }
    legend('top',paste(format(100*orTPR[point],digits=2),'% covered',sep=''),bg = NULL,bty = 'n')

    return(point)



}

OT15_summaries<-function(depMat,essentialGeneSet,tradeoffnc=21){

    ncellLines<-ncol(depMat)

    OU<-matrix(NA,4,ncellLines,
               dimnames = list(c('pred core-fitness in BAGEL ref','Other CoreFitness',
                                 'others in BAGEL ref', 'context specific'),
                               colnames(depMat)))

    coreFitness<-names(which(rowSums(depMat)>=tradeoffnc))

    coreFitnessKnown<-intersect(coreFitness,essentialGeneSet)
    coreFitnessOther<-setdiff(coreFitness,essentialGeneSet)
    KnownOther<-setdiff(essentialGeneSet,coreFitness)


    for (i in 1:ncellLines){
        depleted<-names(which(depMat[,i]>0))

        OU[1,i]<-length(intersect(depleted,coreFitnessKnown))
        OU[2,i]<-length(intersect(depleted,coreFitnessOther))
        OU[3,i]<-length(intersect(depleted,KnownOther))
        OU[4,i]<-length(setdiff(depleted,union(coreFitness,KnownOther)))

        if (i==1){
            totalContextSpec<-setdiff(depleted,union(coreFitness,KnownOther))
        }else{

            totalContextSpec<-union(totalContextSpec,
                                    setdiff(depleted,union(coreFitness,KnownOther)))
        }
    }

    par(mar=c(5,8,1,2))
    barplot(OU[,order(colSums(OU))],horiz=TRUE,las=2,xlim=c(0,2000),border=NULL,
            xlab='n. depleted genes - FDR (1 - precision) < 5%')
    #save(totalContextSpec,file='../../RESULTS/BAGEL-R_output/_predictedContextSpec.rdata')
}

OT15_adaptiveDaisyModel<-function(){

    source('Libraries/OT15.priorKnown_essentiality.R')
    load('../../DATA/R/annotations/TissueColors.Rdata')

    load('../../RESULTS/CoreFitnessGenes/breast_coreFitnessGenes.Rdata')
    load('../../RESULTS/CoreFitnessGenes/breast_min.N.cells.Rdata')

    brca_cf<-
        coreFitnessGenes
    brca_min<-
        minimalN.cells


    load('../../RESULTS/CoreFitnessGenes/large_intestine_coreFitnessGenes.Rdata')
    load('../../RESULTS/CoreFitnessGenes/large_intestine_min.N.cells.Rdata')
    coread_cf<-
        coreFitnessGenes
    coread_min<-
        minimalN.cells

    load('../../RESULTS/CoreFitnessGenes/melanoma_coreFitnessGenes.Rdata')
    load('../../RESULTS/CoreFitnessGenes/melanoma_min.N.cells.Rdata')
    skcm_cf<-
        coreFitnessGenes
    skcm_min<-
        minimalN.cells

    load('../../RESULTS/BAGEL-R_output/ES8_GeneLevelBF.rdata')
    load('../../RESULTS/BAGEL-R_MAGECK_comparison_summaries/BAGEL_BEST_PPV_th.Rdata')

    ews_cf<-
        names(which(GeneLevelBFs[,1]>BAGEL_BEST_PPV_th['ES8']))
    ews_min<-1

    load('../../RESULTS/BAGEL-R_output/NCI-H3122_GeneLevelBF.rdata')
    load('../../RESULTS/BAGEL-R_MAGECK_comparison_summaries/BAGEL_BEST_PPV_th.Rdata')

    luad_cf<-
        names(which(GeneLevelBFs[,1]>BAGEL_BEST_PPV_th['NCI-H3122']))
    luad_min<-1

    ESSgenes<-OT15_assembleEssGenes()$essential

    coveredEssGenes<-c(length(intersect(coread_cf,ESSgenes)),
                       length(intersect(brca_cf,ESSgenes)),
                       length(intersect(skcm_cf,ESSgenes)),
                       length(intersect(ews_cf,ESSgenes)),
                       length(intersect(luad_cf,ESSgenes)))

    coread_new_cf<-setdiff(coread_cf,ESSgenes)
    brca_new_cf<-setdiff(brca_cf,ESSgenes)
    skcm_new_cf<-setdiff(skcm_cf,ESSgenes)
    ews_new_cf<-setdiff(ews_cf,ESSgenes)
    luad_new_cf<-setdiff(luad_cf,ESSgenes)

    otherEssGenes<-c(length(coread_new_cf),
                     length(brca_new_cf),
                     length(skcm_new_cf),
                     length(ews_new_cf),
                     length(luad_new_cf))


    toPlot<-c(length(coread_cf),length(brca_cf),length(skcm_cf),length(ews_cf),length(luad_cf))

    tosPlot<-rbind(coveredEssGenes,otherEssGenes)

    par(mfrow=c(2,1))

    tmp<-barplot(toPlot,
            col=c(TissueColors[c('large_intestine','breast','melanoma',"ewings_sarcoma","lung_NSCLC_adenocarcinoma")]),
            main='Core fitness genes per cancer type',ylab='n. genes',
            names.arg = c('COREAD','BRCA','SKCM','EWS','LUAD'),ylim=c(0,1700))

    COLS<-c(TissueColors[c('large_intestine','breast','melanoma',"ewings_sarcoma","lung_NSCLC_adenocarcinoma")])

    tmp<-barplot(tosPlot,
                 main='Core fitness genes per cancer type',ylab='n. genes',
                 names.arg = c('COREAD','BRCA','SKCM','EWS','LUAD'),ylim=c(0,1700))

    text(tmp,toPlot,c(coread_min,brca_min,skcm_min,ews_min,luad_min),pos = 3)

    SETS<-list(COREAD=coread_new_cf,BRCA=brca_new_cf,SKCM=skcm_new_cf,EWS=ews_new_cf,LUAD=luad_new_cf)

    allG<-unique(unlist(SETS))

    membMat<-matrix(0,length(allG),length(SETS),dimnames = list(allG,names(SETS)))

    for (i in 1:length(SETS)){
        membMat[SETS[[i]],names(SETS)[i]]<-1
    }

    panCancer_CFG<-names(which(rowSums(membMat)==ncol(membMat)))
    panCancer_CFG<-FGtable(panEssentialGenes = panCancer_CFG,sortingCriteria = 'all')
    save(panCancer_CFG,file='../../RESULTS/CoreFitnessGenes/panCancer_CFG.Rdata')

    coread_CFG<-setdiff(names(which(membMat[,'COREAD']==1)),rownames(panCancer_CFG))
    coread_CFG<-FGtable(panEssentialGenes = coread_CFG,sortingCriteria = 'large_intestine')
    save(coread_CFG,file='../../RESULTS/CoreFitnessGenes/coread_CFG.Rdata')

    brca_CFG<-setdiff(names(which(membMat[,'BRCA']==1)),rownames(panCancer_CFG))
    brca_CFG<-FGtable(panEssentialGenes = brca_CFG,sortingCriteria = 'breast')
    save(brca_CFG,file='../../RESULTS/CoreFitnessGenes/brca_CFG.Rdata')

    skcm_CFG<-setdiff(names(which(membMat[,'SKCM']==1)),rownames(panCancer_CFG))
    skcm_CFG<-FGtable(panEssentialGenes = skcm_CFG,sortingCriteria = 'melanoma')
    save(skcm_CFG,file='../../RESULTS/CoreFitnessGenes/skcm_CFG.Rdata')

}


FGtable<-function(panEssentialGenes,sortingCriteria='all'){

    cellLines<-
        unique(PS_inventory$Cell_Line_Name)

    depMat<-OT15_assembleDepletionMatrix(cellLines = cellLines,inputFolder = '../../RESULTS/BAGEL-R_output/',BFths = BAGEL_BEST_PPV_th)

    depMat<-depMat[panEssentialGenes,]

    types<-unique(PS_inventory$GDSC_Description_2)

    for (i in 1:length(types)){

        if(!is.element(types[i],c("ewings_sarcoma","lung_NSCLC_adenocarcinoma"))){
            currentId<-rowSums(depMat[,unique(PS_inventory$Cell_Line_Name[PS_inventory$GDSC_Description_2==types[i]])])
        }else{
            currentId<-depMat[,unique(PS_inventory$Cell_Line_Name[PS_inventory$GDSC_Description_2==types[i]])]
        }

        if (i == 1){
            summa<-currentId
        }else{
            summa<-cbind(summa,currentId)
        }
    }
    colnames(summa)<-types

    if(sortingCriteria=='all'){
        totSum<-rowSums(summa)
        summa<-summa[order(totSum,decreasing=TRUE),]
    }else{
        summa<-summa[order(summa[,sortingCriteria],decreasing=TRUE),]
    }

    if(sortingCriteria=='all'){

        par(xpd=TRUE)
        summa<-summa[,c('large_intestine','breast','melanoma','ewings_sarcoma','lung_NSCLC_adenocarcinoma')]
        col<-TissueColors[c('large_intestine','breast','melanoma','ewings_sarcoma','lung_NSCLC_adenocarcinoma')]
        corplot<-barplot(t(summa[min(100,nrow(summa)):1,]),horiz = TRUE,las=2,col = col,
                border=FALSE,cex.names = 0.5)
        tmp<-OT15_retrievegeneInfo(rownames(summa))
        tmp1<-
            OT15_retrieveKinases_iGof_GSKtar(rownames(summa))



    }else{
        summa<-summa[,sortingCriteria]
        col<-TissueColors[sortingCriteria]

        par(mfrow=c(1,2))

        corplot<-barplot(rev(summa[1:min(100,length(summa))]),horiz = TRUE,las=2,col = col,
                border=FALSE,cex.names = 0.5,xlab='n. vulnerable cell lines',
                main='TOP100 cancer type specific essential genes')

        barplot(rev(summa[intersect(names(summa),union(cancerDrivers$ActingDriver_Symbol,InigoList))]),horiz = TRUE,las=2,col = col,
                border=FALSE,cex.names = 0.5,xlab='n. vulnerable cell lines',main='HC Cancer Genes')

        barplot(rev(summa[intersect(names(summa),GSK_epiTargets)]),horiz = TRUE,las=2,col = col,
                border=FALSE,cex.names = 0.5,xlab='n. vulnerable cell lines',main='GSK epigenetic targets')

        barplot(rev(summa[intersect(names(summa),union(kinase_genes_and_syn,humKinome))]),horiz = TRUE,las=2,col = col,
                border=FALSE,cex.names = 0.5,xlab='n. vulnerable cell lines',main='Kinases (Uniprot & HumKinome)')


    }

    ranges<-max(c(summa))

    summa<-data.frame(summa)



    return(summa)
}
