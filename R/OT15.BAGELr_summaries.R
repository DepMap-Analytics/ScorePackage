
retrieveAMPgenesPicnic<-function(cellLine){
    cid<-as.character(MASTER_LIST$COSMIC.ID[MASTER_LIST$Analysis.Set.Name==cellLine])
    genes<-rownames(geneLevCNA$CNV)

    if(length(intersect(cid,colnames(geneLevCNA$CNV)))>0){

        amplifiedGenes<-sort(genes[unique(c(grep(',8,',geneLevCNA$CNV[,cid]),
                                            grep(',9,',geneLevCNA$CNV[,cid]),
                                            grep(',10,',geneLevCNA$CNV[,cid]),
                                            grep(',11,',geneLevCNA$CNV[,cid]),
                                            grep(',12,',geneLevCNA$CNV[,cid]),
                                            grep(',13,',geneLevCNA$CNV[,cid]),
                                            grep(',14,',geneLevCNA$CNV[,cid]),
                                            grep(',15,',geneLevCNA$CNV[,cid]),
                                            grep(',16,',geneLevCNA$CNV[,cid])))])

        return(amplifiedGenes)
    }else{
        return(NULL)
    }
}

retrieveGisticAmp<-function(cellLine){

    cid<-as.character(MASTER_LIST$COSMIC.ID[MASTER_LIST$Analysis.Set.Name==cellLine])

    if(length(cid)>0 & length(intersect(cid,colnames(gisticCNA)))>0){
        gisticProf<-gisticCNA[,cid]


        tmpGe<-rownames(gisticCNA)[which(gisticProf== 2)]

        return(tmpGe)
    }else{
        return(NULL)
    }

}


findZeroBFcoordinate<-function(x,y,val=0){

    return(x[which(y==min(y[which(y>val)]))])

}

scatterhist <- function(x, y,essentialGs,nonEssentialGs,essentialInd,cytotox,CL,BEST_ppv_th_training,BEST_ppv_th_test){
    zones=matrix(c(2,0,1,3), ncol=2, byrow=TRUE)
    layout(zones, widths=c(4/5,1/5), heights=c(2/6,4/6))
    xhist = hist(x, plot=FALSE)
    yhist = hist(abs(y), plot=FALSE)
    top = max(c(xhist$counts, yhist$counts))
    par(mar=c(4,4,1,4))


    cytotox<-intersect(rownames(x),cytotox)
    essentialInd<-intersect(rownames(x),essentialInd)
    essentialGs<-intersect(rownames(x),essentialGs)
    nonEssentialGs<-intersect(rownames(x),nonEssentialGs)
    otherGenes<-setdiff(rownames(x),c(essentialGs,nonEssentialGs,essentialInd,cytotox))

    xother<-x[otherGenes,]
    xessentialInd<-x[essentialInd,]
    xessential<-x[essentialGs,]
    xnonessential<-x[nonEssentialGs,]
    xcytotox<-x[cytotox,]

    yother<-abs(y[otherGenes])
    yessentialInd<-abs(y[essentialInd])
    yessential<-abs(y[essentialGs])
    ynonessential<-abs(y[nonEssentialGs])
    ycytotox<-abs(y[cytotox])

    plot(xother,
         yother,
         col=opacDarkGray,
         pch=16,
         xlab='Corrected log FC',
         ylab='| Bayesian Factor |',
         xlim=c(min(x),max(x)),ylim=c(min(abs(y)),max(abs(y))),cex=1.2)

    points(xnonessential,ynonessential,col=opacDarkBlue,pch=16,cex=1.2)
    points(xessential,yessential,col=opacGreen,pch=16,cex=1.2)
    points(xcytotox,ycytotox,col=opacDarkGreen,pch=16,cex=1.2)
    points(xessentialInd,yessentialInd,bg=opacCornFlowerBlue,col='white',pch=21,cex=1.2)

    o<-order(yother,decreasing=TRUE)[1:50]

    points(xother[o],yother[o],bg=opacRed,col='white',pch=21,cex=1.2)

    zerox<-findZeroBFcoordinate(x,y)
    abline(v=zerox,col='darkgray',lty=3)

    abline(h=0,col='darkgray',lty=3)
    abline(h=1,col='darkgray',lty=1)

    abline(h=BEST_ppv_th_training,col='lightblue',lty=2)
    abline(h=BEST_ppv_th_test,col='blue',lty=2)

    abline(h=2,col='darkgray',lty=1)


    legend('topright',legend=c('top 50 | BF |','MsigDB ess','Ind.ess','Bagel ess','Bagel non-ess'),
           border=c('white','darkgreen','white','green','darkblue'),
           fill = c(opacRed,
                    opacDarkGreen,
                    opacCornFlowerBlue,
                    opacGreen,
                    opacDarkBlue),
           xpd=TRUE)


    par(mar=c(0,4,1,4))

    multDensPlot(TOPLOT = list(density(xother,na.rm = TRUE),
                               density(xessential,na.rm = TRUE),
                               density(xnonessential,na.rm = TRUE),
                               density(xcytotox,na.rm=TRUE),
                               density(xessentialInd,na.rm = TRUE)),
                 COLS = c('darkblue','green','darkgrey','darkgreen','cornflowerblue'),XLAB='',XLIMS=c(min(x),max(x)),XAXT = 'n')

    par(mar=c(4,3,0,1))
    barplot(sort(xother)[seq(50,1,-1)],horiz = TRUE,las=2,cex.names = 0.7,col=opacRed,border = 'white',xlab='|Bayesian Factor|')

    par(oma=c(0,0,3,3))
    mtext(text=CL, side=3, line=1, outer=TRUE, adj=0,cex =2,
          at=.5 * (mean(x) - min(x))/(max(x)-min(x)))


}

multDensPlot<-function(TOPLOT,COLS,XLIMS=NULL,YLIMS=NULL,TITLE='',LEGentries=NULL,XLAB='',XAXT='s',YAXT='s',VERTICAL=FALSE){

    YM<-vector()
    for (i in 1:length(TOPLOT)){
        YM[i]<-max(TOPLOT[[i]]$y,na.rm = TRUE)
    }

    Ymax<-max(YM,na.rm=TRUE)
    Xmax<-max(YM,na.rm=TRUE)


    if(length(XLIMS)){
        plot(0,0,col=NA,ylab='density',xlab=XLAB,frame.plot = FALSE,xaxt=XAXT,yaxt=YAXT,
             xlim=XLIMS,ylim=c(0,Ymax),type='l',main=TITLE)
    }else{
        plot(0,0,col=NA,ylab='density',xlab=XLAB,frame.plot = FALSE,xaxt=XAXT,yaxt=YAXT,
             ylim=YLIMS,xlim=c(0,Xmax),type='l',main=TITLE)
    }

    for (i in 1:length(TOPLOT)){

        if(!VERTICAL){
            cord.x <- c(TOPLOT[[i]]$x)
            cord.y <- c(TOPLOT[[i]]$y)
        }else{
            cord.x <- c(TOPLOT[[i]]$y)
            cord.y <- c(TOPLOT[[i]]$x)
        }

        rgbc<-col2rgb(COLS[i])
        currCol<-rgb(rgbc[1],rgbc[2],rgbc[3],alpha = 100,maxColorValue = 255)
        polygon(cord.x,cord.y,col=currCol,border = NA)

        if(!VERTICAL){
            lines(TOPLOT[[i]],col=COLS[i],lwd=3)
        }

        if (i == 1 & length(LEGentries)>0){
            legend('topleft',legend = LEGentries,col=COLS,lwd=3,bty = 'n')
        }
    }
}

bangmag.geneLevelFCs<-function(dirname,outDir){
    if(!exists(outDir)){
        dir.create(outDir)
    }

    fc<-dir(dirname)
    fc<-grep('_foldChanges.Rdata',fc,value=TRUE)

    nf<-length(fc)

    for (i in 1:nf){

        load(paste(dirname,fc[i],sep=''))

        genes<-unique(foldchanges$gene)
        ng<-length(genes)

        GL_FoldChanges<-matrix(NA,nrow = ng,ncol = 1,dimnames = list(genes,'GL_fc'))

        CL<-str_split(fc[i],'_')[[1]][1]
        ssam<-ncol(foldchanges)
        for (j in 1:ng){
            print(c(i,CL,j,ng))
            id<-which(foldchanges$gene==genes[j])
            GL_FoldChanges[genes[j],1]<-mean(unlist(foldchanges[id,3:ssam]))
        }

        save(GL_FoldChanges,file=paste(outDir,CL,'_GLFC.Rdata',sep=''))
    }
}

bangmag.geneLevel_correctedFCs<-function(dirname,outDir){
    if(!exists(outDir)){
        dir.create(outDir)
    }

    fc<-dir(dirname)
    fc<-grep('_00_correctedFCs.RData',fc,value=TRUE)

    nf<-length(fc)

    for (i in 1:nf){

        load(paste(dirname,fc[i],sep=''))

        genes<-unique(correctedFCs$genes)
        ng<-length(genes)

        GL_CorrectedFoldChanges<-matrix(NA,nrow = ng,ncol = 1,dimnames = list(genes,'GL_Cfc'))

        CL<-str_split(fc[i],'_')[[1]][1]
        ssam<-ncol(correctedFCs)
        for (j in 1:ng){
            print(c(i,CL,j,ng))
            id<-which(correctedFCs$genes==genes[j])
            GL_CorrectedFoldChanges[genes[j],1]<-mean(unlist(correctedFCs[id,'newFC']))
        }

        save(GL_CorrectedFoldChanges,file=paste(outDir,CL,'_GLcFC.Rdata',sep=''))
    }
}


bangmag.roc<-function(BFs,mgkdata_pval,mgkdata_fdr,ess_genes,non_ess_genes,CL,indSet=FALSE,th){
    Genes<-rownames(BFs)

    ess_genes<-intersect(ess_genes,Genes)
    non_ess_genes<-intersect(non_ess_genes,Genes)

    allg<-c(ess_genes,non_ess_genes)

    essentiality<-rep(0,length(allg))
    names(essentiality)<-allg

    essentiality[ess_genes]<-1
    predBFs<- BFs[allg,1]

    BF_roc<-roc(essentiality,predBFs-min(predBFs,na.rm=TRUE))

    if(indSet){
        TITLE<-'Independent gene sets'
    }else{
        TITLE<-'BAGEL reference gene sets'
    }

    plot(BF_roc,col='blue',lwd=3,main=TITLE)

    legendEntry<-c(paste('AUC = ',format(BF_roc$auc,digits=3),sep=''))
    legend('bottomright',legend=legendEntry,inset = c(0.1,0.1),bty = 'n',y.intersp = 2)

    bf_coords<-coords(BF_roc,'all',ret=c('threshold','sensitivity','specificity','ppv'),transpose=TRUE)
    bf_coords['threshold',]<-bf_coords['threshold',]+min(predBFs,na.rm=TRUE)

    id<-min(which(bf_coords['threshold',] > 0))
    bf_0<-c(0,bf_coords['specificity',id],bf_coords['sensitivity',id])
    points(bf_coords['specificity',id],bf_coords['sensitivity',id],col='blue',lwd=3,cex=1.5)

    id<-min(which(bf_coords['threshold',] > 1))
    bf_1<-c(0,bf_coords['specificity',id],bf_coords['sensitivity',id])
    points(bf_coords['specificity',id],bf_coords['sensitivity',id],col='blue',lwd=3,cex=1.5,pch=16)

    id<-min(which(bf_coords['ppv',]>(1-th)))
    if(id=="Inf"){
    id<-min(which(bf_coords['ppv',]>=(1-th)))}
    if(id=="Inf"){
      id<-max(which(round(bf_coords['ppv',])>=(1-th)))
    }
    bf_best_prec<-c(bf_coords['threshold',id],bf_coords['specificity',id],bf_coords['sensitivity',id],bf_coords['ppv',id])
    abline(h=bf_coords['sensitivity',id],col='blue',lwd=1)
        text(0.8,pos = 4,bf_coords['sensitivity',id]+0.015,
              paste('PPV > 95% (BF > ',format(bf_coords['threshold',id],digits=3),')',sep=''),col='blue',cex = 0.9)

    bestPrecisionTh<-rbind(bf_best_prec)
    colnames(bestPrecisionTh)<-c('thresholds','specificity','sensitivity','ppv')
    rownames(bestPrecisionTh)<-c('Ban_BF')

    return(list(bfAUC=BF_roc$auc,bestPrecisionTh))

}

bangmag.roc2<-function(BFs,mgkdata_pval,mgkdata_fdr,ess_genes,non_ess_genes,CL,indSet=FALSE,th){
  Genes<-rownames(BFs)

  ess_genes<-intersect(ess_genes,Genes)
  non_ess_genes<-intersect(non_ess_genes,Genes)

  allg<-c(ess_genes,non_ess_genes)

  essentiality<-rep(0,length(allg))
  names(essentiality)<-allg

  essentiality[ess_genes]<-1
  predBFs<- BFs[allg,1]

  BF_roc<-roc(essentiality,predBFs-min(predBFs,na.rm=TRUE))


  bf_coords<-coords(BF_roc,'all',ret=c('threshold','sensitivity','specificity','ppv'),transpose=TRUE)
  bf_coords['threshold',]<-bf_coords['threshold',]+min(predBFs,na.rm=TRUE)


  id<-min(which(bf_coords['ppv',]>(1-th)))
  if(id=="Inf"){
    id<-min(which(bf_coords['ppv',]>=(1-th)))}
  if(id=="Inf"){
    id<-max(which(round(bf_coords['ppv',])>=(1-th)))
  }
  bf_best_prec<-c(bf_coords['threshold',id],bf_coords['specificity',id],bf_coords['sensitivity',id],bf_coords['ppv',id])

  bestPrecisionTh<-rbind(bf_best_prec)
  colnames(bestPrecisionTh)<-c('thresholds','specificity','sensitivity','ppv')
  rownames(bestPrecisionTh)<-c('Ban_BF')

  return(list(bfAUC=BF_roc$auc,bestPrecisionTh))

}


bangmag.ind_test_sets<-function(){




    aguirreTest<-setdiff(Aguirre_essential_genes,union(BAGEL_essential,BAGEL_nonEssential))
    ribosomal<-setdiff(EssGenes.ribosomalProteins,union(BAGEL_essential,BAGEL_nonEssential))
    splice<-setdiff(EssGenes.SPLICEOSOME,union(BAGEL_essential,BAGEL_nonEssential))
    dnarep<-setdiff(EssGenes.DNA_REPLICATION,union(BAGEL_essential,BAGEL_nonEssential))

    ess_ind_test<-sort(unique(c(aguirreTest,ribosomal,splice,dnarep)))

    return(ess_ind_test)
}
bangmag.ind_test_sets_non_ess<-function(cellLine){

    cosmicID<-unique(as.character(PS_inventory$COSMIC_ID[PS_inventory$CellLineName==cellLine]))

    if(is.element(cosmicID,colnames(EXPpc))){
        nonExpGenes<-names(which(EXPpc[,cosmicID]<0.05))
    }else{
        nonExpGenes<-NULL
    }

    if(is.element(cosmicID,colnames(geneLevCNA$CNV))){
        ntmp<-geneLevCNA$CNV[,cosmicID]
        ntmp<-unlist(str_split(ntmp,','))
        max<-as.numeric(ntmp[seq(1,length(ntmp),4)])
        min<-as.numeric(ntmp[seq(2,length(ntmp),4)])
        HH<-ntmp[seq(3,length(ntmp),4)]
        deletedGenes<-rownames(geneLevCNA$CNV)[which(min==0 | max==0)]
    }else{
        deletedGenes<-NULL
    }

    non_ess_ind_test<-union(nonExpGenes,deletedGenes)
    return(non_ess_ind_test)
}

