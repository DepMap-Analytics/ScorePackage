##### Miscellaneous
plotPCnnmd<-function(inputdata,essential,nonessential,pcno=1,colvec,shapevec,colcode,shapecode,legval=""){
  nnmds<-nnmd(inputdata,essential,nonessential)
  PCval<-ExtractPCs(inputdata,pcno)
  plot(PCval,nnmds,col=colvec,pch=shapevec)
  legend("topleft",legend=legval,col=colcode,pch=shapecode)
  print(paste("Cor:",cor(PCval,nnmds,method="spearman")))
}
BoxplotQC<-function(inputdata,essential,nonessential){
  boxplot(inputdata[intersect(rownames(inputdata),essential),])
  boxplot(inputdata[intersect(rownames(inputdata),nonessential),])
}
SLAPE.heuristic_mut_ex_sorting<-function (mutPatterns)
{
  mutPatterns <- sign(mutPatterns)
  if (dim(mutPatterns)[1] == 1) {
    mutPatterns <- matrix(c(mutPatterns[, order(mutPatterns,
                                                decreasing = TRUE)]), 1, ncol(mutPatterns), dimnames = list(rownames(mutPatterns),
                                                                                                            colnames(mutPatterns)))
    return(mutPatterns)
  }
  if (dim(mutPatterns)[2] == 1) {
    mutPatterns <- matrix(c(mutPatterns[order(mutPatterns,
                                              decreasing = TRUE), ]), nrow(mutPatterns), 1, dimnames = list(rownames(mutPatters),
                                                                                                            colnames(mutPatterns)))
    return(mutPatterns)
  }
  nsamples <- ncol(mutPatterns)
  coveredGenes <- NA
  uncoveredGenes <- rownames(mutPatterns)
  if (length(uncoveredGenes) > 1) {
    idNull <- which(colSums(mutPatterns) == 0)
    nullCol <- matrix(c(mutPatterns[, idNull]), nrow(mutPatterns),
                      length(idNull), dimnames = list(rownames(mutPatterns),
                                                      colnames(mutPatterns)[idNull]))
    idNonNull <- which(colSums(mutPatterns) > 0)
    mutPatterns <- matrix(c(mutPatterns[, idNonNull]), nrow(mutPatterns),
                          length(idNonNull), dimnames = list(rownames(mutPatterns),
                                                             colnames(mutPatterns)[idNonNull]))
    coveredSamples <- NA
    uncoveredSamples <- colnames(mutPatterns)
    BS <- NA
    while (length(uncoveredGenes) > 0 & length(uncoveredSamples) >
           0) {
      patterns <- matrix(c(mutPatterns[uncoveredGenes,
                                       uncoveredSamples]), nrow = length(uncoveredGenes),
                         ncol = length(uncoveredSamples), dimnames = list(uncoveredGenes,
                                                                          uncoveredSamples))
      if (length(uncoveredGenes) > 1) {
        bestInClass <- SLE.findBestInClass(patterns)
      }
      else {
        bestInClass <- uncoveredGenes
      }
      if (is.na(BS[1])) {
        BS <- bestInClass
      }
      else {
        BS <- c(BS, bestInClass)
      }
      if (is.na(coveredGenes[1])) {
        coveredGenes <- bestInClass
      }
      else {
        coveredGenes <- c(coveredGenes, bestInClass)
      }
      uncoveredGenes <- setdiff(uncoveredGenes, coveredGenes)
      toCheck <- matrix(c(patterns[bestInClass, uncoveredSamples]),
                        nrow = 1, ncol = ncol(patterns), dimnames = list(bestInClass,
                                                                         uncoveredSamples))
      if (length(coveredGenes) == 1) {
        coveredSamples <- names(which(colSums(toCheck) >
                                        0))
      }
      else {
        coveredSamples <- c(coveredSamples, names(which(colSums(toCheck) >
                                                          0)))
      }
      uncoveredSamples <- setdiff(uncoveredSamples, coveredSamples)
    }
    BS <- c(BS, uncoveredGenes)
    CID <- SLE.rearrangeMatrix(mutPatterns, BS)
    FINALMAT <- mutPatterns[BS, CID]
    FINALMAT <- cbind(FINALMAT, nullCol[rownames(FINALMAT),
                                        ])
    return(FINALMAT)
  }
}

my.boxplot<-function(toPlot,MAIN,NAMES,COLORS=NULL,MANIFEST=NULL,YLAB='',YLIM=NULL){

    toPlot<-lapply(toPlot,function(x){x[!is.na(x)]})

    if(length(COLORS)>0){

        pwc<-list()
        for (i in 1:length(toPlot)){
            pwc[[i]]<-COLORS[CELL_LINE_INVENTORY$Tissue[match(names(toPlot[[i]]),CELL_LINE_INVENTORY$Name)]]
        }

    }

    if(length(YLIM)==0){
        YLIM<-range(unlist(toPlot))
    }
        beeswarm(toPlot,pch=16,pwcol = pwc,corral = 'wrap',
                 ylab=YLAB,labels = NAMES,
                 ylim=range(unlist(toPlot)))
        par(new=TRUE)
        boxplot(toPlot,main=MAIN,names = NAMES,ylim=range(unlist(toPlot)),outline=FALSE,frame.plot=FALSE,xaxt='n',yaxt='n',
                main=MAIN,col=NA,lwd=2,
                ylim=range(unlist(toPlot)))

    MEANS<-unlist(lapply(toPlot,FUN = 'mean',na.rm=TRUE))
    SD<-unlist(lapply(toPlot,FUN = 'sd',na.rm=TRUE))
    for (i in 1:length(MEANS)){
        lines(x=c(i-0.3,i+0.3),y=c(MEANS[i],MEANS[i])+SD[i],col='purple',lwd=3)
        lines(x=c(i-0.2,i+0.2),y=c(MEANS[i],MEANS[i]),col='red',lwd=5)
        lines(x=c(i-0.3,i+0.3),y=c(MEANS[i],MEANS[i])-SD[i],col='purple',lwd=3)
    }

    if(length(toPlot)==2){
        print(t.test(toPlot[[1]],toPlot[[2]]))
    }

}


col2hex <- function(cname)
{
    colMat <- col2rgb(cname)
    rgb(
        red=colMat[1,]/255,
        green=colMat[2,]/255,
        blue=colMat[3,]/255
    )
}

makeTransparent<-function(someColor, alpha=100)
{
    newColor<-col2rgb(someColor)
    apply(newColor, 2, function(curcoldata){rgb(red=curcoldata[1], green=curcoldata[2],
                                                blue=curcoldata[3],alpha=alpha, maxColorValue=255)})
}


makePieChart<-function(percPres,colorP,border){

    df <- data.frame(
        group = c("present", "absent"),
        value = c(percPres, 100-percPres)
    )

    blank_theme <- theme_minimal()+
        theme(
            axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            panel.border = element_blank(),
            panel.grid=element_blank(),
            axis.ticks = element_blank(),
            plot.title=element_text(size=14, face="bold"),
            legend.position="none"
        )

    bp<- ggplot(df, aes(x="", y=value, fill=group))+
        geom_bar(width = 5, stat = "identity",color=border,size=1)+coord_polar("y", start=0) + blank_theme +
        theme(axis.text.x=element_blank())+scale_fill_manual(values=colorP)
    return(bp)
}

convert_CCLE_CNAdata<-function(){

    fc<-read.table('../../DATA/raw/absoluteCNA/data_CNA.txt',sep='\t',header=TRUE)

    rn<-fc[,1]
    fc<-fc[,4:ncol(fc)]

    rownames(fc)<-rn

    commonC<-read.xlsx('../../DATA/raw/absoluteCNA/CCLE_GDSC_cellLineMappoing.xlsx',sheetIndex = 1)

    nn<-colnames(fc)
    id<-match(nn,commonC$CCLE.name)

    fc<-fc[,!is.na(id)]
    id<-id[!is.na(id)]

    cids<-
        commonC$GDSC1000.cosmic.id[id]

    colnames(fc)<-cids


}

panel.cor <- function(x, y, digits=2, prefix="", cex.cor, ...)
{

    usr <- par("usr"); on.exit(par(usr))
    par(usr = c(0, 1, 0, 1))
    r <- abs(cor(x, y))
    txt <- format(c(r, 0.123456789), digits=digits)[1]
    txt <- paste(prefix, txt, sep="")
    if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)

    COL<-rep('darkgreen',length(r))
    COL[which(r<MAXm)]<-'gold1'
    COL[which(r<minM)]<-'red'

    text(0.5, 0.5, txt, cex = cex.cor * r,col=COL)
}

perc.rank <- function(x, xo)  length(x[x >= xo])/length(x)
generateCombosIdxs<-function(n){

    combos<-list()

    flag<-1

    for (i in 1:n){

        CB<-combn(1:n,i)

        for (j in 1:ncol(CB)){
            combos[[flag]]<-CB[,j]
            flag<-flag+1
        }
    }

    return(combos)
}
my.hypTest<-function(x,k,n,N){


    PVALS<-phyper(x-1,n,N-n,k,lower.tail=FALSE)

    return(PVALS)
}
my.compress_identical_patterns<-function(DATA){

    flag<-1

    while(flag){
        JD<-as.matrix(dist(DATA,method='binary'))
        JD[lower.tri(JD,diag=TRUE)]<-Inf

        toMix<-which(JD==0,arr.ind=TRUE)

        if (length(toMix)>0){
            toMix<-toMix[1,]
            tm1<-row.names(JD)[toMix[1]]
            tm2<-row.names(JD)[toMix[2]]
            cat(paste('\n\t\t\t\t\t\t\tMerging',tm1,'and',tm2))
            #added 29.6.19 to allow for BEM matrix collapsed to single row.
            if(is.matrix(DATA)){
              rownames(DATA)[toMix[2]]<-paste(tm1,tm2,sep=', ')
              DATA<-DATA[-toMix[1],]}else{
              names(DATA)<-paste(tm1,tm2,sep=", ")
            }
        }else{
            flag<-0
        }

    }

    print('Done!')
    return(DATA)
}

rearrangeMatrix<-function(patterns,GENES){

    remainingSamples<-colnames(patterns)

    toAdd<-NULL

    for (g in GENES){
        remainingGenes<-setdiff(GENES,g)

        P1<-matrix(c(patterns[g,remainingSamples]),length(g),length(remainingSamples),dimnames = list(g,remainingSamples))
        P2<-matrix(c(patterns[remainingGenes,remainingSamples]),length(remainingGenes),length(remainingSamples),
                   dimnames=list(remainingGenes,remainingSamples))

        if(length(remainingGenes)>1){
            DD<-colnames(P1)[order(P1-colSums(P2),decreasing=TRUE)]
        }else{
            DD<-colnames(P1)[order(P1-P2,decreasing=TRUE)]
        }

        toAdd<-c(toAdd,names(which(patterns[g,DD]>0)))
        remainingSamples<-setdiff(remainingSamples,toAdd)
        if(length(remainingSamples)==0){
            break
        }
    }

    toAdd<-c(toAdd,remainingSamples)

    return(toAdd)
}

findBestInClass<-function(patterns){

    if(nrow(patterns)==1){
        return(rownames(patterns))
    }

    if(ncol(patterns)==1){
        return(rownames(patterns)[1])
    }

    genes<-rownames(patterns)

    exclCov<-rep(NA,length(genes))
    names(exclCov)<-genes
    for (g in genes){
        residGenes<-setdiff(genes,g)
        if (length(residGenes)>1){
            exclCov[g]<-sum(patterns[g,]-colSums(patterns[residGenes,]))
        }else{
            exclCov[g]<-sum(patterns[g,]-patterns[residGenes,])
        }
    }

    return(names(sort(exclCov,decreasing=TRUE))[1])
}

HeuristicMutExSorting<-function(mutPatterns){

    mutPatterns<-sign(mutPatterns)

    if(dim(mutPatterns)[1]==1){
        mutPatterns<-matrix(c(mutPatterns[,order(mutPatterns,decreasing=TRUE)]),
                            1,ncol(mutPatterns),
                            dimnames = list(rownames(mutPatterns),colnames(mutPatterns)))

        return(mutPatterns)
    }

    if(dim(mutPatterns)[2]==1){
        mutPatterns<-matrix(c(mutPatterns[order(mutPatterns,decreasing=TRUE),]),
                            nrow(mutPatterns),1,
                            dimnames = list(rownames(mutPatters),colnames(mutPatterns)))
        return(mutPatterns)
    }

    nsamples<-ncol(mutPatterns)

    coveredGenes<-NA
    uncoveredGenes<-rownames(mutPatterns)

    if (length(uncoveredGenes)>1){

        idNull<-which(colSums(mutPatterns)==0)
        nullCol<-matrix(c(mutPatterns[,idNull]),nrow(mutPatterns),length(idNull),dimnames = list(rownames(mutPatterns),colnames(mutPatterns)[idNull]))

        idNonNull<-which(colSums(mutPatterns)>0)
        mutPatterns<-matrix(c(mutPatterns[,idNonNull]),nrow(mutPatterns),length(idNonNull),dimnames=list(rownames(mutPatterns),colnames(mutPatterns)[idNonNull]))

        coveredSamples<-NA
        uncoveredSamples<-colnames(mutPatterns)
        BS<-NA

        while(length(uncoveredGenes)>0 & length(uncoveredSamples)>0){

            patterns<-matrix(c(mutPatterns[uncoveredGenes,uncoveredSamples]),
                             nrow = length(uncoveredGenes),
                             ncol = length(uncoveredSamples),
                             dimnames = list(uncoveredGenes,uncoveredSamples))

            if(length(uncoveredGenes)>1){
                bestInClass<-findBestInClass(patterns)
            }else{
                bestInClass<-uncoveredGenes
            }

            if(is.na(BS[1])){
                BS<-bestInClass
            }else{
                BS<-c(BS,bestInClass)
            }

            if(is.na(coveredGenes[1])){
                coveredGenes<-bestInClass
            }else{
                coveredGenes<-c(coveredGenes,bestInClass)
            }

            uncoveredGenes<-setdiff(uncoveredGenes,coveredGenes)
            toCheck<-matrix(c(patterns[bestInClass,uncoveredSamples]),nrow = 1,ncol=ncol(patterns),dimnames = list(bestInClass,uncoveredSamples))

            if (length(coveredGenes)==1){
                coveredSamples<-names(which(colSums(toCheck)>0))
            }else{
                coveredSamples<-c(coveredSamples,names(which(colSums(toCheck)>0)))
            }

            uncoveredSamples<-setdiff(uncoveredSamples,coveredSamples)

        }

        BS<-c(BS,uncoveredGenes)

        CID<-rearrangeMatrix(mutPatterns,BS)

        FINALMAT<-mutPatterns[BS,CID]

        FINALMAT<-cbind(FINALMAT,nullCol[rownames(FINALMAT),])

        return(FINALMAT)
    }

}



my.quick_GSEA<-function(RANKING,VALUES,SIGNATURE,objectsType='',valueType='',signature_name='',
                        p=1,display=TRUE,returnRS=FALSE){
    layout(mat = c(1,2),heights = c(7,4))
    par(mar=c(4,4,2,4))
    SIGNATURE<-unique(intersect(SIGNATURE,RANKING))

    HITS<-is.element(RANKING,SIGNATURE)+0
    R<-VALUES[RANKING]*HITS

    hitCases<-cumsum(abs(R)^p)
    NR<-sum(abs(R)^p)

    missCases<-cumsum(1-HITS)

    N<-length(RANKING)
    N_Nh<-length(RANKING)-length(SIGNATURE)

    Phit<-hitCases/NR
    Pmiss<-missCases/N_Nh

    m<-max(abs(Phit-Pmiss))
    t<-which(abs(Phit-Pmiss)==m)

    if (length(t)>1){t<-t[1]}

    peak<-t
    ES<-Phit[t]-Pmiss[t]
    RS<-Phit-Pmiss

    names(ES)<-NULL

    if (display){
        if (ES>=0){
            c<-"darkgreen"
            YLIM=c(0,abs(ES)+0.10*(abs(ES)))
        }else{
            c<-"darkgreen"
            YLIM=c(-(abs(ES)+0.10*(abs(ES))),0)
        }

        plot(0:N,c(0,Phit-Pmiss),col=c,type="l",lwd=5,xlim=c(0,N),
             ylim=YLIM,xaxs="i",bty="l",axes=FALSE,
             xlab="",ylab="Enrichment score (ES)",
             main=paste('Enrichment plot:\n',signature_name))
        par(new=TRUE)
        plot(0:N,rep(0,N+1),col='gray',type="l",new=FALSE,xlab="",ylab="",ylim=YLIM,xlim=c(0,N),xaxt='n')

        #     axis(side=2)
        #
        col=colorRampPalette(c('blue','white','red'))(100)

        tmp<-VALUES[RANKING]
        col<-col[ceiling(99*(tmp-min(tmp))/(max(tmp)-min(tmp)))+1]

        plot(0:N,rep(0,N+1),col='gray',type="l",new=FALSE,xlab="Rank in Odered Dataset",ylab="",ylim=c(min(VALUES),max(VALUES)),axes=FALSE)
        abline(v=which(HITS>0))
        par(new=TRUE)
        plot(VALUES[RANKING],ylab=valueType,pch=16,col=col,cex=1,xlab='')


    }


    P<-NA

    if (returnRS){
        POSITIONS<-which(HITS==1)
        names(POSITIONS)<-ranking[which(HITS==1)]

        POSITIONS<-POSITIONS[order(names(POSITIONS))]
        names(POSITIONS)<-names(POSITIONS)[order(names(POSITIONS))]

        return(list(ES=ES,RS=RS,POSITIONS=POSITIONS,PEAK=t))
    } else {return(ES)}
}

my.GSEA<-function(RANKING,VALUES,SIGNATURES,objectsType='',valueType='',signature_name='',
                  p=1,display=TRUE,nperm=1000){

    Ng<-length(SIGNATURES)


    phi <- matrix(nrow = Ng, ncol = nperm)
    phi.norm <- matrix(nrow = Ng, ncol = nperm)
    obs.phi <- matrix(nrow = Ng, ncol = nperm)

    OES<-vector()

    for (i in 1:Ng){

        OES[i]<-my.quick_GSEA(RANKING,VALUES,SIGNATURES[[i]],objectsType,valueType,signature_name,p,display)

        for (j in 1:nperm){
            print(j)
            reshuffledList<-RANKING[sample(1:length(RANKING))]

            phi[i,j]<-my.quick_GSEA(reshuffledList,VALUES,SIGNATURES[[i]],objectsType,valueType,signature_name,p,display=FALSE)
        }
    }

    print("Computing nominal p-values...")

    p.vals <- matrix(0, nrow = Ng, ncol = 2)

    for (i in 1:Ng) {

        pos.phi<-phi[i,which(phi[i,]>=0)]
        neg.phi<-phi[i,which(phi[i,]<0)]

        ES.value <- OES[i]
        if (ES.value >= 0) {

            NULLES <- fitdist(pos.phi,"gamma")
            pval<-1-pgamma(ES.value,shape=NULLES$estimate[1],rate=NULLES$estimate[2])

            #
            #      p.vals[i, 1] <- signif(sum(pos.phi >= ES.value)/length(pos.phi), digits=5)
        } else {
            NULLES <- fitdist(-neg.phi,"gamma")
            pval<-1-pgamma(-ES.value,shape=NULLES$estimate[1],rate=NULLES$estimate[2])

        }
    }

    return(list(ES=ES.value,pval=pval))

}


ELnormaliser<-function(signal){

    PDF<-density(signal,width=sd(signal)/4)
    PDF<-approxfun(PDF$x,PDF$y,yleft=0,yright=0)

    nsamples<-length(signal)

    normsignal<-rep(NA,length(signal))

    for (j in 1:nsamples){
        print(j)
        PLL<-quad(PDF,0,signal[j])

        CLL<-quad(PDF,signal[j],max(signal)+1)

        if (PLL>1){PLL<-1}
        if (CLL==0){CLL<-.Machine$double.eps}

        normsignal[j]<-log(PLL/CLL)
    }

    return(normsignal)

}


multDensPlot<-function(TOPLOT, COLS, XLIMS, TITLE, LEGentries, XLAB = "",legendPosition='topright')
{
    YM <- vector()
    for (i in 1:length(TOPLOT)) {
        YM[i] <- max(TOPLOT[[i]]$y, na.rm = TRUE)
    }
    Ymax <- max(YM, na.rm = TRUE)
    plot(0, 0, col = NA, ylab = "density", xlab = XLAB, xlim = XLIMS,
         ylim = c(0, Ymax), type = "l", main = TITLE)
    for (i in 1:length(TOPLOT)) {
        cord.x <- c(TOPLOT[[i]]$x)
        cord.y <- c(TOPLOT[[i]]$y)
        rgbc <- col2rgb(COLS[i])
        currCol <- rgb(rgbc[1], rgbc[2], rgbc[3], alpha = 100,
                       maxColorValue = 255)
        polygon(cord.x, cord.y, col = currCol, border = NA)
        lines(TOPLOT[[i]], col = COLS[i], lwd = 3)
        if (i == 1) {
            legend(legendPosition, legend = LEGentries, col = COLS,
                   lwd = 3, bty = "n")
        }
    }
}
plotCase<-function(genomicFeature,DepletedGene,assocId=NULL,MSI=FALSE){

  xlim<-c(min(ESSprofiles[,DepletedGene],na.rm = TRUE),
          max(ESSprofiles[,DepletedGene],na.rm = TRUE))

  if (MSI){
    featPattern<-InputFeatures$MSI_VARIABLE
  }else{
    featPattern<-InputFeatures$BEM[genomicFeature,]
  }


  #featPattern<-InputFeatures$MSI_VARIABLE
  cids<-manifest$Cell_Line_Name[match(rownames(ESSprofiles),manifest$COSMIC_ID)]

  tissues<-manifest$Cancer.Type[match(rownames(ESSprofiles),manifest$COSMIC_ID)]
  cols<-TissueColors[tissues]

  sbfPattern<- -sBFs[DepletedGene,cids]
  qnCorrectedlogFC<- ESSprofiles[,DepletedGene]

  ylim<-c(min(sbfPattern,na.rm = TRUE),
          max(sbfPattern,na.rm = TRUE))

  PCH<-rep(21,length(sbfPattern))
  PCH[which(sbfPattern<0)]<-23

  beeswarm(qnCorrectedlogFC~featPattern,pwbg  = makeTransparent(cols,200),
           pwpch = PCH,
           corral = "wrap",
           ylab='QN corrected log FC',labels=c('wt','alt'),xlab=genomicFeature,ylim=xlim,
           main=paste(assocId,': ',DepletedGene,'depletion'),cex=1.2)

  par(new=TRUE)
  boxplot(qnCorrectedlogFC~featPattern,ylim=xlim,frame.plot=FALSE,xaxt='n',yaxt='n',col=NA,lwd=2,outline=FALSE,
          border=makeTransparent('black',200))



}

globalVolcano<-function(TOTRES){

  y<- -log10(as.numeric(TOTRES[,"FEATURE_ANOVA_pval"]))
  x<- sign(as.numeric(TOTRES[,"FEATURE_deltaMEAN_ESS"]))*as.numeric(TOTRES[,"FEATURE_ESS_effect_size"])

  smoothScatter(x,y,nrpoints = 0,frame.plot=FALSE,xlab='signed effect size',ylab='-log p',bandwidth = c(1,1),
                colramp = colorRampPalette(c("white", 'black')),nbin = c(400,400))

  id<-which(as.numeric(TOTRES[,"ANOVA FEATURE FDR %"])<25 & as.numeric(TOTRES[,"FEATURE_ESS_effect_size"])>1)

  COL<-rep('darkgreen',length(id))
  COL[which(x[id]>0)]<-'red'

  br<-rep('darkgreen',length(id))
  br[which(x[id]>0)]<-'darkred'

  abline(v= -1,col='gray')
  abline(v= 1,col='gray')

  abline(h= -log10(max(as.numeric(TOTRES[which(as.numeric(TOTRES[,"ANOVA FEATURE FDR %"])<25),"FEATURE_ANOVA_pval"]))),
         lty=4,col='gray')

  abline(h= -log10(max(as.numeric(TOTRES[which(as.numeric(TOTRES[,"ANOVA FEATURE FDR %"])<10),"FEATURE_ANOVA_pval"]))),
         lty=2,col='gray')

  abline(h= -log10(max(as.numeric(TOTRES[which(as.numeric(TOTRES[,"ANOVA FEATURE FDR %"])<5),"FEATURE_ANOVA_pval"]))),
         lty=1,col='gray')

  points(x[id],y[id],bg=makeTransparent(COL),pch=21,col=br,cex=1.2)

  par(xpd=TRUE)
  identify(x[id],y[id],paste(TOTRES[id,"FEATURE"],':',TOTRES[id,"Depleted Gene"],'_dep',sep=''),cex=0.5)
}
makeNameMap<-function(allgenes,allSymbol){
  symbolS1<-allSymbol[allSymbol$Match.type=="Approved symbol",]
  symbolS2<-allSymbol[allSymbol$Match.type=="Previous symbol",]
  symbolS3<-allSymbol[allSymbol$Match.type=="Synonyms",]
  allMap<-symbolS1[match(allgenes,symbolS1$Input),"Approved.symbol"]
  names(allMap)<-allgenes
  m1<-names(allMap)[is.na(allMap)]
  m1M<-symbolS2[match(m1,symbolS2$Input),"Approved.symbol"]
  names(m1M)<-m1
  m2<-names(m1M)[is.na(m1M)]
  m2M<-symbolS3[match(m2,symbolS3$Input),"Approved.symbol"]
  names(m2M)<-m2

  allMap[names(m1M)]<-m1M
  allMap[names(m2M)]<-m2M


  allMap[is.na(allMap)]<-names(allMap)[is.na(allMap)]
  return(allMap)
}
updateNames<-function(inputdata,Map,uset=TRUE){

  newnames<-Map[inputdata]
  names(newnames)<-inputdata
  newnames[is.na(newnames)]<-names(newnames)[is.na(newnames)]
  if(uset){
    newnames<-newnames[newnames%in%names(table(newnames))[table(newnames)==1]]
  }
  return(newnames)
}
updateRownames<-function(inputdata,Map,allSymbol=NULL){

  newnames<-Map[rownames(inputdata)]
  names(newnames)<-rownames(inputdata)
  newnames[is.na(newnames)]<-names(newnames)[is.na(newnames)]
  if(!is.null(allSymbol)){
    set1<-newnames[newnames%in%names(table(newnames))[table(newnames)==1]]
    set2<-newnames[newnames%in%names(table(newnames))[table(newnames)>1]]
    usenames<-names(set1)
    if(length(set2)>0){
      uset<-unique(set2)
      for(i in uset){
        check<-set2[set2==i]
        res<-allSymbol[allSymbol$Input%in%names(check),]
        sel<-res[res$Match.type=="Approved symbol","Input"]
        usenames<-c(usenames,sel)
      }
    }
    newnames<-newnames[usenames]
  }else{
    newnames<-newnames[newnames%in%names(table(newnames))[table(newnames)==1]]
  }
  inputdata<-inputdata[names(newnames),]
  rownames(inputdata)<-newnames
  return(inputdata)
}
CLnameMapping<-function(inputdata,refdata,compare=c("col","row"),annotation=NULL){
  #assume both using same ID info. For cross-matching need a reference column in CMP model list
  #annotation is model_list_latest from CMP
  #refdata is the one that will have its column names matched to the colnames of the input data
  compare<-match.arg(compare)
  if(compare=="col"){
    colnames(inputdata)<-make.names(colnames(inputdata))
    colnames(refdata)<-make.names(colnames(refdata))
    inref<-intersect(colnames(inputdata),colnames(refdata))
    if(length(inref)==ncol(inputdata)){
      return(list(inputdata=inputdata,refdata=refdata[,colnames(inputdata)]))
    }else{
      missing<-setdiff(colnames(inputdata),colnames(refdata))
      #if annotation is NULL don't do anything as have already used 'bridging' e.g cosmic id.
      if(!is.null(annotation)){
        missingsanger<-missing[grep("SIDM",missing)]
        missingbroad<-missing[grep("ACH",missing)]
        annotation$BROAD_ID<-make.names(annotation$BROAD_ID)
        matchmatrix<-NULL
        if(length(missingsanger)>0){
          matchmatrix<-cbind(missingsanger,annotation[match(missingsanger,annotation$model_id),c("BROAD_ID")])

        }
        if(length(missingbroad)>0){
          temp<-cbind(missingbroad,annotation[match(missingbroad,annotation$BROAD_ID),c("model_id")])
          matchmatrix<-rbind(matchmatrix,temp)
        }
        #update refdata to the matchmatrix values then redo intersection between input and ref data
        orignames<-colnames(refdata)
        colnames(refdata)<-matchmatrix[match(colnames(refdata),matchmatrix[,2]),1]
        colnames(refdata)[is.na(colnames(refdata))]<-orignames[is.na(colnames(refdata))]
        inref<-intersect(colnames(inputdata),colnames(refdata))
        inputdata<-inputdata[,inref]
        refdata<-refdata[,inref]
        #check for second column in the refdata
        matchmatrix<-matchmatrix[matchmatrix[,1]%in%colnames(refdata),]

      }
    }
  }
  if(compare=="row"){

  }
  return(list(inputdata=inputdata,refdata=refdata,matchmatrix=matchmatrix))
}
# SIGNATURE<-
#     rownames(pvals_ach)[sample(nrow(pvals_ach))[1:100]]
#
# # # ####! DEMO !
# VALUES1<- -log10(pvals_ach[,'A375_1234'])
# VALUES1<-VALUES1[unique(names(VALUES1))]
# RANKING1<-names(sort(VALUES1,decreasing=TRUE))
# # ## SIGNATURE should be a list containing th cosmic Ids of the cell lines you want to query
# #
# my.quick_GSEA(RANKING1,VALUES1,SIGNATURE,
#               objectsType = 'Genes',valueType = 'Mageck -log10 p',signature_name = 'Broad A375')
#
# VALUES2<- -log10(pvals_sanger[,'A375_123'])
# VALUES2<-VALUES2[unique(names(VALUES2))]
# RANKING2<-names(sort(VALUES2,decreasing=TRUE))
# # ## SIGNATURE should be a list containing th cosmic Ids of the cell lines you want to query
# # SIGNATURE<-
# #     ribosomalProteins
# #
# my.quick_GSEA(RANKING2,VALUES2,SIGNATURE,
#               objectsType = 'Genes',valueType = 'Mageck -log10 p',signature_name = 'Sanger A375')
#
#
