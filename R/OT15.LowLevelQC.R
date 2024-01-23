bigcor <- function(
  x, 
  y = NULL,
  fun = "cor", 
  size = 50, 
  verbose = FALSE, 
  ...)
{
  fun <- match.arg(fun)
  if (fun == "cor") FUN <- cor else FUN <- cov
  if (fun == "cor") STR <- "Correlation" else STR <- "Covariance" 
  if (!is.null(y) & NROW(x) != NROW(y)) stop("'x' and 'y' must have compatible dimensions!")
  
  NCOL <- ncol(x)
  if (!is.null(y)) YCOL <- NCOL(y)
  
  ## calculate remainder, largest 'size'-divisible integer and block size
  REST <- NCOL %% size
  LARGE <- NCOL - REST  
  NBLOCKS <- NCOL %/% size
  
  ## preallocate square matrix of dimension
  ## ncol(x) in 'ff' single format
  if (is.null(y)) resMAT <- ff(vmode = "double", dim = c(NCOL, NCOL))  
  else resMAT <- ff(vmode = "double", dim = c(NCOL, YCOL))
  
  ## split column numbers into 'nblocks' groups + remaining block
  GROUP <- rep(1:NBLOCKS, each = size)
  if (REST > 0) GROUP <- c(GROUP, rep(NBLOCKS + 1, REST))
  SPLIT <- split(1:NCOL, GROUP)
  
  ## create all unique combinations of blocks
  COMBS <- expand.grid(1:length(SPLIT), 1:length(SPLIT))
  COMBS <- t(apply(COMBS, 1, sort))
  COMBS <- unique(COMBS)  
  if (!is.null(y)) COMBS <- cbind(1:length(SPLIT), rep(1, length(SPLIT)))
  
  ## initiate time counter
  #timeINIT <- proc.time() 
  
  ## iterate through each block combination, calculate correlation matrix
  ## between blocks and store them in the preallocated matrix on both
  ## symmetric sides of the diagonal
  for(i in 1:nrow(COMBS)) {
    COMB <- COMBS[i, ]    
    G1 <- SPLIT[[COMB[1]]]
    G2 <- SPLIT[[COMB[2]]]    
    
    ## if y = NULL
    if (is.null(y)) {
      if (verbose) cat(sprintf("#%d: %s of Block %s and Block %s (%s x %s) ... ", i, STR,  COMB[1],
                               COMB[2], length(G1),  length(G2)))      
      RES <- FUN(x[, G1], x[, G2], ...)
      resMAT[G1, G2] <- RES
      resMAT[G2, G1] <- t(RES) 
    } else ## if y = smaller matrix or vector  
    {
      if (verbose) cat(sprintf("#%d: %s of Block %s and 'y' (%s x %s) ... ", i, STR,  COMB[1],
                               length(G1),  YCOL))    
      RES <- FUN(x[, G1], y, ...)
      resMAT[G1, ] <- RES             
    }
    
    if (verbose) {
      timeNOW <- proc.time() - timeINIT
      cat(timeNOW[3], "s\n")
    }
    
    gc()
  } 
  
  return(resMAT)
}


bigcorPar <- function(x, nblocks = 10, verbose = TRUE, ncore="all", ...){


  NCOL <- ncol(x)
  
  ## test if ncol(x) %% nblocks gives remainder 0
  if (NCOL %% nblocks != 0){stop("Choose different 'nblocks' so that ncol(x) %% nblocks = 0!")}
  
  ## preallocate square matrix of dimension
  ## ncol(x) in 'ff' single format
  corMAT <- ff(vmode = "single", dim = c(NCOL, NCOL))
  
  ## split column numbers into 'nblocks' groups
  SPLIT <- split(1:NCOL, rep(1:nblocks, each = NCOL/nblocks))
  
  ## create all unique combinations of blocks
  COMBS <- expand.grid(1:length(SPLIT), 1:length(SPLIT))
  COMBS <- t(apply(COMBS, 1, sort))
  COMBS <- unique(COMBS)
  
  ## iterate through each block combination, calculate correlation matrix
  ## between blocks and store them in the preallocated matrix on both
  ## symmetric sides of the diagonal
  results <- foreach(i = 1:nrow(COMBS)) %dopar% {
    COMB <- COMBS[i, ]
    G1 <- SPLIT[[COMB[1]]]
    G2 <- SPLIT[[COMB[2]]]
    if (verbose) cat("Block", COMB[1], "with Block", COMB[2], "\n")
    flush.console()
    COR <- cor(x[, G1], x[, G2], ...)
    corMAT[G1, G2] <- COR
    corMAT[G2, G1] <- t(COR)
    COR <- NULL
  }
  
  gc()
  return(corMAT)
}
getCorRep<-function(guides,rawDataset){
  comparisons<-as.dist(cor(t(rawDataset[guides,])))
  comparisons<-tril(as.matrix(comparisons))
  mean(comparisons)
}

selectReproducibleGuides<-function(rawDataset,libAnnotation,th=0.6,removeEssential=FALSE,essentials=NULL,incgeneset=NULL,subsample=FALSE,numberguides=NA){
    genes<-unique(libAnnotation$GENES)
    if(removeEssential){
      genes<-setdiff(genes,essentials)
    }
    ngenes<-length(genes)
    
    reproducibleGuides<-NULL
    if(!is.null(incgeneset)){
      genes<-incgeneset
      ngenes<-length(genes)
    }
    reproducibleGuides<-foreach(i=1:ngenes)%dopar%{
        
        #print(c(i,length(reproducibleGuides)))
        currentGuides<-intersect(rownames(rawDataset),
                                 rownames(libAnnotation)[libAnnotation$GENES==genes[i]])
        
        comparisons<-as.dist(cor(t(rawDataset[currentGuides,])))
        comparisons<-tril(as.matrix(comparisons))
        
        repG<-currentGuides[unique(c(which(comparisons>th,arr.ind = TRUE)))]
        
        if(length(repG)>0){
            repG
        }
        
    }
    if(!is.na(numberguides)){

      corlist<-list()
      corlist<-lapply(reproducibleGuides,function(x) getCorRep(x,rawDataset))
      corlist<-unlist(corlist)
      rorder<-order(corlist,decreasing=TRUE)
      totalg<-0
      idx<-1
      reproducibleGuides<-reproducibleGuides[rorder]
      while(totalg<numberguides){
        totalg<-totalg+length(reproducibleGuides[idx][[1]])
        idx=idx+1
      }
      idx<-idx-1
      reproducibleGuides<-unlist(reproducibleGuides[1:idx])
    }else{
      reproducibleGuides<-unlist(reproducibleGuides)
    }
    if(subsample&length(reproducibleGuides)>1000){
      nguides<-length(reproducibleGuides)
      reproducibleGuidesList<-list()
      reproducibleGuidesList<-sapply(1:1000,sample(1:nguides,800))
      reproducibleGuidesList<-lapply(reproducibleGuidesList,function(x) reproducibleGuides[x])
      return(list(reproducibleGuides,reproducibleGuidesList))
    }else{
    
    return(list(reproducibleGuides,NULL))}
}


OT015_BG_corrEstimates<-function(inventory,Dataset,outpath,fn,xlim=c(-0.2,1)){
    
    #loggedDataset<-log10(Dataset[,3:ncol(Dataset)]+1)
 
    loggedDataset<-Dataset
    
    
    #ii<-which(apply(loggedDataset,MARGIN = 1,FUN='sd')>1)
    #loggedDataset<-loggedDataset[ii,]
    
    PWc<-as.matrix(cor(loggedDataset,use = 'complete.obs'))
    
    PWc[lower.tri(PWc,diag = TRUE)]<-NA
 
    IDX<-which(!is.na(PWc),arr.ind = TRUE)
    PWvalues<-PWc[!is.na(PWc)]
    comparisons<-cbind(colnames(PWc)[IDX[,1]],colnames(PWc)[IDX[,2]])
    
    sameLibrary<-inventory$LIBRARY[match(comparisons[,1],inventory$SAMPLE_ID)]==
                 inventory$LIBRARY[match(comparisons[,2],inventory$SAMPLE_ID)]
    
    sameTissue<-inventory$TISSUE[match(comparisons[,1],inventory$SAMPLE_ID)]==
        inventory$TISSUE[match(comparisons[,2],inventory$SAMPLE_ID)]
    
    sameCellLine<-inventory$CellLineName[match(comparisons[,1],inventory$SAMPLE_ID)]==
        inventory$CellLineName[match(comparisons[,2],inventory$SAMPLE_ID)]
    
    pdf(paste(outpath,fn,'.pdf',sep=''),8.5,11)
    par(mfrow=c(3,1))
    
    ccr.multDensPlot(list(AllComparisons=density(PWvalues),
                          sameLibrary=density(PWvalues[sameLibrary]),
                          differentLibrary=density(PWvalues[!sameLibrary])),
                     COLS = c('gray','red','pink'),
                     XLIMS=xlim,
                     TITLE = TITLE,LEGentries = c('all pw comparisons',
                                                     'same library',
                                                     'different library'))
    
    abline(v=median(PWvalues[sameLibrary]),col='red',lwd=2)
    abline(v=median(PWvalues[!sameLibrary]),col='pink',lwd=2)
    
    ccr.multDensPlot(list(AllComparisons=density(PWvalues),
                          sameTissue=density(PWvalues[sameTissue]),
                          differentTissue=density(PWvalues[!sameTissue])),
                     COLS = c('gray','blue','cyan'),
                     XLIMS=xlim,
                     TITLE = TITLE,LEGentries = c('all pw comparisons',
                                                                                        'same tissue',
                                                                                        'different tissue'))
    abline(v=median(PWvalues[sameTissue]),col='blue',lwd=2)
    abline(v=median(PWvalues[!sameTissue]),col='cyan',lwd=2)
    
    ccr.multDensPlot(list(AllComparisons=density(PWvalues),
                          sameTissue=density(PWvalues[sameCellLine]),
                          differentTissue=density(PWvalues[!sameCellLine])),
                     COLS = c('gray','darkgreen','green'),
                     XLIMS=xlim,
                     TITLE = TITLE,LEGentries = c('all pw comparisons',
                                                                                        'same cell line',
                                                                                        'different cell line'))
    abline(v=median(PWvalues[sameCellLine]),col='darkgreen',lwd=2)
    abline(v=median(PWvalues[!sameCellLine]),col='green',lwd=2)
    dev.off()
      
    pdf(paste(outpath,fn,'_mixed.pdf',sep=''),8.5,11)
    par(mfrow=c(2,1))  
    ccr.multDensPlot(list(AllComparisons=density(PWvalues),
                          dcstsl=density(PWvalues[!sameCellLine & !sameTissue & sameLibrary]),
                          differentTissue=density(PWvalues[sameCellLine])),
                     COLS = c('gray','purple','darkgreen'),
                     XLIMS=xlim,
                     TITLE = TITLE,LEGentries = c('all pw comparisons',
                                                                                        'different cell line different tissue same library',
                                                                                        'same cell line'))
    abline(v=median(PWvalues[!sameCellLine & !sameTissue & sameLibrary]),col='purple',lwd=2)
    abline(v=median(PWvalues[sameCellLine]),col='darkgreen',lwd=2)
    
    ccr.multDensPlot(list(AllComparisons=density(PWvalues),
                          dcstsl=density(PWvalues[!sameCellLine & sameTissue & sameLibrary]),
                          differentTissue=density(PWvalues[sameCellLine])),
                     COLS = c('gray','blue','darkgreen'),
                     XLIMS=xlim,
                     TITLE = TITLE,LEGentries = c('all pw comparisons',
                                                                                        'different cell line same tissue same library',
                                                                                        'same cell line'))
    
    abline(v=median(PWvalues[!sameCellLine & sameTissue & sameLibrary]),col='blue',lwd=2)
    abline(v=median(PWvalues[sameCellLine]),col='darkgreen',lwd=2)
    dev.off()
    
    CELLline1<-manifest$CellLineName[match(comparisons[,1],manifest$SAMPLE_ID)]
    CELLline2<-manifest$CellLineName[match(comparisons[,2],manifest$SAMPLE_ID)]
    
    PWvalues_same<-PWvalues[CELLline1==CELLline2]
    comparisons_same<-comparisons[CELLline1==CELLline2,]
    
    nSamples<-vector()
    nCellLines<-vector()
    flag<-1
    for (i in seq(0.95,0.45,-0.025)){
       
        
        nSamples[flag]<-length(unique(c(comparisons_same[which(PWvalues_same>=i),])))
        
        
        nCellLines[flag]<-length(unique(
            manifest$CellLineName[match(unique(c(comparisons_same[which(PWvalues_same>=i),])),manifest$SAMPLE_ID)]))
        flag<-flag+1
    }
    
    par(mfrow=c(1,2))
    plot(seq(0.95,0.45,-0.025),nSamples,xlab='R threshold',ylab='n.samples',type='b')
    abline(v=0.5,lty=3)
    abline(v=0.625,lty=2)
    plot(seq(0.95,0.45,-0.025),nCellLines,xlab='R threshold',ylab='n.CellLines',type='b')
    abline(v=0.5,lty=3)
    abline(v=0.625,lty=2)
  
    
}

OT015_individual_logFCprofiles<-function(inventory,Dataset,controlnames){
    
    #identify and remove control samples from main matrix
    if(site=="Broad"){
    	removenames<-colnames(Dataset)[grep("pDNA",colnames(Dataset))]
    	DatasetSample<-Dataset[,!is.element(colnames(Dataset),removenames)]
    }else{

    		DatasetSample<-Dataset[,!is.element(colnames(Dataset),controlnames)]
    }

    DatasetSample<-DatasetSample[,is.element(colnames(DatasetSample),inventory$SAMPLE_ID)]
    nsamples<-ncol(DatasetSample)
    

    logFCs<-foreach(i=1:nsamples,.combine=cbind)%dopar%{
        #print(paste('processing sample ', colnames(DatasetSample)[i],', ',i,' out of ',nsamples))
        #change this to control sample matching instead
        #so when have multiple controls, expect these to be listed as individual files in one vector in CONTROLS separated by comma (necessary for Broad data)
        Lib<-unlist(strsplit(as.character(controlnames[i]),":"))
       
        
        controlsamples<-Dataset[,make.names(Lib)]

        if(length(Lib)==1){
        	subDataset<-cbind(controlsamples,DatasetSample[,i])
        }else{
        	control<-apply(controlsamples,1,median)
        	subDataset<-cbind(control,DatasetSample[,i])
        }
 
        numd<-subDataset
        
        numd[which(numd[,1]<30),]<-NA
        normFact<-t(matrix(rep(colSums(numd,na.rm = TRUE),nrow(numd)),2,nrow(numd)))
        
        numd<-numd/normFact*10000000
        
        #logFCs<-cbind(logFCs,log2((numd[,2]+0.5)/(numd[,1]+0.5)))
        t<-log2((numd[,2]+0.5)/(numd[,1]+0.5))
   
        return(t)
        #colnames(logFCs)[ncol(logFCs)]<-colnames(DatasetSample)[i]
    }
    colnames(logFCs)<-as.character(colnames(DatasetSample))
    logFCs<-logFCs[which(rowSums(is.na(logFCs))==0),]

    return(logFCs)
}


OT015_individual_GL_logFCprofiles<-function(Dataset,Guide_library){



    gDataset<-foreach(i=1:ncol(Dataset),.combine=cbind,.packages="CRISPRcleanR")%dopar%{
        #print(paste('processing sample ', colnames(Dataset)[i],', ',i,' out of ',ncol(Dataset),sep=''))
        
        ccr.geneMeanFCs(Dataset[,i],Guide_library)
        
        
    }
    colnames(gDataset)<-colnames(Dataset)
    
    return(gDataset)
    #save(gDataset,file='../../20171207_PAPERfreeze_data/sgRNA_counts/lowLevelDatasets/individualFCs_GeneLevel_20171207.Rdata')
}




OT015_lowLevel_allPW_corPlots<-function(inventory,Dataset,BGcors,
                                        path='../../RESULTS/lowLevelQC/October2016/replicate_correlation_Oct2016/'){
    
    if(!dir.exists(path)){
        dir.create(path)
    }
    
    Dataset<-Dataset[which(rowSums(is.na(Dataset))==0),]
    
    CellLines<-names(which(summary(as.factor(inventory$CellLineName),maxsum = length(unique(inventory$CellLineName)))>1))
    
    totalAVG<-NULL
    
    for (Cell_line in CellLines){
        
        f1<-paste(path,Cell_line,'_clustering.pdf',sep='')
        f2<-paste(path,Cell_line,'_pwcorr.png',sep='')
        f3<-paste(path,Cell_line,'_empP.pdf',sep='')
        
        inv_idxs<-which(inventory$CellLineName==Cell_line)
        
        tmp<-inventory$CGaP_Sample_Name[inv_idxs]
        
        OT015_lowLevel_clustering(inventory = inventory,
                                  figFile=f1,figEight = 105,
                                  MAIN = Cell_line,
                                  Dataset = Dataset,
                                  toHighLight = tmp)
        
        OT015_lowLevel_repCorr(inventory = inventory,
                               Dataset = Dataset,
                               cellLine = Cell_line,
                               figFile=f2)
        
        Distr<-OT015_lowLevel_repVsBackGround_corr(inventory = PS_inventory,
                                                   Dataset = Dataset,
                                                   MAIN='Replicates vs background R',
                                                   toHighLight = Cell_line,
                                                   figFile=f3,BGcors = BGcors)
        
        
        unsamples<-sort(union(colnames(Distr$corValues),rownames(Distr$corValues)))
        
        avgCorr<-vector()
        
        for (us in 1:length(unsamples)){
            ttmp<-
                c(Distr$corValues[unsamples[us],],Distr$corValues[,unsamples[us]])
            avgCorr[us]<-mean(ttmp[!is.na(ttmp)])
        }
        
        names(avgCorr)<-unsamples
        
        if (length(totalAVG)==0){
            totalAVG<-avgCorr
        }else{
            totalAVG<-c(totalAVG,avgCorr)
        }
    }
    
    AVG_rep_CORR<-rep(NA,nrow(inventory))
    invIDXS<-match(names(totalAVG),inventory$CGaP_Sample_Name)
    AVG_rep_CORR[invIDXS]<-totalAVG               
    inventory<-cbind(inventory,AVG_rep_CORR)
    colnames(inventory)[ncol(inventory)]<-'avg_rep_corr'
    
    return(inventory)
    
}

OT015_lowLevel_clustering<-function(inventory,
                                    Dataset,
                                    toHighLight=NULL,
                                    highLibrary=FALSE,
                                    highTissue=FALSE,
                                    highCtype=FALSE,
                                    highReps=FALSE,
                                    CEXLAB=2,
                                    CEXAX=1,
                                    CEXNAMES=1,
                                    CEXLEG=1.5,
                                    CEXMAIN=1,MAIN='',figFile=NULL,figEight=50,figWidth=10,whatToHighLight=NULL){

    Dataset<-
        Dataset[which(rowSums(is.na(Dataset))==0),is.element(colnames(Dataset),inventory$SAMPLE_ID)]
   #print(dim(Dataset))
    PD<-as.dist(1-cor(Dataset))

    hc<-hclust(PD,method = 'complete')
    
    dend <- as.dendrogram(hc)
    
    idxs<-match(labels(dend),inventory$SAMPLE_ID)
    
    if(length(toHighLight)>0){
        labels_colors(dend)<-'black'
        labels_cex(dend)<-1
        labels_colors(dend)[!is.element(labels(dend),toHighLight)]<-'lightgray'
        labels_cex(dend)[!is.element(labels(dend),toHighLight)]<-0.5
    }

    if(highLibrary){
        labels_colors(dend)<-'blue'
        
        toHighLight<-inventory$SAMPLE_ID[which(inventory$LIBRARY=="V1.1")]
        labels_colors(dend)[is.element(labels(dend),toHighLight) | labels(dend)=='Plasmid V1.1']<-'red'
    }
    
    if(highTissue){
        
        tmp<-TissueColors[inventory$TISSUE[match(labels(dend),inventory$SAMPLE_ID)]]
        
        tmp[is.na(tmp)]<-'black'
        
        labels_colors(dend)<-tmp
        
    }
    
    if(highCtype){
        
        tmp<-TissueColors[manifest$CANCER_TYPE[match(labels(dend),inventory$SAMPLE_ID)]]
        
        tmp[is.na(tmp)]<-'black'
        
        labels_colors(dend)<-tmp
        
    }
    
    if(highReps){
        
        ucells<-unique(inventory$CellLineName)
        cellC<-rainbow(length(ucells))
        names(cellC)<-ucells
        
        tmp<-cellC[inventory$CellLineName[match(labels(dend),inventory$SAMPLE_ID)]]
        labels_colors(dend)<-tmp
        
    }
    
    if(length(whatToHighLight)>0){
        
        ufactor<-unique(inventory[,whatToHighLight])
        factorC<-rainbow(length(ufactor))
        names(factorC)<-ufactor
        
        tmp<-factorC[inventory[match(labels(dend),inventory$SAMPLE_ID),whatToHighLight]]
        labels_colors(dend)<-tmp
        
    }
    
    
    if (length(figFile)){
        pdf(figFile,figWidth,figEight)    
    }
    
    par(mar=c(4,4,8,12))
    par(cex=CEXNAMES)
    plot(dend,horiz = TRUE,xlab='distance',cex.axis=CEXAX,cex.lab=CEXLAB,cex.main=CEXMAIN,main=MAIN)
    
    if(highLibrary){
       legend('topleft',legend=c(unique(inventory$LIBRARY)),cex=CEXLEG,title='sgRNA Library',lty=1,lwd=10,col=c('red','blue'))
    }
    
    if(highTissue){
        plot(0,0,col=NA,xaxt='n',yaxt='n',frame.plot = FALSE,xlab='',ylab='')
        legend('topleft',legend=unique(inventory$TISSUE),
               cex=CEXLEG,title='Tissue',lty=1,lwd=10,col=TissueColors[unique(inventory$TISSUE)])
    }
    
    if(highCtype){
        plot(0,0,col=NA,xaxt='n',yaxt='n',frame.plot = FALSE,xlab='',ylab='')
        legend('topleft',legend=unique(inventory$CANCER_TYPE),
               cex=CEXLEG,title='Tissue',lty=1,lwd=10,col=TissueColors[unique(inventory$CANCER_TYPE)])
    }
    
    if(highReps){
        legend('topleft',legend='',
               cex=CEXLEG,title='different colors = different cell lines')
    }
    
    if(length(whatToHighLight)>0){
        
        legend('topleft',legend=ufactor,lty=1,col=factorC,bty='n',
               cex=CEXLEG,title=whatToHighLight)
        # 
        # ufactor<-unique(inventory[,whatToHighLight])
        # factorC<-rainbow(length(ufactor))
        # names(factorC)<-ufactor
        # 
        # tmp<-factorC[inventory[match(labels(dend),inventory$CGaP_Sample_Name),whatToHighLight]]
        # labels_colors(dend)<-tmp
        # 
    }
    
    if (length(figFile)){
        dev.off()    
    }
}

OT015_lowLevel_selectCorrSet<-function(domain,inventory,Dataset,toRemove=NULL){
    
    Dataset<-
        Dataset[rowSums(is.na(Dataset))==0,]
    
    PD<-cor(Dataset)
    
    PD[lower.tri(PD,diag = TRUE)]<-NA
    
    cset<-NULL
    
    if (domain=='all'){
        cset<-c(PD)
        cset<-cset[!is.na(cset)]
    }
    
    
    if (domain=='diff_ctype'){
        nr<-nrow(inventory)
        for (i in 1:nr){
            current_Id<-
                inventory$SAMPLE_ID[i]
            
            current_ctype<-inventory$CANCER_TYPE[i]
            
            other_Id<-
                inventory$SAMPLE_ID[which(inventory$CANCER_TYPE!=current_ctype)]
        
            current_c_set<-
                c(PD[current_Id,other_Id],PD[other_Id,current_Id])
            
            if (i == 1){
                cset<-current_c_set
            }else{
                cset<-c(cset,current_c_set)    
            }
        }
        cset<-cset[!is.na(cset)]
    }
    if (domain=='same_ctype'){
        nr<-nrow(inventory)
        for (i in 1:nr){
            current_Id<-
                inventory$SAMPLE_ID[i]
            
            current_ctype<-inventory$CANCER_TYPE[i]
            
            other_Id<-
                inventory$SAMPLE_ID[which(inventory$CANCER_TYPE==current_ctype)]
            
            current_c_set<-
                c(PD[current_Id,other_Id],PD[other_Id,current_Id])
            
            if (i == 1){
                cset<-current_c_set
            }else{
                cset<-c(cset,current_c_set)    
            }
        }
        cset<-cset[!is.na(cset)]
    }
    
    if (domain=='diff_cell'){
        nr<-nrow(inventory)
        for (i in 1:nr){
            current_Id<-
                inventory$SAMPLE_ID[i]
            
            current_cell<-inventory$INSTITUTE_ID[i]
            
            other_Id<-
                inventory$SAMPLE_ID[which(inventory$INSTITUTE_ID!=current_cell)]
           
            current_c_set<-c(PD[current_Id,other_Id],PD[other_Id,current_Id])
            
            if (i == 1){
                cset<-current_c_set
            }else{
                cset<-c(cset,current_c_set)    
            }
        }
        cset<-cset[!is.na(cset)]
    }
    
    if (domain=='diff_library'){
        nr<-nrow(inventory)
        for (i in 1:nr){
            current_Id<-
                inventory$SAMPLE_ID[i]
            
            current_cell<-inventory$LIBRARY[i]
            
            other_Id<-
                inventory$SAMPLE_ID[which(inventory$LIBRARY!=current_cell)]
            
            current_c_set<-c(PD[current_Id,other_Id],PD[other_Id,current_Id])
            
            if (i == 1){
                cset<-current_c_set
            }else{
                cset<-c(cset,current_c_set)    
            }
        }
        cset<-cset[!is.na(cset)]
    }
    
    if (domain=='same_cell'){
        nr<-nrow(inventory)
        for (i in 1:nr){
            current_Id<-
                inventory$SAMPLE_ID[i]
            
            current_cell<-inventory$INSTITUTE_ID[i]
            
            if (length(toRemove)==0 | !is.element(current_cell,toRemove)){
                other_Id<-
                    inventory$SAMPLE_ID[which(inventory$INSTITUTE_ID==current_cell)]
                
                current_c_set<-
                    c(PD[current_Id,other_Id],PD[other_Id,current_Id])
                
                if (length(cset) == 0){
                    cset<-current_c_set
                }else{
                    cset<-c(cset,current_c_set)    
                }
            }
        }
        cset<-cset[!is.na(cset)]
    }
    
     
    if (domain=='same_library'){
        nr<-nrow(inventory)
        for (i in 1:nr){
            current_Id<-
                inventory$SAMPLE_ID[i]
            
            current_cell<-inventory$LIBRARY[i]
            
            if (length(toRemove)==0 | !is.element(current_cell,toRemove)){
                other_Id<-
                    inventory$SAMPLE_ID[which(inventory$LIBRARY==current_cell)]
                
                current_c_set<-
                    c(PD[current_Id,other_Id],PD[other_Id,current_Id])
                
                if (length(cset) == 0){
                    cset<-current_c_set
                }else{
                    cset<-c(cset,current_c_set)    
                }
            }
        }
        cset<-cset[!is.na(cset)]
    }
    
      
    
    return(cset)
    
}






OT015_lowLevel_repVsBackGround_corr<-function(inventory,Dataset,toHighLight=NULL,MAIN='',figFile=NULL,BGcors){
    
    
    allScores<-OT015_lowLevel_selectCorrSet(domain = 'all',Dataset = Dataset,inventory = inventory)
    diff_cell_corr_scores<-OT015_lowLevel_selectCorrSet(domain = 'diff_cell',Dataset = Dataset,inventory = inventory)
    
    same_cell_corr_scores<-OT015_lowLevel_selectCorrSet(domain = 'same_cell',
                                                        Dataset = Dataset,
                                                        inventory = inventory,
                                                        toRemove=toHighLight)
	compareLibraries<-length(unique(inventory$LIBRARY))>1
	
    if(length(toHighLight)==0){
    	tempinventory<-inventory
    	tempinventory[tempinventory$CANCER_TYPE=="","CANCER_TYPE"]<-NA
    	tC<-table(tempinventory$CANCER_TYPE)
    	removeCt<-names(tempinventory)[tC>3]
		tempinventory[tempinventory$CANCER_TYPE%in%removeCt,"CANCER_TYPE"]<-NA 
		subinventory<-tempinventory[!is.na(tempinventory$CANCER_TYPE),]  
		nCtype<-length(unique(subinventory$CANCER_TYPE)) 
		if(nCtype>2){	
        diff_ctype_corr_scores<-OT015_lowLevel_selectCorrSet(domain = 'diff_ctype',Dataset = Dataset,inventory = subinventory)
        same_ctype_corr_scores<-OT015_lowLevel_selectCorrSet(domain = 'same_ctype',Dataset = Dataset,inventory = subinventory,
                                                        toRemove=toHighLight)}
    if(compareLibraries){
        diff_library_corr_scores<-OT015_lowLevel_selectCorrSet(domain = 'diff_library',Dataset = Dataset,inventory = inventory)}
        same_library_corr_scores<-OT015_lowLevel_selectCorrSet(domain = 'same_library',Dataset = Dataset,inventory = inventory,
                                                         toRemove=toHighLight)
          if(compareLibraries){                                               
        TRES_lib<-t.test(diff_library_corr_scores,same_library_corr_scores)}
        if(nCtype>2){
        TRES_ctype<-t.test(diff_ctype_corr_scores,same_ctype_corr_scores)}
        
    }
    
    TRES_cell<-t.test(diff_cell_corr_scores,same_cell_corr_scores)
    
    if (length(figFile)){
        pdf(figFile,8,10)    
    }
    
    if (length(toHighLight)==0){
        
        par(mfrow=c(1,3))
        par(mar=c(10,4,4,0))
        if(compareLibraries){
        boxplot(list(diff_library_corr_scores,same_library_corr_scores),names = c('Diff Library','Same Library'),col=c('gray','blue',NA),frame.plot=FALSE,
                ylab='R',main=paste(MAIN,'p = ',format(TRES_lib$p.value,scientific=TRUE,digits=3),sep=''),ylim=c(0,1),las=2,outline = FALSE)}
               if(nCtype>2){
        boxplot(list(diff_ctype_corr_scores,same_ctype_corr_scores),names = c('Diff Cancer Type','Same Cancer Type'),col=c('gray','blue',NA),frame.plot=FALSE,
                ylab='R',main=paste(MAIN,'p = ',format(TRES_ctype$p.value,scientific=TRUE,digits=3),sep=''),ylim=c(0,1),las=2,outline=FALSE)}
        boxplot(list(diff_cell_corr_scores,same_cell_corr_scores),names = c('Diff Cell Lines','Same Cell Lines'),col=c('gray','blue',NA),frame.plot=FALSE,
                ylab='R',main=paste(MAIN,'p = ',format(TRES_cell$p.value,scientific=TRUE,digits=3),sep=''),ylim=c(0,1),las=2,outline = FALSE)
        
        if (length(figFile)){
            dev.off() 
        }
        return(list(BGscores=diff_cell_corr_scores,REPscores=same_cell_corr_scores))
    }else{
        
        #PD<-cor(Dataset)
        #PD<-PD[setdiff(rownames(PD),"ERS717283.plasmid"),setdiff(colnames(PD),"ERS717283.plasmid")]
        #PD[lower.tri(PD,diag = TRUE)]<-NA
        
        ids<-inventory$SAMPLE_ID[which(inventory$CellLineName==toHighLight)]
        
        subM<-PD[ids,ids]
        
        nn<-which(!is.na(subM),arr.ind = TRUE)
        toHighLightSet<-subM[nn]
        rowN<-rownames(subM)[nn[,1]]
        colN<-colnames(subM)[nn[,2]]
        
        
        names(toHighLightSet)<-paste(rowN,'vs',colN)
        
        COL<-rep('darkgreen',length(toHighLightSet))
        COL[toHighLightSet<MAXm]<-'gold1'
        COL[toHighLightSet<minM]<-'red'
        
        boxplot(list(diff_cell_corr_scores,same_cell_corr_scores,toHighLightSet),names = c('Bg','Other Reps',paste(toHighLight,'Reps')),
                ylab='R',main=paste(MAIN,'',sep=''),ylim=c(0,1),col=c('gray','blue',NA),outline = FALSE)
        points(rep(3,length(toHighLightSet)),toHighLightSet,ylim=c(0,1),pch=21,cex=2,bg=COL,col='gray')
        
        if (length(figFile)){
            dev.off() 
        }
        
        return(list(BGscores=diff_cell_corr_scores,REPscores=same_cell_corr_scores,corValues=subM))
    }
    
    
    
}


OT015_repQuality<-function(BGcors,valuesTotest){
    
    nval<-length(valuesTotest)
    
    scores<-rep(NA,nval)
    
    for (i in 1:nval){
        scores[i]<- log10(perc.rank(BGcors$REPscores,valuesTotest[i])/perc.rank(BGcors$BGscores,valuesTotest[i]))
    }
    
    return(scores)
}
OT015_lowLevel_repCorr<-function(inventory,Dataset,cellLine,figFile){
    inventory<-
        inventory[which(inventory$fns!='data not available'),]
    
    png(figFile,1182,1131)
    
    idxs<-inventory$CGaP_Sample_Name[which(inventory$CellLineName==cellLine)]
    
    subData<-Dataset[,idxs]
    
    tmp_name<-inventory$CGaP_Sample_Name[which(inventory$CellLineName==cellLine)]
    
    colnames(subData)<-tmp_name
    
    fc<-subData
    
    nfc<-ncol(fc)
    
    pairs(as.formula(paste('~',paste(colnames(fc),collapse='+'),sep='')),
          data=fc,lower.panel=panel.smooth,upper.panel=panel.cor,
          pch=20)
    
    dev.off()
    
}
OT015_lowLevel_corrThresholds<-function(BGcors){
    RANGE<-seq(0,1,0.0001)
    qc<-vector()
    for (i in 1:length(RANGE)){
        
        qc[i]<-OT015_repQuality(BGcors = BGcors,valuesTotest = RANGE[i])    
    }
    
    MAXm<-RANGE[min(which(qc>=1))]
    minM<-0.5
    
    return(list(MAXm=MAXm, minM=minM))
}

OT015_lowLevel_corrSummaries<-function(inventory,lowLevDataset,figFileName,minM=0.5,TITLE){
    
    lowLevDataset<-lowLevDataset[which(rowSums(is.na(lowLevDataset))==0),]
    
    cellLines<-unique(inventory$Cell_Line_Name)
    
    ncellLines<-length(cellLines)
    
    nreps<-vector()
    avgCorr<-vector()
    
    for (i in 1:ncellLines){
        currSname<-
            inventory$sample_id[inventory$Cell_Line_Name==cellLines[i]]
        
        nreps[i]<-length(currSname)
        if(nreps[i]>1){
            cmat<-cor(lowLevDataset[,currSname])
            avgCorr[i]<-mean(cmat[lower.tri(cmat)])    
        }else{
            avgCorr[i]<-0.01
        }
    }
    
    names(nreps)<-cellLines
    names(avgCorr)<-cellLines
    
    nreps<-nreps[order(avgCorr,decreasing=TRUE)]
    avgCorr<-sort(avgCorr,decreasing=TRUE)
    
    pdf(figFileName,width = 35,height = 8)
    par(mar=c(8,4,4,1))
    avgCorr[is.na(avgCorr)]<-0
    coords<-barplot(avgCorr,col=TissueColors[inventory$GDSC_Description_2[match(names(avgCorr),inventory$Cell_Line_Name)]],
            ylim=c(0,1),
            las=2,
            ylab='avg R',
            main=TITLE,
            border = NA)
    
    
    abline(h = minM,col='darkgreen',lwd=2,lty=2)
    text(coords,avgCorr,nreps,pos = 3)
    dev.off()
    
    return(avgCorr)
}