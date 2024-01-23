library(pROC)

bagelR.createAllInputFiles<-function(inventory,outDir,control_header = NULL,forceLibraryV1=FALSE){

    cellLines<-unique(inventory$Cell_Line_Name)
    mcellLines<-length(cellLines)

    for (i in 1:mcellLines){

        print(cellLines[i])
        
        Library<-unique(inventory$Library[which(is.element(inventory$Cell_Line_Name,cellLines[i]))])

        if (length(Library)>1){
            Library<-"Human V1"
        }
        
        if(length(control_header)==0){
            if(Library=="Human V1.1"){
                controlHeader<-'CRISPR_C6596666.sample'
            }else{
                controlHeader<-'ERS717283.plasmid'
            }    
        }else{
            controlHeader<-control_header
        }
        

        if(forceLibraryV1){
            Library<-"Human V1"
        }
        OT015_AllReplicate_MageckInput(inventory = inventory,
                                       executeMageck_and_remove = FALSE,
                                       control_header = controlHeader,
                                       CL = cellLines[i],
                                       Library = Library,
                                       outPath=outDir)
    }
}

bagelR.createAllFCFiles<-function(inputDir,outputDir,MIN_READS=30){
    fn<-dir(inputDir)
    fn<-grep('_mgk_input.tsv',fn,value = TRUE)

    expnames<-str_split(fn,'_mgk_input.tsv')
    expnames<-unlist(expnames)
    expnames<-expnames[seq(1,length(expnames),2)]
    for (i in 1:length(fn)){
        bagelR.foldChanges(FILENAME = paste(inputDir,fn[i],sep=''),
                           MIN_READS=MIN_READS,DISPLAY = TRUE,OUTDIR = outputDir,
                           EXPname = expnames[i])
    }
}

bagelR.round_to_hundredth<-function(x){
    return (round(x*100) / 100.0)
}

bagelR.foldChanges<-function(FILENAME,MIN_READS,DISPLAY,OUTDIR,EXPname){

    pdf(paste(OUTDIR,EXPname,'_prePostNorm.pdf',sep=''),width = 10,height = 5)
    par(mfrow=c(1,2))
    counts<-read.table(FILENAME,sep='\t',header=TRUE,stringsAsFactors = FALSE)

    if(DISPLAY){
        my.boxplot(counts[,3:ncol(counts)],MAIN = paste(EXPname,'Raw sgRNA counts'),NAMES = c('plsmd',paste('r',1:(ncol(counts)-3))))
    }

    numd<-
        counts[,3:ncol(counts)]

    normFact<-t(matrix(rep(colSums(numd),nrow(numd)),ncol(counts)-2,nrow(numd)))

    numd<-numd/normFact*10000000

    normed<-cbind(counts[,1:2],numd)

    if(DISPLAY){
        my.boxplot(numd,MAIN = paste(EXPname,'Raw sgRNA normalised counts'),NAMES = c('plsmd',paste('r',1:(ncol(counts)-3))))
    }

    dev.off()
    IDX<-which(counts[,3]>=MIN_READS)

    normed<-normed[IDX,]

    nsamples<-ncol(counts)-3

    for ( i in 1:nsamples){

        c_foldchanges<-log2((normed[,3+i]+0.5)/(normed[,3]+0.5))
        c_foldchanges<-matrix(c_foldchanges,length(c_foldchanges),1,dimnames = list(normed$sgRNA,colnames(normed)[3+i]))

        if (i == 1){
            foldchanges<-c_foldchanges
        }else{
            foldchanges<-cbind(foldchanges,c_foldchanges )
        }
    }

    pdf(paste(OUTDIR,EXPname,'_fcs.pdf',sep=''),width = 5,height = 5)
    
    if(DISPLAY){
        my.boxplot(foldchanges,MAIN = paste(EXPname,'log Fold Changes'),NAMES = paste('r',1:(ncol(counts)-3)))
    }

    dev.off()
    foldchanges<-cbind(normed[,1:2],foldchanges)

    save(normed,file=paste(OUTDIR,EXPname,'_normCounts.Rdata',sep=''))
    save(foldchanges,file=paste(OUTDIR,EXPname,'_foldChanges.Rdata',sep=''))
}

bagelR.prePostNormPlots<-function(inventory,inputDir,prenormDir){
    
    
    outputDir<-paste(inputDir,'__plots/',sep='')
    
    if(!exists(outputDir)){
        dir.create(outputDir)
    }
    
    cellLines<-unique(inventory$Cell_Line_Name)
 
    ncellLines<-length(cellLines)
    
    for (i in 1:ncellLines){
        par(mfrow=c(2,1))
        par(mar=c(2,4,2,2))
        fc<-read.table(paste(prenormDir,cellLines[i],'_mgk_input.tsv',sep=''),sep='\t',header=TRUE)
        boxplot(fc[,3:ncol(fc)],las=2,names = c('plsmd',paste('r',1:(ncol(fc)-3))),
                main=paste(cellLines[i],'pre/post normalisation'))
        
        load(paste(inputDir,'A375_normCounts.Rdata',sep=''))
        boxplot(normed[,3:ncol(fc)],las=2,names = rep('',ncol(fc)-2))
        
    }
}

bagelR.computeAllGuidesBFs_v2<-function(cellLine,
                                     NUM_BOOTSTRAPS = 1000,
                                     ESSENTIAL_GENES,
                                     NON_ESSENTIAL_GENES,
                                     inputFolder='../../DATA/R/normalisedCounts_and_FCs/',
                                     outputFolder='../../RESULTS/BAGEL-R_output/',
                                     refGuidesLibrary,diagnosticPlots=TRUE,
                                     whatToTest='newFC'){
    
    fn<-paste(inputFolder,cellLine,'_CCRoutput.RData',sep='')
    load(fn)
    
    print('Computing sgRNA boostrapped Bayesian Factors and classification performances for corrected FCs...')
                    
    correctedFCs<-correctedFCs$corrected_logFCs
    
    correctedFCs<-cbind(rownames(correctedFCs),correctedFCs)
    colnames(correctedFCs)[1]<-'sgRNA'
    
    set.seed(0xA5EED)
    
    rep1<-bagelR.sgRNA_bootStrapped_BFactors(FOLD_CHANGES = correctedFCs,
                                             NUM_BOOTSTRAPS = NUM_BOOTSTRAPS,
                                             ESSENTIAL_GENES = ESSENTIAL_GENES,
                                             NON_ESSENTIAL_GENES = NON_ESSENTIAL_GENES,
                                             SAMPLE_TO_TEST = whatToTest,
                                             refGuidesLibrary = refGuidesLibrary,
                                             compPerformances = TRUE)
    sgRNA_BFs<-rep1
    print('DONE')

    save(sgRNA_BFs,file=paste(outputFolder,cellLine,'_sgRNAs_BFs.Rdata',sep=''))
    if (diagnosticPlots){
        bagelR.diagnosticPlots(sgRNA_BFs,cellLine = cellLine,outDir = outputFolder)
    }
    
    return(sgRNA_BFs)
}

bagelR.sgRNA_bootStrapped_BFactors<-function(FOLD_CHANGES,
                                             ESSENTIAL_GENES,
                                             NON_ESSENTIAL_GENES,
                                             NUM_BOOTSTRAPS=1000,
                                             SAMPLE_TO_TEST,
                                             FC_DENS_THRESHOLD=2^-7,
                                             percTrainingSample=80,
                                             refGuidesLibrary,
                                             compPerformances=TRUE){
    
    rownames(FOLD_CHANGES)<-FOLD_CHANGES$sgRNA
    
    nguides<-nrow(FOLD_CHANGES)
    
    genes<-unique(FOLD_CHANGES$gene)
    guides<-FOLD_CHANGES$sgRNA
    
    ngenes<-length(genes)
    
    print(paste('n. of sgRNA =',nguides))
    print(paste('n. of unique genes =',ngenes))
    
    print(paste('n. of considered essential genes = ',length(intersect(ESSENTIAL_GENES,FOLD_CHANGES$gene)),
                ' (out of ',length(ESSENTIAL_GENES),' in reference set)',sep=''))
    print(paste('          targeted by a total number of ',length(which(is.element(FOLD_CHANGES$gene,ESSENTIAL_GENES))),
                ' sgRNAs',sep=''))
    
    
    print(paste('n. of considered non essential genes = ',length(intersect(NON_ESSENTIAL_GENES,FOLD_CHANGES$gene)),
                ' (out of ',length(NON_ESSENTIAL_GENES),' in reference set)',sep=''))
    print(paste('          targeted by a total number of ',length(which(is.element(FOLD_CHANGES$gene,NON_ESSENTIAL_GENES))),
                ' sgRNAs',sep=''))
    
    BFs_across_loops<-matrix(-Inf,nguides,NUM_BOOTSTRAPS,dimnames = list(guides,1:NUM_BOOTSTRAPS))
    BFs_across_loops_including_training<-matrix(-Inf,nguides,NUM_BOOTSTRAPS,dimnames = list(guides,1:NUM_BOOTSTRAPS))
    
    EssentialGuides<-as.character(FOLD_CHANGES$sgRNA[is.element(FOLD_CHANGES$gene,ESSENTIAL_GENES)])
    nonEssentialGuides<-as.character(FOLD_CHANGES$sgRNA[is.element(FOLD_CHANGES$gene,NON_ESSENTIAL_GENES)])
    
    nEssGuides<-length(EssentialGuides)
    nNonEssGuides<-length(nonEssentialGuides)
    
    print('Bootstrap iterations in progress...')
    
    pb <- txtProgressBar(min=1,max=NUM_BOOTSTRAPS,style=3)
    
    for (loop in 1:NUM_BOOTSTRAPS){
        setTxtProgressBar(pb, loop)
        guide_train_ess<-EssentialGuides[sample(nEssGuides)[1:ceiling(nEssGuides*percTrainingSample/100)]]
        guide_train_non_ess<-nonEssentialGuides[sample(nNonEssGuides)[1:ceiling(nNonEssGuides*percTrainingSample/100)]]
        
        guide_training<- union(guide_train_ess,guide_train_non_ess)
        guide_test <- setdiff(guides,guide_training)
        
        ess_train_fc <- as.numeric(FOLD_CHANGES[guide_train_ess,SAMPLE_TO_TEST])
        non_ess_train_fc <- as.numeric(FOLD_CHANGES[guide_train_non_ess,SAMPLE_TO_TEST])
        
        kess<-density(ess_train_fc, kernel = "gaussian")
        knon<-density(non_ess_train_fc, kernel = "gaussian")
        
        x <- seq(-10,2,0.01)
        nonfitx <- approx(knon$x,knon$y,x)$y
        
        f <- which(nonfitx > FC_DENS_THRESHOLD)
        xmin <- bagelR.round_to_hundredth( min(x[f]) )
        
        subx <- seq(xmin,max(x[f]),0.01)
        logratio_sample <- log2( approx(kess$x,kess$y,subx)$y / approx(knon$x,knon$y,subx)$y )
        
        f <- which(logratio_sample == min(logratio_sample,na.rm = TRUE))
        xmax <- bagelR.round_to_hundredth(subx[f])
        
        RANGE<-seq(xmin,xmax+0.01,0.01)
        RANGE<-round(RANGE,digits = 2)
        
        logratio_lookup <- log2(approx(kess$x,kess$y,RANGE)$y / approx(knon$x,knon$y,RANGE)$y)
        names(logratio_lookup)<-RANGE
        
        foldchanges <- FOLD_CHANGES[guide_test,SAMPLE_TO_TEST]
        
        foldchanges[foldchanges>xmax]<-xmax
        foldchanges[foldchanges<xmin]<-xmin
        foldchanges<-round(foldchanges,digits=2)
        
        currentLoop_bf<-rep(NA,length(guide_test))
        
        currentLoop_bf<-logratio_lookup[as.character(foldchanges)]
        
        names(currentLoop_bf)<-guide_test
        BFs_across_loops[names(currentLoop_bf),loop]<-currentLoop_bf
        BFs_across_loops_including_training[names(currentLoop_bf),loop]<-currentLoop_bf
        
        foldchanges <- FOLD_CHANGES[guide_training,SAMPLE_TO_TEST]
        
        foldchanges[foldchanges>xmax]<-xmax
        foldchanges[foldchanges<xmin]<-xmin
        foldchanges<-round(foldchanges,digits=2)
        
        currentLoop_bf<-rep(NA,length(guide_test))
        
        currentLoop_bf<-logratio_lookup[as.character(foldchanges)]
        
        names(currentLoop_bf)<-guide_training
        BFs_across_loops_including_training[names(currentLoop_bf),loop]<-currentLoop_bf
        
    }
    
    print('')
    print('DONE')
    
    Sys.sleep(1)
    close(pb)
    
    if(compPerformances){
        PERFORMANCES<-
            bagelR.sgRNA_bootStrapped_performances(BS_BF = BFs_across_loops,
                                                   refGuidesLibrary = refGuidesLibrary,
                                                   essentialGuides = EssentialGuides,
                                                   nonessentialGuides = nonEssentialGuides)
    }else{
        PERFORMANCES<-NULL
    }
    
    BFs_across_loops[BFs_across_loops==-Inf]<-NA
    
    sgRNA_BFs<-apply(BFs_across_loops,MARGIN = 1,'mean',na.rm=TRUE)
    sgRNA_BFs_sd<-apply(BFs_across_loops,MARGIN = 1,'sd',na.rm=TRUE)
    
    sgRNA_BFs<-cbind(sgRNA_BFs,sgRNA_BFs_sd)
    colnames(sgRNA_BFs)<-c('avg_bootstr_BFs','sd_bootstr_BFs')
    
    
    sgRNA_BFs_inclTr<-apply(BFs_across_loops_including_training,MARGIN = 1,'mean',na.rm=TRUE)
    sgRNA_BFs_sd_inclTr<-apply(BFs_across_loops_including_training,MARGIN = 1,'sd',na.rm=TRUE)
    
    sgRNA_BFs_inclTr<-cbind(sgRNA_BFs_inclTr,sgRNA_BFs_sd_inclTr)
    colnames(sgRNA_BFs_inclTr)<-c('avg_bootstr_BFs','sd_bootstr_BFs')
    
    return(list(sgRNA_BFs=sgRNA_BFs,sgRNA_BFs_inclTr=sgRNA_BFs_inclTr,boostPERF=PERFORMANCES))
}

bagelR.sgRNA_bootStrapped_performances<-function(BS_BF,essentialGuides,
                                                nonessentialGuides,refGuidesLibrary){


    print('Computing Bootstrapped performances across iterations...')

    #essentialGuides<-
       # refGuidesLibrary[match(intersect(ESSENTIAL_GENES,refGuidesLibrary[,2]),refGuidesLibrary[,2]),1]

    #nonessentialGuides<-
        #refGuidesLibrary[match(intersect(NON_ESSENTIAL_GENES,refGuidesLibrary[,2]),refGuidesLibrary[,2]),1]

    tmp<-BS_BF[intersect(rownames(BS_BF),union(essentialGuides,nonessentialGuides)),]
    
    minimum<-min(tmp[tmp>-Inf])
    #changed 20.6.19 to ensure finite upper bound
    #maximum<-max(tmp[tmp>-Inf])
    maximum<-max(tmp[tmp<Inf])

    nboostLoops<-ncol(BS_BF)
   
    range<-seq(minimum,maximum,(maximum-minimum)/999)
    
    SENSITIVITY<-matrix(NA,1000,nboostLoops)
    SPECIFICITY<-matrix(NA,1000,nboostLoops)
    THRESHOLD<-matrix(NA,1000,nboostLoops)
    PPV<-matrix(NA,1000,nboostLoops)
    AUC<-rep(NA,nboostLoops)
    
    pb <- txtProgressBar(min=1,max=nboostLoops,style=3)

    for (i in 1:nboostLoops){
        setTxtProgressBar(pb, i)
        currentTestSet<-names(which(BS_BF[,i]!='-Inf'))
        
        testEssential<-intersect(essentialGuides,currentTestSet)
        testNonessential<-intersect(nonessentialGuides,currentTestSet)

        currentTestSet<-union(testEssential,testNonessential)
        
        predictions<-BS_BF[currentTestSet,i]
        
        essentiality<-is.element(currentTestSet,testEssential)+0

        ROC<-roc(essentiality,predictions)
        
        AUC[i]<-ROC$auc

        COORDS<-coords(ROC,x = range,ret = c('sensitivity','specificity','ppv'),transpose=TRUE)
        
        SENSITIVITY[,i]<-COORDS['sensitivity',]
        SPECIFICITY[,i]<-COORDS['specificity',]
        PPV[,i]<-COORDS['ppv',]
    }

    THRESHOLD<-range

    avgSens<-apply(SENSITIVITY,MARGIN = 1,'mean',na.rm=TRUE)
    avgSpec<-apply(SPECIFICITY,MARGIN = 1,'mean',na.rm=TRUE)
    avgPpv<-apply(PPV,MARGIN=1,'mean',na.rm=TRUE)

    sdSens<-apply(SENSITIVITY,MARGIN = 1,'sd',na.rm=TRUE)
    sdSpec<-apply(SPECIFICITY,MARGIN = 1,'sd',na.rm=TRUE)
    sdPpv<-apply(PPV,MARGIN=1,'sd',na.rm=TRUE)

    print('')
    print('DONE')
    Sys.sleep(1)
    close(pb)


    return(list(th=THRESHOLD,
                ppv=avgPpv,sens=avgSens,spec=avgSpec,
                sd_ppv=sdPpv,sd_sens=sdSens,sd_spec=sdSpec,AUROC=mean(AUC)))

}

bagelR.diagnosticPlots<-function(sgRNA_BFs,cellLine,outDir){
    
        fn<-paste(outDir,cellLine,'.pdf',sep='')
        pdf(fn)
        
        plot(sgRNA_BFs$boostPERF$th,frame.plot = FALSE,
             sgRNA_BFs$boostPERF$ppv,type='l',lwd=2,
             xlab='sgRNA Bayesian Factor',col='red',
             main=paste(cellLine,': avg precision/sensitivity across boostrap iterations',sep=''),
             ylab='')
        
        par(new=TRUE)
        par(mar=c(4,4,4,5))
        plot(sgRNA_BFs$boostPERF$th,frame.plot = FALSE,
             sgRNA_BFs$boostPERF$sens,type='l',lwd=2,
             xaxt='n',yaxt='n',ylab='',xlab='',col='blue')
        axis(4)

        plot(1-sgRNA_BFs$boostPERF$spec,
             sgRNA_BFs$boostPERF$sens,
             type='l',lwd=2,
             xlab='FPR',col='blue',
             main=paste(cellLine,': avg FPR/TPR across boostrap iterations (avg AUC ',
                        format(sgRNA_BFs$boostPERF$AUROC,digits=4),')',sep=''),
             ylab='TPR')
        lines(x = c(0,1),y=c(0,1))
        dev.off()

        }


bagelR.geneLevel_BFs<-function(allguidesBFs,refGuidesLibrary){
    nreplicates<-length(allguidesBFs)
    allBFsmatrix<-matrix(NA,nrow(allguidesBFs$sgRNA_BFs),nreplicates,
                         dimnames = list(rownames(allguidesBFs$sgRNA_BFs),names(allguidesBFs)))

    allBFsmatrix[,1]<-allguidesBFs$sgRNA_BFs_inclTr[,1]

    red_refGuidesLibrary<-
        refGuidesLibrary[which(is.element(rownames(refGuidesLibrary),rownames(allBFsmatrix))),]
    print(paste('dim of red_refGuides: ',dim(red_refGuidesLibrary)))
    genes<-unique(refGuidesLibrary[,"GENES"])
    ngenes<-length(genes)

    GENElevel_BF<-rep(NA,length(genes))
    GENElevel_BF_sd<-rep(NA,length(genes))
    GENElevel_BF<-cbind(GENElevel_BF,GENElevel_BF_sd)
    rownames(GENElevel_BF)<-genes
    colnames(GENElevel_BF)<-c('Avg','SD')

    print('Computing gene-level Bayes Factors...')
    pb <- txtProgressBar(min=1,max=ngenes,style=3)

    for (i in 1:ngenes){
        setTxtProgressBar(pb, i)

        #guidesetids<-intersect(rownames(allBFsmatrix),as.character(KY_Library_v1.0_list[[genes[i]]]))
      ###CHANGED 24.5.19 to work with Avana and more general libraries without hard coding :
        guidesetids<-intersect(rownames(allBFsmatrix),
                               rownames(red_refGuidesLibrary[which(red_refGuidesLibrary[,"GENES"]==genes[i]),]))
        if(i%in%c(1,67,900)){
          print(guidesetids)
        }
        
        GENElevel_BF[i,1]<-mean(c(allBFsmatrix[guidesetids,1]))
        GENElevel_BF[i,2]<-sd(c(allBFsmatrix[guidesetids,1]))
    }
    
    print('')
    print('DONE')

    Sys.sleep(1)
    close(pb)

    return(GENElevel_BF)
}

bagelR.gene_bootStrapped_BFactors<-function(FOLD_CHANGES,
                                             ESSENTIAL_GENES,
                                             NON_ESSENTIAL_GENES,
                                             NUM_BOOTSTRAPS=1000,
                                             SAMPLE_TO_TEST,
                                             FC_DENS_THRESHOLD=2^-7,
                                             percTrainingSample=80,
                                             refGuidesLibrary,
                                             compPerformances=FALSE){
    
    #rownames(FOLD_CHANGES)<-FOLD_CHANGES$sgRNA
    
    ngenes<-nrow(FOLD_CHANGES)
    
    genes<-unique(FOLD_CHANGES$gene)
    rownames(FOLD_CHANGES)<-FOLD_CHANGES$gene
    #genes<-unique(FOLD_CHANGES$gene)
    #guides<-FOLD_CHANGES$sgRNA
    
    ESSENTIAL_GENES<-intersect(ESSENTIAL_GENES,genes)
    NON_ESSENTIAL_GENES<-intersect(NON_ESSENTIAL_GENES,genes)
    
    ESSENTIAL_GENES<-ESSENTIAL_GENES[!is.na(ESSENTIAL_GENES)]
    NON_ESSENTIAL_GENES<-NON_ESSENTIAL_GENES[!is.na(NON_ESSENTIAL_GENES)]
    print(paste('n. of unique genes =',ngenes))
    print(paste('n. of ess genes = ',length(ESSENTIAL_GENES)))
    print(paste('n of non ess genes = ',length(NON_ESSENTIAL_GENES)))

    
    BFs_across_loops<-matrix(-Inf,ngenes,NUM_BOOTSTRAPS,dimnames = list(genes,1:NUM_BOOTSTRAPS))
    BFs_across_loops_including_training<-matrix(-Inf,ngenes,NUM_BOOTSTRAPS,dimnames = list(genes,1:NUM_BOOTSTRAPS))
    
	nEssGenes<-length(ESSENTIAL_GENES)
	nNonEssGenes<-length(NON_ESSENTIAL_GENES)
	
    print('Bootstrap iterations in progress...')
    

    
   for(loop in 1:NUM_BOOTSTRAPS){
      
        gene_train_ess<-ESSENTIAL_GENES[sample(nEssGenes)[1:ceiling(nEssGenes*percTrainingSample/100)]]
        gene_train_non_ess<-NON_ESSENTIAL_GENES[sample(nNonEssGenes)[1:ceiling(nNonEssGenes*percTrainingSample/100)]]
        
        gene_training<- union(gene_train_ess,gene_train_non_ess)
        gene_test <- setdiff(genes,gene_training)
        
        ess_train_fc <- as.numeric(FOLD_CHANGES[gene_train_ess,SAMPLE_TO_TEST])
        non_ess_train_fc <- as.numeric(FOLD_CHANGES[gene_train_non_ess,SAMPLE_TO_TEST])
        
        kess<-density(ess_train_fc, kernel = "gaussian")
        knon<-density(non_ess_train_fc, kernel = "gaussian")
     
        x <- seq(-10,2,0.01)
        nonfitx <- approx(knon$x,knon$y,x)$y
        
        f <- which(nonfitx > FC_DENS_THRESHOLD)
        xmin <- bagelR.round_to_hundredth( min(x[f]) )
        
        subx <- seq(xmin,max(x[f]),0.01)
  
        logratio_sample <- log2( approx(kess$x,kess$y,subx)$y / approx(knon$x,knon$y,subx)$y )
        
        f <- which(logratio_sample == min(logratio_sample,na.rm = TRUE))
        xmax <- bagelR.round_to_hundredth(subx[f])
        
        RANGE<-seq(xmin,xmax+0.01,0.01)
        RANGE<-round(RANGE,digits = 2)
        
        logratio_lookup <- log2(approx(kess$x,kess$y,RANGE)$y / approx(knon$x,knon$y,RANGE)$y)
        names(logratio_lookup)<-RANGE
        
        foldchanges <- FOLD_CHANGES[gene_test,SAMPLE_TO_TEST]
        
        foldchanges[foldchanges>xmax]<-xmax
        foldchanges[foldchanges<xmin]<-xmin
        foldchanges<-round(foldchanges,digits=2)
        
        currentLoop_bf<-rep(NA,length(gene_test))
        
        currentLoop_bf<-logratio_lookup[as.character(foldchanges)]
        
        names(currentLoop_bf)<-gene_test
        BFs_across_loops[names(currentLoop_bf),loop]<-currentLoop_bf
        BFs_across_loops_including_training[names(currentLoop_bf),loop]<-currentLoop_bf
        
        foldchanges <- FOLD_CHANGES[gene_training,SAMPLE_TO_TEST]
        
        foldchanges[foldchanges>xmax]<-xmax
        foldchanges[foldchanges<xmin]<-xmin
        foldchanges<-round(foldchanges,digits=2)
        
        currentLoop_bf<-rep(NA,length(gene_training))
        
        currentLoop_bf<-logratio_lookup[as.character(foldchanges)]
        
        names(currentLoop_bf)<-gene_training
       
        BFs_across_loops_including_training[names(currentLoop_bf),loop]<-currentLoop_bf
        
        
    }
    

    
    if(compPerformances){
        PERFORMANCES<-
            bagelR.sgRNA_bootStrapped_performances(BS_BF = BFs_across_loops,
                                                   refGuidesLibrary = NULL,
                                                   ESSENTIAL_GENES = ESSENTIAL_GENES,
                                                   NON_ESSENTIAL_GENES = NON_ESSENTIAL_GENES)
    }else{
        PERFORMANCES<-NULL
    }
    
    BFs_across_loops[BFs_across_loops==-Inf]<-NA
    
    gene_BFs<-apply(BFs_across_loops,MARGIN = 1,'mean',na.rm=TRUE)
    gene_BFs_sd<-apply(BFs_across_loops,MARGIN = 1,'sd',na.rm=TRUE)
    
    gene_BFs<-cbind(gene_BFs,gene_BFs_sd)
    colnames(gene_BFs)<-c('avg_bootstr_BFs','sd_bootstr_BFs')
    
    
    gene_BFs_inclTr<-apply(BFs_across_loops_including_training,MARGIN = 1,'mean',na.rm=TRUE)
    gene_BFs_sd_inclTr<-apply(BFs_across_loops_including_training,MARGIN = 1,'sd',na.rm=TRUE)
    
    gene_BFs_inclTr<-cbind(gene_BFs_inclTr,gene_BFs_sd_inclTr)
    colnames(gene_BFs_inclTr)<-c('avg_bootstr_BFs','sd_bootstr_BFs')
    
    return(list(gene_BFs=gene_BFs,gene_BFs_inclTr=gene_BFs_inclTr,boostPERF=PERFORMANCES))
}


bagelR.computeAllGeneBFs<-function(cellLine,
                                        NUM_BOOTSTRAPS = 1000,
                                        ESSENTIAL_GENES,
                                        NON_ESSENTIAL_GENES,
                                        inputFolder='../../DATA/R/normalisedCounts_and_FCs/',
                                        outputFolder='../../RESULTS/BAGEL-R_output/',
                                        refGuidesLibrary,diagnosticPlots=TRUE,
                                        whatToTest='newFC',fn=paste0(inputFolder,'/CorrectedData.Rdata'),glb=TRUE){
  
  #fn<-paste(inputFolder,cellLine,'_CCRoutput.RData',sep='')
  

  load(fn)
  
  print('Computing gene boostrapped Bayesian Factors and classification performances for corrected FCs...')
  
  
  correctedFCs<-data.frame(gene=rownames(CorrectedData),correctedFC=CorrectedData[,cellLine])
  #save(correctedFCs,file=paste0(outputFolder,"/correctedFCs.Rdata"))
  #correctedFCs<-cbind.data.frame(rownames(CorrectedData),as.numeric(CorrectedData[,cellLine]),stringsAsFactors=FALSE)
  #colnames(correctedFCs)<-c("gene",'correctedFC')

  set.seed(0xA5EED)

  #rep1<-bagelR.gene_bootStrapped_BFactors(FOLD_CHANGES = correctedFCs,
  #rep1<-bagelRV2.gene_bootStrapped_BFactors(FOLD_CHANGES = correctedFCs,
  #                                         NUM_BOOTSTRAPS = NUM_BOOTSTRAPS,
  #                                         ESSENTIAL_GENES = ESSENTIAL_GENES,
  #                                         NON_ESSENTIAL_GENES = NON_ESSENTIAL_GENES,
  #                                         SAMPLE_TO_TEST = whatToTest,
  #                                         refGuidesLibrary = NULL,
  #                                         compPerformances = FALSE,glb=glb)
  
  GeneLevelBFsV2<-bagelRV2.gene_bootStrapped_BFactors(FOLD_CHANGES=correctedFCs,
                                                      ESSENTIAL_GENES=curated_BAGEL_essential,
                                                      NON_ESSENTIAL_GENES=curated_BAGEL_nonEssential,
                                                      NUM_BOOTSTRAPS=1000,
                                                      SAMPLE_TO_TEST="correctedFC",
                                                      FC_DENS_THRESHOLD=2^-7,
                                                      percTrainingSample=80,
                                                      refGuidesLibrary=Guide_library,
                                                      compPerformances=FALSE,testPlot=FALSE,glb=F)
  
  sgRNA_BFs<-GeneLevelBFsV2
  print('DONE')
  
  save(sgRNA_BFs,file=paste(outputFolder,cellLine,'_gene_BFs.Rdata',sep=''))
  if (diagnosticPlots){
    bagelR.diagnosticPlots(sgRNA_BFs,cellLine = cellLine,outDir = outputFolder)
  }
  
  return(sgRNA_BFs)
}


bagelR.computeAllGeneBFsTCGA<-function(cellLine,
                                   NUM_BOOTSTRAPS = 1000,
                                   ESSENTIAL_GENES,
                                   NON_ESSENTIAL_GENES,
                                   inputFolder='../../DATA/R/normalisedCounts_and_FCs/',
                                   outputFolder='../../RESULTS/BAGEL-R_output/',
                                   refGuidesLibrary,diagnosticPlots=TRUE,
                                   whatToTest='newFC',startrange=seq(-10,2,0.01)){
  
  #fn<-paste(inputFolder,cellLine,'_CCRoutput.RData',sep='')
  #load in the batch corrected FC data This should be called CorrectedData

  load(paste0(inputFolder,"/CellignerData.Rdata"))
#called tCI
  
  print('Computing gene boostrapped Bayesian Factors and classification performances for corrected FCs...')
  
  pcgene<-read.delim("~/CRISPR_Pipelines/ExternalData/protein-coding_gene.txt",header=T,stringsAsFactors = F,sep="\t")
  temprnames<-cbind.data.frame(rownames(tCI),pcgene[match(rownames(tCI),pcgene$ensembl_gene_id),"symbol"])
  temprnames<-temprnames[!duplicated(temprnames[,2]),]
  temprnames<-temprnames[!is.na(temprnames[,2]),]
  tCI<-tCI[temprnames[,1],]
  rownames(tCI)<-temprnames[,2]
  correctedFCs<-cbind.data.frame(rownames(tCI),as.numeric(tCI[,cellLine]),stringsAsFactors=FALSE)
  colnames(correctedFCs)<-c("gene",'correctedFC')
  
  set.seed(0xA5EED)
  
  rep1<-bagelRV2.gene_bootStrapped_BFactors(FOLD_CHANGES = correctedFCs,
                                          NUM_BOOTSTRAPS = NUM_BOOTSTRAPS,
                                          ESSENTIAL_GENES = ESSENTIAL_GENES,
                                          NON_ESSENTIAL_GENES = NON_ESSENTIAL_GENES,
                                          SAMPLE_TO_TEST = whatToTest,
                                          refGuidesLibrary = NULL,
                                          compPerformances = FALSE,startrange=startrange)
  sgRNA_BFs<-rep1
  print('DONE')
  
  save(sgRNA_BFs,file=paste(outputFolder,cellLine,'_gene_BFs.Rdata',sep=''))
  if (diagnosticPlots){
    bagelR.diagnosticPlots(sgRNA_BFs,cellLine = cellLine,outDir = outputFolder)
  }
  return(sgRNA_BFs)
}


bagelR.computeAllGeneBFsTCGAvoom<-function(cellLine,
                                         NUM_BOOTSTRAPS = 1000,
                                         ESSENTIAL_GENES,
                                         NON_ESSENTIAL_GENES,
                                         inputFolder='../../DATA/R/normalisedCounts_and_FCs/',
                                         outputFolder='../../RESULTS/BAGEL-R_output/',
                                         refGuidesLibrary,diagnosticPlots=TRUE,
                                         whatToTest='newFC',subtype=NULL){
    
    #fn<-paste(inputFolder,cellLine,'_CCRoutput.RData',sep='')
    #load in the batch corrected FC data This should be called CorrectedData
    load(paste0(inputFolder,"/",subtype,"_ExpPC.Rdata"))
    TCGAData<-Exp_PC
    #called tCI
    colnames(TCGAData)<-make.names(colnames(TCGAData))
    print('Computing gene boostrapped Bayesian Factors and classification performances for corrected FCs...')
    

    correctedFCs<-cbind.data.frame(rownames(TCGAData),as.numeric(TCGAData[,cellLine]),stringsAsFactors=FALSE)
    colnames(correctedFCs)<-c("gene",'correctedFC')
    save(correctedFCs,file=paste0(outputFolder,"/correctedFCs",cellLine,".Rdata"))
    set.seed(0xA5EED)
    
    rep1<-bagelR.gene_bootStrapped_BFactorsVoom(FOLD_CHANGES = correctedFCs,
                                            NUM_BOOTSTRAPS = NUM_BOOTSTRAPS,
                                            ESSENTIAL_GENES = ESSENTIAL_GENES,
                                            NON_ESSENTIAL_GENES = NON_ESSENTIAL_GENES,
                                            SAMPLE_TO_TEST = whatToTest,
                                            refGuidesLibrary = NULL,
                                            compPerformances = FALSE)
    sgRNA_BFs<-rep1
    print('DONE')
    
    save(sgRNA_BFs,file=paste(outputFolder,cellLine,'_gene_BFs.Rdata',sep=''))
    if (diagnosticPlots){
      bagelR.diagnosticPlots(sgRNA_BFs,cellLine = cellLine,outDir = outputFolder)
    }
  
  return(sgRNA_BFs)
}
  
bagelR.gene_bootStrapped_BFactorsVoom<-function(FOLD_CHANGES,
                                              ESSENTIAL_GENES,
                                              NON_ESSENTIAL_GENES,
                                              NUM_BOOTSTRAPS=1000,
                                              SAMPLE_TO_TEST,
                                              FC_DENS_THRESHOLD=2^-7,
                                              percTrainingSample=80,
                                              refGuidesLibrary,
                                              compPerformances=FALSE){
    
    #rownames(FOLD_CHANGES)<-FOLD_CHANGES$sgRNA
    
    ngenes<-nrow(FOLD_CHANGES)
    
    genes<-unique(FOLD_CHANGES$gene)
    rownames(FOLD_CHANGES)<-FOLD_CHANGES$gene
    #genes<-unique(FOLD_CHANGES$gene)
    #guides<-FOLD_CHANGES$sgRNA
    
    ESSENTIAL_GENES<-intersect(ESSENTIAL_GENES,genes)
    NON_ESSENTIAL_GENES<-intersect(NON_ESSENTIAL_GENES,genes)
    
    ESSENTIAL_GENES<-ESSENTIAL_GENES[!is.na(ESSENTIAL_GENES)]
    NON_ESSENTIAL_GENES<-NON_ESSENTIAL_GENES[!is.na(NON_ESSENTIAL_GENES)]
    print(paste('n. of unique genes =',ngenes))
    print(paste('n. of ess genes = ',length(ESSENTIAL_GENES)))
    print(paste('n of non ess genes = ',length(NON_ESSENTIAL_GENES)))
    
    
    BFs_across_loops<-matrix(-Inf,ngenes,NUM_BOOTSTRAPS,dimnames = list(genes,1:NUM_BOOTSTRAPS))
    BFs_across_loops_including_training<-matrix(-Inf,ngenes,NUM_BOOTSTRAPS,dimnames = list(genes,1:NUM_BOOTSTRAPS))
    
    nEssGenes<-length(ESSENTIAL_GENES)
    nNonEssGenes<-length(NON_ESSENTIAL_GENES)
    
    print('Bootstrap iterations in progress...')
    
    
    
    for(loop in 1:NUM_BOOTSTRAPS){
      
      gene_train_ess<-ESSENTIAL_GENES[sample(nEssGenes)[1:ceiling(nEssGenes*percTrainingSample/100)]]
      gene_train_non_ess<-NON_ESSENTIAL_GENES[sample(nNonEssGenes)[1:ceiling(nNonEssGenes*percTrainingSample/100)]]
      
      gene_training<- union(gene_train_ess,gene_train_non_ess)
      gene_test <- setdiff(genes,gene_training)
      
      ess_train_fc <- as.numeric(FOLD_CHANGES[gene_train_ess,SAMPLE_TO_TEST])
      non_ess_train_fc <- as.numeric(FOLD_CHANGES[gene_train_non_ess,SAMPLE_TO_TEST])
      
      kess<-density(ess_train_fc, kernel = "gaussian")
      knon<-density(non_ess_train_fc, kernel = "gaussian")
 
      x <- seq(min(FOLD_CHANGES[,SAMPLE_TO_TEST]),max(FOLD_CHANGES[,SAMPLE_TO_TEST]),0.01)
      nonfitx <- approx(knon$x,knon$y,x)$y
      
      f <- which(nonfitx > FC_DENS_THRESHOLD)
      xmin <- bagelR.round_to_hundredth( min(x[f]) )
      #this part isn't really necessary:
      essfitx<-approx(kess$x,kess$y,x)$y
      f<-which(essfitx> FC_DENS_THRESHOLD)
      xmax<-bagelR.round_to_hundredth( max(x[f]) )
      
      subx <- seq(xmin,xmax,0.05)
      #print(summary(subx))
      #print(approx(kess$x,kess$y,subx)$y)
      #print(approx(knon$x,knon$y,subx)$y)
      #logratio_sample <- log2( approx(kess$x,kess$y,subx)$y / approx(knon$x,knon$y,subx)$y )
 
      #f <- which(logratio_sample == min(logratio_sample,na.rm = TRUE))
      #xmax <- bagelR.round_to_hundredth(subx[f[1]])
      
      RANGE<-seq(xmin,xmax+0.1,0.01)
      RANGE<-round(RANGE,digits = 2)
    
      #instead of logratio_lookup change to linear approximation as in BAGELv2 26.9.20- See BagelRv2 code
      logratio_lookup <- log2(approx(kess$x,kess$y,RANGE)$y / approx(knon$x,knon$y,RANGE)$y)
      names(logratio_lookup)<-RANGE
      logratio_lookup[is.na(logratio_lookup)&RANGE<min(ess_train_fc)]<- min(logratio_lookup,na.rm=TRUE)-1
      logratio_lookup[is.na(logratio_lookup)&RANGE>max(non_ess_train_fc)]<- max(logratio_lookup,na.rm=TRUE)+1
   
      foldchanges <- FOLD_CHANGES[gene_test,SAMPLE_TO_TEST]
      
      foldchanges[foldchanges>xmax]<-xmax
      foldchanges[foldchanges<xmin]<-xmin
      foldchanges<-round(foldchanges,digits=2)
      
      currentLoop_bf<-rep(NA,length(gene_test))
      
      currentLoop_bf<-logratio_lookup[as.character(foldchanges)]
      
      names(currentLoop_bf)<-gene_test
      BFs_across_loops[names(currentLoop_bf),loop]<-currentLoop_bf
      BFs_across_loops_including_training[names(currentLoop_bf),loop]<-currentLoop_bf
      
      foldchanges <- FOLD_CHANGES[gene_training,SAMPLE_TO_TEST]
      
      foldchanges[foldchanges>xmax]<-xmax
      foldchanges[foldchanges<xmin]<-xmin
      foldchanges<-round(foldchanges,digits=2)
      
      currentLoop_bf<-rep(NA,length(gene_training))
      
      currentLoop_bf<-logratio_lookup[as.character(foldchanges)]
      
      names(currentLoop_bf)<-gene_training
      
      BFs_across_loops_including_training[names(currentLoop_bf),loop]<-currentLoop_bf
      
      
    }
    
    
    
    if(compPerformances){
      PERFORMANCES<-
        bagelR.sgRNA_bootStrapped_performances(BS_BF = BFs_across_loops,
                                               refGuidesLibrary = NULL,
                                               ESSENTIAL_GENES = ESSENTIAL_GENES,
                                               NON_ESSENTIAL_GENES = NON_ESSENTIAL_GENES)
    }else{
      PERFORMANCES<-NULL
    }
    
    BFs_across_loops[BFs_across_loops==-Inf]<- NA
    BFs_across_loops[BFs_across_loops==Inf]<- NA
    gene_BFs<-apply(BFs_across_loops,MARGIN = 1,'mean',na.rm=TRUE)
    gene_BFs_sd<-apply(BFs_across_loops,MARGIN = 1,'sd',na.rm=TRUE)
    
    gene_BFs<-cbind(gene_BFs,gene_BFs_sd)
    colnames(gene_BFs)<-c('avg_bootstr_BFs','sd_bootstr_BFs')
    
    
    gene_BFs_inclTr<-apply(BFs_across_loops_including_training,MARGIN = 1,'mean',na.rm=TRUE)
    gene_BFs_sd_inclTr<-apply(BFs_across_loops_including_training,MARGIN = 1,'sd',na.rm=TRUE)
    
    gene_BFs_inclTr<-cbind(gene_BFs_inclTr,gene_BFs_sd_inclTr)
    colnames(gene_BFs_inclTr)<-c('avg_bootstr_BFs','sd_bootstr_BFs')
    
    return(list(gene_BFs=gene_BFs,gene_BFs_inclTr=gene_BFs_inclTr,boostPERF=PERFORMANCES))
  }
  

bagelRV2.gene_bootStrapped_BFactors<-function(FOLD_CHANGES,
                                                ESSENTIAL_GENES,
                                                NON_ESSENTIAL_GENES,
                                                NUM_BOOTSTRAPS=1000,
                                                SAMPLE_TO_TEST,
                                                FC_DENS_THRESHOLD=2^-7,
                                                percTrainingSample=80,
                                                refGuidesLibrary,
                                                compPerformances=FALSE,testPlot=FALSE,startrange=seq(-10,2,0.01),glb=TRUE){
  
  #rownames(FOLD_CHANGES)<-FOLD_CHANGES$sgRNA
  
  ngenes<-nrow(FOLD_CHANGES)
  
  genes<-unique(FOLD_CHANGES$gene)
  rownames(FOLD_CHANGES)<-FOLD_CHANGES$gene
  #genes<-unique(FOLD_CHANGES$gene)
  #guides<-FOLD_CHANGES$sgRNA
  
  ESSENTIAL_GENES<-intersect(ESSENTIAL_GENES,genes)
  NON_ESSENTIAL_GENES<-intersect(NON_ESSENTIAL_GENES,genes)
  
  ESSENTIAL_GENES<-ESSENTIAL_GENES[!is.na(ESSENTIAL_GENES)]
  NON_ESSENTIAL_GENES<-NON_ESSENTIAL_GENES[!is.na(NON_ESSENTIAL_GENES)]
  print(paste('n. of unique genes =',ngenes))
  print(paste('n. of ess genes = ',length(ESSENTIAL_GENES)))
  print(paste('n of non ess genes = ',length(NON_ESSENTIAL_GENES)))
  
  
  BFs_across_loops<-matrix(-Inf,ngenes,NUM_BOOTSTRAPS,dimnames = list(genes,1:NUM_BOOTSTRAPS))
  BFs_across_loops_including_training<-matrix(-Inf,ngenes,NUM_BOOTSTRAPS,dimnames = list(genes,1:NUM_BOOTSTRAPS))
  
  nEssGenes<-length(ESSENTIAL_GENES)
  nNonEssGenes<-length(NON_ESSENTIAL_GENES)
  
  print('Bootstrap iterations in progress...')
  exponent<-(-1.1535 * log(length(NON_ESSENTIAL_GENES) + 13.324) + 0.7728)
  FC_DENS_THRESHOLD<- 2*(10^(exponent))
  
  for(loop in 1:NUM_BOOTSTRAPS){
    
    gene_train_ess<-ESSENTIAL_GENES[sample(nEssGenes)[1:ceiling(nEssGenes*percTrainingSample/100)]]
    gene_train_non_ess<-NON_ESSENTIAL_GENES[sample(nNonEssGenes)[1:ceiling(nNonEssGenes*percTrainingSample/100)]]
    
    gene_training<- union(gene_train_ess,gene_train_non_ess)
    gene_test <- setdiff(genes,gene_training)
    
    ess_train_fc <- as.numeric(FOLD_CHANGES[gene_train_ess,SAMPLE_TO_TEST])
    non_ess_train_fc <- as.numeric(FOLD_CHANGES[gene_train_non_ess,SAMPLE_TO_TEST])
    
    kess<-density(ess_train_fc, kernel = "gaussian")
    knon<-density(non_ess_train_fc, kernel = "gaussian")
    #change for the gene expression data as logTPM expression not fitness fold change depletions
    #x <- seq(min(FOLD_CHANGES[,SAMPLE_TO_TEST])-1,max(FOLD_CHANGES[,SAMPLE_TO_TEST])+1,0.01)
    x<-startrange
    nonfitx <- approx(knon$x,knon$y,x)$y
    
    f <- which(nonfitx > FC_DENS_THRESHOLD)
    if(glb){
      low<-which(nonfitx < FC_DENS_THRESHOLD+1e-8)
      values<-x[low]
      xmin<-bagelR.round_to_hundredth(max(values[values<0]))
    }else{
      xmin <- bagelR.round_to_hundredth( min(x[f]) )
    }
    

    xmax<-bagelR.round_to_hundredth(max(x[f]))

    subx <- seq(xmin,xmax,0.01)
    
    logratio_sample <- log2( approx(kess$x,kess$y,subx)$y / approx(knon$x,knon$y,subx)$y )
    
    f <- which(logratio_sample == min(logratio_sample,na.rm = TRUE))
    xmax <- bagelR.round_to_hundredth(subx[f])
  
    if(testPlot){
      print({plot(kess)
      lines(knon,col="blue")
      abline(v=xmin)
      abline(v=xmax)})
      
    }
    RANGE<-seq(xmin,xmax+0.01,0.01)
    RANGE<-round(RANGE,digits = 2)
    
    logratio_lookup <- log2(approx(kess$x,kess$y,RANGE)$y / approx(knon$x,knon$y,RANGE)$y)
 
    fitdf<-data.frame(logRatio=logratio_lookup,fc=RANGE)
    if(testPlot){
      print(plot(RANGE,logratio_lookup))
    }
    #then fit a linear model to the log ratio values instead:
    fitdf<-na.omit(fitdf)
    fitdf<-fitdf[fitdf$logRatio!="Inf",]
    
    logratio_lm<-lm(logRatio~fc,data=fitdf)

    foldchanges <- data.frame(fc=FOLD_CHANGES[gene_test,SAMPLE_TO_TEST])

    currentLoop_bf<-rep(NA,length(gene_test))
    
    currentLoop_bf<-predict(logratio_lm,newdata=foldchanges)
    
    names(currentLoop_bf)<-gene_test
    BFs_across_loops[names(currentLoop_bf),loop]<-currentLoop_bf
    BFs_across_loops_including_training[names(currentLoop_bf),loop]<-currentLoop_bf
    
    foldchanges <- data.frame(fc=FOLD_CHANGES[gene_training,SAMPLE_TO_TEST])
    
  
    currentLoop_bf<-rep(NA,length(gene_training))
    
    currentLoop_bf<-predict(logratio_lm,newdata=foldchanges)
    
    names(currentLoop_bf)<-gene_training
    
    BFs_across_loops_including_training[names(currentLoop_bf),loop]<-currentLoop_bf
    
    
  }
  
  
  
  if(compPerformances){
    PERFORMANCES<-
      bagelR.sgRNA_bootStrapped_performances(BS_BF = BFs_across_loops,
                                             refGuidesLibrary = NULL,
                                             ESSENTIAL_GENES = ESSENTIAL_GENES,
                                             NON_ESSENTIAL_GENES = NON_ESSENTIAL_GENES)
  }else{
    PERFORMANCES<-NULL
  }
  
  BFs_across_loops[BFs_across_loops==-Inf]<- NA
  BFs_across_loops[BFs_across_loops==Inf]<- NA
  gene_BFs<-apply(BFs_across_loops,MARGIN = 1,'mean',na.rm=TRUE)
  gene_BFs_sd<-apply(BFs_across_loops,MARGIN = 1,'sd',na.rm=TRUE)
  
  gene_BFs<-cbind(gene_BFs,gene_BFs_sd)
  colnames(gene_BFs)<-c('avg_bootstr_BFs','sd_bootstr_BFs')
  
  
  gene_BFs_inclTr<-apply(BFs_across_loops_including_training,MARGIN = 1,'mean',na.rm=TRUE)
  gene_BFs_sd_inclTr<-apply(BFs_across_loops_including_training,MARGIN = 1,'sd',na.rm=TRUE)
  
  gene_BFs_inclTr<-cbind(gene_BFs_inclTr,gene_BFs_sd_inclTr)
  colnames(gene_BFs_inclTr)<-c('avg_bootstr_BFs','sd_bootstr_BFs')
  
  return(list(gene_BFs=gene_BFs,gene_BFs_inclTr=gene_BFs_inclTr,boostPERF=PERFORMANCES))
}


computeAllGeneGexpTCGA<-function(cellLine = args[1],
                       inputFolder = args[3],fdr=0.1){
  
load(paste0(inputFolder,"/CellignerData.Rdata"))
#called tCI

print('Computing gene boostrapped Bayesian Factors and classification performances for corrected FCs...')

pcgene<-read.delim("~/CRISPR_Pipelines/ExternalData/protein-coding_gene.txt",header=T,stringsAsFactors = F,sep="\t")
temprnames<-cbind.data.frame(rownames(tCI),pcgene[match(rownames(tCI),pcgene$ensembl_gene_id),"symbol"])
temprnames<-temprnames[!duplicated(temprnames[,2]),]
temprnames<-temprnames[!is.na(temprnames[,2]),]
tCI<-tCI[temprnames[,1],]
rownames(tCI)<-temprnames[,2]
correctedFCs<-cbind.data.frame(rownames(tCI),as.numeric(tCI[,cellLine]),stringsAsFactors=FALSE)
colnames(correctedFCs)<-c("gene",'correctedFC')
rownames(correctedFCs)<-correctedFCs$gene
gaussmean<-median(correctedFCs[,2])
gausssd<-mad(correctedFCs[,2])

p_FCs<-2*pnorm(-abs(correctedFCs[,2]),mean=gaussmean,sd=gausssd)
p_FCsadjust<-p.adjust(p_FCs)
LowExpr<-p_FCsadjust<fdr&p_FCs<0+0
HighExpr<-p_FCsadjust<fdr&p_FCs>0+0
names(LowExpr)<-rownames(correctedFCs)
names(HighExpr)<-rownames(correctedFCs)
return(list(LowExpr=LowExpr,HighExpr=HighExpr))


}
