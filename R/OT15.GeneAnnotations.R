buildPathMemb<-function(GENES,pathway_ids,PATH_COLLECTION){
    np<-length(pathway_ids)
    ngenes<-length(GENES)
    MM<-matrix(0,ngenes,np,dimnames = list(GENES,pathway_ids))

    for (i in 1:np){
        MM[,i]<-is.element(GENES,PATH_COLLECTION$HGNC_SYMBOL[[pathway_ids[i]]])+0
    }

    return(MM)
}


CoreComponent<-function (PFP, EM, PATH = "./", PATH_COLLECTION,BGgs)
{
    filename <- PATH
    #EM <- sign(EM)

    NAMES <- PFP$pathway
    np <- length(NAMES)

    print("Producing heatmaps for core-components of enriched pathways...")

    for (i in 1:np) {
        currentGenes <- PATH_COLLECTION$HGNC_SYMBOL[[match(NAMES[[i]],PATH_COLLECTION$PATHWAY)]]
        currentGenes <- intersect(currentGenes, rownames(EM))
        if (i == 1) {
            mutG <- currentGenes
        }
        else {
            mutG <- union(mutG, currentGenes)
        }
    }

    mutG <- sort(mutG)
    #FREQS <- sort(rowSums(EM[mutG, ]==1), decreasing = TRUE)

    ids<-match(PFP$pathway,PATH_COLLECTION$PATHWAY)

    MM <- buildPathMemb(mutG, ids, PATH_COLLECTION)

    cm <- fastgreedy.community(graph.incidence(MM))

    heatmapsMatrices<-list()

    flag<-1
    for (i in 1:length(cm)){
        currentMembers<-cm$names[cm$membership==i]
        currentGenes<-intersect(currentMembers,rownames(EM))
        currentPathways<-setdiff(currentMembers,currentGenes)
        currentSubM<-t(MM[currentGenes,currentPathways])

        if(length(currentPathways)>1 & length(unique(c(currentSubM)))>1){
            pvals<-PFP$pvals[match(PATH_COLLECTION$PATHWAY[as.numeric(rownames(currentSubM))],PFP$pathway)]
            names<-PFP$pathway[match(PATH_COLLECTION$PATHWAY[as.numeric(rownames(currentSubM))],PFP$pathway)]
            CG<-PFP$GENES[match(PATH_COLLECTION$PATHWAY[as.numeric(rownames(currentSubM))],PFP$pathway)]
            CG<-str_split(CG,', ')
            CG<-unlist(lapply(CG,'length'))

            allG<-PATH_COLLECTION$HGNC_SYMBOL[as.numeric(rownames(currentSubM))]
            allG<-unlist(lapply(lapply(allG,intersect,BGgs),length))
            names(allG)<-NULL


            tokened<-str_split(names,'//')
            llen<-unlist(lapply(tokened,'length'))
            tokened<-unlist(lapply(tokened,function(x){x[1]}))
            suffix<-rep('',length(llen))
            suffix[llen>1]<-', ...'

            pvalsPattern<-rep('',length(pvals))
            pvalsPattern[pvals<0.01]<-' *'
            pvalsPattern[pvals<0.001]<-' **'
            pvalsPattern[pvals<0.0001]<-' ***'
            tokened<-paste(CG,'/',allG,' ',tokened,' ',suffix,' ',pvalsPattern,sep='')


            rownames(currentSubM)<-tokened
            currentSubM<-currentSubM[order(rowSums(currentSubM),decreasing=TRUE),order(colSums(currentSubM),decreasing=TRUE)]


            pheatmap(currentSubM,cluster_rows = FALSE,cluster_cols = FALSE,col=c('white','blue'),legend=FALSE)


            heatmapsMatrices[[flag]]<-currentSubM
            flag<-flag+1

            }
        }

    return(heatmapsMatrices)
}

profileCS<-function(genesToprofile,perctgs,nNAs=3){

    load('../../20180226_Paperfreeze_data/other annotations/TissueColors.RData')
    tmpnames<-colnames(perctgs)
    tmpnames<-str_split(tmpnames,'Vuln_CL_perc.')
    tmpnames<-lapply(tmpnames,function(x){x[2]})
    tmpnames<-str_replace_all(tmpnames,'[.]',' ')
    tmpnames[which(tmpnames=='Non Small Cell Lung Carcinoma')]<-'Non-Small Cell Lung Carcinoma'

    colnames(perctgs)<-tmpnames

    tmp<-as.numeric(t(as.matrix(perctgs[GenesToProfile,])))

    ngenes<-length(genesToprofile)
    toPlot<-NULL
    complexCol<-NULL
    for (i in 1:ngenes){
        toPlot<-c(toPlot,tmp[((i-1)*12+1):(i*12)],rep(NA,nNAs))
        complexCol<-c(complexCol,TissueColors[colnames(perctgs)])
    }

    par(xpd=TRUE)
    x<-barplot(toPlot,ylim=c(0,100),border=NA,col=complexCol)

    labsCoords<-x[seq(1,max(x),15)]+17
    text(labsCoords,rep(-10,length(genesToprofile)),genesToprofile,srt=45,pos = 2,
         cex=0.8)

    lines(x = range(x),y=c(50,50))

    }



kinases.fn<-'../../DATA/R/kinase_genes_and_syn.Rdata'
humKin.fn<-'../../DATA/R/humKinome.Rdata'
intoGen_GoF.fn<-'../../DATA/R/intoGen_GoF.Rdata'
GSK_epiTargets.fn<-'../../DATA/R/GSK_epiTargets.Rdata'

OT15_retrievegeneInfo<-function(GENES){

    fc<-read.delim(paste0(dir.ExternalData,'/protein-coding_gene.txt'),sep = '\t',header=TRUE,stringsAsFactors = FALSE)
    rownames(fc)<-fc$symbol

    commong<-intersect(rownames(fc),GENES)

    fc<-fc[commong,]

    description<-rep('',length(GENES))
    names(description)<-GENES

    hgnc_id<-rep('',length(GENES))
    names(hgnc_id)<-GENES

    entrez_id<-rep('',length(GENES))
    names(entrez_id)<-GENES

    ensemble_id<-rep('',length(GENES))
    names(ensemble_id)<-GENES

    location<-rep('',length(GENES))
    names(location)<-GENES

    family<-rep('',length(GENES))
    names(family)<-GENES

    pubmed_id<-rep('',length(GENES))
    names(pubmed_id)<-GENES

    string_id<-rep('',length(GENES))
    names(string_id)<-GENES

    description[commong]<-fc[,'name']
    family[commong]<-fc[,'gene_family']
    hgnc_id[commong]<-fc[,'hgnc_id']
    entrez_id[commong]<-fc[,'entrez_id']
    ensemble_id[commong]<-fc[,'ensembl_gene_id']
    location[commong]<-fc[,'location_sortable']
    pubmed_id[commong]<-fc[,'pubmed_id']
    string_id[commong]<-fc[,'string_id']
    res<-cbind(description,family,hgnc_id,entrez_id,ensemble_id,location,pubmed_id,string_id)

    return(res)
}
OT15_retrieveKinases_iGof_GSKtar<-function(GENES){

    load(kinases.fn)
    genes<-GENES
    uniprotKinase<-is.element(genes,kinase_genes_and_syn)
    names(uniprotKinase)<-GENES

    load(humKin.fn)
    genes<-GENES
    inHumanKinome<-is.element(genes,humKinome)
    names(inHumanKinome)<-GENES

    uniprotKinase<-rep(FALSE,length(GENES))
    names(uniprotKinase)<-GENES
    uniprotKinase<-is.element(names(uniprotKinase),kinase_genes_and_syn)
    names(uniprotKinase)<-GENES

    load(intoGen_GoF.fn)

    igGOF<-rep(FALSE,length(GENES))
    names(igGOF)<-GENES
    igGOF<-is.element(names(igGOF),intoGen_GoF)
    names(igGOF)<-GENES

    load(GSK_epiTargets.fn)
    GSK_EPIG_TAR<-rep(FALSE,length(GENES))
    names(GSK_EPIG_TAR)<-GENES
    GSK_EPIG_TAR<-is.element(names(GSK_EPIG_TAR),GSK_epiTargets)
    names(GSK_EPIG_TAR)<-GENES

    res<-cbind(inHumanKinome,uniprotKinase,igGOF,GSK_EPIG_TAR)
    return(res)
}
enrichedGeneFamilies<-function(geneset,BGgs){
    geneInfos<-OT15_retrievegeneInfo(BGgs)[,2]
    tokened<-str_split(geneInfos,'[|]')
    names(tokened)<-names(geneInfos)

    observed<-tokened[geneset]

    observed<-observed[observed!='']

    k<-length(observed)

    tokened<-tokened[tokened!='']

    N<-length(tokened)

    toTest<-summary(as.factor(unlist(observed)))
    toTest<-sort(toTest[which(toTest>1)],decreasing=TRUE)
    toTest<-names(toTest)

    nt<-length(toTest)


    pvals<-vector()
    GENES<-vector()
    for (i in 1:nt){
        print(i)
        n<-length(which((unlist(lapply(lapply(tokened,'intersect',toTest[i]),'length')))>0))
        currSet<-names(which((unlist(lapply(lapply(observed,'intersect',toTest[i]),'length')))>0))
        x<-length(currSet)

        pvals[i]<-my.hypTest(x,k,n,N)

        GENES[i]<-paste(sort(currSet),collapse=', ')
    }

    names(pvals)<-toTest


    id<-which(pvals<0.05)

    pvals<-pvals[id]
    GENES<-GENES[id]

    GENES<-GENES[order(pvals)]
    pvals<-sort(pvals)

    RES<-cbind(names(pvals),pvals,GENES)
    RES<-as.data.frame(RES)
    RES$pvals<-as.numeric(as.character(RES$pvals))
    colnames(RES)[1]<-'family'
    return(RES)
}


enrichedPathways<-function(geneset,BGgs,ming=2){

    k<-length(geneset)
    N<-length(BGgs)

    pvals<-vector()
    GENES<-vector()

    nt<-length(PATHCOM_HUMAN$PATHWAY)
    flag<-1


    ii<-unlist(lapply(lapply(PATHCOM_HUMAN$HGNC_SYMBOL,is.element,geneset),'sum'))
    names(ii)<-NULL

    toTest<-which(ii>=ming)

    nt<-length(toTest)
    testedP<-NULL
    for (i in 1:nt){
        print(i)
        n<-length(intersect(BGgs,PATHCOM_HUMAN$HGNC_SYMBOL[[toTest[i]]]))

        if (n >=ming){

            currSet<-intersect(geneset,PATHCOM_HUMAN$HGNC_SYMBOL[[toTest[i]]])
            x<-length(currSet)

            pvals[flag]<-my.hypTest(x,k,n,N)

            GENES[flag]<-paste(sort(currSet),collapse=', ')
            flag<-flag+1
            testedP<-c(testedP,as.character(PATHCOM_HUMAN$PATHWAY[[toTest[i]]]))
        }
    }


    id<-which(p.adjust(pvals,'fdr')<0.05)

    names(pvals)<-testedP
    pvals<-pvals[id]
    GENES<-GENES[id]

    GENES<-GENES[order(pvals)]
    pvals<-sort(pvals)

    RES<-cbind(names(pvals),pvals,GENES)
    RES<-data.frame(RES,stringsAsFactors = FALSE)
    RES$pvals<-as.numeric(as.character(RES$pvals))
    colnames(RES)[1]<-'pathway'
    rownames(RES)<-NULL
    return(RES)
}
