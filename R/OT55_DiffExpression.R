GSEAfunction<-function(rankedlist,geneset){
  guse<-intersect(rankedlist,geneset)
  Inset<-(rankedlist%in%guse)*1
  sizeset<-length(guse)
  N_r<-sizeset
  PosScores<-Inset/N_r
  NegScores<-(!rankedlist%in%guse)*1
  N_Nh<-length(rankedlist)-sizeset
  NegScores<-NegScores/N_Nh
  Allscores<-PosScores-NegScores
  RunningSum<-cumsum(Allscores)
  maxdev<-RunningSum[which.max(abs(RunningSum))]
  return(list(ESscore=maxdev,RunningSum=RunningSum))
}

testGmt<-function(gmtlist,DElist,ExprsList){
  gmtExpr<-sum(tolower(gmtlist)%in%tolower(ExprsList))
  over<-sum(tolower(gmtlist)%in%tolower(DElist))
  
  phyper(over-1,gmtExpr,length(ExprsList)-gmtExpr,length(DElist),lower.tail=FALSE)
}

EMscore<-function(FCs,geneexpression,gmtfile,dep,OUTPUT_DIR){

  emout<-tryCatch(normalmixEM(FCs,lambda=0.5,mu=c(0,-1),sigma=c(0.5,1)),error=function(e){NA})
  
  if(!is.na(emout)){
    mus<-emout$mu
    nonsig<-which.max(mus)
    #nnmd_em<-(mus[-nonsig]-mus[nonsig])/emout$sigma[nonsig]
    
    clgroup<-apply(emout$posterior,1,which.max)
    g1size<-sum(clgroup==1)
    g2size<-sum(clgroup==2)
    nnmd_em<-ANOVA_cohens_d(FCs[clgroup==1],FCs[clgroup==2])
    if(g1size>2&g2size>2){
      clgroup<-as.factor(clgroup)
      names(clgroup)<-names(FCs)
      DgeAll<-DGEList(counts=geneexpression,samples=clgroup[colnames(geneexpression)],genes=rownames(geneexpression))
      designmatrix<-model.matrix(~clgroup)
      DgeAll<-estimateDisp(DgeAll,designmatrix)
      fitAll<-glmFit(DgeAll,designmatrix)
      lrtAll<-glmLRT(fitAll,coef=2)
      ttSig<-topTags(lrtAll,n=nrow(lrtAll),p.value=0.05)
      ttAll<-topTags(lrtAll,n=nrow(lrtAll),p.value=1)
      if(nrow(ttSig)>20){
        logFC<-ttAll$table[,"logFC"]
        names(logFC)<-ttAll$table[,'genes']
        genes<-tolower(names(sort(logFC,decreasing=TRUE)))
        gseares<-lapply(gmtfile,function(x) GSEAfunction(genes,tolower(x))$ESscore)
    
        hyperres<-lapply(gmtfile,function(x) testGmt(x,ttSig$table[,'genes'],rownames(geneexpression)))
        gseares<-sort(abs(unlist(gseares)),decreasing=TRUE)
        topgsea<-gseares[1:5]
        adjhyper<-p.adjust(unlist(hyperres))
        tophyperres<-adjhyper[adjhyper<0.05]
        if(length(tophyperres)>5){
          tophyperres<-sort(tophyperres)[1:5]
        }
        #pdf(paste0(OUTPUT_DIR,"/EMsplit_",dep,".pdf"))

        #  cols <- rep("darkgray", length(FCs))
       #   cols[which(clgroup == 1)] <- "darkgreen"
       #   beeswarm(FCs ~ clgroup, corral = "wrap", bg = c(makeTransparent("gray"), 
         #                                                  makeTransparent("darkgreen")), pch = 21, col = c("gray", 
        #                                                                                                    "darkgreen"), cex = 1.5, ylim = range(FCs), las = 2, 
       #            labels = c("g1", "g2"), xlab = dep, ylab = "",axes=F)
        #  par(new = TRUE)
        #  
        #  boxplot(FCs ~ clgroup, col = NA, ylim = range(FCs), frame.plot = FALSE, 
        #          xaxt = "n",  outline = FALSE,tcl=0.5,tck=-0.01)
       # dev.off()
        return(data.frame(nnmd=nnmd_em,gsea=paste0(topgsea,collapse="//"),gname=paste0(names(topgsea),collapse="//"),hyper=paste0(tophyperres,collapse="//"),hypername=paste0(names(tophyperres),collapse = "//"),dep=dep,DEgenes=paste0(ttSig$table[,'genes'],collapse="//")))
      
      }else{
        return(data.frame(nnmd=nnmd_em,gsea="NoDE",gname="NoDE",hyper="NoDE",hypername="NoDE",dep=dep,DEgenes="NoDE"))
      }
    }else{
      return(data.frame(nnmd=nnmd_em,gsea="NoBimodal",gname="NoBimodal",hyper="NoBimodal",hypername="NoBimodal",dep=dep,DEgenes="NoBimodal"))
    }
  }else{
    return(data.frame(nnmd="NoEM",gsea="NoEM",gname="NoEM",hyper="NoEM",hypername="NoEM",dep=dep,DEgenes="NoEM"))
    
  }
}



