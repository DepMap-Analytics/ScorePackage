removeCoreFitness<-function(data,listToRemove){
    for(i in 1:length(listToRemove)){
      data<-data[setdiff(rownames(data),listToRemove[[i]]),]
    }
  return(data)
}
