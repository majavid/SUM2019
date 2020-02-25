`similarity.pat` <- function(source, target)
{ ###Compares a (learned) chain graph pattern to the (supposed) true pattern. The two patterns should
  ###have the same vertex set in order for the function to return a meaningful result.
  vset <- colnames(source)
  #rownames(source)<-vset
  ###true directed edges
  #truearr <- which(source - t(source) == 1)
  ###true bidirected edges
  #truebow<-which(source == 100)
  ###true undirected edges
  #truearc<-which(source == 10)
  target <- target[vset, vset]
  `skelet` <- function(amat)
  {
    0 + (amat + t(amat) > 0)
  }
  ###true learned skelet
  trueskel<-which(skelet(source)+skelet(target)==2)
  truepos<-length(trueskel)/2
  ###sensitivity, recall, hit rate, or true positive rate (TPR): TPR=TP/P=TP/TP+FN
  ####condition positive (P): the number of real positive cases
  #P<-length(truearr)+(length(truebow)/2)+(length(truearc)/2)
  P<-sum(skelet(source))/2
  TPR<-truepos/P
  total<-choose(length(vset),2)
  falsepos<-(length(which(skelet(target)!=0))/2)-truepos
  ####fall-out or false positive rate (FPR): FPR=FP/N=FP/FP+TN
  ####condition negative (N): the number of real negative cases
  N<-total-P
  FPR<-falsepos/N
  ####accuracy (ACC): ACC=TP+TN/P+N
  ####N=TN+FP
  trueneg<-N-falsepos
  ACC<-(truepos+trueneg)/(P+N)
  ####False Negative: P=TP+FN
  falseneg<-P-truepos
  
  t.total<-sum(skelet(target))/2
  TDR<-truepos/t.total
  
  ## computing structural Hamming distance
  #Structural Hamming distance is defined as the total number of operations needed to convert one
  # graph to the other. Each operation must be one of the following: (1) add or delete an
  # edge, or (2) add, remove or reverse an orientation of an edge.
  r=nrow(source)
  shd<-0
  for (i in 1:(r-1)) {
    for (j in (i+1):r) {
      if((source[i,j]!=0 || source[j,i]!=0)&&(target[i,j]==0)&&(target[j,i]==0)){
        ###missing edge in the target: add appropriate edge
        shd<-shd+1
      }
      if((target[i,j]!=0 || target[j,i]!=0)&&(source[i,j]==0)&&(source[j,i]==0)){
        ###extra edge in the target: remove it
        shd<-shd+1
      }
      if((source[i,j]+source[j,i]==1)&&(target[i,j]+target[j,i]>1)){
        ###there is a directed edge in the source graph, but the corresponding edge in target is un/bi-directed: add/remove orientation 
        shd<-shd+1
      }
      if((source[i,j]+source[j,i]==1)&&(target[i,j]+target[j,i]==1)&&(source[i,j]!=target[i,j])){
        ###-->/<-- in the source graph, but <--/--> in the target, respectively: reverse the orientation 
        shd<-shd+1
      }
      if((source[i,j]+source[j,i]==200)&&((target[i,j]+target[j,i]==1)||(target[i,j]+target[j,i]==20))){
        ###there is a bidirected edge in the source, but the corresponding edge in target is directed -->/<-- or undirected: add appropriate orientation
        shd<-shd+1
      }
    }
  }
  return(list(TP = truepos,
              FN = falseneg,
              FP = falsepos,
              TN = trueneg,
              TPR = TPR,
              TDR = TDR,
              FPR = FPR,
              ACC = ACC,
              SHD = shd))
}