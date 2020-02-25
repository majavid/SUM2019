`learn.mvr.pclike`<-function(suffStat, indepTest, alpha, labels, p,
                             fixedGaps = NULL, fixedEdges = NULL, NAdelete = TRUE, m.max = Inf,
                             skel.method = c("stable", "original", "stable.fast"),
                             conservative = FALSE, maj.rule = FALSE,
                             solve.confl = FALSE, numCores = 1, verbose = FALSE,min.bidirected=FALSE){
  ## Purpose: Perform PC-like Algorithm, i.e., estimate essential graph of MVR CGs given data
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## - dm: Data matrix (rows: samples, cols: nodes)
  ## - C: correlation matrix (only for continuous)
  ## - n: sample size
  ## - alpha: Significance level of individual partial correlation tests
  ## - corMethod: "standard" or "Qn" for standard or robust correlation
  ##              estimation
  ## - G: the adjacency matrix of the graph from which the algorithm
  ##      should start (logical)
  ## - datatype: distinguish between discrete and continuous data
  ## - NAdelete: delete edge if pval=NA (for discrete data)
  ## - m.max: maximal size of conditioning set
  ## - gTrue: Graph suffStatect of true DAG
  ## - conservative: If TRUE, conservative PC is done
  ## - numCores: handed to skeleton(), used for parallelization
  ## - min.bidirected: If TRUE, essential graph to min.bidirected MVR CG is done
  ## ----------------------------------------------------------------------
  ## Author: Mohammad Ali Javidian, Date: 23 June 2019; 
  
  ## Initial Checks
  cl <- match.call()
  if(!missing(p)) stopifnot(is.numeric(p), length(p <- as.integer(p)) == 1, p >= 2)
  if(missing(labels)) {
    if(missing(p)) stop("need to specify 'labels' or 'p'")
    labels <- as.character(seq_len(p))
  } else { ## use labels ==> p  from it
    stopifnot(is.character(labels))
    if(missing(p)) {
      p <- length(labels)
    } else if(p != length(labels))
      stop("'p' is not needed when 'labels' is specified, and must match length(labels)")
    else
      message("No need to specify 'p', when 'labels' is given")
  }
  
  if (conservative && maj.rule) stop("Choose either conservative PC or majority rule PC!")
  
  ## Skeleton
  skel <- skeleton(suffStat, indepTest, alpha, labels = labels, method = skel.method,
                   fixedGaps = fixedGaps, fixedEdges = fixedEdges,
                   NAdelete=NAdelete, m.max=m.max, numCores=numCores, verbose=verbose)
  skel@call <- cl # so that makes it into result
  zmat<-as(skel@graph,"matrix")
  ## Orient edges
  if (!conservative && !maj.rule) {
    #amat<-as(skel@graph,"matrix")####gives adjacency matrix of pcalg
    # sepsets<-skel@sepset
    # print(sepsets)
    zmat<- ug2essential(skel, verbose = verbose, solve.confl = solve.confl)
  }
  else { ## conservative _or_ maj.rule
    
    ## version.unf defined per default
    ## Tetrad CPC works with version.unf=c(2,1)
    ## see comment on pc.cons.intern for description of version.unf
    pc. <- pc.cons.intern(skel, suffStat, indepTest, alpha,
                          version.unf = c(2,1), maj.rule = maj.rule, verbose = verbose)
    zmat<-ug2essential(pc.$sk, verbose = verbose,
                     unfVect = pc.$unfTripl, solve.confl = solve.confl)
    # udag2pdagRelaxed(pc.$sk, verbose = verbose,
    #                  unfVect = pc.$unfTripl, solve.confl = solve.confl)
  }
  
  if(min.bidirected){
    ###################################################################################################################################
    ###A graph is triangulated (chordal) if and only if there exists a perfect ordering (perfect elimination sequence ) of its vertices.
    ###Lines 10-13 ####################################################################################################################
    ###################################################################################################################################
    cmat<-zmat
    for(i in 1:(nrow(zmat)-1)){
      for(j in (i+1):(nrow(zmat))){
        if(zmat[i,j]!=1 || zmat[j,i]!=1){
          cmat[i,j]<-cmat[j,i]<-0
        }
      }
    }
    
    if(!all(cmat == 0)){  
      # print(cmat)
      # return(print(maxcard.search(cmat)$perfect.numbering))
      pOrdering<-maxcard.search(cmat)$perfect.numbering
      for(i in 1:(nrow(zmat)-1)){
        for(j in (i+1):(nrow(zmat))){
          ###################################################################################
          ###line15: Orient the undirected edges according to the obtained perfect ordering 
          ###################################################################################
          if(cmat[i,j]==1 && (which(pOrdering==i)<which(pOrdering==j))){
            zmat[i,j]<-1
            zmat[j,i]<-0
          }
          if(cmat[i,j]==1 && (which(pOrdering==i)>which(pOrdering==j))){
            zmat[i,j]<-0
            zmat[j,i]<-1
          }
        }
      }
    }
  }
  ############### Adjacency matrices for mixed graphs(ggm format) ########
  for(i in 1:(nrow(zmat)-1)){
    for(j in (i+1):(nrow(zmat))){
      if(zmat[i,j]==1 && zmat[j,i]==1){
        zmat[i,j]<- zmat[j,i]<- 10 
      }else if(zmat[i,j]==1 && zmat[j,i]==2){
        zmat[i,j]<- 1 
        zmat[j,i]<- 0
      }else if(zmat[i,j]==2 && zmat[j,i]==1){
        zmat[i,j]<- 0 
        zmat[j,i]<- 1
      }else if(zmat[i,j]==2 && zmat[j,i]==2){
        zmat[i,j]<- zmat[j,i]<- 100
      }else{
        zmat[i,j]<- zmat[i,j] 
        zmat[j,i]<- zmat[j,i]
      }
    }
  }
  zmat#Adjacency matrices for mixed graphs(ggm format)
}