`ug2essential` <- function(gInput, verbose = FALSE, unfVect = NULL, solve.confl = FALSE, orientCollider = TRUE, rules = rep(TRUE, 3))
{
  
  ##################################################
  ## Internal functions
  ##################################################

  ## TODO: include correct VERBOSE statements
  rule1 <- function(pdag, solve.confl = FALSE, unfVect = NULL) {
    ## Rule 1: b\in S_ac  a o-> b - c goes to a o-> b -> c
    ## Interpretation: No new collider is introduced
    ## Out: Updated pdag
    search.pdag <- pdag
    ind <- which(pdag != 0 & t(pdag) == 2, arr.ind = TRUE)
    for (i in seq_len(nrow(ind))) {
      a <- ind[i, 1]
      b <- ind[i, 2]
      ## find all undirected neighbours of b not adjacent to a
      isC <- which(search.pdag[b, ] == 1 & search.pdag[, b] == 1 &
                     search.pdag[a, ] == 0 & search.pdag[, a] == 0)
      if (length(isC) > 0) {
        for (ii in seq_along(isC)) {
          c <- isC[ii]
          ##check that b\in S_ac: i.e., b is a non-collider node
          if((b %in% gInput@sepset[[a]][[c]]) ||
              (b %in% gInput@sepset[[c]][[a]])){
            
          
          ## if the edge between b and c has not been oriented previously,
          ## orient it using normal R1
          if (!solve.confl | (pdag[b,c] == 1 & pdag[c,b] == 1) ) { ## no conflict
            ## !! before, we checked search.pdag, not pdag !!
            if (!is.null(unfVect)) { ## deal with unfaithful triples
              if (!any(unfVect == triple2numb(p, a, b, c), na.rm = TRUE) &&
                  !any(unfVect == triple2numb(p, c, b, a), na.rm = TRUE)) {
                ## if unfaithful triple, don't orient
                pdag[c,b] <- 2
              }
            } else {
              ## don't care about unfaithful triples -> just orient
              pdag[c,b] <- 2
              ## cat("Rule 1\n")
            }
            if (verbose)
              cat("\nRule 1':", a, "->", b, " and ",
                  b, "-", c, " where ", a, " and ", c,
                  " not connected and ", a, b, c, " faithful triple: ",
                  b, "->", c, "\n")
          } else if (pdag[b,c] == 2 & pdag[c,b] == 1) {
            ## conflict that must be solved
            ## solve conflict: if the edge is b <- c because of a previous
            ## orientation within for loop then output <->
            if (!is.null(unfVect)) { ## deal with unfaithful triples
              if (!any(unfVect == triple2numb(p, a, b, c), na.rm = TRUE) &&
                  !any(unfVect == triple2numb(p, c, b, a), na.rm = TRUE)) {
                pdag[b, c] <- 2
                pdag[c, b] <- 2
                if (verbose)
                  cat("\nRule 1':", a, "->", b, "<-",
                      c, " but ", b, "->", c, "also possible and",
                      a, b, c, " faithful triple: ", a,"->", b, "<->", c,"\n")
              }
            } else {
              ## don't care about unfaithful triples -> just orient
              pdag[b, c] <- 2
              pdag[c, b] <- 2
              if (verbose)
                cat("\nRule 1':", a, "->", b, "<-",
                    c, " but ", b, "->", c, "also possible and",
                    a, b, c, " faithful triple: ", a,"->", b, "<->", c,"\n")
            } ## unfVect: if else
          } ## conflict: if else
          }##b\in S_ac
        } ## for c
      } ## if length(isC)
      if (!solve.confl) search.pdag <- pdag
    } ## for ind
    pdag
  }
  
  rule2 <- function(pdag, solve.confl = FALSE) {
    ## Rule 2: a o-> c o-> b with a o- b: a o-> b
    ## Interpretation: Avoid cycle
    ## normal version = conservative version
    search.pdag <- pdag
    ind <- which(search.pdag != 0 & t(search.pdag) == 1, arr.ind = TRUE)
    for (i in seq_len(nrow(ind))) {
      a <- ind[i, 1]
      b <- ind[i, 2]
      isC <- which(search.pdag[a, ] != 0 & search.pdag[, a] == 2 &
                     search.pdag[, b] != 0 & search.pdag[b, ] == 2)
      for (ii in seq_along(isC)) {
        c <- isC[ii]
        ## if the edge has not been oriented yet, orient it with R2
        ## always do this if you don't care about conflicts
        if (!solve.confl | (pdag[a, b] != 0 & pdag[b, a] == 1) ) {
          pdag[b, a] <- 2
          if (verbose)
            cat("\nRule 2: Chain ", a, "o->", c,
                "o->", b, ":", a, "o->", b, "\n")
        }
        # ## else if the edge has been oriented as a <- b by a previous R2
        # else if (pdag[a, b] == 0 & pdag[b, a] == 1) {
        #   pdag[a, b] <- 2
        #   pdag[b, a] <- 2
        #   if (verbose)
        #     cat("\nRule 2: Chain ", a, "->", c,
        #         "->", b, ":", a, "<->", b, "\n")
        # }
      }
      if (!solve.confl) search.pdag <- pdag
    }
    pdag
  }
  
  rule3 <- function(pdag, solve.confl = FALSE, unfVect = NULL) {
    ## Rule 3: a-b, a o-o c1, a o-o c2, c1 o->b, c2 o->b but c1 and c2 not connected and a\in S_c1c2;
    ## then a-b => a -> b
    search.pdag <- pdag
    ind <- which(search.pdag == 1 & t(search.pdag) == 1, arr.ind = TRUE)
    for (i in seq_len(nrow(ind))) {
      a <- ind[i, 1]
      b <- ind[i, 2]
      isc1 <- which(search.pdag[a, ] != 0 & search.pdag[, a] != 0 &
                    search.pdag[, b] != 0 & search.pdag[b, ] == 2)
      isc2 <- which(search.pdag[a, ] == 1 & search.pdag[, a] == 1 &
                   search.pdag[, b] != 0 & search.pdag[b, ] == 2)
        for (j in seq_along(isc1)) {
          c1 <- isc1[j]
          for (k in seq_along(isc2)) {
          c2 <- isc2[k]
          # message("c1=", c1, " c2=",c2)
          # message("gInput@sepset[[c1]][[c2]]: ",gInput@sepset[[c1]][[c2]])
          # message("gInput@sepset[[c2]][[c1]]: ",gInput@sepset[[c2]][[c1]])
          if (search.pdag[c1, c2] == 0 && search.pdag[c2,c1] == 0 && ((a %in% gInput@sepset[[c1]][[c2]]) ||
              (a %in% gInput@sepset[[c2]][[c1]]))) {
            
            if (!is.null(unfVect)) {
              if (!any(unfVect == triple2numb(p, c1, a, c2), na.rm = TRUE) &&
                  !any(unfVect == triple2numb(p, c2, a, c1), na.rm = TRUE)) {
                ## if the edge has not been oriented yet, orient it with R3
                if (!solve.confl | (pdag[a, b] == 1 & pdag[b, a] == 1) ) {
                  pdag[b, a] <- 2
                  if (!solve.confl) search.pdag <- pdag
                  if (verbose)
                    cat("\nRule 3':", a, c1, c2, "faithful triple: ",
                        a, "->", b, "\n")
                  break
                }
                ## else if: we care about conflicts and  the edge has been oriented as a <- b by a previous R3
                else if (pdag[a, b] == 2 & pdag[b, a] == 1) {
                  pdag[a, b] <- pdag[b, a] <- 2
                  if (verbose)
                    cat("\nRule 3':", a, c1, c2, "faithful triple: ",
                        a, "<->", b, "\n")
                  break
                } ## if solve conflict
              } ## if unf. triple found
            } else { ## if care about unf. triples; else don't care
              if (!solve.confl | (pdag[a, b] == 1 & pdag[b, a] == 1) ) {
                pdag[b, a] <- 2
                if (!solve.confl) search.pdag <- pdag
                if (verbose)
                  cat("\nRule 3':", a, c1, c2, "faithful triple: ",
                      a, "->", b, "\n")
                break
              }
              ## else if: we care about conflicts and  the edge has been oriented as a <- b by a previous R3
              else if (pdag[a, b] == 2 & pdag[b, a] == 1) {
                pdag[a, b] <- pdag[b, a] <- 2
                if (verbose)
                  cat("\nRule 3':", a, c1, c2, "faithful triple: ",
                      a, "<->", b, "\n")
                break
              } ## if solve conflict
            } ## if care about unf. triples
          } ## if c1 and c2 are not adjecent
        } ## for all of c2's
        } ## for all of c1's  
      #} ## if at least two c's are found
    } ## for all undirected edges
    pdag
  }
  ##################################################
  ## Main
  ##################################################
  
  ## prepare adjacency matrix of skeleton
  if (graph::numEdges(gInput@graph) == 0)
    return(gInput)
  g <- as(gInput@graph, "matrix")
  p <- nrow(g)
  pdag <- g
  
  ## orient collider
  if (orientCollider) {
    ind <- which(g == 1, arr.ind = TRUE)
    for (i in seq_len(nrow(ind))) {
      x <- ind[i, 1]
      y <- ind[i, 2]
      allZ <- setdiff(which(g[y, ] == 1), x) ## x - y - z
      for (z in allZ) {
        ## check collider condition
        if (g[x, z] == 0 &&
            !((y %in% gInput@sepset[[x]][[z]]) ||
              (y %in% gInput@sepset[[z]][[x]]))) {
          if (length(unfVect) == 0) { ## no unfaithful triples
            pdag[y,x] <- pdag[y,z] <- 2
          } else { ## unfaithful triples are present
            if (!any(unfVect == triple2numb(p, x, y, z), na.rm = TRUE) &&
                !any(unfVect == triple2numb(p, z, y, x), na.rm = TRUE)) {
              pdag[y,x] <- pdag[y,z] <- 2
            }
          }
        }
      } ## for z
    } ## for i
    #print(pdag)
  } ## end: Orient collider
  
  ## Rules 1 - 3
  repeat {
    old_pdag <- pdag
    if (rules[1]) {
      pdag <- rule1(pdag, solve.confl = solve.confl, unfVect = unfVect)
    }
    if (rules[2]) {
      pdag <- rule2(pdag, solve.confl = solve.confl)
    }
    if (rules[3]) {
      pdag <- rule3(pdag, solve.confl = solve.confl, unfVect = unfVect)
    }
    if (all(pdag == old_pdag))
      break
  } ## repeat
  
  # gInput@graph <- as(pdag, "graphNEL")
  # gInput
  pdag
  
}