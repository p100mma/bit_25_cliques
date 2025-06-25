
unlapply<- function(...) unlist(lapply(...))

complexity_score<-function (clique_labels) 
{
  clique_labels <- uniqueSingletonLabels(clique_labels)
  tabulation <- fastTable(clique_labels)
  -(lfactorial(length(tabulation$value)) + sum(lfactorial(tabulation$count)))
}

zeroOutSmall<- function (labelV, mS) {
  rment <- labelV
  non0labs <- unique(labelV[labelV != 0])
  non0cnts <- unlapply(non0labs, function(lab) sum(labelV == 
                                                     lab))
  for (l in seq_along(non0labs)) if (non0cnts[[l]] < mS) 
    rment[labelV == non0labs[[l]]] = 0
  return(rment)
}


optimize_hcl<- function(hcl, scorer=complexity_score){
 sapply(hcl$height, function(x) {
    P <- cutree(hcl, h = x)
    scorer(P)
  })-> scores
  which.max(scores)
  
}

mst_thresholds<- function(mst_graph){
  mst_D<- mst_graph
  sort( igraph::E(mst_D)$weight, decreasing=TRUE)-> dMSTsorted
  dMSTsorted<- dMSTsorted[ dMSTsorted!=0 ]
  ddiffs= -1*(diff(dMSTsorted))
  dmeans= (dMSTsorted[ 2: length(dMSTsorted) ] + dMSTsorted[ 1: (length(dMSTsorted)-1) ])/2
  dmeans[order(ddiffs, decreasing = TRUE) ]-> d_t
  d_t
  
}
