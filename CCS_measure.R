
cliqueSimilarity<- function (cl_mem, WorA) 
{
    stopifnot(length(dim(WorA)) == 2)
    stopifnot(nrow(WorA) == ncol(WorA))
    stopifnot(length(cl_mem) == ncol(WorA))
    non0 <- cl_mem[cl_mem != 0]
    labelsOK = all(sort(unique(non0)) == c(1:max(cl_mem)))
    if (!labelsOK) {
        old_mem <- cl_mem
        cl_mem <- tidyUpLabels(cl_mem)
        labelMap = data.frame(new_label = cl_mem, old_label = old_mem)
        labelMap = labelMap[!duplicated(LabelMap), ]
        labelMap = labelMap[order(LabelMap$new_label), ]
    }
    diag(WorA) = 0
    A = (WorA > 0) * 1
    C_S = cs_matrix(cl_mem = cl_mem, A = A)
    if (!labelsOK) 
        attr(C_S, "labelMap") = labelMap
    C_S
}

cs_matrix<- function (cl_mem, A) 
{
    stopifnot(length(unique(cl_mem)) <= ncol(A))
    stopifnot(all(A %in% c(0, 1)))
    cl_mem_uq <- 1:max(cl_mem)
    componentIndicator <- do.call(cbind, lapply(cl_mem_uq, function(k) cl_mem == 
        k))
    n_con <- A %*% componentIndicator
    numer <- t(componentIndicator) %*% n_con
    compSizes <- colSums(componentIndicator)
    denom <- compSizes %o% compSizes
    diag(numer) <- diag(denom) <- 1
    numer/denom
}


expand_CS_toNodes<- function(CS,cl_mem){

if (0 %in% cl_mem) 
	stop("0 labels detected in cl_mem, call uniqueSingletonLabels(cl_mem)")
uq_cl_mem<- sort(unique(cl_mem))
labels_OK<-all(sort(unique(cl_mem))==(1:max(cl_mem)))

if (!labels_OK)
	stop("labels of cliques should be consecutive integers 1:n_cliques")

node_CS<-matrix(nrow=length(cl_mem), ncol=length(cl_mem))
#print(sum(is.na(node_CS)))
clique_placements<- lapply(uq_cl_mem, function(k) cl_mem==k)

for (k_i in seq_along(uq_cl_mem)){
	cl_i<- clique_placements[[k_i]]
#	print(sum(cl_i))
	node_CS[ cl_i, cl_i ] = 1
#	print(sum(is.na(node_CS)))
	for (k_j in seq_along(uq_cl_mem)) {
		if (k_i < k_j){
		cl_j<- clique_placements[[k_j]]
		node_CS[cl_i,cl_j]= CS[ k_i, k_j ]
		node_CS[cl_j,cl_i]= CS[ k_j, k_i ]
		}}
}
node_CS
}

#unlapply<- function(...) unlist(lapply(...))
#
#complexity_score<-function (clique_labels) 
#{
#  clique_labels <- uniqueSingletonLabels(clique_labels)
#  tabulation <- fastTable(clique_labels)
#  -(lfactorial(length(tabulation$value)) + sum(lfactorial(tabulation$count)))
#}
#
#zeroOutSmall<- function (labelV, mS) {
#  rment <- labelV
#  non0labs <- unique(labelV[labelV != 0])
#  non0cnts <- unlapply(non0labs, function(lab) sum(labelV == 
#                                                     lab))
#  for (l in seq_along(non0labs)) if (non0cnts[[l]] < mS) 
#    rment[labelV == non0labs[[l]]] = 0
#  return(rment)
#}
#
#
#
#
#
#clique_consensus_similairty<- function(GE, B, similarity_fun, seed) {
#  set.seed(seed)
#  bootstrap_indexes<- lapply(1:B, function(b) sample(1:nrow(GE),
#                                                     nrow(GE),
#                                                     replace=TRUE) )
#  CCS<- matrix(0, nrow=ncol(GE),ncol=ncol(GE))
#  for (b in seq_along(bootstrap_indexes)){
#    similarity_fun(GE[bootstrap_indexes[[b]], ])-> S_b
#    hclust(as.dist(max(S_b) - S_b ), method="complete")-> hcl
#    optimize_hcl(hcl)-> best_cut_idx
#    P<- cutree(hcl, h= hcl$height[[best_cut_idx]])
#    D_b<- max(S_b) - S_b
#    D_b [ D_b > hcl$height[[best_cut_idx]] ]=0
#    clqS_b<- cliqueSimilarity(cl_mem=P, WorA=D_b)
#    CCS= CCS + expand_CS_toNodes(clqS_b, cl_mem = P)
#  }
#  CCS/B

corsq<- function(x) similarity_matrix(fastPearsonData(x))

complexity<- function(P) {
  labs<- unique(P)
  (sum(c( lfactorial(length(labs) ),
       sapply(labs, function(clq) lfactorial(sum(P==clq)) )
  )))
  
}

objective<- function(x,hcl){
  cutree(hcl,h=x)->P
  f_x<- complexity(P)
#  print(c(x,f_x))
  return(f_x)
}


exhaustively_optimize_hcl<- function(hcl, loss=complexity){
 sapply(hcl$height, function(x) {
    P <- cutree(hcl, h = x)
    loss(P)
  })-> scores
  hcl$height[[which.min(loss)]]
}

bobyqa_optimize_hcl<- function(hcl, I,xtol_abs){
bobyqa(x0 = mean(I), fn = objective, lower = min(I), upper = max(I),
       control = list(xtol_abs=xtol_abs,
                      xtol_rel=0
                      ), hcl=hcl)-> search_res
search_res$par
}


clique_consensus_similarity<- function(GE, bootstrap_indexes, similarity_fun, optimisation_mode="bobyqa" ) {
  B<- length(bootstrap_indexes)
  CCS<- matrix(0, nrow=ncol(GE),ncol=ncol(GE))
  for (b in seq_along(bootstrap_indexes)){
    similarity_fun(GE[bootstrap_indexes[[b]], ])-> S_b
    hclust(as.dist(max(S_b) - S_b ), method="complete")-> hcl
    if (optimisation_mode=="bobyqa")
	h_<- bobyqa_optimize_hcl(hcl=hcl,I= c(min(S_b[S_b>0]),1),nrow(GE)) 
    if (optimisation_mode=="exhaustive")
	h_<- exhaustively_optimize_hcl(hcl=hcl,scorer=complexity)
    P<- cutree(hcl, h= h_)
    D_b<- max(S_b) - S_b
    D_b [ D_b > h_ ]=0
    clqS_b<- cliqueSimilarity(cl_mem=P, WorA=D_b)
    CCS= CCS + expand_CS_toNodes(clqS_b, cl_mem = P)
  }
  CCS/B
}









