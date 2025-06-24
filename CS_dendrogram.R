
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


