library(HCCsim)
library(igraph)
library(matrixStats)
data("brca")
conductance_cluster<- function(A, cl_indicator){
  diag(A)=0
  c_s<- A[cl_indicator, 
          !cl_indicator] |> sum()
  m_s<- A[cl_indicator,
          cl_indicator] |> sum()
  c_s/(2*m_s + c_s)
}

conductance_network<- function(A, members){
  diag(A)=0
  1- (sapply(unique(members), function(x){
    conductance_cluster(A, members==x)
  }) |> mean())
  
}


brca_sub<- brca[, order(colVars(brca), decreasing=TRUE)[1:400] ]

library(apcluster)


G_S<- graph_from_adjacency_matrix(S, mode="max", weighted=TRUE, diag=FALSE)
modularity(G_S, weights=E(G_S)$weight, membership=ap_membership)

n_trials=10

ap_clusters=list()
ap_mods<-vector()
ap_cond<-vector()
for (j in 1:n_trials){
  set.seed(j)
  sample(1:nrow(brca_sub), nrow(brca_sub), replace=TRUE)-> sample_idx
  S<-cor(brca_sub[sample_idx, ], method="pearson") |> abs()
  res<-apcluster(s=S,seed=j)
  labels(res, type="enum")-> ap_membership
  G_S<- graph_from_adjacency_matrix(S, mode="max", weighted=TRUE, diag=FALSE)
  ap_mods[j]<-modularity(G_S, weights=E(G_S)$weight, membership=ap_membership)
  ap_clusters[[j]]<-ap_membership
  ap_cond[j]<-conductance_network(S, ap_membership)
  message(
    sprintf(
      "%d/%d trial(s) completed ",j,n_trials
    )
  )
}

ARI<- matrix(nrow=length(ap_clusters), ncol=length(ap_clusters), dimnames=list(1:length(ap_clusters), 1:length(ap_clusters)))
for (i in 1:length(ap_clusters)){
  for (j in 1:length(ap_clusters)){
    ARI[i,j]<- compare(ap_clusters[[i]], ap_clusters[[j]], "adjusted.rand")
  }
}

hist(ARI[lower.tri(ARI)])


n_trials=10
B=10
ap_mods_CS<-vector()
ap_clusters_CS=list()
ap_cond_CS<-vector()
for (j in 1:n_trials){
  CS= matrix(0, ncol=ncol(brca_sub),nrow=ncol(brca_sub))
  for (b in 1:B){
    set.seed(j*(b-1) + b)
    sample(1:nrow(brca_sub), nrow(brca_sub), replace=TRUE)-> sample_idx
    S<-cor(brca_sub[sample_idx, ], method="pearson") |> abs()
    res<-apcluster(s=S,seed=(j*(b-1)) + b)
    labels(res, type="enum")-> ap_membership
    
    CS<-CS + outer(ap_membership, ap_membership, function(x, y) as.numeric(x == y))
  }
  CS<- CS/B
  res_CS<-apcluster(s=CS,seed=100*B+j)
  labels(res_CS, type="enum")-> ap_membership_CS
  G_S<- graph_from_adjacency_matrix(S, mode="max", weighted=TRUE, diag=FALSE)
  ap_mods_CS[j]<-modularity(G_S, weights=E(G_S)$weight, 
                         membership=ap_membership_CS)
  ap_clusters_CS[[j]]<-ap_membership_CS
  ap_cond_CS[j]<-conductance_network(S, ap_membership_CS)
  message(
    sprintf(
      "%d/%d trial(s) completed ",j,n_trials
    )
  )
}
ARI_CS<- matrix(nrow=length(ap_clusters_CS), ncol=length(ap_clusters_CS), dimnames=list(1:length(ap_clusters_CS), 1:length(ap_clusters_CS)))
for (i in 1:length(ap_clusters_CS)){
  for (j in 1:length(ap_clusters_CS)){
    ARI_CS[i,j]<- compare(ap_clusters_CS[[i]], ap_clusters_CS[[j]], "adjusted.rand")
  }
}



ARI<- matrix(nrow=length(ap_clusters), ncol=length(ap_clusters), dimnames=list(1:length(ap_clusters), 1:length(ap_clusters)))
for (i in 1:length(ap_clusters)){
  for (j in 1:length(ap_clusters)){
    ARI[i,j]<- compare(ap_clusters[[i]], ap_clusters[[j]], "adjusted.rand")
  }
}

hist(ARI[lower.tri(ARI)])

source("cliques_rework/CS_dendrogram.R")
library(cliquePartitioneR)
n_trials=10
B=10
ap_mods_CCS<-vector()
ap_clusters_CCS=list()
ap_cond_CCS<-vector()
for (j in 1:n_trials){
  CCS= matrix(0, ncol=ncol(brca_sub),nrow=ncol(brca_sub))
  for (b in 1:B){
    set.seed(j*(b-1) + b)
    sample(1:nrow(brca_sub), nrow(brca_sub), replace=TRUE)-> sample_idx
    S<-cor(brca_sub[sample_idx, ], method="pearson") |> abs()
    hcl<- hclust(as.dist(1-S), method="complete")
    candidates<- list()
    for (h in hcl$height){
      candidates[[length(candidates)+1]]<- cutree(hcl, h=h)
      D<- 1 -S
      D[D>h]<-0
      stopifnot( 
         areCliques(D, cutree(hcl,h=h))
        )
    }
    
    lapply(candidates,cliqueClusteR:::complexity_clqScore)-> clq_scores
    best_idx<- which.max(clq_scores)
    clq_membership<- candidates[[best_idx]]
    cutHeight= hcl$height[[best_idx]]
    cliqueSimilarity(clq_membership, (S %thr% (1-cutHeight)) )->clqSim
    expand_CS_toNodes(clqSim, clq_membership)-> CCS_b
    CCS<-CCS + CCS_b
  }
  CCS<- CCS/B
  res_CCS<-apcluster(s=CCS,seed=100*B+j)
  labels(res_CS, type="enum")-> ap_membership_CCS
  G_S<- graph_from_adjacency_matrix(S, mode="max", weighted=TRUE, diag=FALSE)
  ap_mods_CCS[j]<-modularity(G_S, weights=E(G_S)$weight, 
                            membership=ap_membership_CCS)
  ap_clusters_CCS[[j]]<-ap_membership_CCS
  ap_cond_CCS[j]<-conductance_network(S, ap_membership_CCS)
  message(
    sprintf(
      "%d/%d trial(s) completed ",j,n_trials
    )
  )
}

ARI_CCS<- matrix(nrow=length(ap_clusters_CCS), 
                 ncol=length(ap_clusters_CCS), 
                 dimnames=list(1:length(ap_clusters_CCS),
                               1:length(ap_clusters_CCS)))
for (i in 1:length(ap_clusters_CCS)){
  for (j in 1:length(ap_clusters_CCS)){
    ARI_CCS[i,j]<- compare(ap_clusters_CCS[[i]], 
                          ap_clusters_CCS[[j]], "adjusted.rand")
  }
}


data.frame(stability= c(mean(ARI[lower.tri(ARI)]), 
                         mean(ARI_CS[lower.tri(ARI_CS)]), 
                         mean(ARI_CCS[lower.tri(ARI_CCS)])), 
           conductance= c(mean(ap_cond), 
                          mean(ap_cond_CS), 
                          mean(ap_cond_CCS)), 
           modularity= c(mean(ap_mods), 
                         mean(ap_mods_CS), 
                         mean(ap_mods_CCS)), 
           method=c("AP", "AP-CS", "AP-CCS"))


diff_methods= list(ap_clusters, 
                   ap_clusters_CS, 
                   ap_clusters_CCS)
ARI_diff= matrix(nrow=length(diff_methods), 
                 ncol=length(diff_methods), 
                 dimnames=list(c("AP", "AP-CS", "AP-CCS"), 
                               c("AP", "AP-CS", "AP-CCS")))

for (k in 1:length(diff_methods)){
  for (l in 1:length(diff_methods)){
    ARI_diff[k,l]<- compare(diff_methods[[k]][[1]], 
                           diff_methods[[l]][[1]], "adjusted.rand")
  }
}

for (j in 2:n_trials){
  for (k in 1:length(diff_methods)){
    for (l in 1:length(diff_methods)){
      ARI_diff[k,l]<- ARI_diff[k,l]+compare(diff_methods[[k]][[j]], 
                              diff_methods[[l]][[j]], "adjusted.rand")
    }
  }
  
  
}


ARI_diff/n_trials


do.call(c,lapply(ap_clusters, table) ) |> hist(xlim=c(0,60)
                                               ,ylim=c(0,200),
                                               main="Cluster sizes of AP, AP-CS and AP-CCS \n across repeats")
do.call(c,lapply(ap_clusters_CS, table) ) |> hist(col=rgb(1,0,0,0.5),add=TRUE)
do.call(c,lapply(ap_clusters_CCS, table) ) |> hist(col=rgb(0,1,0,0.5),add=TRUE)
legend("topright", 
       pch=16, 
       col=c("grey","red","green"), 
       legend=c("AP","AP-CS","AP-CSS"))

