library(dbscan)

library(cliquePartitioneR)


library(FCPS)

n_neigh=15

datasets<-c("TwoDiamonds", "Chainlink", "Atom", "Target", "Lsun3D")

for (ds_name in datasets){

bunch<-get(ds_name)
X<-bunch$Data
y<-bunch$Cls

D<- dist(X, method = "euclidean" ) |> as.matrix()

S<- max(D) + 1e-3 - D



source("thresholding_partitioning.R")


critical_mst_thr(S)->ts

tsp<- as.numeric(ts)

S[ S < tsp ]= 0

critical_mst_thr(S)->ts

S[ (S < ts) & (S >= tsp) ]= 0.5*S[ (S < ts) & (S >= tsp) ]


library(igraph)

P<- greedyCliquePartitioner(W=S, unique_singletonLabels = TRUE)$membership

source("CS_dendrogram.R")

CS_ex<-expand_CS_toNodes(cliqueSimilarity(P,S), P)

hclust(as.dist(max(CS_ex) - CS_ex), method = "single") -> hc

dbscan::extractFOSC(hc, minPts = 15, prune_unstable=TRUE)-> fin_cl

VCS_ex <- CS_ex[ order(fin_cl$cluster, P),
                 order(fin_cl$cluster, P)
                 ]

orig_colors= colors()[ bunch$Cls[order(fin_cl$cluster, P)] ]
pred_colors = orig_colors
pred_colors[ fin_cl$cluster==0 ] = "gray"

fcl_ordered<- fin_cl$cluster[order(fin_cl$cluster, P)]

pred_colors[pred_colors!="gray"]<- colors()[fcl_ordered]
pred_cliques_colors<- colors() [ P[order(fin_cl$cluster, P)] ]
#heatmap(VCS_ex, Colv=NA, Rowv=NA,
#        scale="none", 
#        ColSideColors = orig_colors,
#        RowSideColors = pred_colors,
#          )
par(mfrow=c(1,2))
plot(X, col = fin_cl$cluster, pch = 19, cex = 0.5,
     main = "Two Diamonds Dataset with Clustering Results")

graph_from_adjacency_matrix(CS_ex, mode = "undirected",
                            weighted=TRUE, diag=FALSE) -> g

plot(g, vertex.label=NA, vertex.size=6, vertex.color=bunch$Cls)

#plot(X,col = P, pch = 19, cex = 0.5,
#     main = "Two Diamonds Dataset with Clustering Results")

}

