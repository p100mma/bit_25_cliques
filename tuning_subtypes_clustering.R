library(data.table)
fread("Brain_GSE50161.csv") |> as.data.frame()-> microarray_data
head(microarray_data)[,1:6]
X<- microarray_data[,3:ncol(microarray_data)]
microarray_data$type -> y_char
y<- as.numeric(factor(y_char))
library(matrixStats)
colMedsX<- colMedians(X|> as.matrix())
plot(sort(colMedsX))
abline(h=5)
X<- X[, colMedsX > 5 ]
colVarsX<- colVars(X |> as.matrix())
plot(sort(colVarsX))
abline(h=1.5)
X<- X[, colVarsX >1.5]
X<- X[, order(-colVars(X|> as.matrix()))]
X<-X[, !grepl("_(x|a|g|r)_at$", colnames(X))]
X<-X[, !grepl("^AFFX-", colnames(X))]

PCA=prcomp(X, scale.=TRUE)
plot(PCA$sdev)
abline(v=20)
k=20
PC= PCA$x[,1:k]

library(Rdimtools)
library(umap)

k_2nn<-Rdimtools::est.twonn(PC)
nn_=seq(from=3,to=ncol(PC)-1,by=1)
k_nXmethod<-array(NA, dim=c(length(nn_),3))
set.seed(1244)
for (n in seq_along(nn_)){
  message(n)
  k_nXmethod[n,1]=Rdimtools::est.mle2(PC, k1=nn_[[n]], k2= nn_[[n]]+10)$estdim
  k_nXmethod[n,2]=Rdimtools::est.danco(PC,k=nn_[[n]])$estdim
  k_nXmethod[n,3]=Rdimtools::est.nearneighbor1(PC,K= nn_[[n]] )$estdim
}
range(c(k_2nn,as.vector(k_nXmethod)))-> ylims
ylims[[2]]= 6
plot(nn_, k_nXmethod[,1], type="line", ylim=ylims, xlab="# nearest neighbors")
lines(nn_,k_nXmethod[,2],col="red")
lines(nn_,k_nXmethod[,3],col="blue")
legend("topright", lwd=2, col=c("black","red","blue"), legend=c("mle","danco","NN") )

k_=c(c(5,5,5),rep(4, 2))
nn_=c(3:5, 10,15)
stopifnot(length(k_)==length(nn_))
n_rep=30
min_dist_ = c(0.01, 0.1, 0.5)
d_rho<-d_sp<-array(NA, dim=c(length(nn_), length(min_dist_),
                             n_rep))
seeds<-sample.int(32423234,size = n_rep,replace = FALSE)
D_x<- as.matrix(dist(X))
i=0
for (n in seq_along(nn_))
  for (r in seq_along(seeds))
  for (m in seq_along(min_dist_)){ i=i+1
    message(sprintf("%d/%d", i, length(nn_)*n_rep*length(min_dist_)))
    UM<- umap(PC, n_neighbors = nn_[[n]],  n_components = k_[[n]] , min_dist=min_dist_[[m]],
     random_state= seeds[[r]] )$layout
    D_u<-as.matrix(dist(UM))
    d_rho[n,m,r]=cor(D_u[lower.tri(D_u)], D_x[lower.tri(D_x)])
    d_sp[n,m,r]=cor(D_u[lower.tri(D_u)], D_x[lower.tri(D_x)], method="spearman")
  }

flat_rho<- melt(d_rho)
head(flat_rho)
colnames(flat_rho)<- c("kXn_neigbors","min_dist","rep","value")
flat_sp<-melt(d_sp)
colnames(flat_sp)<- c("kXn_neigbors","min_dist","rep","value")
flat_rho$type="pearson"
flat_sp$type="spearman"
flat_RHO= rbind(flat_rho,flat_sp)
flat_RHO$min_dist<- min_dist_[ flat_RHO$min_dist ] 
flat_RHO$kXn_neigbors<- paste0( nn_,"x",k_)[ flat_RHO$kXn_neigbors ]
library(ggplot2)
ggplot(flat_RHO, aes(x= kXn_neigbors, y= value, 
                     fill=as.factor(min_dist) )) + geom_boxplot() + facet_wrap(type~.)
flat_RHO$type

UM<-umap(PC, seed=1234, min_dist=0.5, n_neighbors=10, n_components=4)
UM_UM<- umap(UM$layout, seed=1234, min_dist=0.5, n_neighbors=10,n_components=2)

plot(UM$layout, pch=16, col= y )
plot(UM_UM$layout, pch=16, col=y)

library(dbscan)

min_cl= seq(from=2, to=floor(nrow(X)/2),length.out=30)

source("clique_utils.R")
source("CS_dendrogram.R")


kNNs<- function(D,mcs){
  apply(D,1,function(x) order(rank(-x),decreasing=TRUE)[1:mcs]  )
}

kNN_consistency<- function( knns){
 
  mutual<-0
  total<-length(as.vector(knns))
  for (i in 1:ncol(knns))
    for (j in knns[,i])
      mutual<-mutual + (j %in% knns[,i])
  mutual/total
}

SND<- function(D, mcs){
  neighborhoods<-kNNs(D,mcs)
  resmat<-D
  for (x in 1:ncol(neighborhoods))
    for (y in 1:ncol(neighborhoods)) 
      resmat[x,y]= sum( neighborhoods[,x]==
                          neighborhoods[,y])
  max(resmat) - resmat
}

mst_thresholds<-function(mst_graph){
  mst_D<- mst_graph
  sort( igraph::E(mst_D)$weight, decreasing=TRUE)-> dMSTsorted
  dMSTsorted<- dMSTsorted[ dMSTsorted!=0 ]
  ddiffs= -1*(diff(dMSTsorted))
  dmeans= (dMSTsorted[ 2: length(dMSTsorted) ] + dMSTsorted[ 1: (length(dMSTsorted)-1) ])/2
  dmeans[order(ddiffs, decreasing = TRUE) ]-> d_t
  list(thrs=d_t, diffs= ddiffs[ order(ddiffs, decreasing = TRUE) ] )
  
}


LDT<- function(D, min_cluster_size, alpha){
  hclust(as.dist(D), method = "complete") -> hcl
  graph_from_adjacency_matrix(D, mode = "undirected", weighted = TRUE, diag=FALSE) -> G_D
  D_mst<- mst(G_D)
  mst_thr_data<-mst_thresholds(D_mst)
  d_t_vector<- mst_thr_data$thrs
  dd_vector<- mst_thr_data$diffs
  if (all(d_t_vector==d_t_vector[[1]])) d_t=0 else d_t=d_t_vector[[1]]
  cutree(hcl, h = d_t) -> P
  D_t<-D
  D_t[D > d_t] = 0
  cliqueSimilarity(P, (D_t>0)*1. )-> clique_simil
  expand_CS_toNodes(clique_simil, P) -> clique_simil_node2node
  hcl_clq<- hclust(as.dist(max(clique_simil_node2node)- clique_simil_node2node), method = "single")
  clique_simil_node2node[
    lower.tri(clique_simil_node2node)]-> cs_ltr
  iter=1
  cs_ltr[cs_ltr==iter] = 1
  cs_ltr[ !(cs_ltr %in% c(1,-1)) ] =0 
  cs_ltr[(cs_ltr==0) ] = -1
  extractFOSC(hcl_clq, minPts = min_cluster_size, prune_unstable = TRUE,
              alpha = alpha, constraints = cs_ltr 
  )$cluster -> flat_partition
  stopifnot(!any(is.na(flat_partition)))
  # for (uq in unique(P))
  #    node_lvl_labels[ P == uq ] <- flat_partition[ uq ]
  return(list(clq=P,
              cl=flat_partition,
              cs_n2n=clique_simil_node2node,
              dtv=d_t_vector,
              dt=d_t,
              dd= dd_vector[[1]],
              ddv= dd_vector
              ))
  
  
}


mst_thr_pv<- function(thr_diff, D){
  matrix(sample( D[lower.tri(D)],100*2,replace=TRUE),ncol=2)-> random_pairs
  pair_diffs<- abs(rowDiffs(random_pairs))
  sum(thr_diff >= pair_diffs)/length(pair_diffs)
}

LDF<- function(D, min_cluster_size,alpha, to_SNN=FALSE,
               n_nn=7 ){
  D_orig<- D
  if (to_SNN) { 
    mcs_ratio= n_nn/ncol(D)
    D<- SND(D, mcs=n_nn )
  }
  base<-LDT(D,min_cluster_size = min_cluster_size, alpha=alpha)
  d=0
  cs_n2n<- base$cs_n2n
  P= base$clq
  CL=base$cl
  P[base$cl==0]=0
  tree_hierarchy<- list("1,1"=base)
  clqs<- unique(P[P!=0])
  #print(table(clqs))
  splittable<- sapply(clqs, function(x) sum(P==x) > min_cluster_size )
  #print(splittable)
  while (any(splittable)){ d= d+1
 # print(d)
  for (j in seq_along(clqs))
  {
    if (splittable[[j]]){
      clique<- P== clqs[[j]]
      D_sub<- D[clique,clique]
      #if (to_SNN) { D_sub<- SND(D= D_sub, mcs = length(clique)*mcs_ratio) }
      LDT_dj=LDT( D_sub,
                  min_cluster_size = min_cluster_size,
                  alpha=alpha)
    diff_thr=LDT_dj$dd
    PV=mst_thr_pv(thr_diff = diff_thr, D = D_sub)
    #print(PV)
    if (PV < 0.05) {
    CL [ clique ] = max(CL) + LDT_dj$cl
    CL [ clique ][ LDT_dj$cl== 0 ] = 0
    #  candid_cl<- CL
    #  candid_cl[ clique ]= max(CL) + LDT_dj$cl
    #  candid_cl[ clique ][ LDT_dj$cl==0 ] = 0
    #  candid_objective<- MinST_DunnIndex(partition=candid_cl,dS=D_orig,dX.Y="hausdorff")
    #  print("####")
    #  print(base_objective)
    #  print(candid_objective)
    #  print("####")
     # if (candid_objective > base_objective) {
     #   message("accepted")
      #  base_objective<-candid_objective
      #  CL<- candid_cl
        cs_n2n[ clique, clique ] = LDT_dj$cs_n2n +  d
        P[ clique ]= max(P) + LDT_dj$clq
        P[ clique ][ LDT_dj$cl==0 ] = 0
        LDT_dj$pv= PV
        LDT_dj$placement = clique
        tree_hierarchy[[paste0(d+1, clqs[[j]])]]=LDT_dj
      } else {splittable[[j]]=FALSE}
    }
  }
  new_clqs<-unique(P[P!=0])  
  new_splittable<- sapply(new_clqs, function(x) sum(P==x) > min_cluster_size)
  new_splittable[ new_clqs %in% clqs[!splittable] ] = FALSE
  clqs<-new_clqs
  splittable<-new_splittable
  }
  
  extractFOSC(hclust(as.dist(max(cs_n2n) - cs_n2n), method = "single"),
              minPts = min_cluster_size)$cluster-> CL2
  fin_CL<-CL2
  base_objective<- MinST_DunnIndex(partition = fin_CL|> uniqueSingletonLabels(),
                                   dS = D_orig,
                                   dX.Y = "hausdorff")
  for (clique_tree in tree_hierarchy){
    candidatus<- fin_CL 
    candidatus[ clique_tree$placement ]= CL[ clique_tree$placement ] 
    candidate_score<- MinST_DunnIndex(partition = candidatus |> uniqueSingletonLabels(),
                                      dS = D_orig,
                                      dX.Y = "hausdorff")
    if (candidate_score>base_objective){
     base_objective<- candidate_score
     fin_CL<-candidatus
    }
  }
  
  list(cl=fin_CL,
       cs_n2n=cs_n2n,
       cl_fine=CL,
       cl_coarse=CL2,
       clq=P,
       clique_forest=tree_hierarchy)
}





X_<- list(original=X, PC=PC)


orig_objective<-vector()
orig_ARI<-vector()
PC_objective<-vector()
PC_ARI<-vector()
for (N in 3:15){
cl_x<-((X_$original) |> dist() |> as.matrix() |> LDF(min_cluster_size = min_cl[[7]],alpha=1,
                                            to_SNN = TRUE, n_nn= N))$cl
cl_pc<-((X_$PC) |> dist() |> as.matrix() |> LDF(min_cluster_size = min_cl[[7]],alpha=1,
                                            to_SNN = TRUE, n_nn= N))$cl
PC_objective[[N-2]]<- MinST_DunnIndex(cl_pc|> uniqueSingletonLabels(), (X_$PC) |> dist() |> as.matrix(),
                                      dX.Y =  "hausdorff")
orig_objective[[N-2]]<- MinST_DunnIndex(cl_x|> uniqueSingletonLabels(), (X_$original) |> dist() |> as.matrix(),
                                      dX.Y =  "hausdorff")
PC_ARI[[N-2]]<- ARI(cl_pc, y)
orig_ARI[[N-2]]<- ARI(cl_x, y)

}

par(mfrow=c(2,1))
plot(PC_objective, ylim=range(PC_objective,PC_ARI),type="line")
lines(PC_ARI, col="red")
plot(orig_objective, ylim=range(orig_objective, orig_ARI),type="line")
lines(orig_ARI,col="red")

which.max(PC_objective);which.max(PC_ARI)

library(igraph)
library(aricode)

ari_methods<- array(NA, dim=c(length(min_cl),
                              2,
                              3),
                    dimnames=list(seq_along(min_cl),
                                  c("original space","PCA"),
                                  c("hdbscan","ldf","ldf-snn")))

hdbs<-ldfs<-lapply(seq_along(min_cl), function(x) 
  list(original=list(), PC=list(), UMAP= list()))
MCS=6
for (m in seq_along(min_cl))
  for (e in 1:2)
  {
    message(m)
    hdbs[[c(m,e)]]<-hdb<- hdbscan(X_[[e]], minPts = min_cl[[m]])
    sclnn<- LDF(D = SND(
      as.matrix(dist(X_[[e]])),
      mcs=MCS),
      to_SNN = FALSE,
      min_cluster_size = min_cl[[m]],alpha =1
      )
    ldfs[[c(m,e)]]<-scl<- LDF(D =
                                  as.matrix(dist(X_[[e]])),
                                 
                              min_cluster_size = min_cl[[m]],alpha =1,
                              to_SNN=FALSE)
    ari_methods[m,e,"hdbscan"]<-ARI(hdb$cluster|>uniqueSingletonLabels(),y)
    ari_methods[m,e,"ldf"]<- ARI(scl$cl|>uniqueSingletonLabels(),y)
    ari_methods[m,e,"ldf-snn"]<- ARI(sclnn$cl|>uniqueSingletonLabels(),y)
  }

melt(ari_methods)-> ari_methods_flat
#head(ari_methods_flat)
colnames(ari_methods_flat)<-c("min_cluster_size","space","algorithm","ARI")
pdf("PCA_orig.pdf", width=8, height=4)

ggplot(ari_methods_flat[,], aes(x=min_cluster_size, y=ARI, color=algorithm,
                                lty=algorithm)) + 
  geom_line(lwd=1) + 
  facet_wrap(space~., nrow=1) + theme_minimal() 
dev.off()
n_rep=30
ari_umap<-array(NA, dim=c(length(min_cl),
                          3,
                          n_rep),
                dimnames=list(seq_along(min_cl),
                              c("hdbscan","ldf","ldf-snn"),
                              seq_len(n_rep)
                              ))
set.seed(214)
seeds_u<- sample.int(3252344, size=n_rep)
MCS=6
for (m in seq_along(min_cl))
  for (r in seq_len(n_rep))
  {
    message(sprintf("%d %d",r,m))
   
    UM_r<-umap(PC, 
               n_neighbors = 15,  
               n_components = 4 , 
               min_dist=0.5, random_state=seeds_u[[r]])$layout
    ldf<-LDF(D =
               as.matrix(dist(UM_r)),
             min_cluster_size = min_cl[[m]],
             alpha =1, to_SNN = FALSE)
    ldf_snn<-LDF(D= SND(D = as.matrix(dist(UM_r)), mcs = MCS),
             min_cluster_size = min_cl[[m]],
             alpha=1, to_SNN=FALSE)
    
    ari_umap[m,"ldf",r]= ARI(ldf$cl, y)
    ari_umap[m,"ldf-snn",r]= ARI(ldf_snn$cl, y)
    hdb<- hdbscan(UM_r, minPts = min_cl[[m]])
    ari_umap[m,"hdbscan",r]= ARI(hdb$cluster,y)
    print(table(ldf$cl_fine,y))
    print(ARI(ldf$cl_fine, y))
    print(table(ldf_snn$cl_fine))
    print(ARI(ldf_snn$cl_fine,y))
    print(table(hdb$cluster))
    print(ARI(hdb$cluster,y))
    
  }

reshape2::melt(ari_umap[,,])->flat_ari_umap
colnames(flat_ari_umap)<- c("minClusterSize","algorithm",
            "rep",
                            "ARI")


library(ggplot2)

#ggplot(flat_ari_umap, aes(x = minClusterSize, y = ARI, color = algorithm,
#                          fill = algorithm,lty=algorithm)) +
#  geom_line(lwd=2)
pdf("umap.pdf",width=5,height=4)
ggplot(flat_ari_umap, aes(x = minClusterSize, y = ARI, color = algorithm, fill = algorithm)) +
  # Show spread (e.g., interquartile range) as a ribbon
  stat_summary(
    fun.min = function(y) quantile(y, 0.25),
    fun.max = function(y) quantile(y, 0.75),
    geom = "ribbon",
    alpha = 0.2,
    color = NA
  ) +
  # Median line
  stat_summary(
    fun = median,
    geom = "line",
    size = 1
  ) +
  # Optional: median points
  stat_summary(
    fun = median,
    geom = "point",
    size = 2
  ) +
  labs(
    x = "minClusterSize",
    y = "Adjusted Rand Index (ARI)",
    title = "UMAP embedding results",
    color = "Algorithm",
    fill = "Algorithm"
  ) +
  theme_minimal()
dev.off()


nndist_heatmap<- function(dist_mat, k_cap=NULL){
 col2use<- if (is.null(k_cap)) ncol(dist_mat)-1 else k_cap
  do.call(rbind,lapply(1:nrow(D), 
                       function(i) sort(D[i,-i])))[,1:col2use]-> d_arr
  colnames(d_arr)<-1:ncol(d_arr)
  d_arr<- d_arr/rowSums(d_arr)
 # d_arr<- log(d_arr + 1e-30)
  #d_arr<- d_arr[,order(colSums(d_arr))]
  d_arr<- d_arr[order(rowSums(d_arr)),]
  rownames(d_arr)<-1:nrow(d_arr)
  reshape2::melt(d_arr)-> d_flat
  #d_flat
  plot(colMeans(d_arr),log="x")
  #ggplot(d_flat, aes(x=Var1,y=Var2,fill=value+1e-10))+geom_tile()+
  #  scale_fill_gradient(low="orange",high="black", trans="sqrt")
  
  
}



library(FCPS)
examples<-c( "Atom", "TwoDiamonds", "Target", "Chainlink")
best_per_ex_clq<-list()
best_per_ex_hdb<-list()
#for (ex in examples){
#  bunch<- get(ex)
#  X<- bunch$Data
#  y<- bunch$Cls
  #plot(X, col=y, main=ex)
  #for (n in 1:30)
#}
  
#nns<- c(30,30,15,20)
#names(nns)<-examples
clq_per_ex<-list()
pdf("synthetics.pdf",width=5,height=4.3)
for (ex in examples){
  
  bunch<- get(ex)
  X<- bunch$Data
  y<- bunch$Cls
  min_cl_sizes= floor(seq(from=3, to= floor(length(y)/2), 
                          length.out=30))
  D<- as.matrix(dist(X, method = "euclidean"))
  lapply(min_cl_sizes, function(mcl) 
    LDF(D,min_cluster_size = mcl,
        alpha=1, to_SNN = FALSE)
  )-> clq_cl
  clq_per_ex[[ex]]<-clq_cl
  ARI_cl<- sapply(clq_cl, function(x) {
    ARI(y, x$cl|> uniqueSingletonLabels())
  })
  
  lapply(min_cl_sizes, function(mcl)
    hdbscan(X, minPts = mcl))-> hdb_cl
  ARI_hdb<- sapply(hdb_cl, function(x) {
    ARI(y, x$cluster|> uniqueSingletonLabels())
  })
  
  best_clq<- clq_cl[[ which.max(ARI_cl)]]
  best_hdb<- hdb_cl[[ which.max(ARI_hdb)]]
  
  par(mfrow=c(1,1))
  plot(min_cl_sizes, ARI_cl, type="l", 
       xlab="Minimum cluster size", ylab="ARI",
       main=paste0("ARI scores \n on ", ex),
       ylim=c(0,1.1))
  lines(min_cl_sizes, ARI_hdb, col="red", 
        xlab="Minimum cluster size", ylab="ARI",
  )
  legend("right", legend = c("LDF","HDBSCAN*"),
         lwd=2,col=c("black","red"), bty = "n")
  best_per_ex_hdb[[ex]]<- best_hdb
  best_per_ex_clq[[ex]]<- best_clq
}


#par(mfcol=c(2,4))

for (ex in examples) {
  best_clq<- best_per_ex_clq[[ex]]
  best_hdb<- best_per_ex_hdb[[ex]]
  bunch<- get(ex)
  X<- bunch$Data
  plot(X, col=best_clq$cl, 
       main=paste0("best Clique-based clustering \n on ", ex)
  )
  
  plot(X, col=best_hdb$cluster,
       main=paste0("best HDBSCAN* clustering \n on ", ex)
  )
}
dev.off()
