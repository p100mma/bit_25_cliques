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
source("corr_utils.R")
source("clique_utils.R")
source("CCS_measure.R")
source("corr_utils.R")
library(cliquePartitioneR)
library(nloptr)
B<-900
set.seed(1324)
BS<- lapply(1:B, function(b) sample(1:nrow(X), nrow(X), replace=TRUE))
Pg_<- list()
Ph_<- list()
t_<-vector()
start_t<-Sys.time()
D<-1 - cor(X )^2
S<- corsq(X|> as.matrix())
hcl<- hclust(as.dist(D), method="complete")
h_<-bobyqa_optimize_hcl(hcl = hcl,I = range(hcl$height),xtol_abs = min(abs(diff(hcl$height))))
t_<- 1 - h_
S_t<- 1- D
S_t[ S_t < t_ ]=0
Ph<- cutree(hcl, h=h_)
Pg<-greedyCliquePartitioner(S_t,expansion_mode = "average", unique_singletonLabels = TRUE)
print(Sys.time() - start_t)
print("threshold from p.value:") 
print(sqrt(min(S[S>0])))
print("threshold from optimisation:")
print(sqrt(min(S_t[S_t > 0])))
for (b in seq_along(BS)){
  if (((b-1)%%20)==0 ) message(b)
start_t<- Sys.time()
D<-1 - corfast(X[ BS[[b]],] |> as.matrix() )^2
hcl<- hclust(as.dist(D), method="complete")
h_opt<-bobyqa_optimize_hcl(hcl = hcl,I = range(hcl$height),xtol_abs = min(abs(diff(hcl$height))))
t_[[b]]<-thr<- 1 - h_opt
S_t<- 1- D
S_t[ S_t < thr ]=0
Ph_[[b]]<- cutree(hcl, h=h_opt)
Pg_[[b]]<-greedyCliquePartitioner(S_t,expansion_mode = "average", unique_singletonLabels = TRUE)
if (((b-1)%%20)==0 ) message(Sys.time() - start_t)
}

Pg<- Pg$membership

Pg_<- lapply(Pg_, function(x) x$membership)

library(aricode)

ARI_h<- sapply(Ph_, function(ph) ARI(ph,Ph))
ARI_g<- sapply(Pg_, function(ph) ARI(ph,Pg))

saveRDS(Pg_,"Pg_.rds")
saveRDS(Ph_,"Ph_.rds")
saveRDS(t_, "t_.rds")
n_rep<-30
library(dynamicTreeCut)
tom<- function(A) {diag(A)=0; deg<-colSums(A);  (A + A%*%A)/( outer(deg,deg,pmin) +1 - A )  }
DTC_groups<- function(sim, mcs=20){
 diag(sim)<-1
 hcl<-hclust(as.dist(1-sim),method="average")
 hcl$height=round(hcl$height,10)
 cutreeDynamic(hcl, minClusterSize = mcs, distM = 1-sim, verbose=0) |> uniqueSingletonLabels()
}
base_S<- corfast(X|>as.matrix())^2
base_dtc<- DTC_groups(base_S,mcs=100)
heatmap_clusters<- function(mat, labs, maxSize=500, zero_thr=1,main=NULL){
  zeroOutSmall(labs,zero_thr+1)-> labs
  subs<- 1:ncol(mat)
  if (length(labs)>maxSize) subs<-sample(1:ncol(mat),maxSize, replace=FALSE)
  mat<- mat[subs,subs]
  labs<-labs[subs]
  mat<- mat[order(labs), order(labs)]
  labs<-labs[order(labs)]
  rainbow(length(unique(labs[labs!=0])))-> clrs
  #print(clrs)
  side_cols<- rep("gray", length(labs)); k=1
  for (uql in unique(labs[labs!=0])){ side_cols[labs==uql]= clrs[[k]];k=k+1}
  heatmap(mat, Rowv=NA,Colv=NA,scale="none",
          ColSideColors=side_cols,
          RowSideColors=side_cols,
          main= main)
}

dtc_par<-list()
for (param in seq(from=10, to=200, by=10))
  dtc_par[[length(dtc_par)+1]]= DTC_groups(base_S,mcs = param)

ARI_mat<- function(cl_list){
  ARImat<-matrix(nrow=length(cl_list),ncol=length(cl_list))
  for (i in seq_along(cl_list))
    for (j in seq_along(cl_list))
      ARImat[i,j]= ARI(cl_list[[i]],cl_list[[j]])
  ARImat
}

library(igraph)
modularity_fromMatrix<- function(W, members){
  G_W<- graph_from_adjacency_matrix(adjmatrix=W, mode="undirected",
                                    weighted=TRUE,diag=FALSE)
  modularity(G_W,membership=members, weights=E(G_W)$weight)
}

# conductance - abstract cluster quanlity function
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


arm<-ARI_mat(dtc_par)
mod_par<- sapply(dtc_par, function(x) modularity_fromMatrix(W = base_S,members = x))
cn_par<-sapply(dtc_par, function(x) conductance_network(A=base_S,members=x) )

library(ggplot2)
library(reshape2)
ggheat<- function(array,line_ys,main=NULL){
  melt(array)-> flat_array
  colors<-rainbow(length(line_ys))
  res<-ggplot(flat_array, aes(x=Var1,y=Var2,fill=value))+geom_tile()+ 
    theme_minimal()#+xlim(min(flat_array$Var1)-0.2*min(flat_array$Var1), max(flat_array$Var1)*1.2 )
  for (k in seq_along(line_ys)){
    Yk<-line_ys[[k]]
    if ((min(Yk)<0)||(max(Yk)>1)) Yk=  (Yk - min(Yk))/(max(Yk)-min(Yk))
  res<-res+geom_line(data=data.frame(x= flat_array$Var1,
                                       y= Yk
                                       ),
                        aes(x=x,
                            y= y*max(flat_array$Var2) ),
              color=colors[[k]],
              lwd=2,
              inherit.aes = FALSE)
  }
  for (k in seq_along(line_ys)){
    lab_x<- sort(unique(flat_array$Var1), decreasing = TRUE)[k:(k+1)] |> mean() 
    print(lab_x)
    res<- res+ annotate("text",x = lab_x, y= min(flat_array$Var2), 
                        label= round(min(line_ys[[k]]),2), color=colors[[k]] )
    res<- res+ annotate("text",x = lab_x, y= max(flat_array$Var2), 
                        label= round(max(line_ys[[k]]),2), color=colors[[k]] )
  }
  if(!is.null(main)) res<- res + ggtitle(main)
  res
}





ggheat(arm, list(mod_par, cn_par))
heatmap_clusters(base_S,dtc_par[[16]])
base_dtc<- dtc_par[[16]]
chosen_mcs_base<- seq(from=10, to=200, by=10)[[16]]
dtc_g<-list()
dtc_h<-list()
tc_cs<-list()
dtc_tom<-list()
CCSg<-matrix(0,nrow=ncol(X),ncol=ncol(X))
CCSh<-matrix(0,nrow=ncol(X),ncol=ncol(X))
CS<- matrix(0, nrow=ncol(X), ncol=ncol(X))
for (b in 1:n_rep){
  message(sprintf("%d/%d",b,n_rep))
  absCor<- abs(corfast(X[ BS[[b]],]|>as.matrix() ))
  S_t<- absCor^2
  TOM<- tom(absCor^2)
  S_t[ S_t < t_[[b]] ] = 0
  csg<- cliqueSimilarity(Pg_[[b]], S_t)
  csh<- cliqueSimilarity(Ph_[[b]], S_t)
  p_dtc<- DTC_groups(absCor^2, mcs = chosen_mcs)
  CS= CS + outer(p_dtc,p_dtc, function(x,y) as.numeric(x==y))
  ecsg<- expand_CS_toNodes(csg, cl_mem=Pg_[[b]])
  #ecsh<- expand_CS_toNodes(csh, cl_mem=Pg_[[b]])
  CCSg<- CCSg + ecsg
 # CCSh<-CCSh + ecsh
  
}
CCSg<- CCSg/n_rep
#CCSh<-CCSh/n_rep
CS<- CS/n_rep

par(mfrow=c(2,2))
hist((CS)[lower.tri(CS)])
hist(TOM[lower.tri(TOM)])
hist(CCSg[lower.tri(CCSg)])
#hist(CCSh[lower.tri(CCSh)])


mcs_range<-seq(10,200,10)
tom_mcs<-cs_mcs<-ccs_mcs<-list()
for (m in seq_along(mcs_range)){i
  tom_mcs[[m]]<-DTC_groups(sim=TOM, mcs = mcs_range[[m]])
  cs_mcs[[m]]<-DTC_groups(sim=CS, mcs = mcs_range[[m]])
  ccs_mcs[[m]]<- DTC_groups(sim=CCSg, mcs=mcs_range[[m]])
}

simXmcs<-list(TO=tom_mcs, CS=cs_mcs, CCS=ccs_mcs)


mcs_metrics_df<- data.frame(mcs=mcs_range)
for (metric in names(simXmcs)){
  #mcs_metrics_df[[paste0(metric,"+","MOD")]]=sapply(simXmcs[[metric]], function(x) 
  #                                                  modularity_fromMatrix(base_S,x))
  #mcs_metrics_df[[paste0(metric,"+","COND")]]=sapply(simXmcs[[metric]], function(x) 
  #                                                    conductance_network(base_S,x))
  mcs_metrics_df[[paste0(metric,"+","cl_frac")]]=sapply(simXmcs[[metric]], function(x)
                                                          sum(zeroOutSmall(x,2)>0)/length(x))
}

ari_mats<-list()
for (metric in names(simXmcs)){
  ari_mats[[metric]]= ARI_mat(simXmcs[[metric]])
                              }
pdf("tuning_params.pdf")
for (metric in names(simXmcs)){
  sub_df<- mcs_metrics_df[, grep(paste0("^",metric),colnames(mcs_metrics_df)) ]
  print(ggheat(array=ari_mats[[metric]], line_ys = sub_df, main=metric))
}
dev.off()
mcs_values<- c(mcs_range[[7]], mcs_range[[16]])


heatmap_clusters(mat=base_S, labs = tom_mcs[[7]],main=paste0("TOM",mcs_values[[1]]))
heatmap_clusters(mat=base_S, labs = cs_mcs[[7]],main=paste0("CS",mcs_values[[1]]))
heatmap_clusters(mat=base_S, labs = ccs_mcs[[7]],main=paste0("CCS",mcs_values[[1]]))
heatmap_clusters(mat=base_S, labs = tom_mcs[[16]],main=paste0("TOM",mcs_values[[2]]))
heatmap_clusters(mat=base_S, labs = cs_mcs[[16]],main=paste0("CS",mcs_values[[2]]))
heatmap_clusters(mat=base_S, labs = ccs_mcs[[16]],main=paste0("CCS",mcs_values[[2]]))





dtc_g[[length(dtc_g)+1]]<- DTC_groups(sim = CCSg, mcs=chosen_mcs)
dtc_h[[length(dtc_h)+1]]<- DTC_groups(sim = CCSh,mcs=chosen_mcs)
dtc_cs[[length(dtc_cs)+1]]<- DTC_groups(sim= CS,mcs=chosen_mcs)

mcsXsimXgene<- array(NA, dim = c(length(mcs_values), length(simXmcs),
                                 B/n_rep,ncol(X) 
                                 ),
                     dimnames= list( mcs_values|> as.character(), 
                                     names(simXmcs),
                          
                                     seq_len(B/n_rep),
                                     colnames(X)
                                     )
                     )

mcsXsimXgene["70","TO",1,]= simXmcs$TO[[7]]
mcsXsimXgene["160","TO",1,]=simXmcs$TO[[16]]
mcsXsimXgene["70","CS",1,]= simXmcs$CS[[7]]
mcsXsimXgene["160","CS",1,]=simXmcs$CS[[16]]
mcsXsimXgene["70","CCS",1,]= simXmcs$CCS[[7]]
mcsXsimXgene["160","CCS",1,]=simXmcs$CCS[[16]]

  


for (r in 2:(B/n_rep)){
  for (m in seq_along(mcs_values)){
    mcs= mcs_values[[m]]
  CCSg<-matrix(0,nrow=ncol(X),ncol=ncol(X))
 # CCSh<-matrix(0,nrow=ncol(X),ncol=ncol(X))
 # CS<- matrix(0, nrow=ncol(X), ncol=ncol(X))
  b_start<- 1 + ( r- 1)*n_rep
  message(b_start)
  start_t<- Sys.time()
 for  (b in c(b_start: (b_start+ (n_rep-1)  ))){
   #message("############")
   absCor<- abs(corfast(X[ BS[[b]],]|>as.matrix() ))
   S_t<- absCor^2
   S_t[ S_t < t_[[b]] ] = 0
   csg<- cliqueSimilarity(Pg_[[b]], S_t)
   #csh<- cliqueSimilarity(Ph_[[b]], S_t)
  # p_dtc<- DTC_groups(absCor^2, mcs = mcs)
  # CS= CS + outer(p_dtc,p_dtc, function(x,y) as.numeric(x==y))
   #message(Sys.time()-start_t)
   ecsg<- expand_CS_toNodes(csg, cl_mem=Pg_[[b]])
   #ecsh<- expand_CS_toNodes(csh, cl_mem=Pg_[[b]])
   #message(Sys.time()-start_t)
   CCSg<- CCSg + ecsg
   #CCSh<-CCSh + ecsh
   #message(Sys.time()-start_t)
   #if (b==b_start) { TOM<- tom(absCor^2)
  #   mcsXsimXgene[ m, "TO", r,]= DTC_groups(TOM, mcs = mcs)
  # }
 }
  CCSg<- CCSg/n_rep
  #CCSh<-CCSh/n_rep
  #CS<- CS/n_rep
  #mcsXsimXgene[m,"CS",r,]= DTC_groups(sim=CS, mcs=mcs)
  mcsXsimXgene[m,"CCS",r,]=DTC_groups(sim=CCSg, mcs=mcs)
  message(Sys.time()-start_t)
  }
  }



save.image(file="after_mcsXsimXgene.RData")
load("after_mcsXsimXgene.RData")

library(apcluster)
ap_grouping<- function(sim, q, seed){
  apcluster(s=sim,seed=seed, p=q) |> labels(type="enum")
}

AP_simXgene<- array(NA, dim=c(3,B/n_rep,ncol(X)),
                    dimnames = list(
                     c("TO","CS","CCS"),
                     seq_len(B/n_rep),
                     colnames(X)
                    )
                      )

library(matrixStats)
for (r in 1:(B/n_rep)){
    CCSg<-matrix(0,nrow=ncol(X),ncol=ncol(X))
    CS<- matrix(0, nrow=ncol(X), ncol=ncol(X))
    b_start<- 1 + ( r- 1)*n_rep
    message(b_start)
    start_t<- Sys.time()
    for  (b in c(b_start: (b_start+ (n_rep-1)  ))){
      gc()
      absCor<- abs(corfast(X[ BS[[b]],]|>as.matrix() ))
      if (b==b_start) { TOM<- tom(absCor^2)
         AP_simXgene["TO", r,]= ap_grouping(TOM, q = quantile(TOM,0.1), seed = b)
         rm(TOM)
       }
      S_t<- absCor^2
      p_ap<-ap_grouping(absCor^2, q= quantile(absCor^2,0.1),seed=b) 
      rm(absCor)
      CS= CS + outer(p_ap,p_ap, function(x,y) as.numeric(x==y))
      S_t[ S_t < t_[[b]] ] = 0
      csg<- cliqueSimilarity(Pg_[[b]], S_t)
      rm(S_t)
      ecsg<- expand_CS_toNodes(csg, cl_mem=Pg_[[b]])
      rm(csg)
      CCSg<- CCSg + ecsg
      rm(ecsg)
    }
    CCSg<- CCSg/n_rep
    CS<- CS/n_rep
    AP_simXgene["CS",r,]= ap_grouping(CS, q = quantile(CS,0.1), seed = b)
    mcsXsimXgene["CCS",r,]=ap_grouping(CCSg, q = quantile(CCSg,0.1), seed = b)
    message(Sys.time()-start_t)
  }






library(aricode)
ARI_matrices<- array(NA, dim=c(2,3,B/n_rep,B/n_rep),
                     dimnames = list(c("70","160"),
                                     names(simXmcs),
                                     seq_len(B/n_rep),
                                     seq_len(B/n_rep)))
for (i_mcs in 1:2)
  for (i_sim in 1:3){
    ARI_matrices[i_mcs,i_sim,,]=ARI_mat(cl_list = as.data.frame(t(mcsXsimXgene[i_mcs,i_sim,,])))
    
  }
library(reshape2)
ari_flat<- melt(ARI_matrices)
colnames(ari_flat)<- c("minClusterSize","similarity","i","j","ARI")
ari_flat<- ari_flat[ ari_flat$i<ari_flat$j,] 
library(ggplot2)
pdf("DTC_results.pdf",height=3,width=4)
ggplot(ari_flat, aes(x=similarity, y=ARI, fill=similarity))+geom_boxplot()+
  facet_wrap(minClusterSize~.) + theme_minimal()+
  ggtitle("Stability of DynamicTreeCut")


performance_array<-array(NA, dim=c(2,3,B/n_rep),
                          dimnames=list(c("70","160"),
                                        names(simXmcs),
                                        seq_len(B/n_rep)))

library(igraph)
for (m in 1:2)
  for (sim in names(simXmcs))
      for (r in seq_len(B/n_rep)){
      performance_array[m,sim,r]= modularity_fromMatrix(base_S, mcsXsimXgene[m,sim,r,] )
      message(sprintf("m=%d;sim=%s,r=%d",m,sim,r))
      }

flat_perf<-melt(performance_array)
colnames(flat_perf)<-c("minClusterSize","similarity","rep","modularity")

ggplot(flat_perf, aes(x=similarity, y=modularity, fill=similarity))+ geom_boxplot()+
  facet_wrap(minClusterSize~., scales="free_y") +
  ggtitle("Abstract cluster quality - DTC") +
  theme_minimal()
dev.off()

cn_ccsh<-cn_ccsg<-cn_cs<-cn_tom<-mod_tom<-mod_cs<-mod_ccsg<-mod_ccsh<-vector()
for (b in seq_along(dtc_tom)){
  message(b)
  mod_tom[[b]]<-modularity_fromMatrix(base_S, dtc_tom[[b]]) 
  mod_cs[[b]]<-modularity_fromMatrix(base_S, dtc_cs[[b]]) 
  mod_ccsg[[b]]<-modularity_fromMatrix(base_S, dtc_g[[b]]) 
  mod_ccsh[[b]]<-modularity_fromMatrix(base_S, dtc_h[[b]])
  
  cn_tom[[b]]<-conductance_network(base_S, dtc_tom[[b]]) 
  cn_cs[[b]]<-conductance_network(base_S, dtc_cs[[b]]) 
  cn_ccsg[[b]]<-conductance_network(base_S, dtc_g[[b]]) 
  cn_ccsh[[b]]<-conductance_network(base_S, dtc_h[[b]])
  
  
}

summary(mod_tom)
summary(mod_cs)
summary(mod_ccsg)
summary(mod_ccsh)

summary(cn_tom)
summary(cn_cs)
summary(cn_ccsg)

heatmap_clusters<- function(mat, labs, maxSize=500, zero_thr=1){
  zeroOutSmall(labs,zero_thr+1)-> labs
  subs<- 1:ncol(mat)
  if (length(labs)>maxSize) subs<-sample(1:ncol(mat),maxSize, replace=FALSE)
  mat<- mat[subs,subs]
  labs<-labs[subs]
  mat<- mat[order(labs), order(labs)]
  labs<-labs[order(labs)]
  rainbow(length(unique(labs[labs!=0])))-> clrs
  #print(clrs)
  side_cols<- rep("gray", length(labs)); k=1
  for (uql in unique(labs[labs!=0])){ side_cols[labs==uql]= clrs[[k]];k=k+1}
  heatmap(mat, Rowv=NA,Colv=NA,scale="none",
          ColSideColors=side_cols,
          RowSideColors=side_cols,
          main= sum(labs==0))
}



ap_grouping<- function(sim, seed){
  apcluster(s=sim,seed=seed) |> labels(type="enum")
}


load("AP_results_aggregated.RData")

pdf("AP_results_fulldata.pdf", width=5, height=5)
AP_ARI_flat<-melt(AP_ARI_matrices)
colnames(AP_ARI_flat)<- c("similarity","i","j","ARI")
AP_ARI_flat<- AP_ARI_flat[ AP_ARI_flat$i< AP_ARI_flat$j, ]
ggplot(AP_ARI_flat, aes(x=similarity, y=ARI, fill=similarity))+geom_boxplot()+
  ggtitle("Stability of Affinity Propagation")+
  theme_minimal()

AP_perf_flat<- melt(AP_performance_array)
head(AP_perf_flat)
colnames(AP_perf_flat)<- c("metric","similarity","r","value")
AP_perf_flat<- AP_perf_flat[ AP_perf_flat$metric=="modularity",] 
ggplot(AP_perf_flat, aes(x=similarity, y=value, fill=similarity))+ geom_boxplot()+
  ggtitle("Abstract clustering quality of Affinity Propagation")+
  facet_wrap(metric~., scales = "free_y")+theme_minimal()

dev.off()


lapply(mcsXsimXgene |> dimnames(), head)
lapply(AP_simXgene |> dimnames(), head)

AP_expanded<- array(AP_simXgene, dim=c(1, dim(AP_simXgene)),
                    dimnames = c(list( algorithm= "AP" ), 
                                 dimnames(AP_simXgene) )
)

DTC_array<- mcsXsimXgene
dimnames(DTC_array)[[1]]<- c("DTC70", "DTC160")

algorithms= c("DTC70","DTC160","AP")
similarities<- c("TO","CS","CCS")

modules_array<- array(NA, dim=c(length(algorithms),
                                length(similarities),
                                B/n_rep,
                                ncol(X)),
                      dimnames = list(algorithms,similarities,
                                      seq_len(B/n_rep),
                                      colnames(X)))

modules_array[1:2, , , ] = DTC_array 
modules_array[3, , , ] =  AP_expanded

unique(y_char)-> types

Y<- lapply(types, function(y) 1*(y_char==y)) |> as.data.frame()
colnames(Y)<- types

algorithmXsimXrepXtraitXmodule<- vector(mode="list", length=3)

names(algorithmXsimXrepXtraitXmodule)<- algorithms
for (al in algorithms) 
  algorithmXsimXrepXtraitXmodule[[al]]= vector(mode="list",length=3)

for (al in algorithms) 
  names(algorithmXsimXrepXtraitXmodule[[al]])<- similarities

for (al in algorithms)
  for (sim in similarities)
    algorithmXsimXrepXtraitXmodule[[c(al,sim)]]=vector(mode="list",
                                                       length=n_rep)

library(MDFS)

for (al in algorithms){
  for (sim in similarities)
    for (r in 1:n_rep)
    {
      P<-zeroOutSmall(modules_array[al, sim, r,],2)
      unique(P[P!=0])-> mods
        eigengenes<- lapply(mods, function(mod)
                                prcomp(X[, P==mod],scale.=TRUE,
                                       rank.=1)$x[,1]) |> as.data.frame()
        PV<-MI<- matrix(nrow=ncol(Y),
                    ncol=ncol(eigengenes))
        rownames(PV)<-rownames(MI)<-types
        for (type in types){
         MDFS(data = eigengenes, decision = Y[,type], divisions = 1, dimensions = 1, 
              discretizations = 30, range=0.3, seed=124)-> mdfs
        MI[type,]=mdfs$statistic
        PV[type,]=mdfs$p.value
        }
        PV[,]= p.adjust(PV)
        MI [ PV >= 0.05 ] = 0
        MI<-MI/rowSums(MI)
       algorithmXsimXrepXtraitXmodule[[c(al,sim,r)]]=MI
    }
       message(sprintf(" %s & %s done",al,sim))
}


n_type_matches<- array(NA, dim=c(3,3,n_rep),
                       dimnames = list(algorithms,
                                       similarities,
                                       1:n_rep))

for (al in algorithms)
  for (sim in similarities)
    for (r in 1:n_rep)
    {
      type_pref<-apply(algorithmXsimXrepXtraitXmodule[[c(al,sim,r)]],1, function(ROW)
        which.max(ROW))
      mod_pref<-apply(algorithmXsimXrepXtraitXmodule[[c(al,sim,r)]],2, function(COL)
        which.max(COL))
      type_matches=0
      for (i_t in seq_along(type_pref)){
        type_matches= type_matches+ 1*( mod_pref[[ type_pref[[i_t]] ]] == i_t )
      }
      n_type_matches[al,sim,r]=type_matches
    }
melt(n_type_matches)->n_match_flat
head(n_match_flat)
colnames(n_match_flat)<-c("algorithm","similarity","repeat","n_matches")
n_match_flat$n_matches<- factor(n_match_flat$n_matches, levels=c(0:5))
n_match_flat$algorithmXsimilarity= interaction(n_match_flat$algorithm,
                                               n_match_flat$similarity)
table(n_match_flat$algorithmXsimilarity,
      n_match_flat$n_matches) |> melt()-> n_match_tabulation

colnames(n_match_tabulation)<- c("algorithmXsimilarity","n_matches","freq")
pdf("tumor_module_matchingFreq.pdf", width=8,height=6)
ggplot(n_match_tabulation[(n_match_tabulation$n_matches!=0)#&
                         # grepl("DTC",n_match_tabulation$algorithmXsimilarity
                          #)
                          ,], aes(x=n_matches, y=algorithmXsimilarity,
                               fill=freq) ) + geom_tile() +
  labs(
    x = "# matches",
    y = "similarity",
    title = "perfect matching frequency between \n tumor subtype and a module for DynamicTreeCut"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 1),
    strip.text = element_text(face = "bold")
  ) + geom_text(aes(label=freq),color="white") +
  theme(legend.postition= "none")
dev.off()


