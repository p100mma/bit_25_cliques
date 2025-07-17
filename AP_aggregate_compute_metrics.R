library(data.table)
fread("Brain_GSE50161.csv") |> as.data.frame()-> microarray_data
X<- microarray_data[,3:ncol(microarray_data)]
microarray_data$type -> y_char
y<- as.numeric(factor(y_char))
library(matrixStats)
colMedsX<- colMedians(X|> as.matrix())
X<- X[, colMedsX > 5 ]
colVarsX<- colVars(X |> as.matrix())
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
base_S<- corfast(X|>as.matrix())^2
library(aricode)

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

n_rep<-30


AP_simXgene<- array(NA, dim=c(3,B/n_rep,ncol(X)),
                    dimnames = list(
                     c("TO","CS","CCS"),
                     seq_len(B/n_rep),
                     colnames(X)
                    )
                      )

for (r in 1:(B/n_rep)){
	AP_r<- readRDS(sprintf("AP_rep_%d.rds",r))
	for (sim in c("TO","CS","CCS"))
		AP_simXgene[sim,r,]=AP_r[[sim]]
  }

#stability estimation

AP_ARI_matrices<- array(NA, dim=c(3,B/n_rep,B/n_rep),
                     dimnames = list(
                     		    c("TO","CS","CCS"),
                                     seq_len(B/n_rep),
                                     seq_len(B/n_rep)))
for (sim in c("TO","CS","CCS"))
    AP_ARI_matrices[sim,,]=ARI_mat(cl_list = as.data.frame(t(AP_simXgene[sim,,] )))

#performance metrics - abstract cluster quality


AP_performance_array<-array(NA, dim=c(2,3,B/n_rep),
                          dimnames=list(c("modularity","conductance"),
                     		    c("TO","CS","CCS"),
                                        seq_len(B/n_rep)))
for (sim in  c("TO","CS","CCS"))
      for (r in seq_len(B/n_rep)){
      AP_performance_array["modularity",sim,r]= modularity_fromMatrix(base_S, AP_simXgene[sim,r,])
      AP_performance_array["conductance",sim,r]= conductance_network(base_S,  AP_simXgene[sim,r,])
      message(sprintf("sim=%s,r=%d",sim,r))
      }

save(AP_simXgene, AP_performance_array,AP_ARI_matrices, file="AP_results_aggregated.RData")

