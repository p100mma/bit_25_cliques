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
set.seed(1324)
BS<- lapply(1:B, function(b) sample(1:nrow(X), nrow(X), replace=TRUE))
Pg_<- readRDS("Pg_.rds")
t_<-readRDS("t_.rds")
tom<- function(A) {diag(A)=0; deg<-colSums(A);  (A + A%*%A)/( outer(deg,deg,pmin) +1 - A )  }
library(apcluster)
ap_grouping<- function(sim, q, seed){
  apcluster(s=sim,seed=seed, p=q) |> labels(type="enum")
}

r=commandArgs(trailingOnly=TRUE)[[1]] |> as.integer()
n_rep=30
stopifnot(r %in% seq_len(B/n_rep))

AP_simXgene<- array(NA, dim=c(3,B/n_rep,ncol(X)),
                    dimnames = list(
                     c("TO","CS","CCS"),
                     seq_len(B/n_rep),
                     colnames(X)
                    )
                      )

    CCSg<-matrix(0,nrow=ncol(X),ncol=ncol(X))
    CS<- matrix(0, nrow=ncol(X), ncol=ncol(X))
    b_start<- 1 + ( r- 1)*n_rep
    message(b_start)
    for  (b in c(b_start: (b_start+ (n_rep-1)  ))){
      start_t<- Sys.time()
      gc()
      absCor<- abs(corfast(X[ BS[[b]],]|>as.matrix() ))
      if (b==b_start) { TOM<- tom(absCor^2)
         AP_TO= ap_grouping(TOM, q = quantile(TOM,0.1), seed = b)
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
      message(sprintf(" %d/ %d",b, b_start + (n_rep-1)))
      message(Sys.time()-start_t)
    }
    CCSg<- CCSg/n_rep
    CS<- CS/n_rep
    AP_CS= ap_grouping(CS, q = quantile(CS,0.1), seed = b)
    AP_CCS=ap_grouping(CCSg, q = quantile(CCSg,0.1), seed = b)


saveRDS( list(TO=AP_TO, CS=AP_CS, CCS=AP_CCS),
sprintf("AP_rep_%d.rds", r))




