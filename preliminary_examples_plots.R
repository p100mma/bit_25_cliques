
library(FCPS)
library(igraph)
library(cliquePartitioneR)
library(magick)
D<-dist(Hepta$Data) |> as.matrix()
dim(D)
set.seed(123)
subsample<- sample(1:nrow(D), 30, replace=FALSE)
D<- D[subsample, subsample]


cmdscale(D, k = 2) -> mds
D<- dist(mds) |> as.matrix()
hcl<- hclust(as.dist(D), method = "complete")

complexity_score<-function (clique_labels) 
{
  clique_labels <- uniqueSingletonLabels(clique_labels)
  tabulation <- fastTable(clique_labels)
  -(lfactorial(length(tabulation$value)) + sum(lfactorial(tabulation$count)))
}

partitions<-lapply(hcl$height, function(h)
      
            cutree(hcl, h = h)
          
  )

scores<- sapply(partitions, function (P) complexity_score(P))


which.max(scores)

hcl$height[which.max(scores)]

examples<- c(3,which.max(scores), 27)
descriptions<- rep(NA, length(scores))
descriptions[examples[[1]]]<- "too high"
descriptions[which.max(scores)]<- "just right"
#descriptions[480]<- "too high"
descriptions[examples[[3]]]<- "too low"

G_D<- graph_from_adjacency_matrix(D, mode = "undirected", 
                                                weighted = TRUE)
set.seed(32)



plot(G_D, layout=mds, vertex.label=NA, vertex.size=4.5)

colors<-rep(NA, length(scores))
colors<- seq(from=1,to=0, length.out=length(scores))
colors<- rgb(colors^0.25,0,0)
#colors[examples]<- c( rgb(1,0,0),rgb(0.5,1,0),rgb(0,0,1))

pdf("complexity_score_vs_cutoff.pdf", width=8, height=6)
par (#mfrow=c(1,1),
  mfrow=c(2,2),
oma = c(1, 1, 0, 0),    # outer margins (bottom, left, top, right)
mar = c(1, 1, 1, 1))    # inner margins for each plot
plot(max(D) - hcl$height, col=colors, pch=19, cex=1,
     scores, type = "b", xlab = "similarity cutoff", ylab = "Complexity Score",
     main = "Complexity Score vs similarity cutoff")

for (i in examples[c(1,3,2)]) {
  D_t<- D
  D_t[D_t > hcl$height[[i]] ]<- 0
  G_Dt<- graph_from_adjacency_matrix(D_t, mode = "undirected", 
                                                weighted = TRUE)
  plot(G_Dt, layout=mds, vertex.label=NA, vertex.size=8,
       vertex.color= partitions[[i]])
  do.call(text, c(as.list(colMeans(mds)-0.5), 
                          list(labels = descriptions[i], cex = 1.5, 
       col = colors[[i]], pos = 3)
  ))
}

dev.off()
library(magick)
img <- image_read_pdf("complexity_score_vs_cutoff.pdf", density = 300)  # higher density = higher DPI
image_write(img, path = "complexity_score_vs_cutoff.png", format = "png")


D<- dist(Chainlink$Data) |> as.matrix()

rotate3D <- function(points, degrees) {
  # Convert degrees to radians
  theta <- degrees * pi / 180
  
  # Rotation matrices around X, Y, Z axes
  Rx <- matrix(c(1, 0, 0,
                 0, cos(theta), -sin(theta),
                 0, sin(theta),  cos(theta)), nrow = 3, byrow = TRUE)
  
  Ry <- matrix(c(cos(theta), 0, sin(theta),
                 0, 1, 0,
                 -sin(theta), 0, cos(theta)), nrow = 3, byrow = TRUE)
  
  Rz <- matrix(c(cos(theta), -sin(theta), 0,
                 sin(theta),  cos(theta), 0,
                 0, 0, 1), nrow = 3, byrow = TRUE)
  
  # Combined rotation: R = Rz * Ry * Rx
  R <- Rz %*% Ry %*% Rx
  
  # Rotate the points: matrix multiplication
  t(R %*% t(points))
}

layout_3D<- rotate3D(Chainlink$Data, -45)
layout_2D<- layout_3D[, 1:2]
ch_col= Chainlink$Cls

col_scale<- function(x) {
  max_val<- max(x);min_val<-min(x)
  (x - min_val) / (max_val - min_val)
}

ch_col= ifelse(ch_col == 1, "red", 
               rgb(0,0,layout_3D[,3] |> col_scale()))


G_D<- graph_from_adjacency_matrix(D, mode = "undirected", 
                                                weighted = TRUE)

mst(G_D)-> mst_D


sort( igraph::E(mst_D)$weight, decreasing=TRUE)-> dMSTsorted
dMSTsorted<- dMSTsorted[ dMSTsorted!=0 ]
ddiffs= -1*(diff(dMSTsorted))
dmeans= (dMSTsorted[ 2: length(dMSTsorted) ] + dMSTsorted[ 1: (length(dMSTsorted)-1) ])/2
dmeans[[ which.max(ddiffs) ]]-> d_t
d_t

E(mst_D)$color="gray"
E(mst_D)$width=0.5

E(mst_D)$color[ order(abs(E(mst_D)$weight - d_t))[1:2] ]<- "red"

#E(mst_D)$color[ order(abs(E(mst_D)$weight - d_t))[1:2] ]<- "red"
E(mst_D)$width[ order(abs(E(mst_D)$weight - d_t))[1:2] ]<- 4



D_t<- D
D_t[ D > d_t] = 0

hcl<-hclust(as.dist(D), method = "complete")

clqs<-cutree(hcl, h = d_t )



G_dt<- graph_from_adjacency_matrix(D_t, mode = "undirected", 
                                                weighted = TRUE,diag=FALSE)

G_dt_mst<- mst(G_dt)
pdf("mst_thr_cliques.pdf")
par(mfrow=c(2,2),
    oma=c(0,0,0,0))

plot(layout_3D[order(layout_3D[,3]),], 
     col= ch_col[order(layout_3D[,3])], pch=16,xlab=NA,
     ylab=NA,
     xaxt="n", yaxt="n", main="vector space")

plot(mst_D, layout=layout_2D, vertex.label=NA, vertex.size=2, 
     edge.color="black", edge.width=0.5, vertex.color= ch_col, 
     main = "MST"
     
)

plot(mst_D, layout=layout_2D, vertex.label=NA, vertex.size=3, 
     vertex.color= ch_col, 
     main = "Critical thr from MST- between \n red distances "
)


plot(G_dt_mst, layout=layout_2D, vertex.label=NA, vertex.size=6, 
     edge.color="gray", edge.width=0.5, vertex.color= clqs, 
     main = "Clique partitioning at \n threshold from MST")

dev.off()


img <- image_read_pdf("mst_thr_cliques.pdf", density = 300)  # higher density = higher DPI
image_write(img, path = "mst_thr_cliques.png", format = "png")
