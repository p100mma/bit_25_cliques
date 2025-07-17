
library(data.table)
fread("Brain_GSE50161.csv") |> as.data.frame()-> microarray_data
head(microarray_data)[,1:6]
X<- microarray_data[,3:ncol(microarray_data)]
microarray_data$type -> y_char
y<- as.numeric(factor(y_char))
library(matrixStats)
colMedsX<- colMedians(X|> as.matrix())
plot(sort(colMedsX), main="median expression values per gene")
abline(h=5)
X<- X[, colMedsX > 5 ]
colVarsX<- colVars(X |> as.matrix())
plot(sort(colVarsX), main="gene variances")
abline(h=1.5)
X<- X[, colVarsX >1.5]
X<- X[, order(-colVars(X|> as.matrix()))]
X<-X[, !grepl("_(x|a|g|r)_at$", colnames(X))]
X<-X[, !grepl("^AFFX-", colnames(X))]
