
a <- read.csv('D:/CodeMesh/ci/data/libsvm/mg_scale.csv', header=F)

summary(lm(V1 ~ ., a))

library(MASS)

lm.ridge(V1 ~ . , a, lambda=10)$xm

lambda = 10
X = cbind(1, as.matrix(a[, -1]))
y = as.matrix(a[, 1])
xx = t(X) %*% X
xy = t(X) %*% y
beta = ginv(xx + lambda*diag(rep(1,dim(X)[2]))) %*% xy

#%%MatrixMarket matrix array real general
#1 7
#0.914397
#0.130746
#0.0880258
#-0.148705
#-0.143011
#0.0471996
#0.0905323


df <- data.frame(x = c(1,2,3,4), y = c(2,3,4,7), z = c(5,6, 9, 1))
f <- formula('z~x*y')
model.frame(f, df)




yX <- read.table('D:/CodeMesh/ci/sdca/sparSDCA/exp/inputs/dataset.csv', sep=' ')
dput(yX)

a <- structure(list(V1 = c(251.3, 251.3, 248.3, 267.5, 273, 276.5, 
                           270.3, 274.9, 285, 290, 297, 302.5, 304.5, 309.3, 321.7, 330.7, 
                           349), V2 = c(41.9, 43.4, 43.9, 44.5, 47.3, 47.5, 47.9, 50.2, 
                                      52.8, 53.2, 56.7, 57, 63.5, 65.3, 71.1, 77, 77.8), 
                    V3 = c(29.1,29.3, 29.5, 29.7, 29.9, 30.3, 30.5, 30.7, 30.8, 30.9, 31.5, 31.7, 
                                                                                                  31.9, 32, 32.1, 32.5, 32.9)), .Names = c("V1", "V2", "V3"), class = "data.frame", row.names = c(NA, 
                                                                                                                                                                                                  -17L))


X <- as.matrix(yX[, 2:3])
y <- as.vector(yX[, 1])

X <- cbind(1, X)
b <- solve(t(X) %*% X) %*% t(X) %*% y

y_hat <- X %*% b

sqrt(mean((y-y_hat)^2))
