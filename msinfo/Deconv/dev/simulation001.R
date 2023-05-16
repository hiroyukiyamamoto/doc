rm(list=ls(all=TRUE))

gauss.density1 <- function(x) 1/sqrt(2*pi)*exp(-x^2/2)  # 標準正規分布の密度
gauss.density2 <- function(x) 1/sqrt(2*pi)*exp(-x^2/2)  # 標準正規分布の密度

x <- seq(-5,5,0.5)
#z1 <- as.matrix(gauss.density1(x)+0.001*rnorm(21),ncol=1)
#z2 <- as.matrix(gauss.density2(x)+0.001*rnorm(21),ncol=1)
#
# z1 <- as.matrix(gauss.density1(x)+0.001*runif(21, min=0, max=1))
# z2 <- as.matrix(gauss.density1(x)+0.001*runif(21, min=0, max=1))

z1 <- as.matrix(gauss.density1(x))
z2 <- as.matrix(0.5*gauss.density2(x))

plot(z1, type="l", ylim=c(0,0.4))
par(new=T)
plot(z2,type="l", col="red", ylim=c(0,0.4))

mz <- c(1:300)
s1 <- rep(0,300)
s1[sample(mz)[1:10]] <- runif(10, min=0, max=1)
s2 <- rep(0,300)
s2[sample(mz)[1:10]] <- runif(10, min=0, max=1)

X <- z1%*%s1+z2%*%s2
#index <- which(apply(X,2,sd)!=0)
#X <- X[,index]

icacom <- ica(t(X),2, method="imax")

par(mfrow = c(1, 2))
plot(icacom$M[,1],type="l")
plot(icacom$M[,2],type="l")

#plot(icacom$S[,1],type="h")
#plot(icacom$S[,2],type="h")

# -------------------------------

cor.test(s1,icacom$S[,1])
cor.test(s1,icacom$S[,2])

cor.test(s2,icacom$S[,1])
cor.test(s2,icacom$S[,2])

# --------------------------

pcacom <- prcomp(X)

plot(pcacom$x[,1],type="l")
plot(pcacom$x[,2],type="l")

cor.test(s1,pcacom$rotation[,1])
cor.test(s1,pcacom$rotation[,2])

cor.test(s2,pcacom$rotation[,1])
cor.test(s2,pcacom$rotation[,2])

# ---------------------------

### 関数化して、suppressする
com <- 2
C <- matrix(NA,com,nrow(X))
for(k in 1:com){
  S <- icacom$M[,k]
  R <- cor(X,S)

  if (R[which.max(abs(R))]<0){
    S <- -S
  }
  S[S<0] <- 0
  C[k,] <- S
}

C <- t(C)

# number of component
lambda <- 0
maxiter <- 5

E <- NULL
for (k in 1:maxiter){
  A <- solve(t(C)%*%C+lambda*diag(1,ncol(C)))%*%t(C)%*%X
  A[A<0] <- 0

  # normalized constraint
  #for (i in 1:com){
  #  A[i,] <- A[i,]/as.numeric(sqrt(t(A[i,])%*%A[i,]))
  #}

  C <- X%*%t(A)%*%solve(A%*%t(A)+lambda*diag(1,ncol(C)))
  C[C<0] <- 0

  # error
  E[k] <- norm(X-C%*%A);
}
