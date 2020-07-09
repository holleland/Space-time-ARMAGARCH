true.par <- c(5.032,.354,-.153,.122,-.051,.513,.142,.171)
par <- true.par[-(1:5)]

library(starmagarch)
W<-create.neighbourhood.array(m = 30, sp = 2, type = "rook", torus = TRUE)
W <- W[,,2]*4
A <- par[2]*W
B <- par[3]*W

k <- seq(1,30^2,31)
k2 <- rep(1,30^2)
k2[k]<-3
k2

EQ<-  kronecker(A,A)%*% diag(k2) + kronecker(A,B)+ kronecker(B,A)+kronecker(B,B)
print(max(s<-Mod(eigen(EQ)$values)))

