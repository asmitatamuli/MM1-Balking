##############################MH algorithm##############################################################
library("stats")
library("GenSA")
g=1000
aa=array(0,dim=c(g,1))
for( j in 1:g) 
{
n=50
x1=array(0,dim=c(n))    
rho=5
p=0.5

target <- function(x, rho) {
  return(rho^(x + 1) / (factorial(x + 1) * (exp(rho)-1)))
}

proposal <- function(x, p) {
  return(dgeom(x, p))
}

metrop <- function(N, y0, rho, p) {
  sample <- 0
  for (i in 1:N) {
    x.star <- rgeom(1, p)
    r <- ((target(x.star, rho))/ (target(y0, rho))) * ((proposal(y0, p))/ (proposal(x.star, p)))
    u <- runif(1)
    if (u < min(1, r)) {
      y1 <- x.star
    } else {
      y1 <- y0
    }
    sample[i] <- y1
    y0 <- y1
  }
  return(sample)
}


entire_sample <- metrop(50000,2, 5, 0.5)
filtered_sample<- entire_sample[49951:50000]
log_likelihood<-function(rho){
-sum(log(rho^(filtered_sample + 1) / (factorial(filtered_sample + 1) * (exp(rho)-1))))
}
global.min <- 0
tol <- 1e-13
lower <- 0
upper <- 1
out <- GenSA(lower = lower, upper = upper, fn = log_likelihood,
control=list(threshold.stop=global.min+tol,verbose=TRUE))
out[c("value","par","counts")]
print(out$par)
bb=out$par
cc <- array(out$par, dim=c(1,1))
for( k in 1:1)
{
aa[j,k]=cc[1,k]
}
print(j)
}
aa
Ls<-((1-exp(aa)+(exp(aa)*aa))/(exp(aa)-1))
Ls
Lq<-((2-(2*exp(aa))+aa+(exp(aa)*aa))/(exp(aa)-1))
Lq
mean(aa[,1])
sd(aa[,1])
mean(Ls[,1])
sd(Ls[,1])
mean(Lq[,1])
sd(Lq[,1])
rhobias<-mean(aa[,1])-5
rhobias
rhovar<-(sd(aa[,1]))^2
rhoMSE<-(rhobias^2)+rhovar
rhoMSE
Lsbias<-mean(Ls[,1])-4.03392
Lsbias
Lsvar<-(sd(Ls[,1]))^2
LsMSE<-(Lsbias^2)+Lsvar
LsMSE
Lqbias<-mean(Lq[,1])-3.06784
Lqbias
Lqvar<-(sd(Lq[,1]))^2
LqMSE<-(Lqbias^2)+Lqvar
LqMSE

#######################################SIR algorithm#############################################################
library(LaplacesDemon)
library(stats)

MM1 <- function(a, b, rho, n, Ls, Lq) {
  start_time <- Sys.time()
  provas <- 200
  P0 <- 0
  
    for (i in 1:provas) {
    P1 <- (1 / factorial(i + 1)) * (rho ^ i)
    P0 <- P1 + P0
  }
  P0 <- 1 / P0    
  Pn <- rep(0, provas)
  N <- rep(0, provas)
  cont <- 1
  
  for (i in 0:provas) {
    P <- (1 / factorial(i + 1)) * (rho ^ i) * P0
    Pn[cont] <- P
    N[cont] <- i
    cont <- cont + 1
  }
  Pn <- Pn / sum(Pn)   
  maxss <- 50
  init <- 0
  final <- maxss
  Tasir <- 5000
  #for beta prior

  Priori <- rinvbeta(Tasir, a, b)
  #for lognormal use rlnorm(Tasir, a, b) and for gamma use rgamma(Tasir, a, b)

  Priori[is.na(Priori)] <- 0  
  
  while (init < n) {
    Ta <- min(final, n) - init
    R <- sample(N, Ta, replace = TRUE, prob = Pn)
    Vero <- rep(0, Tasir)
    
    for (j in 1:Tasir) {
      Pri <- Priori[j]
      P0 <- 0
      
       for (i in 1:provas) {
        P1 <- (1 / factorial(i + 1)) * (Pri ^ i)
        P0 <- P1 + P0
      }
      P0 <- 1 / (P0 + 1)  
      
      E1 <- 1
      for (i in 1:Ta) {
        i1 <- R[i]
        E1a <- ((1 / factorial(i1 + 1)) * (Pri ^ i1)) * P0
        E1 <- E1 * E1a
      }
      Vero[j] <- E1
    }
    
    Vero[is.na(Vero)] <- 0  
    Vero <- Vero / sum(Vero)  
    Posteriori <- sample(Priori, Tasir, replace = TRUE, prob = Vero)
    
    init <- min(final, n)
    final <- final + maxss
    Priori <- Posteriori
  }
  
  rhonew <- mean(Posteriori)
  rhobias<-mean(Posteriori)-rho
  rhoMSE <- (rhobias^2)+((sd(Posteriori))^2)
  Ls1<-((1-exp(Posteriori)+(exp(Posteriori)*Posteriori))/(exp(Posteriori)-1))
  Ls2<-mean(Ls1)
  Lsbias<- Ls2-Ls
  Ls1MSE <- (Lsbias^2)+((sd(Ls1))^2)
  Lq1<-((2-(2*exp(Posteriori))+Posteriori+(exp(Posteriori)*Posteriori))/(exp(Posteriori)-1))
  Lq2<-mean(Lq1)
  Lqbias<- Lq2-Lq
  Lq1MSE <- (Lqbias^2)+((sd(Lq1))^2)

  return(list(rhonew = rhonew, rhoMSE = rhoMSE, Ls2 = Ls2, Ls1MSE =Ls1MSE, Lq2 = Lq2, Lq1MSE =Lq1MSE))
}

MM1(5, 5, 5, 50, 4.03392, 3.06784)
