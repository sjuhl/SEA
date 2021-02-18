##########################
# Spatial Eigenfunction
# Analysis
##########################

# clean environment and set working directory
rm(list=ls())
setwd("~/Dropbox/Spatial Filtering IV/Eigenvectors")

# load necessary packages
library(spdep)
library(spfilteR)
library(foreach)
library(doParallel)
library(doRNG)

############
# SIMULATION
############

# generate W2: 10x10 grid
W2 <- nb2mat(cell2nb(nrow=10,ncol=10,type="rook",torus=T))

# range of rho
c(1/min(eigen(W2)$values),1/max(eigen(W2)$values))

# EXPERIMENT #1
# number of simulation iterations
nsim <- 1000

# parameters
n <- nrow(W2)
p.exp1 <- seq(from=-.9,to=.9,by=.1)

# spatial multipliers
multi1 <- list()
for(i in seq_along(p.exp1)){
  multi1[[i]] <- solve(diag(1,n) - p.exp1[i]*W2)
}

# first input
input1 <- as.data.frame(cbind(rep(p.exp1,each=nsim),rep(1:length(p.exp1),each=nsim)))
colnames(input1) <- c("p","multi_id")
nrow(input1)
unique(input1)


# firstsimulation function
simfunc1 <- function(n,multi_id,alpha=.05,adj=T){
  u <- rnorm(n=n,0,1)
  x <- multi1[[multi_id]] %*% u
  MItwo <- MI.vec(x=x,W=W2,alternative="two.sided")$pI <= alpha
  MIg <- MI.vec(x=x,W=W2,alternative="greater")$pI <= alpha
  MIl <- MI.vec(x=x,W=W2,alternative="lower")$pI <= alpha
  MIlocal <- sum(MI.local(x=x,W=W2,alternative="two.sided")$pIi  <= alpha)
  d <- MI.decomp(x=x,W=W2,nsim=1000)
  adj.alpha <- ifelse(adj,alpha/2,alpha) # bonferroni adjustment
  MIdec <- d["pI+"]  <= adj.alpha || d["pI-"]  <= adj.alpha
  
  # output
  out <- data.frame(MItwo=MItwo,MIg=MIg,MIl=MIl,MIlocal=MIlocal,MIdec=MIdec)
  return(out)
}

# set up parallel computing
# parallel computing
nworkers <- 7
cl <- makeCluster(nworkers)
registerDoParallel(cl)
registerDoRNG(123)

sim_out1 <- foreach(i=1:nrow(input1), .combine=rbind,
                    .packages = c("spfilteR")
                    ) %dopar% {
                      simfunc1(n=n,multi_id=input1$multi_id[i],alpha=.05,adj=T)
}
stopCluster(cl)




# EXPERIMENT #2
p.exp2 <- seq(from=0,to=.9,by=.1)

multi21 <- multi22 <- list()
for(i in seq_along(p.exp2)){
  multi21[[i]] <- solve(diag(1,n) - p.exp2[i]*W2)
  multi22[[i]] <- solve(diag(1,n) - -p.exp2[i]*W2)
}
length(multi21)


# second input
input2 <- as.data.frame(cbind(rep(p.exp2,each=nsim),-rep(p.exp2,each=nsim),rep(1:length(p.exp2),each=nsim)))
colnames(input2) <- c("p1","p2","multi_id")
nrow(input2)
unique(input2)


# simulation function
simfunc2 <- function(n,multi_id,alpha=.05,adj=T){
  u <- rnorm(n=n,0,1)
  v <- rnorm(n=n,0,1)
  x <- multi21[[multi_id]] %*% u + multi22[[multi_id]] %*% v
  MItwo <- MI.vec(x=x,W=W2,alternative="two.sided")$pI <= alpha
  MIg <- MI.vec(x=x,W=W2,alternative="greater")$pI <= alpha
  MIl <- MI.vec(x=x,W=W2,alternative="lower")$pI <= alpha
  MIlocal <- sum(MI.local(x=x,W=W2,alternative="two.sided")$pIi  <= alpha)
  d <- MI.decomp(x=x,W=W2,nsim=1000)
  adj.alpha <- ifelse(adj,alpha/2,alpha)# bonferroni adjustment
  MIdec <- d["pI+"]  <= adj.alpha || d["pI-"]  <= adj.alpha
  
  # output
  out <- data.frame(MItwo=MItwo,MIg=MIg,MIl=MIl,MIlocal=MIlocal,MIdec=MIdec)
  return(out)
}

simfunc2(n=n,multi_id=input2$multi_id[1],alpha=.05,adj=T)


# parallel computing
nworkers <- 7
cl <- makeCluster(nworkers)
registerDoParallel(cl)
registerDoRNG(123)

sim_out2 <- foreach(i=1:nrow(input2), .combine=rbind,
                    .packages = c("spfilteR")
                    ) %dopar% {
                      simfunc2(n=n,multi_id=input2$multi_id[i],alpha=.05,adj=T)
}
stopCluster(cl)




# SAVE
#save(sim_out1,sim_out2,input1,input2,p.exp1,p.exp2,nsim,file="./MCExperiment.RData")
