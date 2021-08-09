##########################
# Spatial Eigenfunction
# Analysis
##########################

# clean environment and set working directory
rm(list=ls())
setwd("~/Dropbox/Spatial Filtering IV/SEA")

# load necessary packages
library(spdep)
library(spatialreg)
library(spfilteR)
library(RColorBrewer)
library(scales)

# read German NUTS 3 regions shapefile
ger3 <- st_read(dsn="./Shapefile/Germany")

# load Eurostat data
load("./Data/eurostat2017.RData")

# specify spatial matrix
nb <- poly2nb(ger3$geometry, queen=TRUE, row.names=ger3$NUTS_ID)
sty <- "W"
W <- nb2mat(nb,style=sty)

# Moran coefficient
MI.vec(x=data[,2:ncol(data)],W=W)

# eigenfunction decomposition
EV <- getEVs(W=W)

all.equal(target=0, current=EV$values[171], tolerance=1e-10)
all.equal(target=0, current=EV$values)

# identify constant eigenvector
which(vapply(EV$values, function(x) isTRUE(all.equal(target=0,current=x,tol=1e-10))
             ,FUN.VALUE=logical(1)))


# define color palette
brewer.pal.info
col <- brewer.pal(n=9, "Purples")


### FIGURE 1
# empty figure
png(paste0("./Figures/empty.png"),width=400, height=450)
par(oma=c(0,0,0,0),mar=c(1,0,0,0))
plot(st_geometry(ger3),col="white",border="white")
dev.off()

# plain map of Germany
png(paste0("./Figures/GER_plain.png"),width=400, height=450)
par(oma=c(0,0,0,0),mar=c(1,0,0,0))
plot(st_geometry(ger3),col="white")
dev.off()

# eigenvectors to plot
evs <- c(1,25,50,75,170,401)

# color codes
colcode <- matrix(NA,ncol=length(evs),nrow=nrow(W))
for(i in evs){
  seq <- quantile(EV$vectors[,i],probs=c(0,1/9,2/9,3/9,4/9,5/9,6/9,7/9,8/9))
  colcode[,match(i,evs)] <- sapply(EV$vectors[,i],function(x) sum(x >=seq))
}

# plots
for(i in seq_along(evs)){
  png(paste0("./Figures/EV",i,".png"),width=400, height=450)
  par(oma=c(0,0,0,0),mar=c(1,0,0,0))
  plot(st_geometry(ger3),col=col[colcode[,i]])
  mtext(paste("MC =",round(EV$moran[evs[i]],3)),side=1,line=-.5,cex=1.5)
  dev.off()
}

length(EV$values[EV$moran>0])
length(EV$values[EV$moran<0])
EV$vectors[,171]
all(vapply(EV$vectors[,171], function(x) isTRUE(all.equal(target=1/sqrt(nrow(W))
                                                          ,current=x,tol=1e-10))
           ,FUN.VALUE=logical(1)))
all.equal(0, diff(range(EV$vectors[,171])))



# MC of all eigenvectors
par(oma=c(0,.2,0,0))#,mgp=c(3, 2, 0)
png("./Figures/EV_MC.png",width=400, height=450)
plot(0,ylim=c(min(EV$moran),max(EV$moran)),xlim=c(1,length(EV$moran))
     ,axes=F,ann=F,type="n",las=1)
# EVs
points(y=EV$moran,x=1:length(EV$moran),pch=18,cex=.8,col=alpha("black",.8))
abline(h=0,lty=2,cex=.5)
# axes & labels
axis(1,cex.axis=1.7)
axis(2,cex.axis=1.7)
mtext("Eigenvectors",side=1,line=2.5,outer=F,cex=2)
mtext("MC",side=2,line=2.5,outer=F,cex=2)
dev.off()



### FIGURE 2
# source functions to calculate local MC
#source("../../spfilteR Package/spfilteR/R/MI.local.R")
#source("../../spfilteR Package/spfilteR/R/utils.R")

# eigenvectors to plot
EV$moran[80]
EV$moran[350]
evs <- c(80,350)

EVs <- cbind(EV$vectors[,evs],apply(EV$vectors[,evs],1,sum))
MI.vec(x=EVs,W=W)

# color codes
colcode <- matrix(NA,ncol=ncol(EVs),nrow=nrow(W))
for(i in seq_len(ncol(EVs))){
  seq <- quantile(EVs[,i],probs=c(0,1/9,2/9,3/9,4/9,5/9,6/9,7/9,8/9))
  colcode[,i] <- sapply(EVs[,i],function(x) sum(x >=seq))
}

# plots
for(i in seq_len(ncol(EVs))){
  # map
  png(paste0("./Figures/decomp",i,".png"),width=400, height=450)
  par(oma=c(0,0,0,0),mar=c(1,0,0,0))
  plot(st_geometry(ger3),col=col[colcode[,i]])
  EI <- MI.vec(x=EVs[,i],W=W)$EI
  mc <- round(MI.vec(x=EVs[,i],W=W)$I,3)
  pval <- round(MI.vec(x=EVs[,i],W=W
                       ,alternative="two.sided")$pI
                ,3)
  mtext(bquote("global MC ="~.(format(mc,nsmall=3))~~~
                 "("*italic("p")~"="~.(format(pval,nsmall=3))*")")
        ,side=1,line=-.5,cex=1.5)
  dev.off()
  # histogram
  #localMCs <- localmoran(x=EVs[,i],listw=mat2listw(W))[,"Z.Ii"]
  localMCs <- MI.local(x=EVs[,i],W=W)[,"zIi"]
  upper <- qnorm(p=.05/2, lower.tail=F)
  lower <- -upper
  freq_upper <- sum(localMCs>upper)
  freq_lower <- sum(localMCs<lower)
  png(paste0("./Figures/hist",i,".png"),width=400, height=250)
  par(oma=c(1.5,1.5,0,0),mar=c(2,3,0,0))
  hist(localMCs,ylim=c(0,250),xlim=c(-9,9),breaks=50,axes=F,ann=F)
  axis(1,at=c(-9,-6,-3,0,3,6,9),cex=1.5)
  axis(2,at=c(0,50,150,250),las=1,cex=1.5)
  mtext("Frequency",side=2,line=3,cex=2)
  mtext("Standardized local MCs",side=1, line=2.5,cex=2)
  abline(v=lower,lty=5)
  abline(v=upper,lty=5)
  text(x=6,y=150,freq_upper,cex=1.5)
  text(x=6,y=110
       ,bquote("("*.(round(freq_upper/length(localMCs)*100,2))*"%)")
       ,cex=1.5)
  text(x=-6,y=150,freq_lower,cex=1.5)
  text(x=-6,y=110
       ,bquote("("*.(round(freq_lower/length(localMCs)*100,2))*"%)")
       ,cex=1.5)
  dev.off()
}

# MC decomposition
MI.vec(x=EVs,W=W,alternative="two.sided")
set.seed(123)
dec <- MI.decomp(x=EVs,W=W,nsim=1000)
round(dec[,-c(4,8,10)],3)



### FIGURE 4
filter <- lmFilter(y=data$medianage,W=W,objfn="MI")
toplot <- cbind(data$medianage,filter$other$sf,filter$residuals[,"Filtered"])
MCs <- c(filter$moran["Initial","Observed"],NA,filter$moran["Filtered","Observed"])
pvals <- c(filter$moran["Initial","p-value"],NA,filter$moran["Filtered","p-value"])
text <- c("Median Age","Spatial Filter", "Filtered Residuals")
# plot
for(i in seq_len(ncol(toplot))){
  seq <- quantile(toplot[,i],probs=c(0,1/9,2/9,3/9,4/9,5/9,6/9,7/9,8/9))
  colcode <- sapply(toplot[,i],function(x) sum(x >=seq))
  # map
  png(paste0("./Figures/filter",i,".png"),width=400, height=450)
  par(oma=c(0,0,0,0),mar=c(1,0,0,0))
  plot(st_geometry(ger3),col=col[colcode])
  if(i!=2){
    mtext(bquote("MC ="~.(format(MCs[i],digits=3))~~~
                   "("*italic("p")~"="~.(format(round(pvals[i],10),nsmall=3,digits=3))*")")
          ,side=1,line=-.5,cex=1.5)
  }
  dev.off()
}

filter$other$nev

filter.R2 <- lmFilter(y=data$medianage,W=W,objfn="R2")
filter.p <- lmFilter(y=data$medianage,W=W,objfn="p",bonferroni=F)
filter.pMI <- lmFilter(y=data$medianage,W=W,objfn="pMI",bonferroni=F)

filter.R2$other$nev
filter.p$other$nev
filter.pMI$other$nev

round(filter$fit["Filtered"],3)
round(filter.R2$fit["Filtered"],3)
round(filter.p$fit["Filtered"],3)
round(filter.pMI$fit["Filtered"],3)

MI.resid(resid=filter$residuals[,"Filtered"],x=filter$selvecs,W=W,alternative="two.sided")
MI.resid(resid=filter.R2$residuals[,"Filtered"],x=filter.R2$selvecs,W=W,alternative="two.sided")
MI.resid(resid=filter.p$residuals[,"Filtered"],x=filter.p$selvecs,W=W,alternative="two.sided")
MI.resid(resid=filter.pMI$residuals[,"Filtered"],x=filter.pMI$selvecs,W=W,alternative="two.sided")

round(filter$moran["Filtered",],3)
round(filter.R2$moran["Filtered",],3)
round(filter.p$moran["Filtered",],3)
round(filter.pMI$moran["Filtered",],3)

round(cor(x=filter$other$sf,data$medianage),3)
round(cor(x=filter.R2$other$sf,data$medianage),3)
round(cor(x=filter.p$other$sf,data$medianage),3)
round(cor(x=filter.pMI$other$sf,data$medianage),3)



# variation partitioning
hist(data$gdp_mio)
hist(log(data$gdp_mio))
MI.vec(x=log(data$gdp_mio),W=W)
y <- log(data$gdp_mio)
X <- cbind(data$medianage)
# use spatial filtering to select relevant eigenvectors
fil <- lmFilter(y=y,x=X,W=W,objfn="MI")
summary(fil)
# variation partitioning and Moran spectral randomization
set.seed(123)
(var <- vp(y,x=X,evecs=fil$selvecs, msr=1000))
round(cbind(var$R2,var$adjR2),3)



# simulating spatial patterns
sar <- lagsarlm(data$medianage~1,listw=mat2listw(W,style=sty))
filter <- lmFilter(y=data$medianage,W=W,objfn="MI")
sar.rho <- sar$rho
sar.con <- sar$coefficients
multiplier <- solve(diag(1,nrow(W))-sar.rho*W)

filter.con <- filter$estimates[,"Estimate"]
filter.sf <- EV$vectors[,filter$other$sel_id] %*% filter$EV[,"Estimate"]

sd(data$medianage)


nsim <- 10000
sar.out <- filter.out <- matrix(NA,nrow=nrow(W),ncol=nsim)
set.seed(123)
for(i in seq_len(nsim)){
  e <- rnorm(n=nrow(W),mean=0,sd=sd(data$medianage))
  sar.out[,i] <- multiplier %*% (sar.con + e)
  filter.out[,i] <- filter.con + filter.sf + e
}

range(data$medianage)
range(apply(sar.out,1,mean))
range(apply(filter.out,1,mean))


### FIGURE 6
seq <- quantile(data$medianage,probs=c(0,1/9,2/9,3/9,4/9,5/9,6/9,7/9,8/9))
# mpas
png(paste0("./Figures/simsar.png"),width=450, height=450)
par(oma=c(0,0,0,0),mar=c(1,0,0,0))
colcode <- sapply(apply(sar.out,1,mean),function(x) sum(x >=seq))
plot(st_geometry(ger3),col=col[colcode])
range <- range(apply(sar.out,1,mean))
mtext(bquote("Range = ["*.(format(range[1],digits=3,nsmall=3))*";"~
               .(format(range[2],digits=3,nsmall=3))*"]")
      ,side=1,line=-.5,cex=1.5)
dev.off()
png(paste0("./Figures/simfilter.png"),width=450, height=450)
par(oma=c(0,0,0,0),mar=c(1,0,0,0))
colcode <- sapply(apply(filter.out,1,mean),function(x) sum(x >=seq))
plot(st_geometry(ger3),col=col[colcode])
range <- range(apply(filter.out,1,mean))
mtext(bquote("Range = ["*.(format(range[1],digits=3,nsmall=3))*";"~
               .(format(range[2],digits=3,nsmall=3))*"]")
      ,side=1,line=-.5,cex=1.5)
dev.off()


# F-test of equal variances
var.filter <- apply(filter.out,1,var)
range(var.filter)
min1 <- min(var.filter)
max1 <- max(var.filter)
2*pf(max1/min1,df1=400,df2=400,lower=F)

var.sar <- apply(sar.out,1,var)
range(var.sar)
min2 <- min(var.sar)
max2 <- max(var.sar)
2*pf(max2/min2,df1=400,df2=400,lower=F)



# correlation between maps
cor.sar <- cor(sar.out)
diag(cor.sar) <- NA
range(cor.sar,na.rm=T)


############
# SIMULATION
############

# load results
load("MCExperiment.RData")


# PLOT EXPERIMENT 1
power <- matrix(NA,ncol=(ncol(sim_out1)-1),nrow=length(p.exp1))
colnames(power) <- colnames(sim_out1)[colnames(sim_out1)!="MIlocal"]
for(i in seq_along(p.exp1)){
  power[i,] <- apply(sim_out1[input1$p==p.exp1[i],colnames(sim_out1)!=("MIlocal")],2,sum)/nsim
}

# PLOT
png("./Figures/power1.png",width=450, height=400)
par(mar=c(4,4,1,1),oma=c(0,0,0,0))
plot(0,xlim=c(1,nrow(power)),ylim=c(0,1)
     ,type="n",axes=F,ann=F)
lines(x=1:nrow(power),y=power[,1],type="c",lty=2)
points(x=1:nrow(power),y=power[,1],pch=1)
lines(x=1:nrow(power),y=power[,4],type="c",lty=1)
points(x=1:nrow(power),y=power[,4],pch=19)
axis(2,las=1)
mtext("Power",2,line=2.5)
lab <- c(1,4,7,10,13,16,19)
axis(1,at=lab
     ,labels=p.exp1[lab])
mtext(expression(italic(rho)),1,line=2.5)
legend(x=1,y=.15,legend=c("Global MC","Decomposed MC"),pch=c(1,19))
dev.off()

# PLOT EXPERIMENT 2
power <- matrix(NA,ncol=(ncol(sim_out2)-1),nrow=length(p.exp2))
colnames(power) <- colnames(sim_out2)[colnames(sim_out2)!="MIlocal"]
for(i in seq_along(p.exp2)){
  power[i,] <- apply(sim_out2[input2$p1==p.exp2[i],colnames(sim_out2)!=("MIlocal")],2,sum)/nsim
}

# PLOT
png("./Figures/power2.png",width=450, height=400)
par(mar=c(4,4,1,1),oma=c(0,0,0,0))
plot(0,xlim=c(1,nrow(power)),ylim=c(0,1)
     ,type="n",axes=F,ann=F)
lines(x=1:nrow(power),y=power[,1],type="c",lty=2)
points(x=1:nrow(power),y=power[,1],pch=1)
lines(x=1:nrow(power),y=power[,4],type="c",lty=1)
points(x=1:nrow(power),y=power[,4],pch=19)
axis(2,las=1)
mtext("Power",2,line=2.5)
axis(1,at=1:10
     ,labels=p.exp2)
mtext(bquote("|"*italic("p")*"|"),1,line=2.5)
legend(x=1,y=.95,legend=c("Global MC","Decomposed MC"),pch=c(1,19))
dev.off()

power[round(p.exp2,10)==0.7,c(1,4)]








seq <- quantile(data$medianage,probs=c(0,1/9,2/9,3/9,4/9,5/9,6/9,7/9,8/9))
# mpas
png(paste0("./Figures/simtrue.png"),width=450, height=450)
par(oma=c(0,0,0,0),mar=c(1,0,0,0))
colcode <- sapply(data$medianage,function(x) sum(x >=seq))
plot(st_geometry(ger3),col=col[colcode])
range <- range(data$medianage)
mtext(bquote("Range = ["*.(format(range[1],digits=3,nsmall=3))*";"~
               .(format(range[2],digits=3,nsmall=3))*"]")
      ,side=1,line=-.5,cex=1.5)
dev.off()

