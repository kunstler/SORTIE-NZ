#################################################
#################################################
#################################################
###### ADULT GROWTH ESTIMATION FOR SORTIE NZ
#################################################
#################################################
#################################################

#############
### READ DATA
Dave.m.ss.nor <- read.csv("data.adult.csv")

Dave.dec <- Dave.m.ss.nor[!is.na(Dave.m.ss.nor$PLHDAltitude) &
                          Dave.m.ss.nor$PLHDAltitude<601,]
Dave.dec <- Dave.dec[Dave.dec$dbh1>00,]

Dave.mdd2<- Dave.dec[Dave.dec$DMTRSpeciesCode %in% c("DACCUP","METUMB","NOTCLI",
                                                     "NOTMEN","PODHAL","PRUFER",
                                                     "WEIRAC"),]  
### keep only the 7 species used in SORTIE NZ T
Dave.mdd2$DMTRSpeciesCode  <- factor(Dave.mdd2$DMTRSpeciesCode)
Dave.md2 <- Dave.mdd2[(Dave.mdd2$PLHDPlotAreaRadius)==4000,]
### select only the plot with the same size

###### calcule données BA greater than the target tree
vecbasup <- rep(0,length(Dave.md2$dbh1))

for ( i in 1:10000)
{
   vec.tempp <- Dave.dec$dbh1[Dave.md2$pol.idi.y[i]==Dave.dec$pol.idi.y]
   vec.tempp <- as.vector(na.exclude(vec.tempp))
   vecbasup[i] <- sum(pi*(vec.tempp[vec.tempp>Dave.md2$dbh1[i]]/20)^2)
}

for ( i in 10001:length(Dave.md2$dbh1))
{
   vec.tempp <- Dave.dec$dbh1[Dave.md2$pol.idi.y[i]==Dave.dec$pol.idi.y]
   vec.tempp <- as.vector(na.exclude(vec.tempp))
   vecbasup[i] <- sum(pi*(vec.tempp[vec.tempp>Dave.md2$dbh1[i]]/20)^2)
}

## compute ba tot per plot
basal.plot <- tapply(pi*(Dave.dec$dbh1/20)^2, INDEX= Dave.dec$pol.idi.y,FUN=sum)
## basal area plot total
basal.plot.10 <- tapply(pi*(Dave.dec$dbh1[Dave.dec$dbh1>100]/20)^2,
                        INDEX= Dave.dec$pol.idi.y[Dave.dec$dbh1>100],FUN=sum)
### basala area total for tree greater than 10cm of DBH
data.basal.plot <- data.frame(pol.idi.y = names(basal.plot),
                              BA.p= basal.plot/1000,
                              BA.p.10 = basal.plot.10/1000 )
## data basal area per plot
## add BA sup to data set
Dave.md2T <- data.frame(Dave.md2,bas=vecbasup/1000)
### merge data
Dave.md2T2 <-  merge(Dave.md2T, data.basal.plot,by="pol.idi.y")
## compute BA tot - BA of the target tree
Dave.md2T2$BA.p.10.B <- Dave.md2T2$BA.p.10 - pi*(Dave.md2T2$dbh1/20)^2/1000

### data for estimation
DATA.GG <- Dave.md2T2[Dave.md2T2$dbh2 != 0,]   ## keep only tree alive
DATA.GG$pol.idi.y <- factor(DATA.GG$pol.idi.y)

### delete 3.620183 % of the data  with a minimum growth of 0.5 mm per year
DATA.GG <- data.frame(DATA.GG,g = (DATA.GG$dbh2-DATA.GG$dbh1)/(DATA.GG$year2-DATA.GG$year1))
DATA.GG <- DATA.GG[DATA.GG$g > -0.05 & DATA.GG$g < 40 ,]
DATA.GG <- DATA.GG[DATA.GG$dbh1>100 ,]    ## keep only tree with dbh > 10 cm
split.gg <-  split(DATA.GG ,DATA.GG$DMTRSpeciesCode)

####################
####################
####################
###### Analysis
####################
####################
####################



library( neighparam )


###################################
###################################
## model with competition  BAsup
###################################
###################################

### function
  model3b <- function( max.g,x.0,x.b,c,DD,diam2,BAsup)
{

 BAC <- (BAsup)
 rr <- max.g*exp(-1/2*(log(diam2/((x.0)))/((x.b)))^2 )*exp(-c*(BAC)^DD)
 return( rr)
}

# estimation
for (i in 1:7)

{
#windows()
datatemmp <- split.gg[[i]][(split.gg[[i]]$dbh1/10)>10,]
growth<-(datatemmp$g)
diam2 <- datatemmp$dbh1/10
BAsup <- datatemmp$bas
datatemp <- data.frame(growth=growth,diam2=diam2,BAsup=BAsup)

par <- list(max.g = (max(datatemp$growth)), x.0 = 250,x.b=2,c=0.001,
            DD=2, diam2 = "diam2", BAsup = "BAsup")
par_lo <- list(max.g =0.0001 , x.0 = 5 ,x.b =0.0001,c=0,DD=0,sd=0)
par_hi <- list(max.g = 20 , x.0 = 600 ,x.b= 50,c=10,DD=5,sd=1000)
par_step <- list(max.g =1, x.0 = 3 ,x.b= 2,c=2,DD=1,sd=1)


par$x <- "growth"
## Mean in normal PDF
par$mean <- "predicted"
par$sd <- 2
## Have it calculate log likelihood
par$log <- TRUE
list.Adult.G.3b.t <- neighanneal(model = model3b, par, source_data= datatemp,
                                 par_lo,par_hi, par_step, pdf=dnorm,
                                 dep_var="growth", max_iter= 50000)
save( list.Adult.G.3b.t ,file=paste("list.Adult.G.3b.t" ,i,sep=""))

print(i)
}

list.Adult.G.3bL <- vector("list",8)
for (i in 1:7)
{
load(file=paste("list.Adult.G.3b.t" ,i,sep=""))
list.Adult.G.3bL[[i]] <- list.Adult.G.3b.t
}

################################
################################
### MODEL WITH NO COMPETITION
#### estimation
################################
################################

## function
 model1 <- function( max.g,x.0,x.b,diam2)
{
 pre <- max.g*exp(-1/2*(log(diam2/x.0)/x.b)^2 )
 return( pre)
}


list.Adult.G.1 <- vector("list",8)
for (i in 1:7)

{
datatemmp <- split.gg[[i]][(split.gg[[i]]$dbh1/10)>10,]
growth<-(datatemmp$g)
diam2 <- datatemmp$dbh1/10
datatemp <- data.frame(growth=growth,diam2=diam2)

par <- list(max.g = (max(datatemp$growth)), x.0 = 250,x.b=10, diam2 = "diam2")
par_lo <- list(max.g =0.0001 , x.0 = 5 ,x.b =0.0001,sd=0)
par_hi <- list(max.g = 30 , x.0 = 500 ,x.b= 50,sd=1000)
par_step <- list(max.g =1, x.0 = 3 ,x.b= 2,sd=1)
par$x <- "growth"
## Mean in normal PDF
par$mean <- "predicted"
par$sd <- 2
## Have it calculate log likelihood
par$log <- TRUE
list.Adult.G.1[[i]] <- neighanneal (model = model1, par, source_data= datatemp,
par_lo, par_hi, par_step, pdf=dnorm, dep_var="growth", max_iter= 50000)
title(main= names(split.gg)[i])
}


###############
###############
### AICc TABLE

###################
####################
######### AIC table

res.sim.an <- matrix(0, nrow=7,ncol=3)
rownames(res.sim.an) <- names(split.gg)
colnames(res.sim.an) <- c("num","AIC.No.BA","AIC.BAsup")
sum.like.1 <- 0
sum.like.2 <- 0
NUMO <- 0
for (i in 1:7)

{
sum.like.1 <- sum.like.1+list.Adult.G.1[[i]]$max_likeli
sum.like.2 <- sum.like.2+list.Adult.G.3bL[[i]]$max_likeli
NUMO      <- NUMO+length(split.gg[[i]]$g)
res.sim.an[i,1] <- length(split.gg[[i]]$g)
res.sim.an[i,2] <- list.Adult.G.1[[i]]$aic_corr
res.sim.an[i,3] <- list.Adult.G.3bL[[i]]$aic_corr
}

########################
########################
########################
##### MATRIX PARAM


#### matrix parameters
Matrix.Param.G.m3b <- matrix(0,nrow=7,ncol=8)
rownames(Matrix.Param.G.m3b) <- names(split.gg)
colnames(Matrix.Param.G.m3b) <-  c("max.g","x.0","x.b","c","D","sd","R²","Repl.")

for (i in 1:7)
{
 Matrix.Param.G.m3b[i,1] <-  list.Adult.G.3bL[[i]]$best_pars[[1]]
 Matrix.Param.G.m3b[i,2] <-  list.Adult.G.3bL[[i]]$best_pars[[2]]
 Matrix.Param.G.m3b[i,3] <-  list.Adult.G.3bL[[i]]$best_pars[[3]]
 Matrix.Param.G.m3b[i,4] <-  list.Adult.G.3bL[[i]]$best_pars[[4]]
 Matrix.Param.G.m3b[i,5] <-  list.Adult.G.3bL[[i]]$best_pars[[5]]
 Matrix.Param.G.m3b[i,6] <-  list.Adult.G.3bL[[i]]$best_pars$sd
 Matrix.Param.G.m3b[i,7] <- list.Adult.G.3bL[[i]]$R2
 Matrix.Param.G.m3b[i,8] <- length(list.Adult.G.3bL[[i]]$source_data[,1])
}

write.table(cbind(Matrix.Param.G.m3b,res.sim.an[,2:3]),
            file="matrix.param.growth.def.txt")

############
## plots 

Matrix.Param.G.m3b <- read.table("matrix.param.growth.def.txt")

x11()
vec.lwd <- c(1,2,1,2,1,2,1)
vec.lty <- c(1:5,1,6)

par (mfrow=c(1,2))
   i <- 1
   pred.H <- model3b( max.g= Matrix.Param.G.m3b[i,1],x.0= Matrix.Param.G.m3b[i,2],x.b= Matrix.Param.G.m3b[i,3]
   ,c= Matrix.Param.G.m3b[i,4],DD= Matrix.Param.G.m3b[i,5],diam2=10:150,BAsup=2)
  
   plot(10:150,pred.H,xlab="DBH (cm)",ylab="Growth (mm/yr)",ylim=c(0,4),type="l")
   i <- 3
   pred.H <- model3b( max.g= Matrix.Param.G.m3b[i,1],x.0= Matrix.Param.G.m3b[i,2],x.b= Matrix.Param.G.m3b[i,3]
   ,c= Matrix.Param.G.m3b[i,4],DD= Matrix.Param.G.m3b[i,5],diam2=10:150,BAsup=2)
  
  lines(10:150,pred.H,lty=2,lwd=2)
   i <- 4
   pred.H <- model3b( max.g= Matrix.Param.G.m3b[i,1],x.0= Matrix.Param.G.m3b[i,2],x.b= Matrix.Param.G.m3b[i,3]
   ,c= Matrix.Param.G.m3b[i,4],DD= Matrix.Param.G.m3b[i,5],diam2=10:150,BAsup=2)
  
  lines(10:150,pred.H,lty=3)
   i <- 5
   pred.H <- model3b( max.g= Matrix.Param.G.m3b[i,1],x.0= Matrix.Param.G.m3b[i,2],x.b= Matrix.Param.G.m3b[i,3]
   ,c= Matrix.Param.G.m3b[i,4],DD= Matrix.Param.G.m3b[i,5],diam2=10:150,BAsup=2)
  
  lines(10:150,pred.H,lty=4,lwd=2)
   i <- 6
   pred.H <- model3b( max.g= Matrix.Param.G.m3b[i,1],x.0= Matrix.Param.G.m3b[i,2],x.b= Matrix.Param.G.m3b[i,3]
   ,c= Matrix.Param.G.m3b[i,4],DD= Matrix.Param.G.m3b[i,5],diam2=10:150,BAsup=2)
  
  lines(10:150,pred.H,lty=5,lwd=1)
   i <- 7
   pred.H <- model3b( max.g= Matrix.Param.G.m3b[i,1],x.0= Matrix.Param.G.m3b[i,2],x.b= Matrix.Param.G.m3b[i,3]
   ,c= Matrix.Param.G.m3b[i,4],DD= Matrix.Param.G.m3b[i,5],diam2=10:150,BAsup=2)
  
  lines(10:150,pred.H,lty=1,lwd=2)

     i <- 2
   pred.H <- model3b( max.g= Matrix.Param.G.m3b[i,1],x.0= Matrix.Param.G.m3b[i,2],x.b= Matrix.Param.G.m3b[i,3]
   ,c= Matrix.Param.G.m3b[i,4],DD= Matrix.Param.G.m3b[i,5],diam2=10:150,BAsup=2)
  
  lines(10:150,pred.H,lty=6,lwd=1)

text(locator(1),labels= "(a)")

 i <- 1
 pred.H <- model3b( max.g= Matrix.Param.G.m3b[i,1],x.0= Matrix.Param.G.m3b[i,2],x.b= Matrix.Param.G.m3b[i,3]
 ,c= Matrix.Param.G.m3b[i,4],DD= Matrix.Param.G.m3b[i,5],diam2=10:150,BAsup=100)

 plot(10:150,pred.H,xlab="DBH (cm)",ylab="Growth (mm/yr)",ylim=c(0,4),type="l")
 i <- 3
 pred.H <- model3b( max.g= Matrix.Param.G.m3b[i,1],x.0= Matrix.Param.G.m3b[i,2],x.b= Matrix.Param.G.m3b[i,3]
 ,c= Matrix.Param.G.m3b[i,4],DD= Matrix.Param.G.m3b[i,5],diam2=10:150,BAsup=100)

lines(10:150,pred.H,lty=2,lwd=2)
 i <- 4
 pred.H <- model3b( max.g= Matrix.Param.G.m3b[i,1],x.0= Matrix.Param.G.m3b[i,2],x.b= Matrix.Param.G.m3b[i,3]
 ,c= Matrix.Param.G.m3b[i,4],DD= Matrix.Param.G.m3b[i,5],diam2=10:150,BAsup=100)

lines(10:150,pred.H,lty=3)
 i <- 5
 pred.H <- model3b( max.g= Matrix.Param.G.m3b[i,1],x.0= Matrix.Param.G.m3b[i,2],x.b= Matrix.Param.G.m3b[i,3]
 ,c= Matrix.Param.G.m3b[i,4],DD= Matrix.Param.G.m3b[i,5],diam2=10:150,BAsup=100)

lines(10:150,pred.H,lty=4,lwd=2)
 i <- 6
 pred.H <- model3b( max.g= Matrix.Param.G.m3b[i,1],x.0= Matrix.Param.G.m3b[i,2],x.b= Matrix.Param.G.m3b[i,3]
 ,c= Matrix.Param.G.m3b[i,4],DD= Matrix.Param.G.m3b[i,5],diam2=10:150,BAsup=100)

lines(10:150,pred.H,lty=5,lwd=1)
 i <- 7
 pred.H <- model3b( max.g= Matrix.Param.G.m3b[i,1],x.0= Matrix.Param.G.m3b[i,2],x.b= Matrix.Param.G.m3b[i,3]
 ,c= Matrix.Param.G.m3b[i,4],DD= Matrix.Param.G.m3b[i,5],diam2=10:150,BAsup=100)

lines(10:150,pred.H,lty=1,lwd=2)

 i <- 2
 pred.H <- model3b( max.g= Matrix.Param.G.m3b[i,1],x.0= Matrix.Param.G.m3b[i,2],x.b= Matrix.Param.G.m3b[i,3]
 ,c= Matrix.Param.G.m3b[i,4],DD= Matrix.Param.G.m3b[i,5],diam2=10:150,BAsup=100)

lines(10:150,pred.H,lty=6,lwd=1)
text(locator(1),labels= "(b)")
legend(locator(1),legend=c("DACCUP" ,"NOTCLI" ,"NOTMEN" ,"PODHAL" ,"PRUFER" ,"WEIRAC","METUMB")
,lty=vec.lty,lwd=vec.lwd,bty="n",ncol=2, cex=0.5)





 
x11()
vec.lwd <- c(1,2,1,2,1,2,1)
vec.lty <- c(1:5,1,6)

par (mfrow=c(1,2))
   i <- 1
   pred.H <- model3b( max.g= Matrix.Param.G.m3b[i,1],x.0= Matrix.Param.G.m3b[i,2],x.b= Matrix.Param.G.m3b[i,3]
   ,c= Matrix.Param.G.m3b[i,4],DD= Matrix.Param.G.m3b[i,5],diam2=10:150,BAsup=2)
  
   plot(10:150,pred.H,xlab="DBH (cm)",ylab="Growth (mm/yr)",ylim=c(0,4),type="l")
   i <- 3
   pred.H <- model3b( max.g= Matrix.Param.G.m3b[i,1],x.0= Matrix.Param.G.m3b[i,2],x.b= Matrix.Param.G.m3b[i,3]
   ,c= Matrix.Param.G.m3b[i,4],DD= Matrix.Param.G.m3b[i,5],diam2=10:150,BAsup=2)
  
  lines(10:150,pred.H,lty=2,lwd=2)
   i <- 4
   pred.H <- model3b( max.g= Matrix.Param.G.m3b[i,1],x.0= Matrix.Param.G.m3b[i,2],x.b= Matrix.Param.G.m3b[i,3]
   ,c= Matrix.Param.G.m3b[i,4],DD= Matrix.Param.G.m3b[i,5],diam2=10:150,BAsup=2)
  
  lines(10:150,pred.H,lty=3)
   i <- 5
   pred.H <- model3b( max.g= Matrix.Param.G.m3b[i,1],x.0= Matrix.Param.G.m3b[i,2],x.b= Matrix.Param.G.m3b[i,3]
   ,c= Matrix.Param.G.m3b[i,4],DD= Matrix.Param.G.m3b[i,5],diam2=10:150,BAsup=2)
  
  lines(10:150,pred.H,lty=4,lwd=2)
   i <- 6
   pred.H <- model3b( max.g= Matrix.Param.G.m3b[i,1],x.0= Matrix.Param.G.m3b[i,2],x.b= Matrix.Param.G.m3b[i,3]
   ,c= Matrix.Param.G.m3b[i,4],DD= Matrix.Param.G.m3b[i,5],diam2=10:150,BAsup=2)
  
  lines(10:150,pred.H,lty=5,lwd=1)
   i <- 7
   pred.H <- model3b( max.g= Matrix.Param.G.m3b[i,1],x.0= Matrix.Param.G.m3b[i,2],x.b= Matrix.Param.G.m3b[i,3]
   ,c= Matrix.Param.G.m3b[i,4],DD= Matrix.Param.G.m3b[i,5],diam2=10:150,BAsup=2)
  
  lines(10:150,pred.H,lty=1,lwd=2)

     i <- 2
   pred.H <- model3b( max.g= Matrix.Param.G.m3b[i,1],x.0= Matrix.Param.G.m3b[i,2],x.b= Matrix.Param.G.m3b[i,3]
   ,c= Matrix.Param.G.m3b[i,4],DD= Matrix.Param.G.m3b[i,5],diam2=10:150,BAsup=2)
  
  lines(10:150,pred.H,lty=6,lwd=1)

text(locator(1),labels= "(a)")

 i <- 1
 pred.H <- model3b( max.g= Matrix.Param.G.m3b[i,1],x.0= Matrix.Param.G.m3b[i,2],x.b= Matrix.Param.G.m3b[i,3]
 ,c= Matrix.Param.G.m3b[i,4],DD= Matrix.Param.G.m3b[i,5],diam2=20,BAsup=1:100)

 plot(1:100,pred.H,xlab="DBH (cm)",ylab="Growth (mm/yr)",ylim=c(0,4),type="l")
 i <- 3
 pred.H <- model3b( max.g= Matrix.Param.G.m3b[i,1],x.0= Matrix.Param.G.m3b[i,2],x.b= Matrix.Param.G.m3b[i,3]
 ,c= Matrix.Param.G.m3b[i,4],DD= Matrix.Param.G.m3b[i,5],diam2=20,BAsup=1:100)

lines(1:100,pred.H,lty=2,lwd=2)
 i <- 4
 pred.H <- model3b( max.g= Matrix.Param.G.m3b[i,1],x.0= Matrix.Param.G.m3b[i,2],x.b= Matrix.Param.G.m3b[i,3]
 ,c= Matrix.Param.G.m3b[i,4],DD= Matrix.Param.G.m3b[i,5],diam2=20,BAsup=1:100)

lines(1:100,pred.H,lty=3)
 i <- 5
 pred.H <- model3b( max.g= Matrix.Param.G.m3b[i,1],x.0= Matrix.Param.G.m3b[i,2],x.b= Matrix.Param.G.m3b[i,3]
 ,c= Matrix.Param.G.m3b[i,4],DD= Matrix.Param.G.m3b[i,5],diam2=20,BAsup=1:100)

lines(1:100,pred.H,lty=4,lwd=2)
 i <- 6
 pred.H <- model3b( max.g= Matrix.Param.G.m3b[i,1],x.0= Matrix.Param.G.m3b[i,2],x.b= Matrix.Param.G.m3b[i,3]
 ,c= Matrix.Param.G.m3b[i,4],DD= Matrix.Param.G.m3b[i,5],diam2=20,BAsup=1:100)

lines(1:100,pred.H,lty=5,lwd=1)
 i <- 7
 pred.H <- model3b( max.g= Matrix.Param.G.m3b[i,1],x.0= Matrix.Param.G.m3b[i,2],x.b= Matrix.Param.G.m3b[i,3]
 ,c= Matrix.Param.G.m3b[i,4],DD= Matrix.Param.G.m3b[i,5],diam2=20,BAsup=1:100)

lines(1:100,pred.H,lty=1,lwd=2)

 i <- 2
 pred.H <- model3b( max.g= Matrix.Param.G.m3b[i,1],x.0= Matrix.Param.G.m3b[i,2],x.b= Matrix.Param.G.m3b[i,3]
 ,c= Matrix.Param.G.m3b[i,4],DD= Matrix.Param.G.m3b[i,5],diam2=20,BAsup=1:100)

lines(1:100,pred.H,lty=6,lwd=1)
text(locator(1),labels= "(b)")
legend(locator(1),legend=c("DACCUP" ,"NOTCLI" ,"NOTMEN" ,"PODHAL" ,"PRUFER" ,"WEIRAC","METUMB")
,lty=vec.lty,lwd=vec.lwd,bty="n",ncol=2, cex=0.5)

