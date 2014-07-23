#######################################
#######################################
#### REANALYSE ALL DATA UNDER 600 M OF ELEVATION

DT <- read.csv("data.adult.csv")

###### calcule données BA et basup bainf
DT <- DT[!is.na(DT$PLHDAltitude) &
         DT$PLHDAltitude<601 &
         DT$dbh1>100 &
         !is.na(DT$dbh1) ,]

### keep only the 7 species used in SORTIE NZ T
DT.sel <- DT[DT$DMTRSpeciesCode %in%
         c("DACCUP","METUMB","NOTCLI",
           "NOTMEN","PODHAL","PRUFER",
           "WEIRAC"),]
DT.sel$DMTRSpeciesCode  <- factor(DT.sel$DMTRSpeciesCode)

Dave.dec.m <- DT.sel[,]


####################
###### Analysis
####################

library( likelihood )




##data
surv <- rep(1,length(Dave.dec.m$DMTRSpeciesCode))
surv[Dave.dec.m$dbh2==0] <- 0
Dave.dec.s <- data.frame(Dave.dec.m,surv=surv,year=Dave.dec.m$year2- Dave.dec.m$year1)
split.ns <- split(Dave.dec.s,  Dave.dec.s$DMTRSpeciesCode)

#################
##### ESTIMATION

## function
model.s1 <- function( max.s,x.0,x.b,diam2,year)
{
 pre <- (max.s*exp(-1/2*(log(diam2/x.0)/x.b)^2 ))^year
 return( pre)
}

###
list.Adult.S.1b <- vector("list",7)
for (i in 1:7)

{
print(i)
#windows()
data.temp <- split.ns[[i]][(split.ns[[i]]$dbh1/10)>10,]
surv<-( data.temp$surv)
diam2 <- data.temp$dbh1/10
year <- data.temp$year
datatemp <- data.frame(surv=surv,diam2=diam2,year=year)
var <- list( diam2 = "diam2",year="year")
par <- list(max.s = .9, x.0 = mean(datatemp$diam2),x.b=10)
par_lo <- list(max.s =0.000 , x.0 = 2 ,x.b =0.0001)
par_hi <- list(max.s = 1 , x.0 = 350 ,x.b= 80)
par_step <- list(max.s=1, x.0 = 3 ,x.b= 2)
var$x <- "surv"
## Mean in normal PDF
var$prob <- "predicted"
var$size <- 1
## Have it calculate log likelihood
var$log <- TRUE
list.Adult.S.1b[[i]] <- anneal (model = model.s1, par, var,
                                source_data= datatemp, par_lo,
                                par_hi,  pdf= dbinom, dep_var="surv",
                                max_iter= 50000)
title(main= names(split.ns)[i])
}



  ## param

  Matrix.Param.S.mALLb <- matrix(0,nrow=7,ncol=5)
  rownames(Matrix.Param.S.mALLb) <- names(split.ns)
  colnames(Matrix.Param.S.mALLb) <-  c("max.s","x.0","x.b","R²","num")

  for (i in 1:7)
  {
   Matrix.Param.S.mALLb[i,1] <-  list.Adult.S.1b[[i]]$best_pars[[1]]
   Matrix.Param.S.mALLb[i,2] <-  list.Adult.S.1b[[i]]$best_pars[[2]]
   Matrix.Param.S.mALLb[i,3] <-  list.Adult.S.1b[[i]]$best_pars[[3]]
   Matrix.Param.S.mALLb[i,4] <-  list.Adult.S.1b[[i]]$R2
   Matrix.Param.S.mALLb[i,5] <-  length(split.ns[[i]][(split.ns[[i]]$dbh1/10)>10,1])
  }


write.table( Matrix.Param.S.mALLb,file="param.s.all.NEW.txt")

Matrix.Param.S.mALLb <- read.table("param.s.all.NEW.txt")
                                       
#######################################
###### PLOTS

 x11()
   
   Dd <- 10:150
   i <- 1
    plot(Dd,1-model.s1( max.s= Matrix.Param.S.mALLb[i,1],x.0= Matrix.Param.S.mALLb[i,2],x.b= Matrix.Param.S.mALLb[i,3]
   ,diam2=Dd,year=1),xlab="DBH (cm)",ylab="Annual probability of mortality",ylim=c(0,0.03),type="l")

   i <- 3
  lines(Dd,1-model.s1( max.s= Matrix.Param.S.mALLb[i,1],x.0= Matrix.Param.S.mALLb[i,2],x.b= Matrix.Param.S.mALLb[i,3]
   ,diam2=Dd,year=1),lty=2,lwd=2)
   i <- 4
  lines(Dd,1-model.s1( max.s= Matrix.Param.S.mALLb[i,1],x.0= Matrix.Param.S.mALLb[i,2],x.b= Matrix.Param.S.mALLb[i,3]
   ,diam2=Dd,year=1),lty=3)
   i <- 5
  lines(Dd,1-model.s1( max.s= Matrix.Param.S.mALLb[i,1],x.0= Matrix.Param.S.mALLb[i,2],x.b= Matrix.Param.S.mALLb[i,3]
   ,diam2=Dd,year=1),lty=4,lwd=2)
   i <- 6

  lines(Dd,1-model.s1( max.s= Matrix.Param.S.mALLb[i,1],x.0= Matrix.Param.S.mALLb[i,2],x.b= Matrix.Param.S.mALLb[i,3]
   ,diam2=Dd,year=1),lty=5,lwd=1)
   i <- 7
  lines(Dd,1-model.s1( max.s= Matrix.Param.S.mALLb[i,1],x.0= Matrix.Param.S.mALLb[i,2],x.b= Matrix.Param.S.mALLb[i,3]
   ,diam2=Dd,year=1),lty=1,lwd=2)
      i <- 2
  lines(Dd,1-model.s1( max.s= Matrix.Param.S.mALLb[i,1],x.0= Matrix.Param.S.mALLb[i,2],x.b= Matrix.Param.S.mALLb[i,3]
   ,diam2=Dd,year=1),lty=6,lwd=1)
  
  vec.lwd <- c(1,2,1,2,1,2,1)
  vec.lty <- c(1:5,1,6)
  
  legend(20, 0.029,legend=names(split.ns)[c( 1,3:7,2)],lty=vec.lty,lwd=vec.lwd,bty="n")
  
                                       
                                       
                                       
#######################################
#######################################
#######################################
#### for sapling survival

DT <- read.csv("data.adult.csv")

DT <- DT[!is.na(DT$PLHDAltitude) &
         DT$PLHDAltitude<601 , ]

### keep only the 7 species used in SORTIE NZ T
DT.sel <- DT[DT$DMTRSpeciesCode %in%
         c("DACCUP","METUMB","NOTCLI",
           "NOTMEN","PODHAL","PRUFER",
           "WEIRAC"),]
DT.sel$DMTRSpeciesCode  <- factor(DT.sel$DMTRSpeciesCode)

Dave.dec.m <- DT.sel[,]


####################
###### Analysis
####################

library( likelihood )

##data
surv <- rep(1,length(Dave.dec.m$DMTRSpeciesCode))
surv[Dave.dec.m$dbh2==0] <- 0
Dave.dec.s <- data.frame(Dave.dec.m,surv=surv,year=Dave.dec.m$year2- Dave.dec.m$year1)
split.ns <- split(Dave.dec.s,  Dave.dec.s$DMTRSpeciesCode)

load(file="plot.pb.tot.Rdata")

Dave.dec.s2 <- Dave.dec.s[!(Dave.dec.s$pol.idi.y %in% plot.pb.tot),]
## Exclude plot with to high mortalty rate (100%) probably due to a disturbance.

split.ns2 <- split(Dave.dec.s2,  Dave.dec.s2$DMTRSpeciesCode)




##### annual mortality rate for the small tree

model.s0 <- function( max.s,year)
{
 pre <- (max.s)^year
 return( pre)
}


list.Adult.S.0.inf10 <- vector("list",7)
for (i in 1:7)

{
#x11()
data.temp <- split.ns2[[i]][(split.ns2[[i]]$dbh1/10)<10 & (split.ns2[[i]]$dbh1/10)>5,]

surv<-( data.temp$surv)
diam2 <- data.temp$dbh1/10
year <- data.temp$year
datatemp <- data.frame(surv=surv,diam2=diam2,year=year)

var <-  list(year="year")
par <- list(max.s = .9)
par_lo <- list(max.s =0.000)
par_hi <- list(max.s = 1)
var$x <- "surv"
## Mean in normal PDF
var$prob <- "predicted"
var$size <- 1
## Have it calculate log likelihood
var$log <- TRUE
list.Adult.S.0.inf10[[i]] <- anneal (model = model.s0, par, var,
                                     source_data= datatemp, par_lo, par_hi,
                                     pdf= dbinom, dep_var="surv", max_iter= 50000)
print(i)
}

list.Adult.S.0.inf10b <- vector("list",7)
for (i in 1:7)

{
#x11()
data.temp <- split.ns[[i]][(split.ns[[i]]$dbh1/10)<10 & (split.ns[[i]]$dbh1/10)>5,]

surv<-( data.temp$surv)
diam2 <- data.temp$dbh1/10
year <- data.temp$year
datatemp <- data.frame(surv=surv,diam2=diam2,year=year)

var <-  list(year="year")
par <- list(max.s = .9)
par_lo <- list(max.s =0.000)
par_hi <- list(max.s = 1)
var$x <- "surv"
## Mean in normal PDF
var$prob <- "predicted"
var$size <- 1
## Have it calculate log likelihood
var$log <- TRUE
list.Adult.S.0.inf10b[[i]] <- anneal(model = model.s0, par, var,
                                     source_data= datatemp, par_lo, par_hi,
                                     pdf= dbinom, dep_var="surv", max_iter= 50000)
print(i)
}


Matrix.Param.S.SAP <- matrix(0,nrow=7,ncol=3)
rownames(Matrix.Param.S.SAP) <- names(split.ns)
colnames(Matrix.Param.S.SAP) <-  c("max.s","R²","Repl.")

for (i in 1:7)

{
 Matrix.Param.S.SAP[i,1] <-  1- list.Adult.S.0.inf10[[i]]$best_pars[[1]]
 Matrix.Param.S.SAP[i,2] <- list.Adult.S.0.inf10[[i]]$R2
 Matrix.Param.S.SAP[i,3] <- length(list.Adult.S.0.inf10[[i]]$source_data[,1])
}

 write.table( Matrix.Param.S.SAP,file="param.sap.NEW.txt")

Matrix.Param.S.SAPb <- matrix(0,nrow=7,ncol=3)
rownames(Matrix.Param.S.SAPb) <- names(split.ns)
colnames(Matrix.Param.S.SAPb) <-  c("max.s","R²","Repl.")

for (i in 1:7)

{
 Matrix.Param.S.SAPb[i,1] <-  1- list.Adult.S.0.inf10b[[i]]$best_pars[[1]]
 Matrix.Param.S.SAPb[i,2] <- list.Adult.S.0.inf10b[[i]]$R2
 Matrix.Param.S.SAPb[i,3] <- length(list.Adult.S.0.inf10b[[i]]$source_data[,1])
}


