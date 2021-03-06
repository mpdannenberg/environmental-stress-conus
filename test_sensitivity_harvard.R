## Estimate sensitivity of A. rubrum, Q. rubra, and T. canadensis to DBH uncertainty

## Load necessary functions and parameters
library(dplR)
library(ggplot2)
library(extrafont)

build.dbh.table.lyf = function(rwl, info){
  n = dim(rwl)[2]
  dbh.table = data.frame(matrix(ncol=3, nrow=n))
  names(dbh.table) <- c("Core", "Est", "Obs")
  dbh.table$Core <- names(rwl)
  for(i in 1:n){
    x = rwl[,i]
    est = sum(x[!is.na(x)]) # cm
    name = dbh.table$Core[i]
    site = substr(name, 1, 3)
    tree = as.numeric(substr(name, 4, 5))
    ind = which(info$Site %in% site & info$Tree==tree)
    try({
      obs = info$DBH[ind]
      dbh.table$Est[i] = est
      dbh.table$Obs[i] = obs})
  }
  return(dbh.table)
}

build.dbh.table.tow = function(rwl, info){
  n = dim(rwl)[2]
  dbh.table = data.frame(matrix(ncol=3, nrow=n))
  names(dbh.table) <- c("Core", "Est", "Obs")
  dbh.table$Core <- names(rwl)
  for(i in 1:n){
    x = rwl[,i]
    est = sum(x[!is.na(x)]) # cm
    name = dbh.table$Core[i]
    site = substr(name, 1, 3)
    tree = as.numeric(substr(name, 4, 6))
    ind = which(info$Site %in% site & info$Tree==tree)
    try({
      obs = info$DBH[ind]
      dbh.table$Est[i] = est
      dbh.table$Obs[i] = obs})
  }
  return(dbh.table)
}

comp.dopt = function(rwl, info, pars){
  noObs = info$Core[is.na(info$Obs)]
  rwl = rwl[, !(names(rwl) %in% noObs)]
  info = info[!(info$Core %in% noObs), ]
  
  G = pars$G
  Dmax = pars$Dmax
  Hmax = pars$Hmax
  
  if(is.na(pars$b2)) b2 = 2 * ((Hmax-137)/Dmax)
  else b2 = pars$b2
  
  if(is.na(pars$b3)) b3 = (Hmax-137)/Dmax^2
  else b3 = pars$b3
  
  n = dim(rwl)[2]
  d.naive = data.frame(matrix(ncol=n, nrow=dim(rwl)[1]), row.names=row.names(rwl))
  d.act = data.frame(matrix(ncol=n, nrow=dim(rwl)[1]), row.names=row.names(rwl))
  dopt.naive = data.frame(matrix(ncol=n, nrow=dim(rwl)[1]), row.names=row.names(rwl))
  dopt.act = data.frame(matrix(ncol=n, nrow=dim(rwl)[1]), row.names=row.names(rwl))
  names(d.naive) = names(rwl)
  names(d.act) = names(rwl)
  names(dopt.naive) = names(rwl)
  names(dopt.act) = names(rwl)
  
  for(i in 1:n){
    inds = which(!is.na(rwl[,i]))
    d.naive[inds, i] = cumsum(rwl[inds,i]) - rwl[inds,i]
    d.act[inds, i] = rev(info$Obs[i] - cumsum(rev(rwl[inds,i])))
    if(sum(na.omit(d.act[,i]) < 0) > 0){
      d.act[,i] = max(na.omit(d.act[,i])) * (d.act[,i] - min(na.omit(d.act[,i]))) / (max(na.omit(d.act[,i])) - min(na.omit(d.act[,i])))
    }
    
    D1 = d.naive[inds, i]
    D2 = d.act[inds, i]
    H1 = 137 + b2*D1 - b3*D1^2
    H2 = 137 + b2*D2 - b3*D2^2
    dopt.naive[inds, i] = G*D1*(1 - (D1*H1)/(Dmax*Hmax)) / (274 + 3*b2*D1 - 4*b3*D1^2)
    dopt.act[inds, i] = G*D2*(1 - (D2*H2)/(Dmax*Hmax)) / (274 + 3*b2*D2 - 4*b3*D2^2)
  }
  
  return(list(d.naive, d.act, dopt.naive, dopt.act, rwl))
}

source("StressIndexFunctions.R")

info.lyf <- read.table('./data/lyf.field.data.csv', header=T, sep=',')
info.tow <- read.table('./data/tow.field.data.csv', header=T, sep=',')


##### TSCA #####

## Optimal growth parameters
pars <- vector(mode="list", length=5)
names(pars) <- c("G", "Dmax", "Hmax", "b2", "b3")
pars$G = 82.9
pars$Dmax = 210
pars$Hmax = 5180
pars$b2 = NA
pars$b3 = NA

## Read in tree-ring data
trw.lyf = read.rwl('./data/LF_TSCA.rwl') * 0.1 * 2 #Multiply by 2 (radial growth --> diameter growth) and 0.1 (mm --> cm)
trw.tow = read.rwl('./data/TP_TSCA.rwl') * 0.1 * 2 #Multiply by 2 (radial growth --> diameter growth) and 0.1 (mm --> cm)

## Make data frame of estimated and observed DBH
trw.dbh.lyf <- build.dbh.table.lyf(trw.lyf, info.lyf)
trw.dbh.tow <- build.dbh.table.tow(trw.tow, info.tow)
write.csv(rbind(trw.dbh.lyf, trw.dbh.tow), "./data/DBH_Comparison_TSCA.csv")

## Calculate actual and naive D and Dopt
trw.comp.dopt.lyf = comp.dopt(trw.lyf, trw.dbh.lyf, pars)
temp.comp.dopt.lyf = na.omit(cbind(c(apply(trw.comp.dopt.lyf[[1]], 2, rbind)), 
                                    c(apply(trw.comp.dopt.lyf[[2]], 2, rbind)), 
                                    c(apply(trw.comp.dopt.lyf[[3]], 2, rbind)), 
                                    c(apply(trw.comp.dopt.lyf[[4]], 2, rbind)), 
                                    c(apply(trw.comp.dopt.lyf[[5]], 2, rbind))))
all.dopt.comp.lyf = data.frame(temp.comp.dopt.lyf)
names(all.dopt.comp.lyf) = c("d.naive", "d.act", "dopt.naive", "dopt.act", "dD")
all.dopt.comp.lyf$dopt.dif = all.dopt.comp.lyf$dopt.naive - all.dopt.comp.lyf$dopt.act

trw.comp.dopt.tow = comp.dopt(trw.tow, trw.dbh.tow, pars)
temp.comp.dopt.tow = na.omit(cbind(c(apply(trw.comp.dopt.tow[[1]], 2, rbind)), 
                                    c(apply(trw.comp.dopt.tow[[2]], 2, rbind)), 
                                    c(apply(trw.comp.dopt.tow[[3]], 2, rbind)), 
                                    c(apply(trw.comp.dopt.tow[[4]], 2, rbind)), 
                                    c(apply(trw.comp.dopt.tow[[5]], 2, rbind))))
all.dopt.comp.tow = data.frame(temp.comp.dopt.tow)
names(all.dopt.comp.tow) = c("d.naive", "d.act", "dopt.naive", "dopt.act", "dD")
all.dopt.comp.tow$dopt.dif = all.dopt.comp.tow$dopt.naive - all.dopt.comp.tow$dopt.act

all.dopt.comp <- rbind(all.dopt.comp.lyf, all.dopt.comp.tow)

## Calculate Stress 
dD = all.dopt.comp$dD
dDopt.naive = all.dopt.comp$dopt.naive
dDopt.act = all.dopt.comp$dopt.act
D.naive = all.dopt.comp$d.naive
D.act = all.dopt.comp$d.act

all.dopt.comp$s.naive.Dbased = all.dopt.comp$dD / all.dopt.comp$dopt.naive
all.dopt.comp$s.act.Dbased = all.dopt.comp$dD / all.dopt.comp$dopt.act
all.dopt.comp$s.dif.Dbased = all.dopt.comp$s.naive.Dbased - all.dopt.comp$s.act.Dbased
rm(list = c("dD", "dDopt.naive", "dDopt.act", "D.naive", "D.act"))

## Write CSVs for use in MATLAB
all.dopt.comp$d.naive.bins = cut(all.dopt.comp$d.naive, seq(from=0, to=150, by=10), include.lowest=TRUE)

dopt.summary = boxplot(dopt.dif~d.naive.bins, data=all.dopt.comp, plot=FALSE)
write.csv(cbind(dopt.summary$names, dopt.summary$n, dopt.summary$stats[3,]), "./data/Dopt_Comparison_TSCA.csv")

sD.summary = boxplot(s.dif.Dbased~d.naive.bins, data=all.dopt.comp, plot=FALSE)
write.csv(cbind(sD.summary$names, sD.summary$n, sD.summary$stats[3,]), "./data/Stress_Comparison_TSCA.csv")

all.dopt.comp$d.naive.bins = NULL
write.csv(all.dopt.comp, "./data/Sensitivity_TSCA.csv")

rm(list = c(ls(pattern='^all'), ls(pattern='^temp'), ls(pattern='^trw'), 'dopt.summary','pars','sD.summary'))


##### QURU #####

## Optimal growth parameters
pars <- vector(mode="list", length=5)
names(pars) <- c("G", "Dmax", "Hmax", "b2", "b3")
pars$G = 93.4
pars$Dmax = 240
pars$Hmax = 4270
pars$b2 = NA
pars$b3 = NA

## Read in tree-ring data
trw.lyf = read.rwl('./data/LF_QURU.rwl') * 0.1 * 2 #Multiply by 2 (radial growth --> diameter growth) and 0.1 (mm --> cm)
trw.tow = read.rwl('./data/TP_QURU.rwl') * 0.1 * 2 #Multiply by 2 (radial growth --> diameter growth) and 0.1 (mm --> cm)

## Make data frame of estimated and observed DBH
trw.dbh.lyf <- build.dbh.table.lyf(trw.lyf, info.lyf)
trw.dbh.tow <- build.dbh.table.tow(trw.tow, info.tow)
write.csv(rbind(trw.dbh.lyf, trw.dbh.tow), "./data/DBH_Comparison_QURU.csv")

## Calculate actual and naive D and Dopt
trw.comp.dopt.lyf = comp.dopt(trw.lyf, trw.dbh.lyf, pars)
temp.comp.dopt.lyf = na.omit(cbind(c(apply(trw.comp.dopt.lyf[[1]], 2, rbind)), 
                                   c(apply(trw.comp.dopt.lyf[[2]], 2, rbind)), 
                                   c(apply(trw.comp.dopt.lyf[[3]], 2, rbind)), 
                                   c(apply(trw.comp.dopt.lyf[[4]], 2, rbind)), 
                                   c(apply(trw.comp.dopt.lyf[[5]], 2, rbind))))
all.dopt.comp.lyf = data.frame(temp.comp.dopt.lyf)
names(all.dopt.comp.lyf) = c("d.naive", "d.act", "dopt.naive", "dopt.act", "dD")
all.dopt.comp.lyf$dopt.dif = all.dopt.comp.lyf$dopt.naive - all.dopt.comp.lyf$dopt.act

trw.comp.dopt.tow = comp.dopt(trw.tow, trw.dbh.tow, pars)
temp.comp.dopt.tow = na.omit(cbind(c(apply(trw.comp.dopt.tow[[1]], 2, rbind)), 
                                   c(apply(trw.comp.dopt.tow[[2]], 2, rbind)), 
                                   c(apply(trw.comp.dopt.tow[[3]], 2, rbind)), 
                                   c(apply(trw.comp.dopt.tow[[4]], 2, rbind)), 
                                   c(apply(trw.comp.dopt.tow[[5]], 2, rbind))))
all.dopt.comp.tow = data.frame(temp.comp.dopt.tow)
names(all.dopt.comp.tow) = c("d.naive", "d.act", "dopt.naive", "dopt.act", "dD")
all.dopt.comp.tow$dopt.dif = all.dopt.comp.tow$dopt.naive - all.dopt.comp.tow$dopt.act

all.dopt.comp <- rbind(all.dopt.comp.lyf, all.dopt.comp.tow)

## Calculate Stress 
dD = all.dopt.comp$dD
dDopt.naive = all.dopt.comp$dopt.naive
dDopt.act = all.dopt.comp$dopt.act
D.naive = all.dopt.comp$d.naive
D.act = all.dopt.comp$d.act

all.dopt.comp$s.naive.Dbased = all.dopt.comp$dD / all.dopt.comp$dopt.naive
all.dopt.comp$s.act.Dbased = all.dopt.comp$dD / all.dopt.comp$dopt.act
all.dopt.comp$s.dif.Dbased = all.dopt.comp$s.naive.Dbased - all.dopt.comp$s.act.Dbased
rm(list = c("dD", "dDopt.naive", "dDopt.act", "D.naive", "D.act"))

## Write CSVs for use in MATLAB
all.dopt.comp$d.naive.bins = cut(all.dopt.comp$d.naive, seq(from=0, to=150, by=10), include.lowest=TRUE)

dopt.summary = boxplot(dopt.dif~d.naive.bins, data=all.dopt.comp, plot=FALSE)
write.csv(cbind(dopt.summary$names, dopt.summary$n, dopt.summary$stats[3,]), "./data/Dopt_Comparison_QURU.csv")

sD.summary = boxplot(s.dif.Dbased~d.naive.bins, data=all.dopt.comp, plot=FALSE)
write.csv(cbind(sD.summary$names, sD.summary$n, sD.summary$stats[3,]), "./data/Stress_Comparison_QURU.csv")

all.dopt.comp$d.naive.bins = NULL
write.csv(all.dopt.comp, "./data/Sensitivity_QURU.csv")

rm(list = c(ls(pattern='^all'), ls(pattern='^temp'), ls(pattern='^trw'), 'dopt.summary','pars','sD.summary'))

##### ACRU #####

## Optimal growth parameters
pars <- vector(mode="list", length=5)
names(pars) <- c("G", "Dmax", "Hmax", "b2", "b3")
pars$G = 127.9 
pars$Dmax = 180
pars$Hmax = 4420
pars$b2 = NA
pars$b3 = NA

## Read in tree-ring data
trw.lyf = read.rwl('./data/LF_ACRU.rwl') * 0.1 * 2 #Multiply by 2 (radial growth --> diameter growth) and 0.1 (mm --> cm)
trw.tow = read.rwl('./data/TP_ACRU.rwl') * 0.1 * 2 #Multiply by 2 (radial growth --> diameter growth) and 0.1 (mm --> cm)

## Make data frame of estimated and observed DBH
trw.dbh.lyf <- build.dbh.table.lyf(trw.lyf, info.lyf)
trw.dbh.tow <- build.dbh.table.tow(trw.tow, info.tow)
write.csv(rbind(trw.dbh.lyf, trw.dbh.tow), "./data/DBH_Comparison_ACRU.csv")

## Calculate actual and naive D and Dopt
trw.comp.dopt.lyf = comp.dopt(trw.lyf, trw.dbh.lyf, pars)
temp.comp.dopt.lyf = na.omit(cbind(c(apply(trw.comp.dopt.lyf[[1]], 2, rbind)), 
                                   c(apply(trw.comp.dopt.lyf[[2]], 2, rbind)), 
                                   c(apply(trw.comp.dopt.lyf[[3]], 2, rbind)), 
                                   c(apply(trw.comp.dopt.lyf[[4]], 2, rbind)), 
                                   c(apply(trw.comp.dopt.lyf[[5]], 2, rbind))))
all.dopt.comp.lyf = data.frame(temp.comp.dopt.lyf)
names(all.dopt.comp.lyf) = c("d.naive", "d.act", "dopt.naive", "dopt.act", "dD")
all.dopt.comp.lyf$dopt.dif = all.dopt.comp.lyf$dopt.naive - all.dopt.comp.lyf$dopt.act

trw.comp.dopt.tow = comp.dopt(trw.tow, trw.dbh.tow, pars)
temp.comp.dopt.tow = na.omit(cbind(c(apply(trw.comp.dopt.tow[[1]], 2, rbind)), 
                                   c(apply(trw.comp.dopt.tow[[2]], 2, rbind)), 
                                   c(apply(trw.comp.dopt.tow[[3]], 2, rbind)), 
                                   c(apply(trw.comp.dopt.tow[[4]], 2, rbind)), 
                                   c(apply(trw.comp.dopt.tow[[5]], 2, rbind))))
all.dopt.comp.tow = data.frame(temp.comp.dopt.tow)
names(all.dopt.comp.tow) = c("d.naive", "d.act", "dopt.naive", "dopt.act", "dD")
all.dopt.comp.tow$dopt.dif = all.dopt.comp.tow$dopt.naive - all.dopt.comp.tow$dopt.act

all.dopt.comp <- rbind(all.dopt.comp.lyf, all.dopt.comp.tow)

## Calculate Stress 
dD = all.dopt.comp$dD
dDopt.naive = all.dopt.comp$dopt.naive
dDopt.act = all.dopt.comp$dopt.act
D.naive = all.dopt.comp$d.naive
D.act = all.dopt.comp$d.act

all.dopt.comp$s.naive.Dbased = all.dopt.comp$dD / all.dopt.comp$dopt.naive
all.dopt.comp$s.act.Dbased = all.dopt.comp$dD / all.dopt.comp$dopt.act
all.dopt.comp$s.dif.Dbased = all.dopt.comp$s.naive.Dbased - all.dopt.comp$s.act.Dbased
rm(list = c("dD", "dDopt.naive", "dDopt.act", "D.naive", "D.act"))

## Write CSVs for use in MATLAB
all.dopt.comp$d.naive.bins = cut(all.dopt.comp$d.naive, seq(from=0, to=150, by=10), include.lowest=TRUE)

dopt.summary = boxplot(dopt.dif~d.naive.bins, data=all.dopt.comp, plot=FALSE)
write.csv(cbind(dopt.summary$names, dopt.summary$n, dopt.summary$stats[3,]), "./data/Dopt_Comparison_ACRU.csv")

sD.summary = boxplot(s.dif.Dbased~d.naive.bins, data=all.dopt.comp, plot=FALSE)
write.csv(cbind(sD.summary$names, sD.summary$n, sD.summary$stats[3,]), "./data/Stress_Comparison_ACRU.csv")

all.dopt.comp$d.naive.bins = NULL
write.csv(all.dopt.comp, "./data/Sensitivity_ACRU.csv")

rm(list = c(ls(pattern='^all'), ls(pattern='^temp'), ls(pattern='^trw'), 'dopt.summary','pars','sD.summary'))



