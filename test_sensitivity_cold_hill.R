## Things to do for chapter 1
#     Optimal growth using known diameter from PIPO sites
#     Curve matching approach (Erika's idea) for estimating D bias
#     Compare to flux tower GPP



## Load necessary functions and parameters
library(dplR)
library(ggplot2)
library(extrafont)

build.dbh.table = function(rwl, info){
  n = dim(rwl)[2]
  dbh.table = data.frame(matrix(ncol=3, nrow=n))
  names(dbh.table) <- c("Core", "Est", "Obs")
  dbh.table$Core <- names(rwl)
  for(i in 1:n){
    x = rwl[,i]
    est = sum(x[!is.na(x)]) # cm
    name = dbh.table$Core[i]
    stand = as.numeric(substr(name, 1, 2))
    plot = as.numeric(substr(name, 3, 4))
    tree = as.numeric(substr(name, 5, 6))
    ind = which(info$stand==stand & info$plot==plot & info$tree==tree)
    try({
    obs = info$dbh[ind]
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

source("C:\\Users\\dannenberg\\Documents\\Publications\\Dannenberg_et_al_ForestEcolMgmt\\Sensitivity\\StressIndexFunctions.R")

setwd("C:\\Users\\dannenberg\\Documents\\Data_Analysis\\Pederson_TreeRingData\\ColdHillPlateauKYTreeRingDataJan2015")
info <- read.table('metadata.csv', header=T, sep=',')
info$dbh = info$dbh * 2.54

workDir = 'C:\\Users\\dannenberg\\Documents\\Publications\\Dannenberg_et_al_ForestEcolMgmt\\Sensitivity\\'



##### QUVE #####

## Optimal growth parameters
pars <- vector(mode="list", length=5)
names(pars) <- c("G", "Dmax", "Hmax", "b2", "b3")
pars$G = 199
pars$Dmax = 240
pars$Hmax = 4570
pars$b2 = NA
pars$b3 = NA

## Read in tree-ring data
trw = read.rwl('ColdHillBlackOakQUVE.rw') * 0.1 * 2 #Multiply by 2 (radial growth --> diameter growth) and 0.1 (mm --> cm)

## Make data frame of estimated and observed DBH
trw.dbh = build.dbh.table(trw, info)
write.csv(trw.dbh, paste(c(workDir, "DBH_Comparison_QUVE.csv"), collapse=""))

## Calculate actual and naive D and Dopt
trw.comp.dopt = comp.dopt(trw, trw.dbh, pars)
temp.comp.dopt = na.omit(cbind(c(apply(trw.comp.dopt[[1]], 2, rbind)), 
                         c(apply(trw.comp.dopt[[2]], 2, rbind)), 
                         c(apply(trw.comp.dopt[[3]], 2, rbind)), 
                         c(apply(trw.comp.dopt[[4]], 2, rbind)), 
                         c(apply(trw.comp.dopt[[5]], 2, rbind))))
all.dopt.comp = data.frame(temp.comp.dopt)
names(all.dopt.comp) = c("d.naive", "d.act", "dopt.naive", "dopt.act", "dD")
all.dopt.comp$dopt.dif = all.dopt.comp$dopt.naive - all.dopt.comp$dopt.act

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
write.csv(cbind(dopt.summary$names, dopt.summary$n, dopt.summary$stats[3,]), paste(c(workDir, "Dopt_Comparison_QUVE.csv"), collapse=""))

sD.summary = boxplot(s.dif.Dbased~d.naive.bins, data=all.dopt.comp, plot=FALSE)
write.csv(cbind(sD.summary$names, sD.summary$n, sD.summary$stats[3,]), paste(c(workDir, "Stress_Comparison_QUVE.csv"), collapse=""))

all.dopt.comp$d.naive.bins = NULL
write.csv(all.dopt.comp, paste(c(workDir, "Sensitivity_QUVE.csv"), collapse=""))

rm(list = c('all.dopt.comp', 'temp.comp.dopt', 'trw','trw.dbh','dopt.summary','pars','sD.summary','trw.comp.dopt'))



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
trw = read.rwl('ColdHillEasternHemlockTSCA.rw') * 0.1 * 2 #Multiply by 2 (radial growth --> diameter growth) and 0.1 (mm --> cm)

## Make data frame of estimated and observed DBH
trw.dbh = build.dbh.table(trw, info)
write.csv(trw.dbh, paste(c(workDir, "DBH_Comparison_TSCA.csv"), collapse=""))

## Calculate actual and naive D and Dopt
trw.comp.dopt = comp.dopt(trw, trw.dbh, pars)
temp.comp.dopt = na.omit(cbind(c(apply(trw.comp.dopt[[1]], 2, rbind)), 
                               c(apply(trw.comp.dopt[[2]], 2, rbind)), 
                               c(apply(trw.comp.dopt[[3]], 2, rbind)), 
                               c(apply(trw.comp.dopt[[4]], 2, rbind)), 
                               c(apply(trw.comp.dopt[[5]], 2, rbind))))
all.dopt.comp = data.frame(temp.comp.dopt)
names(all.dopt.comp) = c("d.naive", "d.act", "dopt.naive", "dopt.act", "dD")
all.dopt.comp$dopt.dif = all.dopt.comp$dopt.naive - all.dopt.comp$dopt.act

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
write.csv(cbind(dopt.summary$names, dopt.summary$n, dopt.summary$stats[3,]), paste(c(workDir, "Dopt_Comparison_TSCA.csv"), collapse=""))

sD.summary = boxplot(s.dif.Dbased~d.naive.bins, data=all.dopt.comp, plot=FALSE)
write.csv(cbind(sD.summary$names, sD.summary$n, sD.summary$stats[3,]), paste(c(workDir, "Stress_Comparison_TSCA.csv"), collapse=""))

all.dopt.comp$d.naive.bins = NULL
write.csv(all.dopt.comp, paste(c(workDir, "Sensitivity_TSCA.csv"), collapse=""))

rm(list = c('all.dopt.comp', 'temp.comp.dopt', 'trw','trw.dbh','dopt.summary','pars','sD.summary','trw.comp.dopt'))



##### QUCO #####

## NOTE: somethings up with this one.... summing ring widths is much higher than measured DBH?

## Optimal growth parameters
pars <- vector(mode="list", length=5)
names(pars) <- c("G", "Dmax", "Hmax", "b2", "b3")
pars$G = 69
pars$Dmax = 150
pars$Hmax = 3110
pars$b2 = NA
pars$b3 = NA

## Read in tree-ring data
trw = read.rwl('ColdHillScarletOakQUCO.rw') * 0.1 #Multiply by 2 (radial growth --> diameter growth) and 0.1 (mm --> cm)

## Make data frame of estimated and observed DBH
trw.dbh = build.dbh.table(trw, info)
write.csv(trw.dbh, paste(c(workDir, "DBH_Comparison_QUCO.csv"), collapse=""))

## Calculate actual and naive D and Dopt
trw.comp.dopt = comp.dopt(trw, trw.dbh, pars)
temp.comp.dopt = na.omit(cbind(c(apply(trw.comp.dopt[[1]], 2, rbind)), 
                               c(apply(trw.comp.dopt[[2]], 2, rbind)), 
                               c(apply(trw.comp.dopt[[3]], 2, rbind)), 
                               c(apply(trw.comp.dopt[[4]], 2, rbind)), 
                               c(apply(trw.comp.dopt[[5]], 2, rbind))))
all.dopt.comp = data.frame(temp.comp.dopt)
names(all.dopt.comp) = c("d.naive", "d.act", "dopt.naive", "dopt.act", "dD")
all.dopt.comp$dopt.dif = all.dopt.comp$dopt.naive - all.dopt.comp$dopt.act

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
write.csv(cbind(dopt.summary$names, dopt.summary$n, dopt.summary$stats[3,]), paste(c(workDir, "Dopt_Comparison_QUCO.csv"), collapse=""))

sD.summary = boxplot(s.dif.Dbased~d.naive.bins, data=all.dopt.comp, plot=FALSE)
write.csv(cbind(sD.summary$names, sD.summary$n, sD.summary$stats[3,]), paste(c(workDir, "Stress_Comparison_QUCO.csv"), collapse=""))

all.dopt.comp$d.naive.bins = NULL
write.csv(all.dopt.comp, paste(c(workDir, "Sensitivity_QUCO.csv"), collapse=""))

rm(list = c('all.dopt.comp', 'temp.comp.dopt', 'trw','trw.dbh','dopt.summary','pars','sD.summary','trw.comp.dopt'))



##### LITU #####

## Optimal growth parameters
pars <- vector(mode="list", length=5)
names(pars) <- c("G", "Dmax", "Hmax", "b2", "b3")
pars$G = 82.9
pars$Dmax = 210
pars$Hmax = 5180
pars$b2 = NA
pars$b3 = NA

## Read in tree-ring data
trw = read.rwl('ColdHillTulipPoplarLITU.rw') * 0.1 * 2 #Multiply by 2 (radial growth --> diameter growth) and 0.1 (mm --> cm)

## Make data frame of estimated and observed DBH
trw.dbh = build.dbh.table(trw, info)
write.csv(trw.dbh, paste(c(workDir, "DBH_Comparison_LITU.csv"), collapse=""))

## Calculate actual and naive D and Dopt
trw.comp.dopt = comp.dopt(trw, trw.dbh, pars)
temp.comp.dopt = na.omit(cbind(c(apply(trw.comp.dopt[[1]], 2, rbind)), 
                               c(apply(trw.comp.dopt[[2]], 2, rbind)), 
                               c(apply(trw.comp.dopt[[3]], 2, rbind)), 
                               c(apply(trw.comp.dopt[[4]], 2, rbind)), 
                               c(apply(trw.comp.dopt[[5]], 2, rbind))))
all.dopt.comp = data.frame(temp.comp.dopt)
names(all.dopt.comp) = c("d.naive", "d.act", "dopt.naive", "dopt.act", "dD")
all.dopt.comp$dopt.dif = all.dopt.comp$dopt.naive - all.dopt.comp$dopt.act

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
write.csv(cbind(dopt.summary$names, dopt.summary$n, dopt.summary$stats[3,]), paste(c(workDir, "Dopt_Comparison_LITU.csv"), collapse=""))

sD.summary = boxplot(s.dif.Dbased~d.naive.bins, data=all.dopt.comp, plot=FALSE)
write.csv(cbind(sD.summary$names, sD.summary$n, sD.summary$stats[3,]), paste(c(workDir, "Stress_Comparison_LITU.csv"), collapse=""))

all.dopt.comp$d.naive.bins = NULL
write.csv(all.dopt.comp, paste(c(workDir, "Sensitivity_LITU.csv"), collapse=""))

rm(list = c('all.dopt.comp', 'temp.comp.dopt', 'trw','trw.dbh','dopt.summary','pars','sD.summary','trw.comp.dopt'))



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
trw = read.rwl('ColdHillRedMapleACRU.rw') * 0.1 #Multiply by 2 (radial growth --> diameter growth) and 0.1 (mm --> cm)

## Make data frame of estimated and observed DBH
trw.dbh = build.dbh.table(trw, info)
write.csv(trw.dbh, paste(c(workDir, "DBH_Comparison_ACRU.csv"), collapse=""))

## Calculate actual and naive D and Dopt
trw.comp.dopt = comp.dopt(trw, trw.dbh, pars)
temp.comp.dopt = na.omit(cbind(c(apply(trw.comp.dopt[[1]], 2, rbind)), 
                               c(apply(trw.comp.dopt[[2]], 2, rbind)), 
                               c(apply(trw.comp.dopt[[3]], 2, rbind)), 
                               c(apply(trw.comp.dopt[[4]], 2, rbind)), 
                               c(apply(trw.comp.dopt[[5]], 2, rbind))))
all.dopt.comp = data.frame(temp.comp.dopt)
names(all.dopt.comp) = c("d.naive", "d.act", "dopt.naive", "dopt.act", "dD")
all.dopt.comp$dopt.dif = all.dopt.comp$dopt.naive - all.dopt.comp$dopt.act

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
write.csv(cbind(dopt.summary$names, dopt.summary$n, dopt.summary$stats[3,]), paste(c(workDir, "Dopt_Comparison_ACRU.csv"), collapse=""))

sD.summary = boxplot(s.dif.Dbased~d.naive.bins, data=all.dopt.comp, plot=FALSE)
write.csv(cbind(sD.summary$names, sD.summary$n, sD.summary$stats[3,]), paste(c(workDir, "Stress_Comparison_ACRU.csv"), collapse=""))

all.dopt.comp$d.naive.bins = NULL
write.csv(all.dopt.comp, paste(c(workDir, "Sensitivity_ACRU.csv"), collapse=""))

rm(list = c('all.dopt.comp', 'temp.comp.dopt', 'trw','trw.dbh','dopt.summary','pars','sD.summary','trw.comp.dopt'))












