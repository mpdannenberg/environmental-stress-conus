## Estimate sensitivity of P. ponderosa and P. menziesii to DBH uncertainty

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
    name = substr(name, 1, nchar(name)-1)
    ind = pmatch(name, info$Name)
    obs = info$DBH[ind]*100 # m --> cm
    dbh.table$Est[i] = est
    dbh.table$Obs[i] = obs
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


#####  PIPO  #####

## Optimal growth parameters for PIPO
pars <- vector(mode="list", length=5)
names(pars) <- c("G", "Dmax", "Hmax", "b2", "b3")
pars$G = 66.2
pars$Dmax = 270
pars$Hmax = 7800
pars$b2 = NA
pars$b3 = NA

bcr <- read.rwl("./data/BCR_QCd_040715.txt") * 0.1 * 2 #Multiply by 2 (radial growth --> diameter growth) and 0.1 (mm --> cm)
fod <- read.rwl("./data/FOD_QCd_032315.txt") * 0.1 * 2 #Multiply by 2 (radial growth --> diameter growth) and 0.1 (mm --> cm)
fod$FOD041C1 = NULL
fod$FOD041C2 = NULL
fod$FOD033H1 = NULL
fod$FOD033H2 = NULL
rrr <- read.rwl("./data/RRR_QCd_041315.txt") * 0.1 * 2 #Multiply by 2 (radial growth --> diameter growth) and 0.1 (mm --> cm)
rrr$RRR019E1 = NULL
rrr$RRR019E2 = NULL
snd <- read.rwl("./data/SND_QCd_040915.txt") * 0.1 * 2 #Multiply by 2 (radial growth --> diameter growth) and 0.1 (mm --> cm)
snd$SND033E1 = NULL
snd$SND033E2 = NULL
sug <- read.rwl("./data/SUG_QCd_021615.txt") * 0.1 * 2 #Multiply by 2 (radial growth --> diameter growth) and 0.1 (mm --> cm)
tmr <- read.rwl("./data/TMR_QCd_043015.txt") * 0.1 * 2 #Multiply by 2 (radial growth --> diameter growth) and 0.1 (mm --> cm)
tmr$TMRX01 = NULL
tmr$TMR024D1 = NULL
tmr$TMR024D2 = NULL
tmr$TMR035D1 = NULL
tmr$TMR035D2 = NULL


## Get measured DBH
bcr.info <- read.table("BCR_info.csv", sep=",", header=T, na.strings=-9999)
fod.info <- read.table("FOD_info.csv", sep=",", header=T, na.strings=-9999)
rrr.info <- read.table("RRR_info.csv", sep=",", header=T, na.strings=-9999)
snd.info <- read.table("SND_info.csv", sep=",", header=T, na.strings=-9999)
sug.info <- read.table("SUG_info.csv", sep=",", header=T, na.strings=-9999)
tmr.info <- read.table("TMR_info.csv", sep=",", header=T, na.strings=-9999)


## Make data frame of estimated and observed DBH
bcr.dbh = build.dbh.table(bcr, bcr.info)
fod.dbh = build.dbh.table(fod, fod.info)
rrr.dbh = build.dbh.table(rrr, rrr.info)
snd.dbh = build.dbh.table(snd, snd.info)
sug.dbh = build.dbh.table(sug, sug.info)
tmr.dbh = build.dbh.table(tmr, tmr.info)
all.dbh = rbind(bcr.dbh, fod.dbh, rrr.dbh, snd.dbh, sug.dbh, tmr.dbh)
write.csv(all.dbh, './data/DBH_Comparison_PIPO.csv')


## Calculate actual and naive D and Dopt
bcr.all = comp.dopt(bcr, bcr.dbh, pars)
fod.all = comp.dopt(fod, fod.dbh, pars)
rrr.all = comp.dopt(rrr, rrr.dbh, pars)
snd.all = comp.dopt(snd, snd.dbh, pars)
sug.all = comp.dopt(sug, sug.dbh, pars)
tmr.all = comp.dopt(tmr, tmr.dbh, pars)


## Reshape into long form
temp.bcr = na.omit(cbind(c(apply(bcr.all[[1]], 2, rbind)), c(apply(bcr.all[[2]], 2, rbind)), c(apply(bcr.all[[3]], 2, rbind)), c(apply(bcr.all[[4]], 2, rbind)), c(apply(bcr.all[[5]], 2, rbind))))
temp.fod = na.omit(cbind(c(apply(fod.all[[1]], 2, rbind)), c(apply(fod.all[[2]], 2, rbind)), c(apply(fod.all[[3]], 2, rbind)), c(apply(fod.all[[4]], 2, rbind)), c(apply(fod.all[[5]], 2, rbind))))
temp.rrr = na.omit(cbind(c(apply(rrr.all[[1]], 2, rbind)), c(apply(rrr.all[[2]], 2, rbind)), c(apply(rrr.all[[3]], 2, rbind)), c(apply(rrr.all[[4]], 2, rbind)), c(apply(rrr.all[[5]], 2, rbind))))
temp.snd = na.omit(cbind(c(apply(snd.all[[1]], 2, rbind)), c(apply(snd.all[[2]], 2, rbind)), c(apply(snd.all[[3]], 2, rbind)), c(apply(snd.all[[4]], 2, rbind)), c(apply(snd.all[[5]], 2, rbind))))
temp.sug = na.omit(cbind(c(apply(sug.all[[1]], 2, rbind)), c(apply(sug.all[[2]], 2, rbind)), c(apply(sug.all[[3]], 2, rbind)), c(apply(sug.all[[4]], 2, rbind)), c(apply(sug.all[[5]], 2, rbind))))
temp.tmr = na.omit(cbind(c(apply(tmr.all[[1]], 2, rbind)), c(apply(tmr.all[[2]], 2, rbind)), c(apply(tmr.all[[3]], 2, rbind)), c(apply(tmr.all[[4]], 2, rbind)), c(apply(tmr.all[[5]], 2, rbind))))

all.dopt.comp = data.frame(rbind(temp.bcr, temp.fod, temp.rrr, temp.snd, temp.sug, temp.tmr))
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
write.csv(cbind(dopt.summary$names, dopt.summary$n, dopt.summary$stats[3,]), './data/Dopt_Comparison_PIPO.csv')

sD.summary = boxplot(s.dif.Dbased~d.naive.bins, data=all.dopt.comp, plot=FALSE)
write.csv(cbind(sD.summary$names, sD.summary$n, sD.summary$stats[3,]), './data/Stress_Comparison_PIPO.csv')

all.dopt.comp$d.naive.bins = NULL
write.csv(all.dopt.comp, './data/Sensitivity_PIPO.csv')

rm(list = setdiff(ls(), lsf.str()))


#####  PSME  #####

## Optimal growth parameters for PSME
pars <- vector(mode="list", length=5)
names(pars) <- c("G", "Dmax", "Hmax", "b2", "b3")
pars$G = 56.7
pars$Dmax = 425
pars$Hmax = 8500
pars$b2 = NA
pars$b3 = NA

cwc <- read.rwl("./data/CWC_PSME.txt") * 0.1 * 2 #Multiply by 2 (radial growth --> diameter growth) and 0.1 (mm --> cm)
cwc[, nchar(names(cwc))>6] = NULL # Get rid of multi-part cores
cwc$CWC23A = NULL # no DBH recorded
efv <- read.rwl("./data/EFV_PSME.txt") * 0.1 * 2 #Multiply by 2 (radial growth --> diameter growth) and 0.1 (mm --> cm)
efv[, nchar(names(efv))>6] = NULL # Get rid of multi-part cores
efv$EFV12A = NULL # no DBH recorded
etc <- read.rwl("./data/ETC_PSME.txt") * 0.1 * 2 #Multiply by 2 (radial growth --> diameter growth) and 0.1 (mm --> cm)
etc[, nchar(names(etc))>6] = NULL # Get rid of multi-part cores
etc$ETC01A = NULL
etc$ETC01B = NULL
etc$ETC03A = NULL
etc$ETC03B = NULL
etc$ETC04A = NULL
etc$ETC04B = NULL
etc$ETC06A = NULL
etc$ETC06B = NULL
etc$ETC07A = NULL
etc$ETC07B = NULL
etc$ETC10A = NULL


## Get measured DBH
cwc.info <- read.table("CWC_info.csv", sep=",", header=T, na.strings=-9999)
efv.info <- read.table("EFV_info.csv", sep=",", header=T, na.strings=-9999)
etc.info <- read.table("ETC_info.csv", sep=",", header=T, na.strings=-9999)


## Make data frame of estimated and observed DBH
cwc.dbh = build.dbh.table(cwc, cwc.info)
efv.dbh = build.dbh.table(efv, efv.info)
etc.dbh = build.dbh.table(etc, etc.info)
all.dbh = rbind(cwc.dbh, efv.dbh, etc.dbh)
write.csv(all.dbh, './data/DBH_Comparison_PSME.csv')


## Calculate actual and naive D and Dopt
cwc.all = comp.dopt(cwc, cwc.dbh, pars)
efv.all = comp.dopt(efv, efv.dbh, pars)
etc.all = comp.dopt(etc, etc.dbh, pars)


## Reshape into long form
temp.cwc = na.omit(cbind(c(apply(cwc.all[[1]], 2, rbind)), c(apply(cwc.all[[2]], 2, rbind)), c(apply(cwc.all[[3]], 2, rbind)), c(apply(cwc.all[[4]], 2, rbind)), c(apply(cwc.all[[5]], 2, rbind))))
temp.efv = na.omit(cbind(c(apply(efv.all[[1]], 2, rbind)), c(apply(efv.all[[2]], 2, rbind)), c(apply(efv.all[[3]], 2, rbind)), c(apply(efv.all[[4]], 2, rbind)), c(apply(efv.all[[5]], 2, rbind))))
temp.etc = na.omit(cbind(c(apply(etc.all[[1]], 2, rbind)), c(apply(etc.all[[2]], 2, rbind)), c(apply(etc.all[[3]], 2, rbind)), c(apply(etc.all[[4]], 2, rbind)), c(apply(etc.all[[5]], 2, rbind))))


## Pool sites
all.dopt.comp = data.frame(rbind(temp.cwc, temp.efv, temp.etc))
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
write.csv(cbind(dopt.summary$names, dopt.summary$n, dopt.summary$stats[3,]), './data/Dopt_Comparison_PSME.csv')

sD.summary = boxplot(s.dif.Dbased~d.naive.bins, data=all.dopt.comp, plot=FALSE)
write.csv(cbind(sD.summary$names, sD.summary$n, sD.summary$stats[3,]), './data/Stress_Comparison_PSME.csv')

all.dopt.comp$d.naive.bins = NULL
write.csv(all.dopt.comp, './data/Sensitivity_PSME.csv')

rm(list = setdiff(ls(), lsf.str()))


