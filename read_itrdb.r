# Read in ITRDB data

adjustLW <- function(lw, ew){
  cores = names(lw)
  lwadj = data.frame(matrix(ncol = dim(lw)[2], nrow = dim(lw)[1]))
  names(lwadj) <- names(lw)
  row.names(lwadj) <- row.names(lw)
  for(j in 1:length(cores)){
    if(names(lw)[j] %in% names(ew)){
      idx = match(names(lw)[j], names(ew))
      mdl1 = lm(lw[,j]~ew[,idx])
      mdl2 = lm(lw[,j]~ew[,idx]+I(ew[,idx]^2))
      if(AIC(mdl2)<AIC(mdl1) & mdl2$coefficients[2]>0 & mdl2$coefficients[3]>0){
        lwadj[!is.na(lw[,j])&!is.na(ew[,idx]), j] = mdl2$residuals+1
      } else lwadj[!is.na(lw[,j])&!is.na(ew[,idx]), j] = mdl1$residuals+1
    }
  }
  return(lwadj)
}



##### ITRDB operations #####

setwd("C:\\Users\\dannenberg\\Documents\\Data_Analysis\\ITRDB_update092717") ## Change to ITRDB directory on local machine!!!
library(dplR)

# List text files
files <- list.files()                             # List all files
txtfiles <- files[grep(glob2rx("*.txt"), files)]  # List txt files

n = length(txtfiles)
trdata <- data.frame(matrix(ncol = 10, nrow = n))
names(trdata) <- c("SITE", "START", "END", "PI", "NAME", "LOCATION", "SPECIES", "LAT", "LON", "ELEV")

for(i in 1:n){
  file = txtfiles[i]
  site = unlist(strsplit(file, "[.]"))[1]
  site = unlist(strsplit(site, "[_]"))[1] # some with "_gap"
  fc = file(file)
  mylist <- strsplit(readLines(fc), ": ")
  e = "great!"
  idx = 1
  idxn = length(mylist)
  repeat{
    testline = mylist[[idx]][1]
    test = grepl("Beginning", testline)
    if(idx==idxn){
      e = "crap!"
      break()
    }
    if(test==TRUE) break()
    idx = idx+1
  }
  trdata$SITE[i] = site
  if(e != "crap!"){
    trdata$START[i] = as.numeric(mylist[[idx]][2])
    trdata$END[i] = as.numeric(mylist[[idx+1]][2])
    trdata$PI[i] = mylist[[idx+2]][2]
    trdata$NAME[i] = mylist[[idx+3]][2]
    trdata$LOCATION[i] = mylist[[idx+4]][2]
    trdata$SPECIES[i] = substring(mylist[[idx+5]][2], 1, 4)
    
    if(grepl("S", mylist[[idx+6]][2])==TRUE){ lat = -1*as.numeric(gsub("\\D", "", mylist[[idx+6]][2])) / 100
    } else lat = as.numeric(gsub("\\D", "", mylist[[idx+6]][2])) / 100
    
    if(!is.na(lat)){
      dsign = sign(lat)
      degree = floor(abs(lat))
      minute = round(100*(abs(lat)-degree))
      
      if(minute < 60) trdata$LAT[i] = dsign*degree + dsign*(minute/60)
      else trdata$LAT[i] = lat
    } else trdata$LAT[i] = NA
    
    if(grepl("W", mylist[[idx+7]][2])==TRUE){ lon = -1*as.numeric(gsub("\\D", "", mylist[[idx+7]][2])) / 100
    } else lon = as.numeric(mylist[[idx+7]][2]) / 100
    
    if(!is.na(lon)){
      dsign = sign(lon)
      degree = floor(abs(lon))
      minute = round(100*(abs(lon)-degree))
      
      if(minute < 60) trdata$LON[i] = dsign*degree + dsign*(minute/60)
      else trdata$LON[i] = lon
    } else trdata$LON[i] = NA
    
    trdata$ELEV[i] = as.numeric(gsub("\\D", "", mylist[[idx+8]][2]))
  }
  
  
  ## Insert code to detrend, chronologize, and export to csv
    # Include exceptions for when it can't find the *.rwl file (write names to txt file)
  rwlFile = paste(c(site, ".rwl"), collapse="")
  if(!length(grep(rwlFile, files))){
    txtToWrite = paste(c(site, "\n"), collapse="")
  } else{
    
    scratchRWL <- tryCatch(
      {
        rwl = read.rwl(rwlFile)
        spl = detrend(rwl, method="Spline")
        
        ids = read.ids(spl)
        window.length = min(50, nrow(spl))
        window.overlap = window.length - 10
        stats <- tryCatch(
          {
            stats = rwi.stats.running(spl, ids=ids, window.length=window.length, window.overlap=window.overlap)
            
          }, error=function(cond){
            stats = rwi.stats.running(spl, window.length=window.length, window.overlap=window.overlap)
            
          }, warning=function(cond){
            stats = rwi.stats.running(spl, window.length=window.length, window.overlap=window.overlap)
            
          }, finally = {
            temp = NULL
          }
        )
        
        write.csv(stats, paste(c(site, "w_stats.csv"), collapse=""))
        
        res = chron(spl, prefix="", prewhiten=TRUE)
        write.csv(res, paste(c(site, "w_crn.csv"), collapse=""))
      }, error=function(e){})
    
    
    ## EW
    site = unlist(strsplit(site, "[w]"))[1]
    rwlFile = paste(c(site, "e.rwl"), collapse="")
    scratchEW <- tryCatch(
      {
        
        ew = read.rwl(rwlFile)
        ew_spl = detrend(ew, method="Spline")
        
        ids = read.ids(ew_spl)
        window.length = min(50, nrow(ew_spl))
        window.overlap = window.length - 10
        stats <- tryCatch(
          {
            stats = rwi.stats.running(ew_spl, ids=ids, window.length=window.length, window.overlap=window.overlap)
            
          }, error=function(cond){
            stats = rwi.stats.running(ew_spl, window.length=window.length, window.overlap=window.overlap)
            
          }, warning=function(cond){
            stats = rwi.stats.running(ew_spl, window.length=window.length, window.overlap=window.overlap)
            
          }, finally = {
            temp = NULL
          }
        )
        
        write.csv(stats, paste(c(site, "e_stats.csv"), collapse=""))
        
        ewres = chron(ew_spl, prefix="", prewhiten=TRUE)
        write.csv(ewres, paste(c(site, "e_crn.csv"), collapse=""))
      }, error=function(e){})
    
    
    ## LW
    rwlFile = paste(c(site, "l.rwl"), collapse="")
    scratchLW <- tryCatch(
      {
        
        lw = read.rwl(rwlFile)
        lw_spl = detrend(lw, method="Spline")
        
        ids = read.ids(lw_spl)
        window.length = min(50, nrow(lw_spl))
        window.overlap = window.length - 10
        stats <- tryCatch(
          {
            stats = rwi.stats.running(lw_spl, ids=ids, window.length=window.length, window.overlap=window.overlap)
            
          }, error=function(cond){
            stats = rwi.stats.running(lw_spl, window.length=window.length, window.overlap=window.overlap)
            
          }, warning=function(cond){
            stats = rwi.stats.running(lw_spl, window.length=window.length, window.overlap=window.overlap)
            
          }, finally = {
            temp = NULL
          }
        )
        
        write.csv(stats, paste(c(site, "l_stats.csv"), collapse=""))
        
        lwres = chron(lw_spl, prefix="", prewhiten=TRUE)
        write.csv(lwres, paste(c(site, "l_crn.csv"), collapse=""))
        
        
        ## Adjusted LW (following Griffin et al. 2011, Tree-Ring Res)
        if(dim(lw_spl)[1]==dim(ew_spl)[1] & dim(lw_spl)[2]==dim(ew_spl)[2]){
          lwadj = adjustLW(lw_spl, ew_spl)
          ids = read.ids(lwadj)
          window.length = min(50, nrow(lwadj))
          window.overlap = window.length - 10
          stats <- tryCatch(
            {
              stats = rwi.stats.running(lwadj, ids=ids, window.length=window.length, window.overlap=window.overlap)
              
            }, error=function(cond){
              stats = rwi.stats.running(lwadj, window.length=window.length, window.overlap=window.overlap)
              
            }, warning=function(cond){
              stats = rwi.stats.running(lwadj, window.length=window.length, window.overlap=window.overlap)
              
            }, finally = {
              temp = NULL
            }
          )
          write.csv(stats, paste(c(site, "la_stats.csv"), collapse=""))
          
          lwadjres = chron(lwadj, prefix="", prewhiten=TRUE)
          write.csv(lwadjres, paste(c(site, "la_crn.csv"), collapse=""))
          
        } else {
          commonYrs = intersect(row.names(lwres), row.names(ew))
          lwadjres = data.frame(matrix(nrow = length(commonYrs), ncol=3))
          names(lwadjres) <- names(lwres)
          row.names(lwadjres) <- commonYrs
          lwadjres[,3] = lwres[row.names(lwres) %in% commonYrs,3]
          
          # STD adjustment
          lwtemp = lwres[row.names(lwres) %in% commonYrs,1]
          ewtemp = ewres[row.names(ewres) %in% commonYrs,1]
          mdl = lm(lwtemp~ewtemp)
          lwadjres[,1] = mdl$residuals +1 
          
          # RES adjustment
          lwtemp = lwres[row.names(lwres) %in% commonYrs,2]
          ewtemp = ewres[row.names(ewres) %in% commonYrs,2]
          naidx = !is.na(lwtemp)&!is.na(ewtemp)
          mdl = lm(lwtemp~ewtemp, na.action = na.exclude)
          lwadjres[naidx,2] = mdl$residuals +1 
          
          write.csv(lwadjres, paste(c(site, "la_crn.csv"), collapse=""))
          
        }
        
      }, error=function(e){})
    
    
  }

  
  
  close(fc)
}



