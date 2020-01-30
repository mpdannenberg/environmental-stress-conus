## Apply gap model stress function to ITRDB sites


library(dplR)
source("StressIndexFunctions.R")
gap_pars = read.table("./data/GAP_MODEL_PARS_FINAL.csv", header=T, sep=',')


###### Now calculate stress index ######
setwd("C:\\Users\\Matthew\\Documents\\Data_Analysis\\ITRDB") # Change to local directory!!!

itrdb.info = read.table("ITRDB_NAm_091115.csv", header=T, sep=",", na.strings="NA")

find.old.trees = function(rwl, year){
  end.year = max(year[!is.na(rwl)])
  return(end.year)
}

find.min.year = function(rwl, year){
  min.year = min(year[!is.na(rwl)])
  return(min.year)
}

# List text files
files <- list.files()                             # List all files
txtfiles <- files[grep(glob2rx("*crn.csv"), files)]  # List txt files
n = length(txtfiles)
fConn = file("ITRDB_Site_Locs.txt")

failed = character()
idx = 1 # for when it fails to find species pars
for(i in 1:n){
  setwd("C:\\Users\\Matthew\\Documents\\Data_Analysis\\ITRDB")
  file = txtfiles[i]
  site = unlist(strsplit(file, "[_]"))[1]
  test = grep(site, itrdb.info$SITE)
  species = as.character(itrdb.info$SPECIES)[test]
  spec.ind = grep(species, as.character(gap_pars$CODE))
  fname = paste(c(site, ".rwl"), collapse="")
  if(substr(site, 1, 2)!="ak" & substr(site, 1, 4)!="cana" & substr(site, 1, 4)!="mexi" & !is.na(itrdb.info$LAT[test]) & !is.na(itrdb.info$LON[test])){
  rwl = suppressWarnings(read.rwl(fname)) * 0.1 * 2 # Multiply by 2 (radius --> diameter) and 0.1 (mm --> cm)
  year = as.numeric(row.names(rwl))
  last.year = max(apply(rwl, 2, find.old.trees, year=year))
  if(!length(spec.ind) | last.year<1982){
    failed[idx] = site
    idx = idx+1
  } else{
    # Get parameters from gap model look-up table
    pars <- vector(mode="list", length=5)
    names(pars) <- c("G", "Dmax", "Hmax", "b2", "b3")
    pars$G = gap_pars$G[spec.ind]
    pars$Dmax = gap_pars$Dmax[spec.ind]
    pars$Hmax = gap_pars$Hmax[spec.ind]
    pars$b2 = gap_pars$b2[spec.ind]
    pars$b3 = gap_pars$b3[spec.ind]
    
    # Drop trees that died before 1900
    end.year = apply(rwl, 2, find.old.trees, year=year)
    to.drop = names(rwl)[end.year<1900]
    rwl = rwl[,!names(rwl) %in% to.drop]
    first.year = min(apply(rwl, 2, find.min.year, year=year))
    rwl = rwl[year>=first.year,]
    nrows = dim(rwl)[1]
    ncors = dim(rwl)[2]
    
    S = data.frame(apply(rwl, 2, StressIndex, pars=pars))
    row.names(S) = row.names(rwl)
    
    strToWrite = paste(site,",", itrdb.info$LAT[test],",", itrdb.info$LON[test],"\n")
    cat(strToWrite, file="ITRDB_Site_Locs.txt", append=T)
    
    
    ##### Grab longest sequence of EPS>0.85 #####
    temp = tryCatch({
      ids = autoread.ids(S)
      temp =1
    }, error=function(e){
      temp = 0
      cat("ERROR :",conditionMessage(e), "\n")
      return(temp)
    }
    )
    

    eps.idx = tryCatch({
      if(temp==1){ 
        eps = rwi.stats.running(S, ids=ids)
      } else eps = rwi.stats.running(S)
      eps = na.omit(eps)
      eps.idx = eps$eps>=0.85
    }, error=function(e){
      eps.idx =0
      cat("ERROR :",conditionMessage(e), "\n")
      return(eps.idx)
    }
    )
      
    # test for at least one segment with EPS>0.85
    if(sum(eps.idx>0)){
    starts = eps$start.year[eps.idx]
    ends = eps$end.year[eps.idx]
    
    temp <- cumsum(c(1, diff(starts) - 25))
    temp2 <- rle(temp)
    starts[which(temp == with(temp2, values[which.max(lengths)]))]
    ends[which(temp == with(temp2, values[which.max(lengths)]))]
    
    good.years = seq(from=min(starts), to=max(ends), by=1)
    rm("eps", "eps.idx")
    
    
    # Build chronology and write to csv file
    crn = chron(S)
    
    crn.years = as.numeric(row.names(crn))
    crn.idx = crn.years %in% good.years
    crn = crn[crn.idx,]
    
    write.csv(crn, paste(c(site, "_stress.csv"), collapse=""))
    }
  }
  }
}
close(fConn)


