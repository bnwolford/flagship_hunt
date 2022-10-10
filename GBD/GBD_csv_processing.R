#prepare GBD stats

path<-"/mnt/work/workbench/bwolford/hunt_flagship/GBD/"
mort1<-fread(paste0(path,"IHME-GBD_2019_DATA-538dd722-1.csv"))
mort2<-fread(paste0(path,"IHME-GBD_2019_DATA-538dd722-2.csv"))
mortality<-rbind(mort1,mort2)

