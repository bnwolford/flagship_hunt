#library(readstata13)
library(dplyr)
library(bit64)
library(haven)
library(data.table)

### phenotypes 

hunt<-"/home/bwolford/work/phenotypes/all-in-cvd-main_phenotypes/hunt/2021-03-09_HUNTing-for-genes-affecting-cardiovascular-related-traits-105828/2021-03-09_105828_Data.sav"
data<-read_sav(hunt)
names(data)[430]<-"Seb5CRP@NT3BLM"; names(data)[646]<-"Seb5CRP@NT4BLM"; names(data)[707]<-"Seb5CRP@NT2HPM"; names(data)[708]<-"Seb5CRP@NT2BLM2a"

bridge<-"/mnt/work/bridge/all-in-cvd-main-bridge/PID@105828-PID@105118.csv"
bdf<-fread(bridge)

master_file<-"/mnt/work/master/DATASET_20170512/SAMPLE_QC/Masterkey_DATASET.20170512.txt.gz"
mdf<-fread(master_file)

prs_pipe_path<-"/home/bwolford/scratch/brooke/prspipe/prspipe/custom_input/hunt/phenotypes/"
#read in the config file with info about phenotypes
#config<-fread("/home/bwolford/scratch/brooke/prspipe/prspipe/config/studies_for_methods_comparison.tsv")
quant<-c("Height","HbA1c","BMI","Creatinine_eGFR","HDL_cholesterol") #these come from the config file, properly formatted

#geno data looks like this 3999485001_R08C01@ALT_1051180324263 
#cvd main looks like this 1058280000011


df<-data %>% mutate(`PID@105828`=as.integer64(`PID@105828`)) %>% left_join(bdf,by=c("PID@105828"="PID.105828")) %>%
  right_join(mdf,by=c("PID.105118"="gid.current"))  %>% filter(Ancestry4=="EUR")

#binary phenotype file
#pheno<-fread("/mnt/work/workbench/bwolford/intervene/endpointsPhenoFormatHUNT.csv")
#pdf<-left_join(pheno,df,by=c("ID"="IID"))

#BMI
df %>% mutate(latest=coalesce(`BMI@NT3BLM`,`BMI@NT2BLM`,`BMI@NT1BLM`)) %>% mutate(FID=IID) %>% select(IID,FID,latest) %>% 
  filter(!is.na(latest)) %>%
fwrite(paste0(prs_pipe_path,"BMI.tsv"),quote=FALSE,sep="\t",col.names=FALSE,row.names=FALSE)

#height (gonna take HUNT4 when we have it)
#"Hei@NT3BLM" "Hei@NT1BLM" "Hei@NT4BLM" "Hei@NT2BLM" 
df %>% mutate(latest=coalesce(`Hei@NT4BLM`,  `Hei@NT3BLM`, `Hei@NT2BLM`,  `Hei@NT1BLM`)) %>% mutate(FID=IID) %>% select(IID,FID,latest) %>% filter(!is.na(latest)) %>%
  fwrite(paste0(prs_pipe_path,"Height.tsv"),quote=FALSE,sep="\t",col.names=FALSE,row.names=FALSE)

#HDl
#"SeHDLChol@NT2BLM" "SeHDLChol@NT3BLM" "SeHDLChol@NT4BLM"
df %>% mutate(latest=coalesce(`SeHDLChol@NT4BLM`, `SeHDLChol@NT3BLM`, `SeHDLChol@NT4BLM`)) %>% mutate(FID=IID) %>% select(IID,FID,latest) %>% filter(!is.na(latest)) %>%
  fwrite(paste0(prs_pipe_path,"HDL_cholesterol.tsv"),quote=FALSE,sep="\t",col.names=FALSE,row.names=FALSE)

#hba1c (there was another measuremnet standard but I chose IFCC)
#"BloHbA1cIFCC@NT2BLM" "BloHbA1cIFCC@NT3Dia2M1" "BloHbA1cIFCC@NT4BLM" 
df %>% mutate(latest=coalesce(`BloHbA1cIFCC@NT4BLM`)) %>% mutate(FID=IID) %>% select(IID,FID,latest) %>% filter(!is.na(latest)) %>%
  fwrite(paste0(prs_pipe_path,"HbA1c.tsv"),quote=FALSE,sep="\t",col.names=FALSE,row.names=FALSE)


#creatinine
vars<-c("SeCrea@NT2BLM", "SeCrea@NT3BLM","SeCrea@NT4BLM","SeCrea@NT1Dia2MI1")
ages<-c("PartAg@NT2BLQ1","PartAg@NT3BLQ1","PartAg@NT4BLM","PartAg@NT2BLQ1") #is BLQ1 equiv to BLM and what about PartAge for the DIa2M1???w
for (v in 1:length(vars)){
  if (!is.na(vars[v])){
    val<-vars[v]
    age<-ages[v]
    next
  }
}
#Creatinine to eGFR:
#  cr[,eGFR := 30849 * (Creatinine^-1.154) * (age^-0.203) * ifelse(sex == 1, yes = 1., no = 0.742)]
#this formula can create new outliers, so replace any values above 200 with 200
## Note: we don't have exact ages for the various measurements, estimating with closest
sub<-df %>% mutate(latest=coalesce(`SeCrea@NT4BLM`, `SeCrea@NT3BLM`, `SeCrea@NT2BLM`,`SeCrea@NT1Dia2MI1`)) %>% 
  mutate(age=case_when(latest==`SeCrea@NT4BLM` ~ `PartAg@NT4BLM`,
                       latest==`SeCrea@NT3BLM` ~ `PartAg@NT3BLQ1`,
                       latest==`SeCrea@NT2BLM` ~ `PartAg@NT2BLQ1`,
                       latest==`SeCrea@NT1Dia2MI1` ~ `PartAg@NT1BLQ1`)) %>% 
  mutate(FID=IID) %>% select(IID,FID,latest,age,Sex.y) %>% filter(!is.na(latest)) %>% 
  mutate(sex_coeff=ifelse(Sex.y==1,1,0.742)) %>% 
  mutate(eGFR=30849 * (latest^-1.154) * (age^-0.203) * sex_coeff) %>% mutate(eGFR=ifelse(eGFR>200,200, eGFR)) %>%
  select(-latest) %>% rename(latest=eGFR) #remove the previous latest and replace with new creatinine to eGFR


#"SeCreaRecal@NT2BLM"
#A recalibrated variable from Jaffe to enzymatic method exists in the HUNT Databank
  
fwrite(sub,paste0(prs_pipe_path,"Creatinine_eGFR.tsv"),quote=FALSE,sep="\t",col.names=FALSE,row.names=FALSE)


### HbA1c shows really low performance. I wonder if there might be outliers, or if the data quality is just poor. Could you maybe try to look at the distribution of the data at some point? In the UK biobank, the values make a nice bell curve (https://biobank.ndph.ox.ac.uk/ukb/field.cgi?id=30750).
#Same thing for BMI , perhaps there are outliers or data codings like “-1”?
 

