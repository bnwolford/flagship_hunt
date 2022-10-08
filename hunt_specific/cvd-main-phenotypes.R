library(readstata13)

### phenotypes 

hunt<-"/home/bwolford/work/phenotypes/all-in-cvd-main_phenotypes/hunt/2021-03-09_HUNTing-for-genes-affecting-cardiovascular-related-traits-105828/2021-03-09_105828_Data.sav"
data<-read_sav(hunt)

bridge<-"/mnt/work/bridge/all-in-cvd-main-bridge/PID@105828-PID@105118.csv"
bdf<-fread(bridge)

master_file<-"/mnt/work/master/DATASET_20170512/SAMPLE_QC/Masterkey_DATASET.20170512.txt.gz"
mdf<-fread(master)

#BMI
"BMI@NT1BLM" "BMI@NT2BLM" "BMI@NT3BLM"

#height
"Hei@NT3BLM" "Hei@NT1BLM" "Hei@NT4BLM" "Hei@NT2BLM" 

#HDl
"SeHDLChol@NT2BLM" "SeHDLChol@NT3BLM" "SeHDLChol@NT4BLM"

#hba1c 
"BloHbA1cIFCC@NT2BLM"
"BloHbA1cIFCC@NT3Dia2M1"
"BloHbA1cIFCC@NT4BLM" 

#creatinine
"SeCrea@NT2BLM"      "SeCreaRecal@NT2BLM" "SeCrea@NT3BLM" "SeCrea@NT4BLM"      "SeCrea@NT1Dia2MI1" 


###editing I9_CHD so we include I20.0

fam_file<-"/mnt/scratch/brooke/bcf/all.log.fam"
output_dir="/mnt/work/workbench/bwolford/intervene/"

#files<-list.files("/mnt/work/phenotypes/allin-phecode-2018_41492/kilde/hnt/",full.names=TRUE)
#files<- files[!grepl("etter",files)] #what is issue with that sav file?

