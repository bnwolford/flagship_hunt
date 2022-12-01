### wrapper script to run hazard ratio scripts

script_path="/mnt/work/workbench/bwolford/hunt_flagship/PRS/"
pheno_file="/mnt/work/workbench/bwolford/intervene/2022_10_06/endpointsPhenoFormatHUNT.csv"
prs_path="/home/bwolford/scratch/brooke/2022_10_07_scores/"
ID1<-"ID"
ID2<-"IID"
covariates<-"PC1 + PC2 + PC3 + PC4 + PC5"
#covariates<- "PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10"
biobank_name<-"HUNT"
output_dir<-"/mnt/work/workbench/bwolford/intervene/2022_10_06/"
custom_covar<-c("BATCH")
covariates<-paste(sep=" + ",covariates,custom_covar)
#NOTE, data must already be subsetted by ancestry
#paths must include last backslash
  
#instructions in readme.md https://github.com/intervene-EU-H2020/flagship

#source(paste0(script_path, "HUNT_HazardRatioperStandardDeviation.R"))
print("HazardRatioperStandardDeviation.R finished running")

#source(paste0(script_path,"HUNT_HazardRatio_AgeStratified.R"))
print("HazardRatio_AgeStratified finished running")

#source(paste0(script_path,"HUNT_HazardRatio_AgeandSexStratified.R"))
print("HazardRatio_AgeandSexStratified.R finished running")

#source(paste0(script_path,"HUNT_HazardRatio_FullSample.R"))
print("HUNT_HazardRatio_FullSample.R finished running")

source(paste0(script_path,"HUNT_HazardRatioperSD_AgeandSexStratified.R"))
print("HUNT_HazardRatioperSD_AgeandSexStratified.R finished running")

source(paste0(script_path,"HUNT_HazardRatioperSD_AgeStratified.R"))
print("HUNT_HazardRatioperSD_AgeStratified.R finished running")

source(paste0(script_path, "HUNT_HazardRatio_SexInteraction.R"))
print("HUNT_HazardRatio_SexInteraction.R finished running")

source(paste0(script_path, "HUNT_HazardRatio_SexStratified.R"))
print("HUNT_HazardRatio_SexStratified.R finished running")

