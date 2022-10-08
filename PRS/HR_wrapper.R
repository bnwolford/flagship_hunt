### wrapper script to run hazard ratio scripts

script_path="/mnt/work/workbench/bwolford/hunt_flagship/PRS/"
pheno_file="/mnt/work/workbench/bwolford/intervene/2022_10_06/endpointsPhenoFormatHUNT.csv"
prs_path="/home/bwolford/scratch/brooke/2022_10_07_scores/"
ID1<-"ID"
ID2<-"IID"
covariates<-"PC1 + PC2 + PC3 + PC4 + PC5"
biobank_name<-"HUNT"
output_dirt<-"/mnt/work/workbench/bwolford/intervene/2022_10_06/"
#NOTE, data must already be subsetted by ancestry
#paths must include last backslash
  
source(paste0(script_path,"HUNT_HazardRatio_AgeandSexStratified.R"))
source(paste0(script_path,"HUNT_HazardRatio_AgeStratified.R"))
source(paste0(script_path,"HUNT_HazardRatio_FullSample.R"))
source(paste0(script_path,"HUNT_HazardRatioperSD_AgeandSexStratified.R"))
source(paste0(script_path,"HUNT_HazardRatioperSD_AgeStratified.R"))
source(paste0(script_path, "HUNT_HazardRatioperStandardDeviation.R"))
source(paste0(script_path, "HUNT_HazardRatio_SexInteraction.R"))
source(paste0(script_path, "HUNT_HazardRatio_SexStratified.R"))

