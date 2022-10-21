#!/usr/bin/env Rscript 
library(system)
####wrapper for countries and studies

#bb_gbd_phenos <- c("Total cancers", "Appendicitis", "Asthma", "Atrial fibrillation and flutter", "Breast cancer", "Ischemic heart disease", "Colon and rectum cancer", "Idiopathic epilepsy", "Gout", "Osteoarthritis hip", "Osteoarthritis knee", "Major depressive disorder", "Malignant skin melanoma", "Prostate cancer", "Rheumatoid arthritis", "Subarachnoid hemorrhage", "Diabetes mellitus type 1", "Diabetes mellitus type 2")
#bb_hr_phenos <- c("C3_CANCER", "K11_APPENDACUT", "J10_ASTHMA", "I9_AF", "C3_BREAST", "I9_CHD", "C3_COLORECTAL", "G6_EPLEPSY", "GOUT", "COX_ARTHROSIS", "KNEE_ARTHROSIS", "F5_DEPRESSIO", "C3_MELANOMA_SKIN", "C3_PROSTATE", "RHEUMA_SEROPOS_OTH", "I9_SAH", "T1D", "T2D")

nk<-50

#HUNT
country<-"Norway"
path<-"/mnt/work/workbench/bwolford/flagship/AbsoluteRiskEstimation/" #must have final backslash
full_HR_path<-"/mnt/work/workbench/bwolford/intervene/2022_10_06/HR_FullSample_HUNT.csv" #must go straight to csv
output_dir<-"/mnt/work/workbench/bwolford/intervene/2022_10_06/RiskEstimates/" #must have final backslash
biobank<-"HUNT"
#source(paste0(path,"LifetimeRiskEstimation_FullSample.R"))

system("Rscript /mnt/work/workbench/bwolford/flagship/AbsoluteRiskEstimation/LifetimeRiskEstimation_FullSample_Script.R --country Norway --path /mnt/work/workbench/bwolford/hunt_flagship/AbsoluteRiskEstimation/ --output_dir /mnt/work/workbench/bwolford/intervene/2022_10_06/RiskEstimates/ --biobank HUNT --k 100 --hr /mnt/work/workbench/bwolford/intervene/2022_10_06/HR_FullSample_HUNT.csv")

#BBJ

#ESTBB
country<-"Estonia"
path<-"/mnt/work/workbench/bwolford/flagship/AbsoluteRiskEstimation/" #must have final backslash
full_HR_path<-"/mnt/work/workbench/bwolford/intervene/GoogleDrive/EstBB_HazardRatios/HR_FullSample_EstBB.csv"
df<-fread(full_HR_path)
fwrite(df[,-1],"/mnt/work/workbench/bwolford/intervene/GoogleDrive/EstBB_HazardRatios/HR_FullSample_EstBB_reformat.csv")
full_HR_path<-"/mnt/work/workbench/bwolford/intervene/GoogleDrive/EstBB_HazardRatios/HR_FullSample_EstBB_reformat.csv"
output_dir<-"/mnt/work/workbench/bwolford/intervene/2022_10_06/RiskEstimates/" #must have final backslash
biobank<-"EstBB"
ancestry<-"EUR"
#source(paste0(path,"LifetimeRiskEstimation_FullSample.R"))

country<-"Japan"
path<-"/mnt/work/workbench/bwolford/flagship/AbsoluteRiskEstimation/" #must have final backslash
full_HR_path<-'/mnt/work/workbench/bwolford/intervene/GoogleDrive/Biobank_Japan_HazardRatios/HR_FullSampleBBJ.csv'
df<-fread(full_HR_path)
fwrite(df[,-1],"/mnt/work/workbench/bwolford/intervene/GoogleDrive/Biobank_Japan_HazardRatios/HR_FullSampleBBJ_reformat.csv")
full_HR_path<-"/mnt/work/workbench/bwolford/intervene/GoogleDrive/Biobank_Japan_HazardRatios/HR_FullSampleBBJ_reformat.csv"
output_dir<-"/mnt/work/workbench/bwolford/intervene/2022_10_06/RiskEstimates/" #must have final backslash
biobank<-"BBJ"
ancestry<-"EAS"

system(paste0("Rscript /mnt/work/workbench/bwolford/flagship/AbsoluteRiskEstimation/LifetimeRiskEstimation_FullSample_Script.R --country ", country," --path ",path, " --output_dir ", output_dir, " --biobank ", biobank, " --k ", nk," --hr ", full_HR_path))

country<-"Scotland"
path<-"/mnt/work/workbench/bwolford/flagship/AbsoluteRiskEstimation/"
biobank<-"GS"
ancestry<-"EUR"
output_dir<-"/mnt/work/workbench/bwolford/intervene/2022_10_06/RiskEstimates/"
full_HR_path<-"/mnt/work/workbench/bwolford/intervene/GoogleDrive/GenerationScotland_HazardRatios/HR_FullSampleGS.csv"

country<-"England"
path<-"/mnt/work/workbench/bwolford/flagship/AbsoluteRiskEstimation/"
biobank<-"GNH"
ancestry<-"SAS"
output_dir<-"/mnt/work/workbench/bwolford/intervene/2022_10_06/RiskEstimates/"
full_HR_path<-"/mnt/work/workbench/bwolford/intervene/GoogleDrive/GNH_HazardRatios/HR_FullSample_GNH.csv"


country<-"'United Kingdom'"
path<-"/mnt/work/workbench/bwolford/flagship/AbsoluteRiskEstimation/"
biobank<-"UKB"
output_dir<-"/mnt/work/workbench/bwolford/intervene/2022_10_06/RiskEstimates/"
full_HR_path<-"/mnt/work/workbench/bwolford/intervene/GoogleDrive/UKB_HazardRatios/HR_FullSample_UKB_AllAncestries.csv"
df<-fread(full_HR_path)
ancestries<-unique(df$Ancestry)
for (i in 1:length(ancestries)){
 df%>% filter(str_detect(Ancestry,ancestries[i])) %>% select(-"Ancestry") %>% fwrite(paste0("/mnt/work/workbench/bwolford/intervene/GoogleDrive/UKB_HazardRatios/HR_FullSample_UKB_",ancestries[i],".csv"))
  ancestry<-ancestries[i]
  full_HR_path<-paste0("/mnt/work/workbench/bwolford/intervene/GoogleDrive/UKB_HazardRatios/HR_FullSample_UKB_",ancestries[i],".csv")
  system(paste0("Rscript /mnt/work/workbench/bwolford/flagship/AbsoluteRiskEstimation/LifetimeRiskEstimation_FullSample_Script.R --country=", country," --path=",path, " --output_dir=", output_dir, " --biobank=", biobank, " --k=", nk," --hr=", full_HR_path))
}






