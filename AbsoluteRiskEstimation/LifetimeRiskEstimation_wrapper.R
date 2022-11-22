#!/usr/bin/env Rscript 

print(.libPaths())

library(system)
####wrapper for countries and studies

#bb_gbd_phenos <- c("Total cancers", "Appendicitis", "Asthma", "Atrial fibrillation and flutter", "Breast cancer", "Ischemic heart disease", "Colon and rectum cancer", "Idiopathic epilepsy", "Gout", "Osteoarthritis hip", "Osteoarthritis knee", "Major depressive disorder", "Malignant skin melanoma", "Prostate cancer", "Rheumatoid arthritis", "Subarachnoid hemorrhage", "Diabetes mellitus type 1", "Diabetes mellitus type 2")
#bb_hr_phenos <- c("C3_CANCER", "K11_APPENDACUT", "J10_ASTHMA", "I9_AF", "C3_BREAST", "I9_CHD", "C3_COLORECTAL", "G6_EPLEPSY", "GOUT", "COX_ARTHROSIS", "KNEE_ARTHROSIS", "F5_DEPRESSIO", "C3_MELANOMA_SKIN", "C3_PROSTATE", "RHEUMA_SEROPOS_OTH", "I9_SAH", "T1D", "T2D")

nk<-50

#HUNT
country<-"Norway"
path<-"/mnt/work/workbench/bwolford/hunt_flagship/AbsoluteRiskEstimation/" #must have final backslash
full_HR_path<-"/mnt/work/workbench/bwolford/intervene/2022_10_06/HR_FullSample_HUNT.csv" #must go straight to csv
output_dir<-"/mnt/work/workbench/bwolford/intervene/2022_10_06/RiskEstimates/" #must have final backslash
biobank<-"HUNT"
#source(paste0(path,"LifetimeRiskEstimation_FullSample.R"))

system("Rscript /mnt/work/workbench/bwolford/flagship/AbsoluteRiskEstimation/LifetimeRiskEstimation_FullSample_Script.R --country Norway --path /mnt/work/workbench/bwolford/hunt_flagship/AbsoluteRiskEstimation/ --output_dir /mnt/work/workbench/bwolford/intervene/2022_10_06/RiskEstimates/ --biobank HUNT --k 100 --hr /mnt/work/workbench/bwolford/intervene/2022_10_06/HR_FullSample_HUNT.csv")

#BBJ

#ESTBB
country<-"Estonia"
path<-"/mnt/work/workbench/bwolford/hunt_flagship/AbsoluteRiskEstimation/" #must have final backslash
full_HR_path<-"/mnt/work/workbench/bwolford/intervene/2022_10_06/HR_FullSample_HUNT.csv" #must go straight to csv
output_dir<-"/mnt/work/workbench/bwolford/intervene/2022_10_06/RiskEstimates/" #must have final backslash
biobank<-"EstBB"
#source(paste0(path,"LifetimeRiskEstimation_FullSample.R"))

system(paste0("Rscript /mnt/work/workbench/bwolford/flagship/AbsoluteRiskEstimation/LifetimeRiskEstimation_FullSample_Script.R --country ", country,
" --path ",path, " --output_dir ", output_dir, " --biobank", biobank, " --k", nk," --hr ", full_HR_path))









