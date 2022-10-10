####wrapper for countries and studies

#HUNT
bb_gbd_phenos <- c("Appendicitis", "Asthma", "Atrial fibrillation and flutter", "Breast cancer", "Ischemic heart disease", "Colon and rectum cancer", "Idiopathic epilepsy", "Gout", "Osteoarthritis hip", "Osteoarthritis knee", "Major depressive disorder", "Malignant skin melanoma", "Prostate cancer", "Rheumatoid arthritis", "Diabetes mellitus type 2", "Interstitial lung disease and pulmonary sarcoidosis", "Tracheal, bronchus, and lung cancer")
bb_hr_phenos <- c("K11_APPENDACUT", "J10_ASTHMA", "I9_AF", "C3_BREAST", "I9_CHD", "C3_COLORECTAL", "G6_EPLEPSY", "GOUT", "COX_ARTHROSIS", "KNEE_ARTHROSIS", "F5_DEPRESSIO", "C3_MELANOMA_SKIN", "C3_PROSTATE", "RHEUMA_SEROPOS_OTH", "T2D", "ILD", "C3_BRONCHUS_LUNG")
#remove T1D and all cancers for HUNT
country<-"Norway"
path<-"/mnt/work/workbench/bwolford/hunt_flagship/AbsoluteRiskEstimation/" #must have final backslash
full_HR_path<-"/mnt/work/workbench/bwolford/intervene/2022_10_06/HR_FullSample_HUNT.csv" #must go straight to csv
output_dir<-"/mnt/work/workbench/bwolford/intervene/2022_10_06/RiskEstimates/" #must have final backslash
biobank<-"HUNT"
source(paste0(path,"LifetimeRiskEstimation_FullSample.R"))
#ESTBB