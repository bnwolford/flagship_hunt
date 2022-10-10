#prepare GBD stats

gbd_phenos <- c("Total cancers", "Appendicitis", "Asthma", "Atrial fibrillation and flutter", "Ischemic heart disease", "Colon and rectum cancer", "Idiopathic epilepsy", "Gout", "Osteoarthritis hip", "Osteoarthritis knee", "Major depressive disorder", "Malignant skin melanoma", "Rheumatoid arthritis", "Diabetes mellitus type 1", "Diabetes mellitus type 2", "Interstitial lung disease and pulmonary sarcoidosis", "Tracheal, bronchus, and lung cancer")
gbd_bcpc <-c("Breast cancer","Prostate cancer")
  
path<-"/mnt/work/workbench/bwolford/hunt_flagship/GBD/"
output_dir<-"/mnt/work/workbench/bwolford/hunt_flagship/GBD/"

dat1<-fread(paste0(path,"IHME-GBD_2019_DATA-538dd722-1.csv"))
dat2<-fread(paste0(path,"IHME-GBD_2019_DATA-538dd722-2.csv"))
dat3<-
dat<-rbind(dat1,dat2,dat3)

#### fill in the gaps
agebins <- c("1 to 4","5 to 9","10 to 14","15 to 19","20 to 24","25 to 29","30 to 34","35 to 39","40 to 44","45 to 49","50 to 54","55 to 59","60 to 64","65 to 69","70 to 74","75 to 79")
sex<-c("Both","Female","Male")
####### Mortality

#check for NA
mortality %>% filter(is.na(val))

dat %>% filter(str_detect(metric_name,"Rate")|str_detect(metric_name,"Number")) %>% filter(str_detect(measure_name,"Deaths")) %>%
  filter(cause_name %in% gbd_phenos) %>% filter(str_detect(sex_name,"Both")) %>% select(measure_name,location_name,sex_name,age_name,cause_name,metric_name,year,val,upper,lower) %>%
  fwrite(paste0(output_dir,"GBD_Mortality.csv"))

dat %>% filter(str_detect(metric_name,"Rate")|str_detect(metric_name,"Number")) %>%filter(str_detect(measure_name,"Deaths")) %>%
  filter(cause_name %in% gbd_phenos) %>% filter(!str_detect(sex_name,"Both")) %>% select(measure_name,location_name,sex_name,age_name,cause_name,metric_name,year,val,upper,lower) %>%
  fwrite(paste0(output_dir,"Sex_Stratified_Mortality.csv"))

dat %>% filter(str_detect(metric_name,"Rate")|str_detect(metric_name,"Number")) %>%filter(str_detect(measure_name,"Deaths")) %>%
  filter(cause_name %in% gbd_bcpc) %>% 
  filter(str_detect(cause_name,"Breast cancer") & str_detect(sex_name,"Female") | str_detect(cause_name,"Prostate cancer") & str_detect(sex_name,"Male")) %>%
select(measure_name,location_name,sex_name,age_name,cause_name,metric_name,year,val,upper,lower) %>%
  fwrite(paste0(output_dir,"BreastCancerProstateCancer_Mortality.csv"))


######### Prevalence

dat %>% filter(str_detect(metric_name,"Rate")|str_detect(metric_name,"Number")) %>% filter(str_detect(measure_name,"Prevalence")) %>%
  filter(cause_name %in% gbd_phenos) %>% filter(str_detect(sex_name,"Both")) %>% select(measure_name,location_name,sex_name,age_name,cause_name,metric_name,year,val,upper,lower) %>%
  fwrite(paste0(output_dir,"GBD_Prevalence.csv"))

dat %>% filter(str_detect(metric_name,"Rate")|str_detect(metric_name,"Number")) %>%filter(str_detect(measure_name,"Prevalence")) %>%
  filter(cause_name %in% gbd_phenos) %>% filter(!str_detect(sex_name,"Both")) %>% select(measure_name,location_name,sex_name,age_name,cause_name,metric_name,year,val,upper,lower) %>%
  fwrite(paste0(output_dir,"Sex_Stratified_Prevalence.csv"))

dat %>% filter(str_detect(metric_name,"Rate")|str_detect(metric_name,"Number")) %>%filter(str_detect(measure_name,"Prevalence")) %>%
  filter(cause_name %in% gbd_bcpc) %>% 
  filter(str_detect(cause_name,"Breast cancer") & str_detect(sex_name,"Female") | str_detect(cause_name,"Prostate cancer") & str_detect(sex_name,"Male")) %>%
  select(measure_name,location_name,sex_name,age_name,cause_name,metric_name,year,val,upper,lower) %>%
  fwrite(paste0(output_dir,"BreastCancerProstateCancer_Prevalence.csv"))

##### Incidence

prev %>% filter(str_detect(metric_name,"Rate")|str_detect(metric_name,"Number")) %>%filter(str_detect(measure_name,"Incidence")) %>%
  filter(cause_name %in% gbd_phenos) %>% filter(str_detect(sex_name,"Both")) %>% select(measure_name,location_name,sex_name,age_name,cause_name,metric_name,year,val,upper,lower) %>%
  fwrite(paste0(output_dir,"GBD_Incidence.csv"))

prev %>% filter(str_detect(metric_name,"Rate")|str_detect(metric_name,"Number")) %>%filter(str_detect(measure_name,"Incidence")) %>%
  filter(cause_name %in% gbd_phenos) %>% filter(!str_detect(sex_name,"Both")) %>% select(measure_name,location_name,sex_name,age_name,cause_name,metric_name,year,val,upper,lower) %>%
  fwrite(paste0(output_dir,"Sex_Stratified_Incidence.csv"))

prev %>% filter(str_detect(metric_name,"Rate")|str_detect(metric_name,"Number")) %>%filter(str_detect(measure_name,"Incidence")) %>%
  filter(cause_name %in% gbd_bcpc) %>% 
  filter(str_detect(cause_name,"Breast cancer") & str_detect(sex_name,"Female") | str_detect(cause_name,"Prostate cancer") & str_detect(sex_name,"Male")) %>%
  select(measure_name,location_name,sex_name,age_name,cause_name,metric_name,year,val,upper,lower) %>%
  fwrite(paste0(output_dir,"BreastCancerProstateCancer_Incidence.csv"))
