#Note: this script assumes the names of the phenotypes are consistent with FinnGen: Refer to https://docs.google.com/spreadsheets/d/1DNKd1KzI8WOIfG2klXWskbCSyX6h5gTu/edit#gid=334983519

#Libraries
library(data.table)
library(dplyr)
library(lubridate)
library(survival)
library(cmprsk)

setwd("/gpfs/space/GI/GV/Projects/Intervene_grs/Absrisk_calibration")
#artikli supplement https://docs.google.com/spreadsheets/d/1u31O_aVrB--SnVJn98mVvu9kexdSWbFsTrD4s_Mh5O8/edit?pli=1#gid=0
abs=read.csv("/gpfs/space/GI/GV/Projects/Intervene_grs/Absrisk_calibration/Supplementary_Table_12_INTERVENE.csv",header=T,sep=";")
abs=abs[abs$Country %in% "Estonia",]

phenocols <- c("C3_CANCER", "K11_APPENDACUT", "J10_ASTHMA", "I9_AF", "C3_BREAST", "I9_CHD", "C3_COLORECTAL", "G6_EPLEPSY", "GOUT", "COX_ARTHROSIS", "KNEE_ARTHROSIS", "F5_DEPRESSIO", "C3_MELANOMA_SKIN", "C3_PROSTATE", "RHEUMA_SEROPOS_OTH", "T1D", "T2D", "ILD", "C3_BRONCHUS_LUNG")
prscols <- c("AllCancers", "Appendicitis", "Asthma", "Atrial_Fibrillation", "Breast_Cancer", "CHD", "Colorectal_Cancer", "Epilepsy","Gout", "Hip_Osteoarthritis", "Knee_Osteoarthritis","MDD", "Melanoma", "Prostate_Cancer", "Rheumatoid_Arthritis", "T1D","T2D", "ILD", "Lung_Cancer")

percentiles <- list(c(0,0.2,0.4,0.6,0.8,0.9,0.95,1)) #5%


results <- c()

for(i in 1:length(phenocols)){
  for(p in percentiles){
    
 
 print(phenocols[i])
    print(prscols[i])
    
    #Read in phenotype file
  pheno <- fread(input="/gpfs/space/home/kristil2/INTERVENE/input_for_flagship_290422.txt", select=c("SEX","Person.agreementDate","ANCESTRY","VKOOD1","VKOOD2","eid","DATE_OF_BIRTH","PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10",phenocols[i],paste0(phenocols[i],"_DATE"),"END_OF_FOLLOWUP"), data.table=FALSE)
    
    pheno[,paste0(phenocols[i],"_DATE")] <- as.Date(pheno[,paste0(phenocols[i],"_DATE")], origin = "1970-01-01")
    
    #Read in PRS scores
    PRS <- fread(input=paste0("/gpfs/space/GI/GV/Projects/Intervene_grs/results_HR_v2/",prscols[i],"_PRS.sscore"), data.table=FALSE)
    
    #Subset columns to the IDs and score only. Note: columns FID or IID may be redundant and can be removed if necessary. Kept in to avoid bugs.
    PRS <- PRS[,c("IID","SCORE1_SUM")]
    if ( phenocols[i]=="C3_CANCER"){
	PRS[,2]=-1*PRS[,2] }

    
#This is important, as some sumstats have estbb older batch included, so for rough estimate, first 50k are excluded from some analyses
amc=read.table("/gpfs/space/home/kristil2/INTERVENE/outliers_based_on_1KG_PCA.txt")


#finalvcodes, et saada oige vkood phenole külge, millega genotud
pheno$ID=ifelse(pheno$VKOOD1 %in% PRS[,1],as.character(pheno$VKOOD1),as.character(pheno$VKOOD2))

if (sum(-which(pheno$ID %in% as.character(amc[,1])))!=0){
pheno=pheno[-which(pheno$ID %in% as.character(amc[,1])),] }
pheno$k50=ifelse(pheno$Person.agreementDate>"2017-12-31",0,1)



if (phenocols[i] %in% c("T2D","I9_CHD","I9_HEARTFAIL_NS","N14_CHRONKIDNEYDIS","AUD_SWEDISH","I9_AF")){
pheno=pheno[pheno$k50==0,]}else{
     pheno=pheno
      }


    #Rename ID column to the name of the ID column in the phenotype file
    colnames(PRS) <- c("ID", paste0(prscols[i],"_prs"))
    
    #left_join to the phenotype file
    pheno <- left_join(pheno, PRS)
    
    pheno <- subset(pheno, !(is.na(pheno[[paste0(phenocols[i])]]) | is.na(pheno[[paste0(prscols[i],"_prs")]])))
    
    #Subset to those of european ancestry/those that have principal components calculated for EUROPEAN ancestry, i.e. within ancestry principal components, not global genetic principal components.
    #As we have been unable to use the standardised method for computing ancestry, if you have this information available from your centralised QC please use this. 
    #Feel free to subset using your own code: only provided as a reminder.
    pheno <- subset(pheno, ANCESTRY=='EUR')
    
    #Assign PRS into percentiles
    q <- quantile(pheno[[paste0(prscols[i],"_prs")]], probs=p)
    
    pheno[[paste0(prscols[i],"_group")]] <- cut(pheno[[paste0(prscols[i],"_prs")]], q, include.lowest=TRUE,
                                                labels=paste("Group",1:(length(p)-1)))
    
    #Make all necessary variables factors
    pheno[[paste0(prscols[i],"_group")]] <- as.factor(pheno[[paste0(prscols[i],"_group")]])
    pheno[[paste0(prscols[i],"_group")]] <- relevel(pheno[[paste0(prscols[i],"_group")]], ref=paste("Group",3))
    
    #Specify age as either the Age at Onset or End of Follow-up (if not a case)
    pheno$AGE <- ifelse(pheno[[phenocols[i]]]==1, time_length(difftime(pheno[[paste0(phenocols[i],"_DATE")]], pheno$DATE_OF_BIRTH), 'years'), time_length(difftime(pheno$END_OF_FOLLOWUP, pheno$DATE_OF_BIRTH), 'years'))
    
    #Adjust to censor at age 80
    pheno[[paste0(phenocols[i])]] <- ifelse(pheno[[paste0(phenocols[i])]]==1 & pheno$AGE > 80, 0, pheno[[paste0(phenocols[i])]])
    pheno$AGE <- ifelse(pheno$AGE > 80, 80, pheno$AGE)


   #competing event:
   pheno$surm=ifelse( pheno$END_OF_FOLLOWUP!="2020-09-14",1,0) #as last date if linking for death registry
   pheno$pheno_comp=ifelse(pheno[,phenocols[i]]==0,2*pheno$surm,pheno[,phenocols[i]])

    
    #Perform survival analysis
    survival <- coxph(as.formula(paste0("Surv(AGE,",phenocols[i],") ~ ",prscols[i],"_group + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10")), data=pheno, na.action=na.exclude)
	# cuminc    
	pheno$gr=pheno[,paste0(prscols[i],"_group")]

if(phenocols[i] %in% c("C3_PROSTATE","C3_BREAST")) {
	inc_full=timepoints(with(pheno,cuminc(AGE, pheno_comp, gr, cencode=0, na.action=na.omit)),seq(0,80,5))$est[1:7,] }else{

	inc_female=timepoints(with(pheno[pheno$SEX=="female",],cuminc(AGE, pheno_comp, gr, cencode=0, na.action=na.omit)),seq(0,80,5))$est[1:7,]
	inc_male=timepoints(with(pheno[pheno$SEX=="male",],cuminc(AGE, pheno_comp, gr, cencode=0, na.action=na.omit)),seq(0,80,5))$est[1:7,]
	inc_full=timepoints(with(pheno,cuminc(AGE, pheno_comp, gr, cencode=0, na.action=na.omit)),seq(0,80,5))$est[1:7,] }


    #Extract hazard ratios, betas, standard errors and p-vals - in the first instance extract all results, for the latter just take the 
      phenotype <- rep(phenocols[i],7)
      prs <- rep(prscols[i],7)
      sub <-c("full","female","male")
      groups <- (c("40-60%","< 20%","20-40%","60-80%","80-90%","90-95%","> 95%"))

if (phenocols[i] %in% c("C3_PROSTATE","C3_BREAST")) {
      result <- cbind(sub[1],phenotype, prs,groups, inc_full)
      results <- rbind(results, result)}
else{

      result <- cbind(sub[1],phenotype, prs,groups, inc_full)
      result2 <- cbind(sub[2],phenotype, prs,groups,inc_female)
      result3 <- cbind(sub[3],phenotype, prs,groups,inc_male)

      results <- rbind(results, result,result2,result3) }
        
  }
}

write.csv(results, file="Inc_FullSample_EstBB211122.csv")
