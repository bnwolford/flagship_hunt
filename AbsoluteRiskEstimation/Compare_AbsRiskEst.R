library(ggplot2)
library(data.table)
library(dplyr)
library(rcartocolor)

hr_phenos <- c("C3_PROSTATE","C3_BREAST", "K11_APPENDACUT", "J10_ASTHMA", "I9_AF", "C3_BREAST", "I9_CHD", "C3_COLORECTAL", "G6_EPLEPSY", "GOUT", "COX_ARTHROSIS", "KNEE_ARTHROSIS", "F5_DEPRESSIO", "C3_MELANOMA_SKIN", "C3_PROSTATE", "RHEUMA_SEROPOS_OTH", "T1D", "T2D")
gbd_phenos <- c("Prostate cancer", "Breast cancer", "Appendicitis", "Asthma", "Atrial fibrillation and flutter", "Breast cancer", "Ischemic heart disease", "Colon and rectum cancer", "Idiopathic epilepsy", "Gout", "Osteoarthritis hip", "Osteoarthritis knee", "Major depressive disorder", "Malignant skin melanoma", "Prostate cancer", "Rheumatoid arthritis", "Diabetes mellitus type 1", "Diabetes mellitus type 2")

#for now remove total cancer because estbb has weird group1 and group 2 lifetime risk flip
biobanks<-c("EstBB","HUNT","UKBiobank","FinnGen")

rel_differences<-data.frame(NULL)
differences<-data.frame(NULL)
for(j in 1:length(hr_phenos)){
  for (b in 1:length(biobanks)){
  file<-paste0("/mnt/work/workbench/bwolford/intervene/results/",hr_phenos[j],"_LifetimeRisk_BootstrappedConfidenceIntervals_",biobanks[b],".csv")
  if (file.exists(file)){
    df<-fread(file)
    
    #Average PRS cumulative risk at age 80
    avg<-df %>% filter(Group %in% c("Group6") & str_detect(Age,"75 to 79")) %>% select(CIneg,CIpos,LifetimeRisk)
    avg$biobank<-biobanks[b]
    avg$pheno<-gbd_phenos[j]
    avg$comparison<-"40-60%"
    
    top<-df %>% filter(Group %in% c("Group11") & str_detect(Age,"75 to 79")) %>% select(CIneg,CIpos,LifetimeRisk)
    top$biobank<-biobanks[b]
    top$pheno<-gbd_phenos[j]
    top$comparison<-">99%"
    
    bottom<-df %>% filter(Group %in% c("Group1") & str_detect(Age,"75 to 79")) %>% select(CIneg,CIpos,LifetimeRisk)
    bottom$biobank<-biobanks[b]
    bottom$pheno<-gbd_phenos[j]
    bottom$comparison<-"<1%"
    
    if (nrow(bottom)==1 & nrow(top)==1){
     diff<-data.frame(biobank=biobanks[b], pheno=gbd_phenos[j],comparison="top vs bottom", CIneg=(top$CIneg-bottom$CIneg), CIpos=(top$CIpos-bottom$CIpos), LifetimeRisk=(top$LifetimeRisk-bottom$LifetimeRisk))
    }
    ####need to handle when Group10 is the top percentile, maybe take max of the Group 10 and Group 11 if both exist?
  
    rel_differences<-rbind(rel_differences,avg,top,bottom)
    differences<-rbind(differences,diff)
  }}}

rel_differences$comparison<-factor(rel_differences$comparison,levels=c("<1%","40-60%",">99%"))
png(file="/mnt/work/workbench/bwolford/intervene/2022_10_06/RiskEstimates/cumulative_risk_differences.png",units="in",res=300,height=6,width=10)
ggplot(rel_differences,aes(x=LifetimeRisk,y=pheno,color=biobank)) + geom_point(size=3,alpha=0.8) +
  theme_bw() + geom_errorbarh(aes(xmin=CIneg,xmax=CIpos),height=0.4) + facet_wrap(~comparison) +
  labs(x="Cumulative Risk at age 80",y='Trait') +
  theme(title = element_text(size = 22),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 18),
        axis.title.x = element_text(size = 18),
        axis.text.x = element_text(size = 12),
        axis.title.y = element_text(size = 18),
        axis.text.y = element_text(size = 16)) +
  scale_color_manual(values=carto_pal(n=length(biobanks), name="Safe"))
  dev.off()
  
  png(file="/mnt/work/workbench/bwolford/intervene/2022_10_06/RiskEstimates/cumulative_risk_top_vs_bottom.png",units="in",res=300,height=6,width=10)
  ggplot(differences,aes(x=LifetimeRisk,y=pheno,color=biobank)) + geom_point(size=3,alpha=0.8) +
    theme_bw() + geom_errorbarh(aes(xmin=CIneg,xmax=CIpos),height=0.4) + 
    labs(x="Cumulative Risk at age 80",y='Trait') +
    theme(title = element_text(size = 22),
          legend.text = element_text(size = 16),
          legend.title = element_text(size = 18),
          axis.title.x = element_text(size = 18),
          axis.text.x = element_text(size = 12),
          axis.title.y = element_text(size = 18),
          axis.text.y = element_text(size = 16)) +
    scale_color_manual(values=carto_pal(n=length(biobanks), name="Safe"))
  dev.off()
  
  
  
  