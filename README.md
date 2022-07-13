# flagship
## Pipeline for the flagship project (WP2)

### Step 0: Requirements and Set-up

* These scripts assume you have plink and R (libraries: data.table, tidyverse/dplyr) installed on your biobanks system. 

* In each step, I have highlighted what parts of the scripts need to be adjusted.

* If there are errors in the scripts please let me (Bradley Jermy) know on slack and I can amend.

### Step 1: Define phenotypes 

#### Step 1a: Run UKBPhenotyper.R
* Run UKBBPhenotyper.R (can be found within Phenotyping folder) to define the phenotypes required for this analysis. Note this will require you to also download the file UKBB_definitions_demo_TEST.csv. The script assumes ICD code individual level data is in long format. 

* Please check your ICD codes are in the correct format within the individual level data before applying the code. 

* Phenotypes of interest after running this script are (you can subset the file to just these):
    1. C3_CANCER (All cancers)
    2. K11_APPENDACUT (Appendicitis)
    3. J10_ASTHMA (Asthma)
    4. I9_AF (Atrial fibrillation)
    5. C3_BREAST (Breast cancer
    6. I9_CHD (Coronary heart disease)
    7. C3_COLORECTAL (Colorectal cancer)
    8. G6_EPLEPSY (Epilepsy)
    9. GOUT (Gout)
    10. COX_ARTHROSIS (Hip osteoarthritis) 
    11. KNEE_ARTHROSIS (Knee osteoarthritis)
    12. F5_DEPRESSIO (Major depressive disorder)
    13. C3_MELANOMA_SKIN (Malignant skin melanoma)
    14. C3_PROSTATE (Prostate cancer)
    15. RHEUMA_SEROPOS_OTH (Rheumatoid arthritis)
    16. T1D (Type 1 diabetes)
    17. T2D (Type 2 diabetes)
    18. ILD (Interstitial lung disease)
    19. C3_BRONCHUS_LUNG (Lung cancer)

* To make sure the script can separate between ICD 10 and ICD 9 codes, it places a '10x' or a '9x' at the start of the code according to whether the code is ICD9 or 10. 

* The regex pattern then searches for strings starting with 10x or 9x before identifying the codes themselves. 

#### Step 1b: Construct phenotype file
* The phenotype file will include all variables required to construct hazard ratios, i.e. genetic principal components, ancestry allocation, variables that allow you to control for technical artefacts (genotype batch, cohort etc.). Please use whatever variables currently exist within your biobank. These do not need to be harmonised.

* To learn more about the variables that are required, please read: https://docs.google.com/document/d/1GbZszpPeyf-hyb0V_YDx828YbM7woh8OBJhvzkEwo2g/edit Note: you can currently ignore the requirement to define smoking and education. 

### Step 2: Download Adjusted Summary Statistics 

* All pre-adjusted summary statistics can be found [here](https://figshare.com/account/home#/projects/131369).

* Download the pre-adjusted summary statistics created using MegaPRS that correspond to the build of the genome for your biobank.

* hg19 files found [here](https://figshare.com/account/projects/131369/articles/19409501) / hg38 files found [here](https://figshare.com/account/projects/131369/articles/19093313). Hg19 contains rsids whereas hg38 contains variant IDs in the format CHR_POS_REF_ALT. 

### Step 3: Merge variant IDs to match those within the .bim file. 

#### Note - only to be performed if the genome build is hg38 and the variant ID structure is CHR_POS_REF_ALT.

* Run script hg38_biobankadjustments.R

* This script is required for all biobanks which have a genome build of hg38/where the variant ID structure within the bim file is in the form CHR_POS_REF_ALT. 

* Due to the nature of the conversion to hg38 within the pre-adjusted summary statistics, we cannot be sure what the REF or ALT allele is within your bim file as we have relied on those within the summary statistics themselves. As such, both combinations are currently saved within the sum stats. 

* The script reads in your bim file and identifies which variant ID is contained within your bim file and the summary statistics and creates a final snp list for the adjusted summary statistics. 

* For this script to work you will have to:
    1. Line 9 - Adjust the path within the 'bim' variable to read in your biobanks bim file.
    2. Line 19 - Adjust the path to the location of the adjusted summary statistics downloaded during step 1. Note: I have assumed you have mainted the same filename structure as when you downloaded the files. 
    3. Line 40 - Adjust the path to the location where you wish to save the amended summary statistics. Recommended to just overwrite the original summary statistics so keep the same path as specified in line 19. You may want to test the code before saving initially just to make sure it is behaving as expected. 

### Step 4: Compute PRS using Plink

* **Note: if you use a job scheduler which allows multiple jobs to be submitted, please adjust the script below and create a single job for each phenotype. This script assumes an interactive job and loops over phenotypes which will be much much slower than running in parallel.**

#### If your plink files are not split by chromosome
* Run script GeneratePRS.sh

* For this script to work you will have to:
    1. Line 8 - Specify the path where PRS are to be saved.
    2. Line 9 - Specify the path of the downloaded summary statistics from Step 1.
    3. Line 10 - Specify the path to the directory where the allele frequencys are stored.
    4. Line 11 - Optional: If you want to speed up computation and calculate PRS on a subset of SNPs, i.e. hapmap, you can specify the list here. 
    5. Line 12 - Specify the path to the directory where the genotypes are stored (it is likely line 8 and line 10 will have the same path). 
    6. Line 19 - Specify the path to plink. If you do not have plink2 installed, also change to the version used by your biobank.
    7. Line 20 - Change 'genotype_plink_files' to the name of your genotype files.
    8. Line 21 - Optional: Name of your SNP list file. 
    9. Line 22 - Change 'frequency_file' to the name of your file containing allele frequencies.
    10. Line 23 - Select 19 or 38 depending on the build of your biobanks genome. Also remove square brackets surrounding the number.
    
* *Note: If you do not have a file containing allele frequencies, we recommend producing one before computing PRS as otherwise plink will redo this step for every phenotype.*

#### If your plink files are split by chromosome


### Step 5: Calculate hazard ratios between PRS and Phenotype - survival analysis

* **Note: this script has four assumptions.** 
    1) You have a phenotype file with case control assignments and dates of first diagnosis for each phenotype. More generally, you have the same column names as specified in the phenotype file specifications. 
    2) You have been able to allocate genetic ancestry for the participants within your biobank. This does not have to be harmonized, any approach taken by your biobank will suffice for this analysis.
    3) You have kept the same short hand names for the phenotypes as within [FinnGen](https://docs.google.com/spreadsheets/d/1DNKd1KzI8WOIfG2klXWskbCSyX6h5gTu/edit#gid=334983519) (column B).
    4) You have kept the same naming structure for the PRS files as when you downloaded them.

* Run script HazardRatioCalculation.R - Can be found within the PRS folder. 

* For this script to work you will have to:
    1. Amend each mention of "path/to/pheno_file" to the path of your phenotype file. *See assumption 1.*
    2. Amend each mention of "path/to/PRS" to the path of your score file.
    3. Lines 31,113,220,305 - Change 'ENTER_ID' to correspond to the ID column within your biobank, i.e. FINNGENID. *Note: there may be redundant info taken by choosing the two columns. If this is the case, just select one of the IDs.*
    5. Lines 41,123,230,315 - If each section does not subset individuals to those of european ancestries, remove and subset to european ancestry yourself.
    6. Amend each mention of "file/path/to/output_ to your chosen output. 

**Send results to bradley.jermy@helsinki.fi :)**

--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
### UKBBPhenotyper.R
v1.0 under development :warning:  
Script by Bradley Jermy can work across biobanks as long as the ICD code individual level data is in long format.  
It extracts information from UKBB_definitions_demo_TEST.csv to extract all the endpoints for the flagship project.  
dependencies: data.table, tidyverse. 

Please check your ICD codes are in the correct format within the individual level data before applying the code. 

To make sure the script can separate between ICD 10 and ICD 9 codes, it places a '10x' or a '9x' at the start of the code according to whether the code is ICD9 or 10. 

The regex pattern then searches for strings starting with 10x or 9x before identifying the codes themselves. 

### biobank_summary_plot.R
Script by Brooke Wolford to analyze endpoint metrics collected from biobanks. 
Reads in data from google drive excel file. 
