# hunt_flagship
Documentation for scripts for hunt biobank specific scripts from intervene flagship project


#hunt_specific
Sandbox for script development up to October 2022

HUNT_Phenotyper.R
HUNT_baseline_stats.R

#PRS
HUNT versions of the scripts in flagship/PRS


UKBBPhenotyper.R
v1.0 under development
Script by Bradley Jermy can work across biobanks as long as the ICD code individual level data is in long format.
It extracts information from UKBB_definitions_demo_TEST.csv to extract all the endpoints for the flagship project.
dependencies: data.table, tidyverse.

Please check your ICD codes are in the correct format within the individual level data before applying the code.

To make sure the script can separate between ICD 10 and ICD 9 codes, it places a '10x' or a '9x' at the start of the code according to whether the code is ICD9 or 10.

The regex pattern then searches for strings starting with 10x or 9x before identifying the codes themselves.

biobank_summary_plot.R
Script by Brooke Wolford to analyze endpoint metrics collected from biobanks. Reads in data from google drive excel file.
