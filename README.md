# de-Groot-et-al-2025_CJFAS
Data and code accompanying de Groot et al. (2025) CJFAS


Description of the data and file structure:

The Atlantic cod population data for this study were collected from the age-0 2019 cohort found in Newman Sound, Bonavista Bay, Newfoundland in Terra Nova National Park. A regular monitoring program to assess cod recruitment has been conducted here since 1996.


Files and variables:

  File: Chung_Cresp_Oikos.csv
  
  Description: 
  
  Variables
  
  HCID: 
  
  d13C:
  
  d18O:
  
  PP_ID: Population ID
  
  location: 
  
  year:
  
  Catch.Date:
  
  length:
  
  Weight_g:
  
  Diet.mean:
  
  Diet.sd:
  
  DIC.mean:
  
  DIC.sd:
  
  Salinity.mean:
  
  Salinity.sd:
  
  Growth Ring:
  
  PredCresp:
  
  PredCresp.sd:
  
  FMR:
  
  FMR.sd:
  
  Temp:
  
  PredT.sd:
  

  File: deGroot_Raw_NFD.csv
  
  Description: 
  
  Variables
  
  HCID:
  
  d13Coto:
  
  d18Ooto:
  
  PP_ID: Population ID
  
  Location:
  
  Sample date:
  
  Month:
  
  Year:
  
  Julian date:
  
  SL_mm:
  
  Weight_g:
  
  N%:
  
  C%:
  
  d15N:
  
  d13C:
  
  CN_weight:
  
  d13C_diet:
  
  d13C_DIC:
  
  DIC.sd:



Code/software

RStudio was used to run the code for the analyses associated with this manuscript. The below packages are needed to run the files. The code for the analyses has been submitted in .R format and is available at https://github.com/ValescadG/de-Groot-et-al-2025_CJFAS_Submission/blob/5e110e1232ba77bc9cd370a4a05d4d76bff72ca4/Valesca_ACcode.R.



library(ggplot2)

library(nlme)

library(AICcmodavg) #for the predictSE.gls function

library(emmeans)

library(dplyr)


