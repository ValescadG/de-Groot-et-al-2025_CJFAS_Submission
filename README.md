# de-Groot-et-al-2025_CJFAS
Data and code accompanying de Groot et al. (2025) CJFAS


Description of the data and file structure:

The Atlantic cod population data for this study were collected from the age-0 2019 cohort found in Newman Sound, Bonavista Bay, Newfoundland in Terra Nova National Park. A regular monitoring program to assess cod recruitment has been conducted here since 1996.


Files and variables:

  File: Chung_Cresp_Oikos.csv
  
  Variables
  
  HCID: Capture ID
  
  d13C: Otolith Carbon Stable Isotope Ratio (δ¹³C) — carbon isotope composition of the otolith carbonate relative to the VPDB standard.
  
  d18O: Otolith Oxygen Stable Isotope Ratio (δ¹⁸O) — oxygen isotope composition of the otolith carbonate relative to the VPDB standard.
  
  PP_ID: Population Identifier — code indicating the population or sampling group to which the individual fish belongs.
  
  location: Sampling Location — geographic site where the fish was captured.
  
  year: Sampling Year — year in which the fish was collected.
  
  Catch.Date: Capture Date — specific date the fish was sampled.
  
  length: Fish Standard Length — standard body length of the fish.
  
  Weight_g: Fish Wet Weight — whole body mass of the fish in grams.
  
  Diet.mean:
  
  Diet.sd:
  
  DIC.mean: Mean Dissolved Inorganic Carbon δ¹³C — average carbon isotope value of dissolved inorganic carbon in the surrounding water.
  
  DIC.sd: Standard Deviation of Dissolved Inorganic Carbon δ¹³C.
  
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
    
  Variables
  
  HCID: Capture ID
  
  d13Coto: Otolith Carbon Stable Isotope Ratio (δ¹³C) — carbon isotope composition of the otolith carbonate relative to the VPDB standard.
  
  d18Ooto: Otolith Oxygen Stable Isotope Ratio (δ¹⁸O) — oxygen isotope composition of the otolith carbonate relative to the VPDB standard.
  
  PP_ID: Population Identifier — code indicating the population or sampling group to which the individual fish belongs.
  
  Location: Sampling Location — geographic site where the fish was captured.
  
  Sample date:
  
  Month:
  
  Year: Sampling Year — year in which the fish was collected.
  
  Julian date:
  
  SL_mm: Fish Standard Length — standard body length of the fish.
  
  Weight_g: Fish Wet Weight — whole body mass of the fish in grams.
  
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


