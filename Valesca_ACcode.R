setwd("/Users/trueman/Desktop/R scripts/OtoMet/cod met data")

library(ggplot2)
library(nlme)
library(AICcmodavg) #for the predictSE.gls function
library(emmeans)
library(dplyr)

View(Valesca_Raw_NFD)

#load databases
NFD_data<- deGroot_Raw_NFD
Skg_data<-Chung_Cresp_Oikos
print(Skg_data)


d18Ow<--2.25 # estimated water d18O value

NFD_data<-NFD_data%>%
  mutate(Cresp = (d13Coto-d13C_DIC) / (d13C-d13C_DIC)) %>% # Calculate Cresp
  mutate(Temp =  ((d18Ooto-d18Ow)-3.90)/-0.2)#  Temp equation for cod from Hoie et al)

# calibration parameters from Chung et al 2019
slopeK <- 0.0088
C <-  0.4



NFD_data<-NFD_data%>%
 mutate (FMR = (log(1-(Cresp/C)))/-slopeK) # Cresp to oxygen consumption


Select_NFD <- NFD_data%>%
  select(HCID, PP_ID, Weight_g, Temp, FMR)

Select_Skg <- Skg_data%>%
  select(HCID, PP_ID, Weight_g, Temp, FMR)

# Select and combine datafiles (ID, PPID, weight, temp, FMR)
Combined_data <-rbind(Select_NFD, Select_Skg)



# z-score function 
# allows to take data points drawn from populations with different means and 
# standard deviations and place them on a common scale

Z.score <- function(data, colname){
  # Get mean and SD
  mu <- mean(data, na.rm = TRUE)
  st_dev <- sd(data, na.rm = TRUE)
  Z_score <- data.frame()
  for(i in 1:length(data)){
    z <- (data[i] - mu) / st_dev
    Z_score <- rbind(Z_score, z)
  }
  colnames(Z_score) <- paste("z", colname, sep = "_")
  return(Z_score)
}


K<-8.62*10^-5 #Boltzman constant

Combined_data<-Combined_data%>%
  mutate(InvTemp=1/((K)*(Temp+273)))%>% # inverse temperature for modelling
  mutate(lnFMR=log(FMR))%>% # natural lof FMR for modelling
  mutate(lnM=log(Weight_g)) # natural lof of mass



# Add mass scaled FMR with a common allometric scaling of between -0.1:-0.2 (0.9 to 0.8) for whole organism scaling - standardise mass to common 10g
alpha=-0.13 # mass scaling estimate used in paper
Combined_data<-Combined_data%>%
  mutate(FMR_MS=FMR*((10/Weight_g)^alpha) )%>%
  mutate(lnFMR_MS=log(FMR_MS))

# Quick check of non mass scaled data with plots
test_Temp <- ggplot(data = Combined_data, aes(x = InvTemp, y = lnFMR, group = PP_ID)) +
  geom_smooth(method = "lm", se = TRUE) +
  geom_point(aes(col = as.factor(PP_ID)))

test_Mass <- ggplot(data = Combined_data, aes(x = lnM, y = lnFMR, group = PP_ID)) +
  geom_smooth(method = "lm", se = TRUE) +
  geom_point(aes(col = as.factor(PP_ID)))

test_Mass_Temp <- ggplot(data = Combined_data, aes(x = lnM, y = Temp, group = PP_ID)) +
  geom_smooth(method = "lm", se = TRUE) +
  geom_point(aes(col = as.factor(PP_ID)))



# Quick check of data with plots_Mass scaled
test_Temp_MS <- ggplot(data = Combined_data, aes(x = InvTemp, y = lnFMR_MS, group = PP_ID)) +
  geom_smooth(method = "lm", se = TRUE) +
  geom_point(aes(col = as.factor(PP_ID)))

test_Mass_MS <- ggplot(data = Combined_data, aes(x = lnM, y = lnFMR_MS, group = PP_ID)) +
  geom_smooth(method = "lm", se = TRUE) +
  geom_point(aes(col = as.factor(PP_ID)))



# Print the plot
print(test_Temp)
print(test_Mass)
print(test_Mass_Temp)# mass-temp covariance


print(test_Temp_MS)
print(test_Mass_MS)

#set up z_score & arrhenius equation values in the database
K<-8.62*10^-5 # for inverse temp calculation

#z.score mass and inverse temp for modelling
Combined_data$z_inv_temp <- Z.score(Combined_data$InvTemp, "inv_temp")$z_inv_temp
Combined_data$z_ln_mass <- Z.score(Combined_data$lnM,"ln_mass")$z_ln_mass

#collect mean and std dev of ln inv temp and ln mass for back-correction
mean_Inv.Temp<-mean(Combined_data$InvTemp, na.rm=TRUE)
mean_lnMass<-mean(Combined_data$lnM, na.rm=TRUE)
sd_Inv.Temp<-sd(Combined_data$InvTemp, na.rm=TRUE)
sd_lnMass<-sd(Combined_data$lnM, na.rm=TRUE)


Combined_data<-na.omit(Combined_data)# remove IDs with NA
##Q1: Find best fitting model for the data and determine whether there are sig
## diffs between mean FMR (intercept) and/or thermal sensitivity of FMR (slope)
## across PPIDs

GLS_1<- gls(lnFMR ~ z_ln_mass * PP_ID + z_inv_temp*PP_ID, method="ML", data = Combined_data) #hypothesis model, is there a relationship between mass/temp/PPID, have more variable data with FJ & NS
GLS_2<- gls(lnFMR ~ z_ln_mass * PP_ID + z_inv_temp*PP_ID, weights=varIdent(form=~1|PP_ID), method="ML",  data = Combined_data) #pick the best out of 1 & 2
anova(GLS_1, GLS_2) 
summary(GLS_2)
# returns a +ve effect of mass on (mass-specific) FMR which is unlikely - probably due to covarianvce between mass and temperature

#variance of NS & FJ are x4 that of P1 and P3, the variances were not equal between groups,
#We model the variance structure of the different groups w varIdent and include the groups as a fixed effect
#GLS_2 has a lower AIC and is a significantly better model (p<0.0001)

GLS_3<- gls(lnFMR ~ z_ln_mass * PP_ID + z_inv_temp, weights=varIdent(form=~1|PP_ID), method="ML",  data = Combined_data) 
anova(GLS_2, GLS_3) 
# GLS_3 has a lower AIC but is not a significantly better model (p=0.2223)
# We choose GLS_3 because it is a simpler model and because the additional 
# PPID * temp interaction adds an additional 2 units to the AIC 

# What does the fit look like
summary(GLS_3)
qqnorm(GLS_3) #nice fit 

#calculate estimated marginal means to check 
emm <- emmeans(GLS_3, ~ PP_ID)
print(emm)

# Perform pairwise comparisons (tukey adjustment)
pairwise_comparisons <- pairs(emm)
print(pairwise_comparisons)


##A1: There are sig diffs in mean FMR between certain groups -
# P1-FJ (0.0003), P1-NS (0.0003), P1-P3(0.0049) only, 
# There is a sig effect of temperature (p=0) and mass (p=0.0005) on FMR, 
# There is a sig diff in mass between P3 and other groups
# the interaction term of PPID on temperature (GLS_2) is not included
# in the best fitting model, therefore there are no significant 
# differences between the slopes (FMR thermal sensitivity) across PPIDs


# Same model with mass scaled to 10g  -i.e. no mass in the regression

GLS_MS_1<- gls(lnFMR_MS ~ z_inv_temp*PP_ID, method="ML", data = Combined_data) #hypothesis model, is there a relationship between mass/temp/PPID, have more variable data with FJ & NS
GLS_MS_2<- gls(lnFMR_MS ~ z_inv_temp*PP_ID, weights=varIdent(form=~1|PP_ID), method="ML",  data = Combined_data) #pick the best out of 1 & 2
anova(GLS_MS_1, GLS_MS_2) 
summary(GLS_MS_2)
#variance of NS & FJ are x4 that of P1 and P3, the variances were not equal between groups,
#We model the variance structure of the different groups w varIdent and include the groups as a fixed effect
#GLS_2 has a lower AIC and is a significantly better model (p<0.0001)

GLS_MS_3<- gls(lnFMR_MS ~ z_inv_temp, weights=varIdent(form=~1|PP_ID), method="ML",  data = Combined_data) 
anova(GLS_MS_2, GLS_MS_3) 
# GLS_2 has a  lower AIC and is a significantly better model (p=<0.001)


# What does the fit look like
summary(GLS_MS_2)
qqnorm(GLS_MS_2) #nice fit 

#calculate estimated marginal means to check 
emm <- emmeans(GLS_MS_2, ~ PP_ID)
print(emm)

# Perform pairwise comparisons (tukey adjustment)
pairwise_comparisons <- pairs(emm)
print(pairwise_comparisons)

##A1 (Mass Scaled): There are sig diffs in mean FMR between certain groups -
# P3-FJ (0.0006)
# There is a sig effect of temperature (p=0)on FMR, 
# There is a sig diff in mass between P3 and other groups
# the interaction term of PPID on temperature (GLS_MS_2) is  included
# in the best fitting model, temp slope for Fj and NS are different





##Q2: Is temperature experience significantly different between pulses (PPID: P1, P3)?
# Checked for this in a boxplot/regression in previous Rscript, 
# found significant differences (p=0.0027) 
# is there a way to check this again using the above model

# Plot experienced temperatures using this dataset, check visually 
tempgg<-ggplot(Combined_data, aes(x=PP_ID, y=Temp, fill=PP_ID)) + geom_boxplot()+
  theme_classic(base_size=19) + geom_jitter(color="black",width =0.05, height = 0.05, size=0.4, alpha=0.9) +
  labs(title = ,
       x = "Pulse",
       y = "Experienced Temperature (°C)") + theme(legend.position="right", legend.spacing.y = unit(0, "mm"), 
                                                   panel.border = element_rect(colour = "black", fill=NA), aspect.ratio = 1, axis.text = element_text(colour = 1, size = 12),
                                                   legend.background = element_blank(), legend.box.background = element_rect(colour = "black")) 
plot(tempgg)


#Check significance between groups in both datasets:
#using this dataset

testtemp<-lm(Temp~PP_ID, data=Combined_data)
summary(testtemp)
#p=0.234 between P3 & P1, not significant
#Adjusted R-squared of 0.6973

#using dataset from old code
old<-lm(Temp~Pulse, data=TempLen)
TempLen$Pulse <- as.factor((TempLen$Pulse))
summary(old)
#p=0.00288 between P3 & P1, significant
#Adjusted R-squared of 0.1475 - very low, new code/dataset a more accurate estimation

##A2 - there are no significant differences in temperature experience between pulses in Newman Sound




##Q3 Predict FMR values to compare to measured FMR and plot thermal sensitivity 
# as log2 FMR vs absolute temp in *C 
# we use log2 FMR, which means that each unit increase is a doubling of metabolic rate 
# e.g. Q10 of 2 would mean a slope of 1/10

#make sure headings are still intact 
Combined_data$z_inv_temp <- Z.score(Combined_data$InvTemp, "inv_temp")$z_inv_temp
Combined_data$z_ln_mass <- Z.score(Combined_data$lnM,"ln_mass")$z_ln_mass
head(Combined_data)
#predict new lnFMR values using chosen model
GLS_3<- gls(lnFMR ~ z_ln_mass * PP_ID + z_inv_temp, weights=varIdent(form=~1|PP_ID), method="ML",  data = Combined_data) 
FMR_predicted<-data.frame(z_ln_mass=Combined_data$z_ln_mass, 
                          z_inv_temp=Combined_data$z_inv_temp, 
                          PP_ID=Combined_data$PP_ID)
predictedlnFMR <- predict(GLS_3, newdata = FMR_predicted)
# Add the predicted values to the GLS_data dataframe
Combined_data$predictedlnFMR <- predictedlnFMR
# Print the first few rows of the updated dataframe to check the results
head(Combined_data)

#convert all values for the log2 FMR vs Temp plot
Combined_data<-Combined_data%>%
  mutate(pred_FMR= exp(predictedlnFMR))%>%
  mutate(log2_FMR = (log(FMR,2)))%>%
  mutate(log2_pred_FMR = log(pred_FMR, 2))%>%
  mutate(Temp=((1/InvTemp)/K)-273       
  )

# plot predicted & measured, add se to geom_smooth (predictedlogFMR)

#choose cold colours for NL pulses and warm colours for Norway ecotypes
colors <- c('P1'="#1F77B4", 'Fjord'="#FF7F0E",'North Sea'= "#D62728",'P3'= "#9467BD")

# define the labels for the legend
custom_labels <- c('P1' = "Pulse 1 - \nNewman Sound", 'Fjord' = "Fjord - \nSkagerrak", 'North Sea' = "North Sea - \nSkagerrak", 'P3' = "Pulse 3 - \nNewman Sound")

#add in predicted log2 FMR at common mass & temp
star_data <- data.frame(
  x = c(7.5, 7.5, 7.5, 7.5),
  y = c(6.72, 6.10, 5.94, 6.09),
  color = c('P1', 'FJ', 'NS', 'P3')  # Assign a color category to each star
)

TempvsFMR <- ggplot(aes(x = Temp, y = log2_pred_FMR), data = Combined_data) +
  geom_point(aes(x = Temp, y = log2_FMR, group = PP_ID, col = as.factor(PP_ID)))+
  geom_smooth(aes(x = Temp, y = log2_pred_FMR, col = as.factor(PP_ID)), method = "lm", se = TRUE)+
  geom_point(data = star_data, aes(x = x, y = y, color = color), shape = 18, size = 5) +  # Add colored stars
  geom_abline(slope=0.1, intercept=5.5, lty=2, color="gray37")+
  geom_vline(xintercept = 7.5, lty = 2, color = "gray37") +  # Add vertical line at x = 7.5
  scale_color_manual(values = colors, labels = custom_labels) +
  scale_fill_manual(values = colors, labels = custom_labels) +
  theme(plot.background = element_rect(fill = "transparent", color = NA)) +
  theme(axis.text.y = element_text(size = 14)) +
  theme(axis.text.x = element_text(size = 14)) +
  theme(axis.title.x = element_text(size = 18)) +
  theme(axis.title.y = element_text(size = 18)) +
  theme(panel.background = element_rect(fill = "white"))  +
  theme(panel.background = element_rect(colour = "black"))  +
  theme(strip.background = element_rect(fill = "white"))   +
  theme(legend.position = "top",
    panel.border = element_rect(colour = "black", fill = NA),
    axis.text = element_text(colour = "black", size = 9, angle = 0, hjust = 0.5),
    legend.background = element_blank(),
    legend.title = element_text(hjust = 0.5, size= 10),
    legend.text = element_text(hjust=0.5, size= 9),
    legend.box.background = element_rect(colour = "black")) +
  labs(x = "Temperature (°C)", y = " log2 O2 consumption (mgO2/kg/hr)",
       color = "Within and Across \nPopulation Classification", fill = "Within and Across \nPopulation Classification")

print(TempvsFMR)

# Plot with stars including pink border effect
TempvsFMR <- ggplot(aes(x = Temp, y = log2_pred_FMR), data = GLS_data) +
  geom_point(aes(x = Temp, y = log2_FMR, group = PPID, col = as.factor(PPID))) +
  geom_smooth(aes(x = Temp, y = log2_pred_FMR, col = as.factor(PPID)), method = "lm", se = TRUE) +
  
  # Add an outline layer with larger, empty shapes
  geom_point(data = star_data, aes(x = x, y = y), shape = 18, size = 5, color = "aquamarine", fill = NA) +
  
  # Add the filled stars on top
  geom_point(data = star_data, aes(x = x, y = y, color = color), shape = 18, size = 3) +
  
  geom_abline(slope = 0.1, intercept = 5.5, lty = 2, color = "gray37") +
  geom_vline(xintercept = 7.5, lty = 2, color = "gray37") +  # Add vertical line at x = 7.5
  
  scale_color_manual(values = colors, labels = custom_labels) +
  scale_fill_manual(values = colors, labels = custom_labels) +
  
  theme(plot.background = element_rect(fill = "transparent", color = NA)) +
  theme(axis.text.y = element_text(size = 14)) +
  theme(axis.text.x = element_text(size = 14)) +
  theme(axis.title.x = element_text(size = 18)) +
  theme(axis.title.y = element_text(size = 18)) +
  theme(panel.background = element_rect(fill = "white")) +
  theme(panel.background = element_rect(colour = "black")) +
  theme(strip.background = element_rect(fill = "white")) +
  theme(legend.position = "top",
        panel.border = element_rect(colour = "black", fill = NA),
        axis.text = element_text(colour = "black", size = 9, angle = 0, hjust = 0.5),
        legend.background = element_blank(),
        legend.title = element_text(hjust = 0.5, size = 10),
        legend.text = element_text(hjust = 0.5, size = 9),
        legend.box.background = element_rect(colour = "black")) +
  labs(x = "Temperature (°C)", y = "log2 O2 consumption (mgO2/kg/hr)",
       color = "Within and Across \nPopulation Classification", fill = "Within and Across \nPopulation Classification")

print(TempvsFMR)










##Q3_Mass scaled Predict FMR values to compare to measured FMR and plot thermal sensitivity 
# as log2 FMR vs absolute temp in *C 
# we use log2 FMR, which means that each unit increase is a doubling of metabolic rate 
# e.g. Q10 of 2 would mean a slope of 1/10

#make sure headings are still intact 
Combined_data$z_inv_temp <- Z.score(Combined_data$InvTemp, "inv_temp")$z_inv_temp
Combined_data$z_ln_mass <- Z.score(Combined_data$lnM,"ln_mass")$z_ln_mass
head(Combined_data)
#predict new lnFMR values using chosen model
GLS_MS_2<- gls(lnFMR_MS ~ z_inv_temp*PP_ID, weights=varIdent(form=~1|PP_ID), method="ML",  data = Combined_data)
FMR_MS_predicted<-data.frame(z_inv_temp=Combined_data$z_inv_temp, 
                          PP_ID=Combined_data$PP_ID)
predictedlnFMR_MS <- predict(GLS_MS_2, newdata = FMR_MS_predicted)
# Add the predicted values to the GLS_data dataframe
Combined_data$predictedlnFMR_MS <- predictedlnFMR_MS
# Print the first few rows of the updated dataframe to check the results
head(Combined_data)

#convert all values for the log2 FMR vs Temp plot
Combined_data<-Combined_data%>%
  mutate(pred_FMR_MS= exp(predictedlnFMR_MS))%>%
  mutate(log2_FMR = (log(FMR,2)))%>%
  mutate(log2_FMR_MS=log(FMR_MS,2))%>%
  mutate(log2_pred_FMR_MS = log(pred_FMR_MS, 2))
  

# plot predicted & measured, add se to geom_smooth (predictedlogFMR)

#choose cold colours for NL pulses and warm colours for Norway ecotypes
colors <- c('P1'="#1F77B4", 'Fjord'="#FF7F0E",'North Sea'= "#D62728",'P3'= "#9467BD")

# define the labels for the legend
custom_labels <- c('P1' = "Pulse 1 - \nNewman Sound", 'Fjord' = "Fjord - \nSkagerrak", 'North Sea' = "North Sea - \nSkagerrak", 'P3' = "Pulse 3 - \nNewman Sound")

#add in predicted log2 FMR at common mass & temp
star_data <- data.frame(
  x = c(7.5, 7.5, 7.5, 7.5),
  y = c(6.3385, 5.5938, 5.8489, 6.7533),
  color = c('P1', 'Fjord', 'North Sea', 'P3')  # Assign a color category to each star
)

TempvsFMR_MS <- ggplot(aes(x = Temp, y = log2_pred_FMR_MS), data = Combined_data) +
  geom_point(aes(x = Temp, y = log2_FMR_MS, group = PP_ID, col = as.factor(PP_ID)))+
  geom_smooth(aes(x = Temp, y = log2_pred_FMR_MS, col = as.factor(PP_ID)), method = "lm", se = TRUE)+
  geom_point(data = star_data, aes(x = x, y = y, color = color), shape = 18, size = 5) +  # Add colored stars
  geom_abline(slope=0.1, intercept=5.5, lty=2, color="gray37")+
  geom_vline(xintercept = 7.5, lty = 2, color = "gray37") +  # Add vertical line at x = 7.5
  scale_color_manual(values = colors, labels = custom_labels) +
  scale_fill_manual(values = colors, labels = custom_labels) +
  theme(plot.background = element_rect(fill = "transparent", color = NA)) +
  theme(axis.text.y = element_text(size = 14)) +
  theme(axis.text.x = element_text(size = 14)) +
  theme(axis.title.x = element_text(size = 18)) +
  theme(axis.title.y = element_text(size = 18)) +
  theme(panel.background = element_rect(fill = "white"))  +
  theme(panel.background = element_rect(colour = "black"))  +
  theme(strip.background = element_rect(fill = "white"))   +
  theme(legend.position = "top",
        panel.border = element_rect(colour = "black", fill = NA),
        axis.text = element_text(colour = "black", size = 9, angle = 0, hjust = 0.5),
        legend.background = element_blank(),
        legend.title = element_text(hjust = 0.5, size= 10),
        legend.text = element_text(hjust=0.5, size= 9),
        legend.box.background = element_rect(colour = "black")) +
  labs(x = "Temperature (°C)", y = " log2 O2 consumption (mgO2/kg/hr) at 10g mass",
       color = "Within and Across \nPopulation Classification", fill = "Within and Across \nPopulation Classification")

print(TempvsFMR_MS)

# Plot with stars including pink border effect
TempvsFMR_MS <- ggplot(aes(x = Temp, y = log2_pred_FMR_MS), data = Combined_data) +
  geom_point(aes(x = Temp, y = log2_FMR_MS, group = PP_ID, col = as.factor(PP_ID))) +
  geom_smooth(aes(x = Temp, y = log2_pred_FMR_MS, col = as.factor(PP_ID)), method = "lm", se = TRUE) +
  
  # Add an outline layer with larger, empty shapes
  geom_point(data = star_data, aes(x = x, y = y), shape = 18, size = 5, color = "aquamarine", fill = NA) +
  
  # Add the filled stars on top
  geom_point(data = star_data, aes(x = x, y = y, color = color), shape = 18, size = 3) +
  
  geom_abline(slope = 0.1, intercept = 5.5, lty = 2, color = "gray37") +
  geom_vline(xintercept = 7.5, lty = 2, color = "gray37") +  # Add vertical line at x = 7.5
  
  scale_color_manual(values = colors, labels = custom_labels) +
  scale_fill_manual(values = colors, labels = custom_labels) +
  
  theme(plot.background = element_rect(fill = "transparent", color = NA)) +
  theme(axis.text.y = element_text(size = 14)) +
  theme(axis.text.x = element_text(size = 14)) +
  theme(axis.title.x = element_text(size = 18)) +
  theme(axis.title.y = element_text(size = 18)) +
  theme(panel.background = element_rect(fill = "white")) +
  theme(panel.background = element_rect(colour = "black")) +
  theme(strip.background = element_rect(fill = "white")) +
  theme(legend.position = "top",
        panel.border = element_rect(colour = "black", fill = NA),
        axis.text = element_text(colour = "black", size = 9, angle = 0, hjust = 0.5),
        legend.background = element_blank(),
        legend.title = element_text(hjust = 0.5, size = 10),
        legend.text = element_text(hjust = 0.5, size = 9),
        legend.box.background = element_rect(colour = "black")) +
  labs(x = "Temperature (°C)", y = "log2 O2 consumption at 10g (mgO2/kg/hr)",
       color = "Within and Across \nPopulation Classification", fill = "Within and Across \nPopulation Classification")

print(TempvsFMR_MS)



##Q4: Are there significant differences in mean FMR across groups at common mass(10g)/temp(7.5)

target_inv_temp<-1/((K)*(7.5+273))

z_common_inv_temp<-(target_inv_temp-mean_Inv.Temp)/sd_Inv.Temp

#predict value for each PPID based on common mass and temp
unique_PP_ID <- unique(Combined_data$PP_ID)
FMR_common<-data.frame(z_inv_temp=z_common_inv_temp, PP_ID=unique_PP_ID)
TargetlnFMR<-predictSE.gls(GLS_MS_2, newdata=FMR_common, se.fit=T, print.matrix=FALSE)
print(TargetlnFMR)

#predicted lnFMR for each PPID is: 4.654700 (P1), 4.218790 (P3), 4.223290 (FJ), 4.117086 (NS)
#se.fit for each PPID is: 0.03455208 (P1), 0.11963300 (P3), 0.09335259 (FJ), 0.12642476 (NS)
#The se.fit values are the standard errors of the predicted values
# calculate the confidence intervals 

# Define a named vector or lookup table to map PPID numbers to names
ppid_names <- c("P1" = 1, "P3" = 2, "FJ" = 3, "NS" = 4)

# Create a reverse lookup from numbers to names
unique_PP_ID <- c(1, 2, 3, 4)  # Replace with actual unique PPID numbers
ppid_names_reverse <- names(ppid_names)[match(unique_PP_ID, ppid_names)]

# Calculate confidence intervals for each
# Extract the fitted values and standard errors
fitted_values <- TargetlnFMR$fit
standard_errors <- TargetlnFMR$se.fit

Ln2_fitted_common_values<-log(exp(fitted_values),2)

# Step 4: Calculate the 95% confidence intervals
lower_bound <- fitted_values - 1.96 * standard_errors
upper_bound <- fitted_values + 1.96 * standard_errors

print(lower_bound)
print(upper_bound)

# Combine the results into a dataframe with correct PPID names
confidence_intervals <- data.frame(
  PP_ID = ppid_names_reverse,
  Fitted_Values = fitted_values,
  Lower_CI = lower_bound,
  Upper_CI = upper_bound
)

# Print the confidence intervals
print(confidence_intervals)

#visually inspect if confidence intervals overlap
for (i in 1:(nrow(confidence_intervals) - 1)) {
  for (j in (i + 1):nrow(confidence_intervals)) {
    if (confidence_intervals$Upper_CI[i] < confidence_intervals$Lower_CI[j] ||
        confidence_intervals$Lower_CI[i] > confidence_intervals$Upper_CI[j]) {
      cat("PPID", confidence_intervals$PPID[i], "and PPID", confidence_intervals$PPID[j], "are significantly different.\n")
    } else {
      cat("PPID", confidence_intervals$PPID[i], "and PPID", confidence_intervals$PPID[j], "are not significantly different.\n")
    }
  }
}



## A4 We compare confidence intervals of different PPIDs to check sig diffs
## At common mass (5g) and temperature (7.5) sig differences in mean FMR differ - 
#same sig diff results as across the range of temperature & masses

#Within pops:
#PPID FJ and PPID NS are  significantly different.
#PPID P1 and PPID P3 are significantly different.

#Across pops:
#PPID P1 and PPID FJ are not significantly different.
#PPID P1 and PPID NS are significantly different.
#PPID P3 and PPID FJ are  significantly different.
#PPID P3 and PPID NS are not significantly different.

#predicted lnFMR for each PPID is: 4.654700 (P1), 4.218790 (P3), 4.223290 (FJ), 4.117086 (NS)
#se.fit for each PPID is: 0.03455208 (P1), 0.11963300 (P3), 0.09335259 (FJ), 0.12642476 (NS)

