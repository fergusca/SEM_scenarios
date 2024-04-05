##########################
## SCENARIOS AND O/E PREDICTIONS
## XERIC wadeable
## POPULATION WEIGHTS
##
## 10/4/2023
##########################
remove(list=ls())

library(tidyverse)
library(lavaan)
library(dplyr)
library(ggplot2)
library(ggpmisc)
library(ggpubr)
library(broom)
library(Metrics) # for rmse

library(grid)
library(gridExtra)

library(kableExtra)
library(knitr)
library(stringr)

library(spsurvey)

devtools::load_all()
library(SEMNRSA)

#########################
# LOAD DATA
##############
## READ PROCESSED NRSA DATA USED IN MODEL
# PROCESSED - all 0809 and only new sites from later surveys
#  VISITS 1 and 2 n = 4578 w/ 321 vars
dat_org<-read.csv("C:/Users/EFergus/OneDrive - Environmental Protection Agency (EPA)/a_NLA_OE_project/Project_repository/SEMNRSA/data_processed/Compiled//nrsa081318_nonresampled_VISIT_12.csv")
#dat_org <- read.csv("data_processed/Compiled/NRSA_compile_proc.csv")#C:/Users/EFergus/OneDrive - Environmental Protection Agency (EPA)/a_NLA_OE_project/Project_repository/SEMNRSA/data_processed/Compiled/nrsa081318_nonresampled_VISIT_12.csv")

# PROCESSED DATA VISIT_NO=1 ONLY n = 4389
dat_proc<-dat_org%>%
  filter(VISIT_NO==1)

###############
## PROCESS DATA DROPPING MISSING PROTOCOL
dat_proc$PROTOCOL<-as.factor(dat_proc$PROTOCOL)

# n = 4371
dat_proc<- dat_proc%>%
  drop_na(PROTOCOL)%>%
  filter(PROTOCOL=="BOATABLE"|PROTOCOL=="WADEABLE")

# DROP NOPHAB class from REALM
dat_proc$PROTOCOL<-droplevels(dat_proc$PROTOCOL)

#########################
## SUBSET BY MANUAL ECOREGION
# XER n = 342
xer<-dat_proc%>%
  filter(AG_ECO9=="XER")%>%
  drop_na(LOE_QLow_cl)

# SCALE CUMULATIVE PRECIPITATION in XER
xer$PSUMPY_SY_WS_sc<-scale(xer$PSUMPY_SY_WS)

# WADEABLE n = 200
xer_w <- xer %>%
  filter(PROTOCOL=="WADEABLE")

########
## TN CONDITION CLASS BASED ON NRSA 1314 thresholds for XER
xer_w_mod<-xer_w%>%
  filter(!is.na(WGT_TP))%>%
  mutate(TN_cond = case_when(
    L_NTL < -0.55 ~ "Good",
    L_NTL > -0.55 & L_NTL < -0.29 ~ "Fair",
    L_NTL >= -0.29 ~ "Poor"))
table(xer_w_mod$TN_cond)


#############
# Convert XER_w dataframe to a simple features (sf) object
# https://stackoverflow.com/questions/29736577/how-to-convert-data-frame-to-spatial-coordinates/45484314#45484314
xer_w_sf <- st_as_sf(x=xer_w_mod,
                           coords=c("LON_DD83","LAT_DD83"),
                           crs=4269) # Albers equal area EPSG = 5070
#################
## ESTIMATE POPULATION OF XER STREAMS BY TN CONDITION CLASS
cat_ests <- cat_analysis(
  xer_w_sf,
  siteID = "SITE_ID",
  vars = "TN_cond",
  weight = "WGT_TP"
)
cat_ests
cat_ests_red<-cat_ests%>%
  filter(Category %in% c("Fair","Good","Poor"))

cat_ests_red$Category<-ordered(cat_ests_red$Category, levels=c("Good","Fair","Poor"))

col_palette <-c("#91bfdb","#ffffbf","#fc8d59")
#Set font
windowsFonts(AR=windowsFont("Arial"))

XER_TN_class<-ggplot(cat_ests_red, aes(x=Category, y=Estimate.P,fill=Category))+
  geom_bar(stat="identity",position = position_dodge())+
  geom_errorbar(aes(ymax=Estimate.P + MarginofError.P,ymin=Estimate.P -MarginofError.P ),
                width=0.2)+
  scale_fill_manual(values=col_palette)+
  geom_text(aes(label=paste0(round(Estimate.P,0),"%")),vjust=-0.2)+
  #theme_bw(base_size=12)+
  theme(plot.title = element_text(family = "AR",face="plain",size=14, hjust=0.5),
        axis.text.x = element_text(family = "AR", angle=45, hjust=1,size=12),#, colour=c(rep("#8c510a",3),rep("#bf812d",1),rep("#dfc27d",3),rep("#c7eae5",2),rep("#80cdc1",1),rep("#c7eae5",1),rep("#35978f",2),rep("#01665e",1))),
        axis.text.y = element_text(family = "AR", size=12),
        axis.title.y = element_text(family="AR"), #element_blank(),#
        strip.text.x = element_text(family="AR", size=12),
        panel.grid.major =  element_line(colour = NA),
        panel.grid.minor=element_line(colour = NA))+
  xlab("")+
  ylab("% of Stream length")+
  ggtitle("Xeric TN class") +
  coord_flip()
XER_TN_class

tiff(filename="C:/Users/efergus/SEMNRSA/inst/analysis_prediction/Routput/Figures/XER_TN_Class_Pop.tiff",
     width=5, height=5, units="in", res=200)
XER_TN_class
dev.off()


#############
## READ SEM OUTPUT PROCESSED W/COEFFICIENT LABELS UNSTANDARDIZED NON-BAYESIAN
coef_proc<-read.csv("inst/analysis_prediction/Routput/XERw_m3_OE_coef_unstd.csv")

# Grab label and coefficient estimate to call in equations
coef_labels<-coef_proc[,3:4]
# Widen table - make coefficient labels columns and estimates the row
coef_use<-pivot_wider(coef_labels,names_from=coeff_name,values_from = est)#


###########################
# CREATE HYPOTHETICAL TN
# Decrease TN on subset of sites where TN is elevated above regional reference 75th percentile values (i.e., TN concentrations where 75% of reference sites are below that value L_NTL= -0.55)
## USE NRSA 1314 Technical report threshold for Good/Fair TN class = 285 ug/L or -0.55

#  Add L_NTL scenario values
xer_mod<-xer_w%>%
  #filter(L_NTL> -0.55)%>%
  # Decrease to 75th percentile
  mutate(LNTL_s2=ifelse(L_NTL> -0.55,-0.55,
                        L_NTL))%>%
  mutate(NTL_s2 = log10(LNTL_s2)-0.01)


#########################
# Set scenario values and predict bug O/E
#Bug O/E ~ 1.18 + 0.06(LRBS_use) + -0.14(L_NTL) + -0.08(L_SULF) + 0.05(LOE_Qbkf_cl) + 0.04(dam) + 0.12(Specific strm pwr)

# SCENARIO 1 - OBSERVED TN
L_NTL.predict.1 <- xer_mod$L_NTL
tn.predict.1.raw<-xer_mod$NTL_RESULT

#SCENARIO 2 - TN decrease to 75th
L_NTL.predict.2 <- xer_mod$LNTL_s2
tn.predict.2.raw<-10^L_NTL.predict.2 -0.01

# PREDICTED DISTRIBUTIONS FOR EXOGENOUS DRIVERS
# VALUES IN DATASET
urb.pred.1 <- xer_mod$asin_PCTURB_CAT
dam.pred.1 <- xer_mod$L_NABD_NrmStorWs_ratio
hag.pred.1 <- xer_mod$W1_HAG
nohag.pred.1 <- xer_mod$W1_HNOAG
rnat.pred.1 <- xer_mod$asin_PCTFORGRS_CATRP100
precip.pred.1 <-xer_mod$PSUMPY_SY_WS_sc
phdi.pred.1<-xer_mod$drought_mean

# ENDOGENOUS DRIVERS BUT NOT ALTERED BY LFLOW
strmpwr.pred.1 <- xer_mod$L_STRM_POWER
xcmgw.pred.1 <- xer_mod$Lpt01_XCMGW
loebflow.pred.1 <- xer_mod$LOE_Qbkf_cl
loelflow.pred.1 <- xer_mod$LOE_QLow_cl
rbs.pred.1 <- xer_mod$LRBS_use
sulf.pred.1 <- xer_mod$L_SULF


############
# PREDICTED DISTRIBUTIONS FOR BENTHIC O/E
# SCENARIO #1
oe.predict.1 <- coef_use$b9.1 + coef_use$b9.2*rbs.pred.1 +
  coef_use$b9.3*L_NTL.predict.1 + coef_use$b9.4*sulf.pred.1 +
  coef_use$b9.5*loebflow.pred.1  + coef_use$b9.6*dam.pred.1 +
  coef_use$b9.7*strmpwr.pred.1

# SCENARIO #2
oe.predict.2 <- coef_use$b9.1 + coef_use$b9.2*rbs.pred.1 +
  coef_use$b9.3*L_NTL.predict.2 + coef_use$b9.4*sulf.pred.1 +
  coef_use$b9.5*loebflow.pred.1  + coef_use$b9.6*dam.pred.1 +
  coef_use$b9.7*strmpwr.pred.1

#######################
## CREATE DATAFRAMES OF PREDICTED VALUES
scenario.1<-data.frame(L_NTL.predict.1,tn.predict.1.raw,
                       oe.predict.1) #xfcnat.predict.1.raw,

scenario.2<-data.frame(L_NTL.predict.2,tn.predict.2.raw,
                       oe.predict.2)

# ADD COLUMNS TO INDIVIDUAL SCENARIO DATASETS TO INDICATE WHICH SCENARIO THEY BELONG TO
scenario.1<-scenario.1%>%
  mutate(scenario=1)
scenario.2<-scenario.2%>%
  mutate(scenario=2)

# change column names
cols_names<-c("LNTL.predict",
              "tn.predict.raw","oe.predict","scenario")
names(scenario.1)<- cols_names
names(scenario.2)<- cols_names


# SCENARIOS GRP 1 - Decrease to reference percentile
scenario_all<-rbind(scenario.1,scenario.2) #,scenario.4
scenario_all$scenario<-as.factor(scenario_all$scenario)

########################
## COMBINE PREDICTED WITH OBSERVED
################
## Combine scenario predictions with NRSA dataset
xer_mod_red<- xer_mod%>%
  select(SITE_ID,VISIT_NO,YEAR,UNIQUE_ID,LAT_DD83,LON_DD83,RT_MASTER,
         OE_SCORE,LOE_QLow_cl,LRBS_use,L_NTL,L_SULF,L_STRM_POWER,WGT_TP)
xer_pred1<-cbind(xer_mod_red,scenario.1)
xer_pred2<-cbind(xer_mod_red,scenario.2)

xer_pred1<- xer_pred1%>%
  mutate(OE_SCORE_diff=oe.predict-OE_SCORE)%>%
  filter(!is.na(OE_SCORE))

xer_pred2<- xer_pred2%>%
  mutate(OE_SCORE_diff=oe.predict-OE_SCORE)%>%
  filter(!is.na(OE_SCORE))

xer_prediction<-rbind(xer_pred1,xer_pred2)
xer_prediction$scenario<- as.factor(xer_prediction$scenario)

# Remove NA OE_SCORE
xer_prediction<-xer_prediction%>%
  filter(!is.na(OE_SCORE_diff))

# Plot differences in observed O/E and predicted across scenarios
OEdiff_box <- ggplot(xer_prediction, aes(scenario,OE_SCORE_diff))+ #OE_bin
  geom_boxplot()+
  geom_jitter(alpha=0.2)+
  #facet_wrap(~RT_MASTER)+
  theme(axis.text.x = element_text(angle=65, hjust=1,size=10))+
  xlab("Scenarios")+
  ggtitle("Difference in predicted from observed macroinvertebrate O/E")
OEdiff_box

# RMSE
rmse(xer_pred1$oe.predict,xer_pred1$OE_SCORE)
# 0.2420565
rmse(xer_pred2$oe.predict,xer_pred2$OE_SCORE)
# 0.2506521
###########
# Predicted O/E and observed has a high degree of error as seen with the rmse in scenario 1 where TN is at the default setting
#  Predicted O/E could be off by 0.24 units which is high when the mean value is 0.79

#####################
## ESTIMATING POPULATION VALUE O/E using spsurvey package
#https://cran.r-project.org/web/packages/spsurvey/vignettes/analysis.html#5_Variance_Estimation

# NOTE THERE ARE 56 sites in XER that are missing WGT_TP values - not sure why...
# NEED TO DROP them to estimate values for population
# SELECTING SITES WITH TN THAT HAVE BEEN ALTERED (sites below FAIR/GOOD) >-0.55 L_NTL

## Combine scenario predictions with NRSA dataset reduced to sites
xer_prediction<-xer_prediction%>%
  mutate(TN_class = case_when(
    L_NTL > -0.55 ~ "hi TN",
    L_NTL <= -0.55 ~ "other"))

table(xer_prediction$scenario, xer_prediction$TN_class)


########
## TN CONDITION CLASS BASED ON NRSA 1314 thresholds for XER
#xer_prediction<-xer_prediction%>%
##  mutate(TN_cond = case_when(
##    L_NTL < -0.55 ~ "good",
#    L_NTL > -0.55 & L_NTL < -0.29 ~ "fair",
#    L_NTL >= -0.29 ~ "poor"))
#table(xer_prediction$TN_cond)

# SCENARIO 1 n=135 - sites with TN greater than GOod/Fair condition
xer_prediction_s1<-xer_prediction%>%
  filter(scenario=="1")%>%
  filter(!is.na(WGT_TP))%>%
  filter(L_NTL>-0.55)

# Convert dataframe to a simple features (sf) object
# https://stackoverflow.com/questions/29736577/how-to-convert-data-frame-to-spatial-coordinates/45484314#45484314
xer_pred_s1_sf <- st_as_sf(x=xer_prediction_s1,
                           coords=c("LON_DD83","LAT_DD83"),
                           crs=4269) # Albers equal area EPSG = 5070

cont_est_s1 <-cont_analysis(
  xer_pred_s1_sf,
  siteID = "SITE_ID",
  vars="oe.predict",
  weight="WGT_TP"
)

# Mean estimates of Scenario 1 predicted O/E
cont_est_s1$Mean
#       Type Subpopulation  Indicator nResp  Estimate   StdError MarginofError  LCB95Pct  UCB95Pct
#1 All_Sites     All Sites oe.predict   135 0.7188487 0.01415532    0.02774391 0.6911048 0.7465927

# CDF estimates
plot(cont_est_s1$CDF)


# SCENARIO 2 - sites NOT in Good/Fair condition
xer_prediction_s2<-xer_prediction%>%
  filter(scenario=="2")%>%
  filter(!is.na(WGT_TP))%>%
  filter(L_NTL>-0.55)
  #filter(L_NTL>-0.28)

# Convert dataframe to a simple features (sf) object
# https://stackoverflow.com/questions/29736577/how-to-convert-data-frame-to-spatial-coordinates/45484314#45484314
xer_pred_s2_sf <- st_as_sf(x=xer_prediction_s2,
                           coords=c("LON_DD83","LAT_DD83"),
                           crs=4269) # Albers equal area EPSG = 5070

cont_est_s2 <-cont_analysis(
  xer_pred_s2_sf,
  siteID = "SITE_ID",
  vars="oe.predict",
  weight="WGT_TP"
)

# Mean estimates of Scenario 2 predicted O/E
cont_est_s2$Mean
#        Type Subpopulation  Indicator nResp  Estimate   StdError MarginofError  LCB95Pct  UCB95Pct
# 1 All_Sites     All Sites oe.predict   135 0.7674417 0.01201239    0.02354385 0.7438978 0.7909855

# CDF estimates
plot(cont_est_s2$CDF)

#################
## COMBINE CDF dataframes
s1_cdf<-cont_est_s1$CDF%>%
  mutate(scenario="scenario_1")

s2_cdf<-cont_est_s2$CDF%>%
  mutate(scenario="scenario_2")

cdf_all_highTN <- rbind(s1_cdf,s2_cdf)

# CDF for WMT & XER
oe_cdf<-ggplot(cdf_all_highTN, aes(x=Value))+
  stat_ecdf(aes(color=scenario, linetype=scenario),
            geom="step", size=1.5)+
  theme(legend.position="bottom")+
  xlab("O/E predicted")+
  ggtitle("XER population O/E \nscenario reduce TN to Good/Fair threshold")

oe_cdf

# PRINT CDF
tiff(filename="C:/Users/efergus/SEMNRSA/inst/analysis_prediction/Routput/Figures/Population_OE_scenario_v3_2023_1005.tiff",
     width=5, height=4, units="in", res=200)
oe_cdf
dev.off()
