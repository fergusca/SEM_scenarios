#####################
## REVISED PATH ANALYSIS MODEL FOR PREDICTION
##  XERIC wadeable
##  Model v3
##  Removed agriculture Ws because very low numbers and getting spurious relationships
##  Removed wetlands,
##  Impervious surface in Ws 100 m buffer
##   AND catchment scale riparian cover (only forest + grassland - dropped shrub)
##
##  9/12/2023
######################

remove(list=ls())

library(tidyverse)
library(lavaan)
library(dplyr)
library(ggplot2)
library(ggpmisc)

#######################
# Need to treat SEMNRSA like a package to be able to run function relabeling predictors
devtools::load_all()
library(SEMNRSA)

###########################
## READ PROCESSED DATA
# COMPILED NRSA SURVEYS WITH SUBSET OF VARIABLES
#  INCLUDES ALL RESAMPLED SITES AND VISITS 1 & 2

# PROCESSED - all 0809 and only new sites from later surveys
#  VISITS 1 and 2 n = 4578 w/ 316 vars
dat_org <- read.csv("data_processed/Compiled/nrsa081318_nonresampled_VISIT_12.csv")
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

# TRANSFORM IMPERVIOUS SURFACE (NLCD
dat_proc<-dat_proc%>%
  mutate(asin_PCTIMP_WS = asin(sqrt(PCTIMP_WS/100)),
         asin_PCTIMP_WsRp100 = asin(sqrt(PCTIMP_WsRp100/100)),
         asin_PCTIMP_CAT = asin(sqrt(PCTIMP_CAT/100)),
         asin_PCTIMP_CATRP100 = asin(sqrt(PCTIMP_CATRP100/100)))

summary(dat_proc$asin_PCTIMP_WsRp100)

# ADD SOME VARIABLES TO TEST
#dat_proc<- dat_proc%>%
#  mutate(PCTURB_CAT_mod = floor(PCTURB_CAT),
#         asin_PCTURB_CAT = asin(sqrt(PCTURB_CAT_mod/100)),
#         asin_PCTAGR_CAT = asin(sqrt(PCTAGR_CAT/100)),
#         PCTURB_CATRP100_mod = floor(PCTURB_CATRP100),
#         asin_PCTURB_CATRP100 = asin(sqrt(PCTURB_CATRP100_mod/100)),
#         asin_PCTAGR_CAT = asin(sqrt(PCTAGR_CATRP100/100))
#         )

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

###########################
# SCATTERPLOT NLCD Developed vs Impervious
urb_imp<-ggplot(xer_w, aes(x=PCTURB_WsRp100, y=PCTIMP_WsRp100))+ #PCTFOR_WsRp100
  stat_poly_line()+
  stat_poly_eq(use_label(c("eq","R2")))+
  geom_point()#+
#facet_wrap(~RT_MASTER)
urb_imp
# Highly correlated but smaller values for impervious surface %

# SCATTERPLOT Watershed area ~~ Catchment area in XER
WS_CAT<-ggplot(xer_w, aes(x=log10(WsAreaSqKm), y=log10(CatAreaSqKm)))+ #PCTFOR_WsRp100
  stat_poly_line()+
  stat_poly_eq(use_label(c("eq","R2")))+
  geom_point()
WS_CAT


###########################################
## PATH MODEL v3 O/E
mymodel_v3_OE <- '
Lpt01_XCMGW ~ asin_PCTIMP_WsRp100 + W1_HAG + W1_HNOAG + L_NABD_NrmStorWs_ratio + asin_PCTFORGRS_CATRP100 + LOE_QLow_cl+ LOE_Qbkf_cl + L_STRM_POWER + PSUMPY_SY_WS_sc + drought_mean
L_STRM_POWER ~PSUMPY_SY_WS_sc + drought_mean + L_NABD_NrmStorWs_ratio + W1_HAG + asin_PCTFORGRS_CATRP100
LOE_QLow_cl~ asin_PCTIMP_WsRp100 + PSUMPY_SY_WS_sc + L_NABD_NrmStorWs_ratio + asin_PCTFORGRS_CATRP100 + W1_HAG + W1_HNOAG + drought_mean
LOE_Qbkf_cl ~ asin_PCTIMP_WsRp100 + PSUMPY_SY_WS_sc + L_NABD_NrmStorWs_ratio + asin_PCTFORGRS_CATRP100 + W1_HAG + W1_HNOAG + drought_mean
evap_index_sc ~ LOE_QLow_cl + LOE_Qbkf_cl + PSUMPY_SY_WS_sc + L_NABD_NrmStorWs_ratio + asin_PCTFORGRS_CATRP100 + drought_mean + L_STRM_POWER
LRBS_use ~ asin_PCTIMP_WsRp100 + L_NABD_NrmStorWs_ratio + Lpt01_XCMGW + LOE_QLow_cl + LOE_Qbkf_cl + W1_HAG + W1_HNOAG + asin_PCTFORGRS_CATRP100 + PSUMPY_SY_WS_sc + L_STRM_POWER + drought_mean
L_PTL ~ asin_PCTIMP_WsRp100 + Lpt01_XCMGW + evap_index_sc + W1_HAG + W1_HNOAG + asin_PCTFORGRS_CATRP100 + LOE_QLow_cl + LOE_Qbkf_cl + PSUMPY_SY_WS_sc + drought_mean
L_NTL ~ asin_PCTIMP_WsRp100 + Lpt01_XCMGW + evap_index_sc + W1_HAG + W1_HNOAG + asin_PCTFORGRS_CATRP100 + LOE_QLow_cl + LOE_Qbkf_cl + PSUMPY_SY_WS_sc + drought_mean
L_CHLR ~ asin_PCTIMP_WsRp100 + Lpt01_XCMGW + evap_index_sc + W1_HAG + W1_HNOAG + asin_PCTFORGRS_CATRP100 + LOE_QLow_cl + LOE_Qbkf_cl + PSUMPY_SY_WS_sc + drought_mean
L_SULF ~  asin_PCTIMP_WsRp100 + Lpt01_XCMGW + evap_index_sc + W1_HAG + W1_HNOAG + asin_PCTFORGRS_CATRP100 + LOE_QLow_cl + LOE_Qbkf_cl + PSUMPY_SY_WS_sc + drought_mean
L_TURB ~ asin_PCTIMP_WsRp100 + Lpt01_XCMGW + evap_index_sc + W1_HAG + W1_HNOAG + asin_PCTFORGRS_CATRP100 + LOE_QLow_cl + LOE_Qbkf_cl + PSUMPY_SY_WS_sc + drought_mean

OE_SCORE ~ Lpt01_XCMGW + LRBS_use + L_PTL + L_NTL + L_CHLR + L_SULF + L_TURB + LOE_QLow_cl + LOE_Qbkf_cl + asin_PCTIMP_WsRp100 + L_NABD_NrmStorWs_ratio + asin_PCTFORGRS_CATRP100 + W1_HAG + W1_HNOAG + L_STRM_POWER + drought_mean

# Covariance
LOE_QLow_cl~~LOE_Qbkf_cl
L_STRM_POWER~~LOE_QLow_cl
L_STRM_POWER~~LOE_Qbkf_cl
L_PTL ~~ L_NTL
L_CHLR ~~ L_SULF
L_PTL ~~ L_TURB

'

##############
## XER O/E
# ROBUST ESTIMATION MAX LIKELIHOOD METHOD
fit_v3_XERw_OE_robust.est<- sem(mymodel_v3_OE, data=xer_w,
                                estimator="MLM")#estimator="MLM")

summary(fit_v3_XERw_OE_robust.est, standardized=TRUE, fit.measures=TRUE, modindices=F,rsquare=TRUE)#, modindices=T, rsquare=TRUE) #

# request modification indices greater than 3.0 - from Grace USGS materials
mi_min <-modindices(fit_v3_XERw_OE_robust.est)
print(mi_min[mi_min$mi >3.0,])

##########################
#####################
## REVISED - Added total effects
mymodel_v3_OE_rev <- '
Lpt01_XCMGW ~ 1 + ha3*W1_HAG + na2*W1_HNOAG + d2*L_NABD_NrmStorWs_ratio + sp3*L_STRM_POWER + p3*PSUMPY_SY_WS_sc
L_STRM_POWER ~ 1 + d3*L_NABD_NrmStorWs_ratio + ha6*W1_HAG + rn4*asin_PCTFORGRS_CATRP100
LOE_QLow_cl~ 1 + i5*asin_PCTIMP_WsRp100 + p8*PSUMPY_SY_WS_sc + d5*L_NABD_NrmStorWs_ratio + rn5*asin_PCTFORGRS_CATRP100 + ha5*W1_HAG + na4*W1_HNOAG
LOE_Qbkf_cl ~ 1 + p7*PSUMPY_SY_WS_sc + d4*L_NABD_NrmStorWs_ratio + ha4*W1_HAG + ph3*drought_mean
evap_index_sc ~ 1 + l5*LOE_QLow_cl + b3*LOE_Qbkf_cl + p6*PSUMPY_SY_WS_sc + d6*L_NABD_NrmStorWs_ratio + sp4*L_STRM_POWER
LRBS_use ~ 1 + i2*asin_PCTIMP_WsRp100 + l2*LOE_QLow_cl + b2*LOE_Qbkf_cl + ha2*W1_HAG + rn2*asin_PCTFORGRS_CATRP100 + p2*PSUMPY_SY_WS_sc + sp2*L_STRM_POWER
L_NTL ~ 1 + i3*asin_PCTIMP_WsRp100 + x2*Lpt01_XCMGW + l3*LOE_QLow_cl + p4*PSUMPY_SY_WS_sc
L_SULF ~  1 + i4*asin_PCTIMP_WsRp100 + e2*evap_index_sc + na3*W1_HNOAG + rn3*asin_PCTFORGRS_CATRP100 + l4*LOE_QLow_cl + p5*PSUMPY_SY_WS_sc + ph2*drought_mean

OE_SCORE ~ 1 + r1*LRBS_use + n1*L_NTL + su1*L_SULF + boe1*LOE_Qbkf_cl + d1*L_NABD_NrmStorWs_ratio + sp1*L_STRM_POWER

# INDIRECT EFFECTS ON OE
imp_oe:= i2*r1 + i3*n1 + i4*su1 + i5*l2*r1 + i5*l3*n1 + i5*l4*su1 + i5*l5*e2*su1
dam_oe:= d2*x2*n1 + d3*sp1 + d3*sp2*r1 +d3*sp3*x2*n1 + d3*sp4*e2*su1 + d4*boe1 + d4*b2*r1 + d4*b3*e2*su1 +d5*l2*r1 + d5*l3*n1 + d5*l4*su1 + d5*l5*e2*su1 + d6*e2*su1
precip_oe:= p2*r1 + p3*x2*n1 + p4*n1 + p5*su1 + p6*e2*su1 + p7*boe1 + p7*b2*r1 + p7*b3*e2*su1 + p8*l2*r1 + p8*l3*n1 + p8*l4*su1 + p8*l5*e2*su1
phdi_oe:= ph2*su1 + ph3*boe1 + ph3*b2*r1 + ph3*b3*e2*su1
strpwr_oe:= sp2*r1 + sp3*x2*n1 + sp4*e2*su1
hag_oe:= ha2*r1 + ha3*x2*n1 + ha4*boe1 + ha4*b2*r1 + ha4*b3*e2*su1 + ha5*l2*r1 + ha5*l3*n1 + ha5*l4*su1 + ha5*l5*e2*su1 + ha6*sp1 + ha6*sp2*r1 + ha6*sp3*x2*n1 + ha6*sp4*e2*su1
nhag_oe:= na2*x2*n1 + na3*su1 + na4*l2*r1 + na4*l3*n1 + na4*l4*su1 + na4*l5*e2*su1
rnat_oe:= rn2*r1 + rn3*su1 + rn4*sp1 + rn4*sp2*r1 + rn4*sp3*x2*n1 + rn4*sp4*e2*su1
oelflow_oe:= l2*r1 + l3*n1 + l4*su1 + l5*e2*su1
oebflow_oe:= b2*r1 + b3*e2*su1
dexcess_oe:= e2*su1
xcmgw_oe:=x2*n1

# TOTAL EFFECTS ON OE
imp_tot:= imp_oe
dam_tot:= d1 + dam_oe
precip_tot:= precip_oe
phdi_tot:= phdi_oe
strpwr_tot:= sp1 + strpwr_oe
hag_tot:= hag_oe
nhag_tot:= nhag_oe
rnat_tot:= rnat_oe
oelflow_tot:= oelflow_oe
oebflow_tot:= boe1 + oebflow_oe
dexcess_tot:= dexcess_oe
xcmgw_tot:= xcmgw_oe
rbs_tot:= r1
tn_tot:= n1
sulf_tot:= su1

# TOTAL EFFECTS RBS
imp_rbstot:= i2+i5*l2
dam_rbstot:= d3*sp2 + d5*l2 + d4*b2
precip_rbstot:= p2+ p8*l2 + p7*b2
phdi_rbstot:= ph3*b2
strpwr_rbstot:= sp2
oelflow_rbstot:= l2
oebflow_rbstot:= b2
hag_rbstot:= ha2 + ha4*b2 + ha5*l2 + ha6*sp2
nhag_rbstot:= na4*l2

# TOTAL EFFECTS LOELflow
imp_lftot:= i5
dam_lftot:= d5
precip_lftot:= p8
hag_lftot:= ha5
nhag_lftot:= na4
rnat_lftot:= rn5

# TOTAL EFFECTS LOEBflow
dam_bftot:= d4
precip_bftot:= p7
phdi_bftot:= ph3
hag_bftot:= ha4

# TOTAL EFFECTS SP POWER
dam_sptot:= d3
rnat_sptot:= rn4
hag_sptot:= ha6

# TOTAL EFFECTS XCMGW
#imp_xcmgwtot:=
dam_xcmgwtot:= d2 + d3*sp3
precip_xcmgwtot:= p3
strpwr_xcmgwtot:= sp3
rnat_xcmgwtot:= rn4*sp3
hag_xcmgwtot:= ha3 + ha6*sp3
nhag_xcmgwtot:= na2

# TOTAL EFFECTS dexcess
imp_evtot:= i5*l5
dam_evtot:= d6 + d3*sp4 + d4*b3 + d5*l5
precip_evtot:= p6 + p8*l5 + p7*b3
phdi_evtot:= ph3*b3
strpwr_evtot:= sp4
oelflow_evtot:= l5
oebflow_evtot:= b3
rnat_evtot:= rn5*l5
hag_evtot:= ha4*b3 + ha5*l5 + ha6*sp4
nhag_evtot:= na4*l5

# TOTAL EFFECTS TN
imp_tntot:= i3 + i5*l3
dam_tntot:= d2*x2 + d5*l3 + d3*sp3*x2
precip_tntot:= p4 + p3*x2 + p8*l3
strpwr_tntot:= sp3*x2
oelflow_tntot:= l3
rnat_tntot:=rn4* sp3*x2 + rn5*l3
hag_tntot:= ha3*x2 + ha6*sp3*x2 + ha5*l3
nhag_tntot:= na2*x2 + na4*l3
xcmgw_tntot:= x2

# TOTAL EFFECTS SULFATE
imp_sulftot:= i4 + i5*l4 + i5*l5*e2
dam_sulftot:= d6*e2 + d4*b3*e2 + d5*l4 + d5*l5*e2 + d3*sp4*e2
precip_sulftot:= p5 + p8*l4 + p8*l5*e2 + p6*e2
phdi_sulftot:= ph2 + ph3*b3*e2
strpwr_sulftot:= sp4*e2
oelflow_sulftot:= l4 + l5*e2
oebflow_sulftot:= b3*e2
dexcess_sulftot:= e2
rnat_sulftot:= rn3 + rn4*sp4*e2 + rn5*l4 + rn5*l5*e2
hag_sulftot:= ha6*sp4*e2 + ha5*l4 + ha5*l5*e2 + ha4*b3*e2
nhag_sulftot:= na3 + na4*l4 + na4*l5*e2
#xcmgw_sulftot:=


# Covariance
LOE_QLow_cl~~LOE_Qbkf_cl
L_STRM_POWER~~LOE_QLow_cl
L_SULF~~L_NTL
#L_STRM_POWER~~LOE_Qbkf_cl

'

# ROBUST ESTIMATION MAX LIKELIHOOD METHOD
fit_v3_XERw_OE_rev_robust.est<- sem(mymodel_v3_OE_rev, data=xer_w,
                                    estimator="MLM")

summary(fit_v3_XERw_OE_rev_robust.est, standardized=TRUE, fit.measures=TRUE, modindices=F,rsquare=TRUE)#, modindices=T, rsquare=TRUE) #

# request modification indices greater than 3.0 - from Grace USGS materials
mi_min <-modindices(fit_v3_XERw_OE_rev_robust.est)
print(mi_min[mi_min$mi >3.0,])


###########
# SEM: PREDICT Y-VALUES

# PREDICT Using lavPredictY - specify variables to include
xvars<-c("PSUMPY_SY_WS_sc","drought_mean","LRBS_use","L_NTL","L_SULF","LOE_Qbkf_cl","LOE_QLow_cl","W1_HAG","W1_HNOAG","asin_PCTIMP_WsRp100","L_NABD_NrmStorWs_ratio","L_STRM_POWER","evap_index_sc","asin_PCTFORGRS_CATRP100","Lpt01_XCMGW")

pred.sem<-lavPredictY(fit_v3_XERw_OE_rev_robust.est,ynames="OE_SCORE",xnames=xvars,newdata = xer_w)
pred.sem


######################
## Bollen.stine bootstrap to estimate parameters -OE
fit_v3_XERw_OE_bootstrap_rev  <- sem(mymodel_v3_OE_rev, data=xer_w,
                                     #group = "ECOREG_rev",
                                     #missing="ML",
                                     test="bollen.stine", se="boot",bootstrap=1000)

summary(fit_v3_XERw_OE_bootstrap_rev, standardized=TRUE, fit.measures=TRUE, modindices=F,rsquare=TRUE)#, modindices=T, rsquare=TRUE) #

# AIC = 3997
#############
# Export R output - XERIC OE MODEL
#https://www.r-bloggers.com/export-r-output-to-a-file/
out_fit_v3_XERw_OE_rev<- capture.output(summary(fit_v3_XERw_OE_bootstrap_rev, standardized=FALSE, fit.measures=TRUE,rsquare=TRUE)) #, modindices=T
write.csv(out_fit_v3_XERw_OE_rev, "C:/Users/EFergus/OneDrive - Environmental Protection Agency (EPA)/a_NLA_OE_project/Project_repository/Routput/Scenario_modeling/SEM_output/XERw_m3_OE_rev.csv" , #"inst/analysis_prediction/Routput/XERw_m3_OE_rev.csv",
          row.names=FALSE)

# Standardized estimates of bootstrap model
std_parameter_se_bootstrap_min<- standardizedSolution(fit_v3_XERw_OE_bootstrap_rev)
write.csv(std_parameter_se_bootstrap_min, "C:/Users/EFergus/OneDrive - Environmental Protection Agency (EPA)/a_NLA_OE_project/Project_repository/Routput/Scenario_modeling/SEM_output/XERw_m3_OE_CI.csv" , #"inst/analysis_prediction/Routput/XERw_m3_OE_CI.csv",
          row.names = FALSE)

# Save model fit parameters
## GRAB R2
r2_XER_OE<-as.data.frame(lavInspect(fit_v3_XERw_OE_bootstrap_rev, "rsquare"))

r2_XER_OE<-cbind(variable=rownames(r2_XER_OE),r2_XER_OE)
rownames(r2_XER_OE) <- 1:nrow(r2_XER_OE) # To remove index column

colnames(r2_XER_OE) <- c("Variable","R2")
write.csv(r2_XER_OE, "C:/Users/EFergus/OneDrive - Environmental Protection Agency (EPA)/a_NLA_OE_project/Project_repository/Routput/Scenario_modeling/SEM_output/XERw_m3_OE_R2.csv", #"inst/analysis_prediction/Routput/XERw_m3_OE_R2.csv",
          row.names=FALSE)

#################
# GRAB R2 to be able to add to fit table
r2_red_OE <- r2_XER_OE%>%
  filter(Variable=="OE_SCORE")%>%
  select(R2)%>%
  mutate(across(where(is.numeric), round,2))

# TABLE OF FIT INDICES comparing 3 RESPONSES
table_fit_oe_xer <- matrix(NA, nrow=2, ncol=11)
colnames(table_fit_oe_xer) = c("Model","Estimation","X2", "df","CFI","TLI", "RMSEA","SRMR","AIC","n","npar")

table_fit_oe_xer[1,]<-c("OE","Standard",round(fitmeasures(fit_v3_XERw_OE_rev_robust.est,
                                                          c("chisq","df","cfi","tli",
                                                            "rmsea","srmr","aic","ntotal","npar")),2))
table_fit_oe_xer[2,]<-c("OE","Robust",round(fitmeasures(fit_v3_XERw_OE_rev_robust.est,
                                                        c("chisq.scaled","df.scaled","cfi.scaled","tli.scaled",
                                                          "rmsea.robust","srmr","aic","ntotal","npar")),2))
# Add R2 for biotic response
table_fit_oe_xer<- data.frame(table_fit_oe_xer)%>%
  mutate(R2 = r2_red_OE$R2)

write.csv(table_fit_oe_xer, "C:/Users/EFergus/OneDrive - Environmental Protection Agency (EPA)/a_NLA_OE_project/Project_repository/Routput/Scenario_modeling/SEM_output/XERw_m3_OE_fit.csv", #"inst/analysis_prediction/Routput/XERw_m3_OE_fit.csv",
          row.names = FALSE)

############################
# OUTPUT DATAFRAME OF MODEL PARAMETERS _ UNSTANDARDIZED
# n = 130 obs 15 vars
coef<-parTable(fit_v3_XERw_OE_bootstrap_rev)

################
# PARAMETER LABELS AND VALUES FOR PREDICTION MODELS
# APPLY FUNCTION TO PROCESS PARAMETERS AND GIVE COEFFICIENTS LABELS
coef_proc<-unstd_coeff(coef)

head(coef_proc)
############
## WRITE UNSTANDARIZED COEFFICIENTS - FOR PREDICTIONS
write.csv(coef_proc, "C:/Users/EFergus/OneDrive - Environmental Protection Agency (EPA)/a_NLA_OE_project/Project_repository/Routput/Scenario_modeling/SEM_output/XERw_m3_OE_coef_unstd.csv",#"inst/analysis_prediction/Routput/XERw_m3_OE_coef_unstd.csv",
          row.names=FALSE)


###########################################
## PATH MODEL v3 MMI
#   Model to use for prediction - removed PCTAGR_WS because low values and may be correlated w/ unaccounted for driver
mymodel_v3_MMI <- '
Lpt01_XCMGW ~ asin_PCTIMP_WsRp100 + W1_HAG + W1_HNOAG + L_NABD_NrmStorWs_ratio + asin_PCTFORGRS_CATRP100 + LOE_QLow_cl+ LOE_Qbkf_cl + L_STRM_POWER + PSUMPY_SY_WS_sc + drought_mean
L_STRM_POWER ~PSUMPY_SY_WS_sc + drought_mean + L_NABD_NrmStorWs_ratio + W1_HAG + asin_PCTFORGRS_CATRP100
LOE_QLow_cl~ asin_PCTIMP_WsRp100 + PSUMPY_SY_WS_sc + L_NABD_NrmStorWs_ratio + asin_PCTFORGRS_CATRP100 + W1_HAG + W1_HNOAG + drought_mean
LOE_Qbkf_cl ~ asin_PCTIMP_WsRp100 + PSUMPY_SY_WS_sc + L_NABD_NrmStorWs_ratio + asin_PCTFORGRS_CATRP100 + W1_HAG + W1_HNOAG + drought_mean
evap_index_sc ~ LOE_QLow_cl + LOE_Qbkf_cl + PSUMPY_SY_WS_sc + L_NABD_NrmStorWs_ratio + asin_PCTFORGRS_CATRP100 + drought_mean + L_STRM_POWER
LRBS_use ~ asin_PCTIMP_WsRp100 + L_NABD_NrmStorWs_ratio + Lpt01_XCMGW + LOE_QLow_cl + LOE_Qbkf_cl + W1_HAG + W1_HNOAG + asin_PCTFORGRS_CATRP100 + PSUMPY_SY_WS_sc + L_STRM_POWER + drought_mean
L_PTL ~ asin_PCTIMP_WsRp100 + Lpt01_XCMGW + evap_index_sc + W1_HAG + W1_HNOAG + asin_PCTFORGRS_CATRP100 + LOE_QLow_cl + LOE_Qbkf_cl + PSUMPY_SY_WS_sc + drought_mean
L_NTL ~ asin_PCTIMP_WsRp100 + Lpt01_XCMGW + evap_index_sc + W1_HAG + W1_HNOAG + asin_PCTFORGRS_CATRP100 + LOE_QLow_cl + LOE_Qbkf_cl + PSUMPY_SY_WS_sc + drought_mean
L_CHLR ~ asin_PCTIMP_WsRp100 + Lpt01_XCMGW + evap_index_sc + W1_HAG + W1_HNOAG + asin_PCTFORGRS_CATRP100 + LOE_QLow_cl + LOE_Qbkf_cl + PSUMPY_SY_WS_sc + drought_mean
L_SULF ~  asin_PCTIMP_WsRp100 + Lpt01_XCMGW + evap_index_sc + W1_HAG + W1_HNOAG + asin_PCTFORGRS_CATRP100 + LOE_QLow_cl + LOE_Qbkf_cl + PSUMPY_SY_WS_sc + drought_mean
L_TURB ~ asin_PCTIMP_WsRp100 + Lpt01_XCMGW + evap_index_sc + W1_HAG + W1_HNOAG + asin_PCTFORGRS_CATRP100 + LOE_QLow_cl + LOE_Qbkf_cl + PSUMPY_SY_WS_sc + drought_mean

MMI_BENT_sc ~ Lpt01_XCMGW + LRBS_use + L_PTL + L_NTL + L_CHLR + L_SULF + L_TURB + LOE_QLow_cl + LOE_Qbkf_cl + asin_PCTIMP_WsRp100 + L_NABD_NrmStorWs_ratio + asin_PCTFORGRS_CATRP100 + W1_HAG + W1_HNOAG + L_STRM_POWER + drought_mean

# Covariance
LOE_QLow_cl~~LOE_Qbkf_cl
L_STRM_POWER~~LOE_QLow_cl
L_STRM_POWER~~LOE_Qbkf_cl
L_PTL ~~ L_NTL
L_CHLR ~~ L_SULF
L_PTL ~~ L_TURB

'

##############
## XER MMI
# ROBUST ESTIMATION MAX LIKELIHOOD METHOD
fit_v3_XERw_MMI_robust.est<- sem(mymodel_v3_MMI, data=xer_w,
                                estimator="MLM")#estimator="MLM")

summary(fit_v3_XERw_MMI_robust.est, standardized=TRUE, fit.measures=TRUE, modindices=F,rsquare=TRUE)#, modindices=T, rsquare=TRUE) #

# request modification indices greater than 3.0 - from Grace USGS materials
mi_min <-modindices(fit_v3_XERw_MMI_robust.est)
print(mi_min[mi_min$mi >3.0,])


##########################
#####################
## MMI REVISED
mymodel_v3_MMI_rev <- '
Lpt01_XCMGW ~ 1 + ha3*W1_HAG + na2*W1_HNOAG + d2*L_NABD_NrmStorWs_ratio + sp3*L_STRM_POWER + p6*PSUMPY_SY_WS_sc
L_STRM_POWER ~ 1 + d3*L_NABD_NrmStorWs_ratio + rn4*asin_PCTFORGRS_CATRP100
LOE_QLow_cl~ 1 + i6*asin_PCTIMP_WsRp100 + p8*PSUMPY_SY_WS_sc + d5*L_NABD_NrmStorWs_ratio + rn5*asin_PCTFORGRS_CATRP100 + na4*W1_HNOAG
LOE_Qbkf_cl ~ 1 + i5*asin_PCTIMP_WsRp100 + p7*PSUMPY_SY_WS_sc + d4*L_NABD_NrmStorWs_ratio + ph3*drought_mean
evap_index_sc ~ 1 + l5*LOE_QLow_cl + b3*LOE_Qbkf_cl + p5*PSUMPY_SY_WS_sc + d6*L_NABD_NrmStorWs_ratio + sp4*L_STRM_POWER
LRBS_use ~ 1 + i2*asin_PCTIMP_WsRp100 + l2*LOE_QLow_cl + b2*LOE_Qbkf_cl + ha2*W1_HAG + rn2*asin_PCTFORGRS_CATRP100 + p2*PSUMPY_SY_WS_sc + sp2*L_STRM_POWER
L_NTL ~ 1 + i3*asin_PCTIMP_WsRp100 + x2*Lpt01_XCMGW + l3*LOE_QLow_cl + p3*PSUMPY_SY_WS_sc
L_SULF ~  1 + i4*asin_PCTIMP_WsRp100 + e2*evap_index_sc + na3*W1_HNOAG + rn3*asin_PCTFORGRS_CATRP100 + l4*LOE_QLow_cl + p4*PSUMPY_SY_WS_sc + ph2*drought_mean

MMI_BENT_sc ~ 1 + r1*LRBS_use + n1*L_NTL + su1*L_SULF + boe1*LOE_Qbkf_cl + d1*L_NABD_NrmStorWs_ratio + na1*W1_HNOAG + sp1*L_STRM_POWER + ph1*drought_mean

# INDIRECT EFFECTS ON MMI
imp_mmi:= i2*r1 + i3*n1 + i4*su1 + i5*boe1 + i5*b2*r1 + i5*b3*e2*su1 + i6*l2*r1 + i6*l3*n1 + i6*l4*su1 + i6*l5*e2*su1
dam_mmi:= d2*x2*n1 + d3*sp1 + d3*sp2*r1 + d3*sp3*x2*n1 + d3*sp4*e2*su1 + d4*boe1 + d4*b2*r1 + d4*b3*e2*su1 + d5*l2*r1 + d5*l3*n1 + d5*l4*su1 + d5*l5*e2*su1 + d6*e2*su1
precip_mmi:= p2*r1 + p3*n1 + p4*su1 + p5*e2*su1 + p6*x2*n1 + p7*boe1 + p7*b2*r1 + p7*b3*e2*su1 + p8*l2*r1 + p8*l3*n1 + p8*l4*su1 + p8*l5*e2*su1
phdi_mmi:= ph2*su1 + ph3*boe1 + ph3*b2*r1 + ph3*b3*e2*su1
strpwr_mmi:= sp2*r1 + sp3*x2*n1 + sp4*e2*su1
hag_mmi:= ha2*r1 + ha3*x2*n1
nhag_mmi:= na2*x2*n1 + na3*su1 + na4*l2*r1 + na4*l3*n1 + na4*l4*su1 + na4*l5*e2*su1
rnat_mmi:= rn2*r1 + rn3*su1 + rn4*sp1 + rn4*sp2*r1 + rn4*sp3*x2*n1 + rn4*sp4*e2*su1 + rn5*l2*r1 + rn5*l3*n1 + rn5*l4*su1 + rn5*l5*e2*su1
oelflow_mmi:= l2*r1 + l3*n1 + l4*su1 + l5*e2*su1
oebflow_mmi:= b2*r1 + b3*e2*su1
dexcess_mmi:= e2*su1
xcmgw_mmi:= x2*n1

# TOTAL EFFECTS ON MMI
imp_tot:= imp_mmi
dam_tot:= d1 + dam_mmi
precip_tot:= precip_mmi
phdi_tot:= ph1 + phdi_mmi
strpwr_tot:= sp1 + strpwr_mmi
hag_tot:= hag_mmi
nhag_tot:= na1 + nhag_mmi
rnat_tot:= rnat_mmi
oelflow_tot:= oelflow_mmi
oebflow_tot:= boe1 + oebflow_mmi
dexcess_tot:= dexcess_mmi
xcmgw_tot:= xcmgw_mmi
rbs_tot:= r1
tn_tot:= n1
sulf_tot:= su1

# Covariance
LOE_QLow_cl~~LOE_Qbkf_cl
L_STRM_POWER~~LOE_QLow_cl
L_STRM_POWER~~LOE_Qbkf_cl
L_SULF~~L_NTL

'

# ROBUST ESTIMATION MAX LIKELIHOOD METHOD
fit_v3_XERw_MMI_rev_robust.est<- sem(mymodel_v3_MMI_rev, data=xer_w,
                                    estimator="MLM")

summary(fit_v3_XERw_MMI_rev_robust.est, standardized=TRUE, fit.measures=TRUE, modindices=F,rsquare=TRUE)#, modindices=T, rsquare=TRUE) #

# request modification indices greater than 3.0 - from Grace USGS materials
mi_min <-modindices(fit_v3_XERw_MMI_rev_robust.est)
print(mi_min[mi_min$mi >3.0,])

######################
## Bollen.stine bootstrap to estimate parameters -MMI
fit_v3_XERw_MMI_bootstrap_rev  <- sem(mymodel_v3_MMI_rev, data=xer_w,
                                     #group = "ECOREG_rev",
                                     #missing="ML",
                                     test="bollen.stine", se="boot",bootstrap=1000)

summary(fit_v3_XERw_MMI_bootstrap_rev, standardized=TRUE, fit.measures=TRUE, modindices=F,rsquare=TRUE)#, modindices=T, rsquare=TRUE) #

#############
# Export R output - XERIC MMI MODEL
#https://www.r-bloggers.com/export-r-output-to-a-file/
out_fit_v3_XERw_MMI_rev<- capture.output(summary(fit_v3_XERw_MMI_bootstrap_rev, standardized=FALSE, fit.measures=TRUE,rsquare=TRUE)) #, modindices=T
write.csv(out_fit_v3_XERw_MMI_rev,"C:/Users/EFergus/OneDrive - Environmental Protection Agency (EPA)/a_NLA_OE_project/Project_repository/Routput/Scenario_modeling/SEM_output/XERw_m3_MMI_rev.csv" , #"inst/analysis_prediction/Routput/XERw_m3_MMI_rev.csv",
          row.names=FALSE)

# Standardized estimates of bootstrap model
std_parameter_se_bootstrap_min<- standardizedSolution(fit_v3_XERw_MMI_bootstrap_rev)
write.csv(std_parameter_se_bootstrap_min, 'C:/Users/EFergus/OneDrive - Environmental Protection Agency (EPA)/a_NLA_OE_project/Project_repository/Routput/Scenario_modeling/SEM_output/XERw_m3_MMI_CI.csv', #"inst/analysis_prediction/Routput/XERw_m3_MMI_CI.csv",
          row.names = FALSE)

# Save model fit parameters
## GRAB R2
r2_XER_MMI<-as.data.frame(lavInspect(fit_v3_XERw_MMI_bootstrap_rev, "rsquare"))

r2_XER_MMI<-cbind(variable=rownames(r2_XER_MMI),r2_XER_MMI)
rownames(r2_XER_MMI) <- 1:nrow(r2_XER_MMI) # To remove index column

colnames(r2_XER_MMI) <- c("Variable","R2")
write.csv(r2_XER_MMI,"C:/Users/EFergus/OneDrive - Environmental Protection Agency (EPA)/a_NLA_OE_project/Project_repository/Routput/Scenario_modeling/SEM_output/XERw_m3_MMI_R2.csv", #"inst/analysis_prediction/Routput/XERw_m3_MMI_R2.csv",
          row.names=FALSE)

#################
# GRAB R2 to be able to add to fit table
r2_red_MMI <- r2_XER_MMI%>%
  filter(Variable=="MMI_BENT_sc")%>%
  select(R2)%>%
  mutate(across(where(is.numeric), round,2))

# TABLE OF FIT INDICES comparing 3 RESPONSES
table_fit_mmi_xer <- matrix(NA, nrow=2, ncol=11)
colnames(table_fit_mmi_xer) = c("Model","Estimation","X2", "df","CFI","TLI", "RMSEA","SRMR","AIC","n","npar")

table_fit_mmi_xer[1,]<-c("MMI","Standard",round(fitmeasures(fit_v3_XERw_MMI_rev_robust.est,
                                                          c("chisq","df","cfi","tli",
                                                            "rmsea","srmr","aic","ntotal","npar")),2))
table_fit_mmi_xer[2,]<-c("MMI","Robust",round(fitmeasures(fit_v3_XERw_MMI_rev_robust.est,
                                                        c("chisq.scaled","df.scaled","cfi.scaled","tli.scaled",
                                                          "rmsea.robust","srmr","aic","ntotal","npar")),2))
# Add R2 for biotic response
table_fit_mmi_xer<- data.frame(table_fit_mmi_xer)%>%
  mutate(R2 = r2_red_MMI$R2)

write.csv(table_fit_mmi_xer, "C:/Users/EFergus/OneDrive - Environmental Protection Agency (EPA)/a_NLA_OE_project/Project_repository/Routput/Scenario_modeling/SEM_output/XERw_m3_MMI_fit.csv", #"inst/analysis_prediction/Routput/XERw_m3_MMI_fit.csv",
          row.names = FALSE)

############################
# OUTPUT DATAFRAME OF MODEL PARAMETERS _ UNSTANDARDIZED
# n = 130 obs 15 vars
coef<-parTable(fit_v3_XERw_MMI_bootstrap_rev)

################
# PARAMETER LABELS AND VALUES FOR PREDICTION MODELS
# APPLY FUNCTION TO PROCESS PARAMETERS AND GIVE COEFFICIENTS LABELS
coef_proc<-unstd_coeff(coef)

head(coef_proc)
############
## WRITE UNSTANDARIZED COEFFICIENTS - FOR PREDICTIONS
write.csv(coef_proc, "C:/Users/EFergus/OneDrive - Environmental Protection Agency (EPA)/a_NLA_OE_project/Project_repository/Routput/Scenario_modeling/SEM_output/XERw_m3_MMI_coef_unstd.csv", #"inst/analysis_prediction/Routput/XERw_m3_MMI_coef_unstd.csv",
          row.names=FALSE)


#############################################
###########################################
## PATH MODEL v3 EPT
#   Model to use for prediction - removed PCTAGR_WS because low values and may be correlated w/ unaccounted for driver
mymodel_v3_EPT <- '
Lpt01_XCMGW ~ asin_PCTIMP_WsRp100 + W1_HAG + W1_HNOAG + L_NABD_NrmStorWs_ratio + asin_PCTFORGRS_CATRP100 + LOE_QLow_cl+ LOE_Qbkf_cl + L_STRM_POWER + PSUMPY_SY_WS_sc + drought_mean
L_STRM_POWER ~PSUMPY_SY_WS_sc + drought_mean + L_NABD_NrmStorWs_ratio + W1_HAG + asin_PCTFORGRS_CATRP100
LOE_QLow_cl~ asin_PCTIMP_WsRp100 + PSUMPY_SY_WS_sc + L_NABD_NrmStorWs_ratio + asin_PCTFORGRS_CATRP100 + W1_HAG + W1_HNOAG + drought_mean
LOE_Qbkf_cl ~ asin_PCTIMP_WsRp100 + PSUMPY_SY_WS_sc + L_NABD_NrmStorWs_ratio + asin_PCTFORGRS_CATRP100 + W1_HAG + W1_HNOAG + drought_mean
evap_index_sc ~ LOE_QLow_cl + LOE_Qbkf_cl + PSUMPY_SY_WS_sc + L_NABD_NrmStorWs_ratio + asin_PCTFORGRS_CATRP100 + drought_mean + L_STRM_POWER
LRBS_use ~ asin_PCTIMP_WsRp100 + L_NABD_NrmStorWs_ratio + Lpt01_XCMGW + LOE_QLow_cl + LOE_Qbkf_cl + W1_HAG + W1_HNOAG + asin_PCTFORGRS_CATRP100 + PSUMPY_SY_WS_sc + L_STRM_POWER + drought_mean
L_PTL ~ asin_PCTIMP_WsRp100 + Lpt01_XCMGW + evap_index_sc + W1_HAG + W1_HNOAG + asin_PCTFORGRS_CATRP100 + LOE_QLow_cl + LOE_Qbkf_cl + PSUMPY_SY_WS_sc + drought_mean
L_NTL ~ asin_PCTIMP_WsRp100 + Lpt01_XCMGW + evap_index_sc + W1_HAG + W1_HNOAG + asin_PCTFORGRS_CATRP100 + LOE_QLow_cl + LOE_Qbkf_cl + PSUMPY_SY_WS_sc + drought_mean
L_CHLR ~ asin_PCTIMP_WsRp100 + Lpt01_XCMGW + evap_index_sc + W1_HAG + W1_HNOAG + asin_PCTFORGRS_CATRP100 + LOE_QLow_cl + LOE_Qbkf_cl + PSUMPY_SY_WS_sc + drought_mean
L_SULF ~  asin_PCTIMP_WsRp100 + Lpt01_XCMGW + evap_index_sc + W1_HAG + W1_HNOAG + asin_PCTFORGRS_CATRP100 + LOE_QLow_cl + LOE_Qbkf_cl + PSUMPY_SY_WS_sc + drought_mean
L_TURB ~ asin_PCTIMP_WsRp100 + Lpt01_XCMGW + evap_index_sc + W1_HAG + W1_HNOAG + asin_PCTFORGRS_CATRP100 + LOE_QLow_cl + LOE_Qbkf_cl + PSUMPY_SY_WS_sc + drought_mean

EPT_RICH_sc ~ Lpt01_XCMGW + LRBS_use + L_PTL + L_NTL + L_CHLR + L_SULF + L_TURB + LOE_QLow_cl + LOE_Qbkf_cl + asin_PCTIMP_WsRp100 + L_NABD_NrmStorWs_ratio + asin_PCTFORGRS_CATRP100 + W1_HAG + W1_HNOAG + L_STRM_POWER + drought_mean

# Covariance
LOE_QLow_cl~~LOE_Qbkf_cl
L_STRM_POWER~~LOE_QLow_cl
L_STRM_POWER~~LOE_Qbkf_cl
L_PTL ~~ L_NTL
L_CHLR ~~ L_SULF
L_PTL ~~ L_TURB

'

##############
## XER EPT
# ROBUST ESTIMATION MAX LIKELIHOOD METHOD
fit_v3_XERw_EPT_robust.est<- sem(mymodel_v3_EPT, data=xer_w,
                                 estimator="MLM")#estimator="MLM")

summary(fit_v3_XERw_EPT_robust.est, standardized=TRUE, fit.measures=TRUE, modindices=F,rsquare=TRUE)#, modindices=T, rsquare=TRUE) #

# request modification indices greater than 3.0 - from Grace USGS materials
mi_min <-modindices(fit_v3_XERw_EPT_robust.est)
print(mi_min[mi_min$mi >3.0,])


##########################
#####################
## EPT REVISED
mymodel_v3_EPT_rev <- '
Lpt01_XCMGW ~ 1 + ha3*W1_HAG + na3*W1_HNOAG + d3*L_NABD_NrmStorWs_ratio + sp3*L_STRM_POWER + p5*PSUMPY_SY_WS_sc
L_STRM_POWER ~ 1 + d2*L_NABD_NrmStorWs_ratio + rn4*asin_PCTFORGRS_CATRP100
LOE_QLow_cl~ 1 + i5*asin_PCTIMP_WsRp100 + p8*PSUMPY_SY_WS_sc + d6*L_NABD_NrmStorWs_ratio + rn5*asin_PCTFORGRS_CATRP100 + ha5*W1_HAG + na4*W1_HNOAG
LOE_Qbkf_cl ~ 1 + p7*PSUMPY_SY_WS_sc + d5*L_NABD_NrmStorWs_ratio + ha4*W1_HAG + ph3*drought_mean
evap_index_sc ~ 1 + l5*LOE_QLow_cl + b3*LOE_Qbkf_cl + p6*PSUMPY_SY_WS_sc + d4*L_NABD_NrmStorWs_ratio + sp4*L_STRM_POWER
LRBS_use ~ 1 + i2*asin_PCTIMP_WsRp100 + l2*LOE_QLow_cl + b2*LOE_Qbkf_cl + ha2*W1_HAG + rn2*asin_PCTFORGRS_CATRP100 + p2*PSUMPY_SY_WS_sc + sp2*L_STRM_POWER
L_NTL ~ 1 + i3*asin_PCTIMP_WsRp100 + x2*Lpt01_XCMGW + l3*LOE_QLow_cl + p3*PSUMPY_SY_WS_sc
L_SULF ~  1 + i4*asin_PCTIMP_WsRp100 + e2*evap_index_sc + na2*W1_HNOAG + rn3*asin_PCTFORGRS_CATRP100 + l4*LOE_QLow_cl + p4*PSUMPY_SY_WS_sc + ph2*drought_mean

EPT_RICH_sc ~ 1 + x1*Lpt01_XCMGW + r1*LRBS_use + n1*L_NTL + su1*L_SULF + boe1*LOE_Qbkf_cl + d1*L_NABD_NrmStorWs_ratio + na1*W1_HNOAG + ph1*drought_mean

# INDIRECT EFFECTS ON EPT
imp_ept:= i2*r1 + i3*n1 + i4*su1 + i5*l2*r1 + i5*l3*n1 + i5*l4*su1 + i5*l5*e2*su1
dam_ept:= d2*sp2*r1 + d2*sp3*x1 + d2*sp3*x2*n1 + d2*sp4*e2*su1 + d3*x1 + d3*x2*n1 + d4*e2*su1 + d5*boe1 + d5*b2*r1 + d5*b3*e2*su1 + d6*l2*r1 + d6*l3*n1 + d6*l4*su1 + d6*l5*e2*su1
precip_ept:= p2*r1 + p3*n1 + p4*su1 + p5*x1 + p5*x2*n1 + p6*e2*su1 + p7*boe1 + p7*b2*r1 + p7*b3*e2*su1 + p8*l2*r1 + p8*l3*n1 + p8*l4*su1 + p8*l5*e2*su1
phdi_ept:= ph2*su1 + ph3*boe1 + ph3*b2*r1 + ph3*b3*e2*su1
strpwr_ept:= sp2*r1 + sp3*x1 + sp3*x2*n1 + sp4*e2*su1
hag_ept:= ha2*r1 + ha3*x1 + ha3*x2*n1 + ha4*boe1 + ha4*b2*r1 + ha4*b3*e2*su1 + ha5*l2*r1 + ha5*l3*n1 + ha5*l4*su1 + ha5*l5*e2*su1
nhag_ept:= na2*su1 + na3*x1 + na3*x2*n1 + na4*l2*r1 + na4*l3*n1 + na4*l4*su1 + na4*l5*e2*su1
rnat_ept:= rn2*r1 + rn3*su1 + rn4*sp2*r1 + rn4*sp3*x1 + rn4*sp3*x2*n1 + rn5*l2*r1 + rn5*l3*n1 + rn5*l4*su1 + rn5*l5*e2*su1
oelflow_ept:= l2*r1 + l3*n1 + l4*su1 + l5*e2*su1
oebflow_ept:= b2*r1 + b3*e2*su1
dexcess_ept:= e2*su1
xcmgw_ept:= x2*n1

# TOTAL EFFECTS ON EPT
imp_tot:= imp_ept
dam_tot:= d1 + dam_ept
precip_tot:= precip_ept
phdi_tot:= ph1 + phdi_ept
strpwr_tot:= strpwr_ept
hag_tot:= hag_ept
nhag_tot:= na1 + nhag_ept
rnat_tot:= rnat_ept
oelflow_tot:= oelflow_ept
oebflow_tot:= boe1 + oebflow_ept
dexcess_tot:= dexcess_ept
xcmgw_tot:= x1 + xcmgw_ept
rbs_tot:= r1
tn_tot:= n1
sulf_tot:= su1

# Covariance
LOE_QLow_cl~~LOE_Qbkf_cl
L_STRM_POWER~~LOE_QLow_cl
#L_STRM_POWER~~LOE_Qbkf_cl
L_SULF~~L_NTL

'


# ROBUST ESTIMATION MAX LIKELIHOOD METHOD
fit_v3_XERw_EPT_rev_robust.est<- sem(mymodel_v3_EPT_rev, data=xer_w,
                                    estimator="MLM")

summary(fit_v3_XERw_EPT_rev_robust.est, standardized=TRUE, fit.measures=TRUE, modindices=F,rsquare=TRUE)#, modindices=T, rsquare=TRUE) #

# request modification indices greater than 3.0 - from Grace USGS materials
mi_min <-modindices(fit_v3_XERw_EPT_rev_robust.est)
print(mi_min[mi_min$mi >3.0,])

######################
## Bollen.stine bootstrap to estimate parameters -EPT
fit_v3_XERw_EPT_bootstrap_rev  <- sem(mymodel_v3_EPT_rev, data=xer_w,
                                     #group = "ECOREG_rev",
                                     #missing="ML",
                                     test="bollen.stine", se="boot",bootstrap=1000)

summary(fit_v3_XERw_EPT_bootstrap_rev, standardized=TRUE, fit.measures=TRUE, modindices=F,rsquare=TRUE)#, modindices=T, rsquare=TRUE) #

#############
# Export R output - XERIC EPT MODEL
#https://www.r-bloggers.com/export-r-output-to-a-file/
out_fit_v3_XERw_EPT_rev<- capture.output(summary(fit_v3_XERw_EPT_bootstrap_rev, standardized=FALSE, fit.measures=TRUE,rsquare=TRUE)) #, modindices=T
write.csv(out_fit_v3_XERw_EPT_rev,"C:/Users/EFergus/OneDrive - Environmental Protection Agency (EPA)/a_NLA_OE_project/Project_repository/Routput/Scenario_modeling/SEM_output/XERw_m3_EPT_rev.csv", #"inst/analysis_prediction/Routput/XERw_m3_EPT_rev.csv",
          row.names=FALSE)

# Standardized estimates of bootstrap model
std_parameter_se_bootstrap_min<- standardizedSolution(fit_v3_XERw_EPT_bootstrap_rev)
write.csv(std_parameter_se_bootstrap_min, "C:/Users/EFergus/OneDrive - Environmental Protection Agency (EPA)/a_NLA_OE_project/Project_repository/Routput/Scenario_modeling/SEM_output/XERw_m3_EPT_CI.csv", #"inst/analysis_prediction/Routput/XERw_m3_EPT_CI.csv",
          row.names = FALSE)

# Save model fit parameters
## GRAB R2
r2_XER_EPT<-as.data.frame(lavInspect(fit_v3_XERw_EPT_bootstrap_rev, "rsquare"))

r2_XER_EPT<-cbind(variable=rownames(r2_XER_EPT),r2_XER_EPT)
rownames(r2_XER_EPT) <- 1:nrow(r2_XER_EPT) # To remove index column

colnames(r2_XER_EPT) <- c("Variable","R2")
write.csv(r2_XER_EPT,"C:/Users/EFergus/OneDrive - Environmental Protection Agency (EPA)/a_NLA_OE_project/Project_repository/Routput/Scenario_modeling/SEM_output/XERw_m3_EPT_R2.csv",#"inst/analysis_prediction/Routput/XERw_m3_EPT_R2.csv",
          row.names=FALSE)

#################
# GRAB R2 to be able to add to fit table
r2_red_EPT <- r2_XER_EPT%>%
  filter(Variable=="EPT_RICH_sc")%>%
  select(R2)%>%
  mutate(across(where(is.numeric), round,2))

# TABLE OF FIT INDICES comparing 3 RESPONSES
table_fit_ept_xer <- matrix(NA, nrow=2, ncol=11)
colnames(table_fit_ept_xer) = c("Model","Estimation","X2", "df","CFI","TLI", "RMSEA","SRMR","AIC","n","npar")

table_fit_ept_xer[1,]<-c("EPT","Standard",round(fitmeasures(fit_v3_XERw_EPT_rev_robust.est,
                                                          c("chisq","df","cfi","tli",
                                                            "rmsea","srmr","aic","ntotal","npar")),2))
table_fit_ept_xer[2,]<-c("EPT","Robust",round(fitmeasures(fit_v3_XERw_EPT_rev_robust.est,
                                                        c("chisq.scaled","df.scaled","cfi.scaled","tli.scaled",
                                                          "rmsea.robust","srmr","aic","ntotal","npar")),2))
# Add R2 for biotic response
table_fit_ept_xer<- data.frame(table_fit_ept_xer)%>%
  mutate(R2 = r2_red_EPT$R2)

write.csv(table_fit_ept_xer, "C:/Users/EFergus/OneDrive - Environmental Protection Agency (EPA)/a_NLA_OE_project/Project_repository/Routput/Scenario_modeling/SEM_output/XERw_m3_EPT_fit.csv",#"inst/analysis_prediction/Routput/XERw_m3_EPT_fit.csv",
          row.names = FALSE)

############################
# OUTPUT DATAFRAME OF MODEL PARAMETERS _ UNSTANDARDIZED
# n = 130 obs 15 vars
coef<-parTable(fit_v3_XERw_EPT_bootstrap_rev)

################
# PARAMETER LABELS AND VALUES FOR PREDICTION MODELS
# APPLY FUNCTION TO PROCESS PARAMETERS AND GIVE COEFFICIENTS LABELS
coef_proc<-unstd_coeff(coef)

head(coef_proc)
############
## WRITE UNSTANDARIZED COEFFICIENTS - FOR PREDICTIONS
write.csv(coef_proc,"C:/Users/EFergus/OneDrive - Environmental Protection Agency (EPA)/a_NLA_OE_project/Project_repository/Routput/Scenario_modeling/SEM_output/XERw_m3_EPT_coef_unstd.csv" ,#"inst/analysis_prediction/Routput/",
          row.names=FALSE)



############################################
###########################################
## PATH MODEL v3b O/E - Developed land in Ws 100 buffer
#  TO COMPARE TO MODEL WITH IMPERVIOUS - Larger AIC and developed land in ws 100 buffer is maybe harder to characterize bc includes impervious and other development activities
#
mymodel_v3b_OE <- '
Lpt01_XCMGW ~ asin_PCTURB_WsRp100 + W1_HAG + W1_HNOAG + L_NABD_NrmStorWs_ratio + asin_PCTFORGRS_CATRP100 + LOE_QLow_cl+ LOE_Qbkf_cl + L_STRM_POWER + PSUMPY_SY_WS_sc + drought_mean
L_STRM_POWER ~PSUMPY_SY_WS_sc + drought_mean + L_NABD_NrmStorWs_ratio + W1_HAG + asin_PCTFORGRS_CATRP100
LOE_QLow_cl~ asin_PCTURB_WsRp100 + PSUMPY_SY_WS_sc + L_NABD_NrmStorWs_ratio + asin_PCTFORGRS_CATRP100 + W1_HAG + W1_HNOAG + drought_mean
LOE_Qbkf_cl ~ asin_PCTURB_WsRp100 + PSUMPY_SY_WS_sc + L_NABD_NrmStorWs_ratio + asin_PCTFORGRS_CATRP100 + W1_HAG + W1_HNOAG + drought_mean
evap_index_sc ~ LOE_QLow_cl + LOE_Qbkf_cl + PSUMPY_SY_WS_sc + L_NABD_NrmStorWs_ratio + asin_PCTFORGRS_CATRP100 + drought_mean + L_STRM_POWER
LRBS_use ~ asin_PCTURB_WsRp100 + L_NABD_NrmStorWs_ratio + Lpt01_XCMGW + LOE_QLow_cl + LOE_Qbkf_cl + W1_HAG + W1_HNOAG + asin_PCTFORGRS_CATRP100 + PSUMPY_SY_WS_sc + L_STRM_POWER + drought_mean
L_PTL ~ asin_PCTURB_WsRp100 + Lpt01_XCMGW + evap_index_sc + W1_HAG + W1_HNOAG + asin_PCTFORGRS_CATRP100 + LOE_QLow_cl + LOE_Qbkf_cl + PSUMPY_SY_WS_sc + drought_mean
L_NTL ~ asin_PCTURB_WsRp100 + Lpt01_XCMGW + evap_index_sc + W1_HAG + W1_HNOAG + asin_PCTFORGRS_CATRP100 + LOE_QLow_cl + LOE_Qbkf_cl + PSUMPY_SY_WS_sc + drought_mean
L_CHLR ~ asin_PCTURB_WsRp100 + Lpt01_XCMGW + evap_index_sc + W1_HAG + W1_HNOAG + asin_PCTFORGRS_CATRP100 + LOE_QLow_cl + LOE_Qbkf_cl + PSUMPY_SY_WS_sc + drought_mean
L_SULF ~  asin_PCTURB_WsRp100 + Lpt01_XCMGW + evap_index_sc + W1_HAG + W1_HNOAG + asin_PCTFORGRS_CATRP100 + LOE_QLow_cl + LOE_Qbkf_cl + PSUMPY_SY_WS_sc + drought_mean
L_TURB ~ asin_PCTURB_WsRp100 + Lpt01_XCMGW + evap_index_sc + W1_HAG + W1_HNOAG + asin_PCTFORGRS_CATRP100 + LOE_QLow_cl + LOE_Qbkf_cl + PSUMPY_SY_WS_sc + drought_mean

OE_SCORE ~ Lpt01_XCMGW + LRBS_use + L_PTL + L_NTL + L_CHLR + L_SULF + L_TURB + LOE_QLow_cl + LOE_Qbkf_cl + asin_PCTURB_WsRp100 + L_NABD_NrmStorWs_ratio + asin_PCTFORGRS_CATRP100 + W1_HAG + W1_HNOAG + L_STRM_POWER + drought_mean

# Covariance
LOE_QLow_cl~~LOE_Qbkf_cl
L_STRM_POWER~~LOE_QLow_cl
L_STRM_POWER~~LOE_Qbkf_cl
L_PTL ~~ L_NTL
L_CHLR ~~ L_SULF
L_PTL ~~ L_TURB

'

##############
## XER O/E
# ROBUST ESTIMATION MAX LIKELIHOOD METHOD
fit_v3b_XERw_OE_robust.est<- sem(mymodel_v3b_OE, data=xer_w,
                                estimator="MLM")#estimator="MLM")

summary(fit_v3b_XERw_OE_robust.est, standardized=TRUE, fit.measures=TRUE, modindices=F,rsquare=TRUE)#, modindices=T, rsquare=TRUE) #

# request modification indices greater than 3.0 - from Grace USGS materials
mi_min <-modindices(fit_v3b_XERw_OE_robust.est)
print(mi_min[mi_min$mi >3.0,])


##########################
#####################
## REVISED - Added total effects
mymodel_v3b_OE_rev <- '
Lpt01_XCMGW ~ W1_HAG + W1_HNOAG + L_NABD_NrmStorWs_ratio + L_STRM_POWER + PSUMPY_SY_WS_sc
L_STRM_POWER ~ L_NABD_NrmStorWs_ratio + asin_PCTFORGRS_CATRP100
LOE_QLow_cl~ PSUMPY_SY_WS_sc + L_NABD_NrmStorWs_ratio + W1_HAG + W1_HNOAG
LOE_Qbkf_cl ~ PSUMPY_SY_WS_sc + L_NABD_NrmStorWs_ratio + W1_HAG
evap_index_sc ~ LOE_QLow_cl + LOE_Qbkf_cl + PSUMPY_SY_WS_sc + L_NABD_NrmStorWs_ratio + L_STRM_POWER
LRBS_use ~ asin_PCTURB_WsRp100 + LOE_QLow_cl + LOE_Qbkf_cl + W1_HAG + asin_PCTFORGRS_CATRP100 + PSUMPY_SY_WS_sc + L_STRM_POWER
L_NTL ~ asin_PCTURB_WsRp100 + Lpt01_XCMGW + W1_HNOAG + LOE_QLow_cl + PSUMPY_SY_WS_sc
L_SULF ~  asin_PCTURB_WsRp100 + evap_index_sc + W1_HNOAG + asin_PCTFORGRS_CATRP100 + LOE_QLow_cl + PSUMPY_SY_WS_sc + drought_mean

OE_SCORE ~ LRBS_use + L_NTL + L_SULF + LOE_Qbkf_cl + L_NABD_NrmStorWs_ratio + L_STRM_POWER

# Covariance
LOE_QLow_cl~~LOE_Qbkf_cl
L_STRM_POWER~~LOE_QLow_cl
#L_STRM_POWER~~LOE_Qbkf_cl
L_SULF~~L_NTL

'

# ROBUST ESTIMATION MAX LIKELIHOOD METHOD
fit_v3b_XERw_OE_rev_robust.est<- sem(mymodel_v3b_OE_rev, data=xer_w,
                                    estimator="MLM")

summary(fit_v3b_XERw_OE_rev_robust.est, standardized=TRUE, fit.measures=TRUE, modindices=F,rsquare=TRUE)#, modindices=T, rsquare=TRUE) #

# request modification indices greater than 3.0 - from Grace USGS materials
mi_min <-modindices(fit_v3b_XERw_OE_rev_robust.est)
print(mi_min[mi_min$mi >3.0,])

######################
## Bollen.stine bootstrap to estimate parameters -OE
fit_v3b_XERw_OE_bootstrap_rev  <- sem(mymodel_v3b_OE_rev, data=xer_w,
                                     #group = "ECOREG_rev",
                                     #missing="ML",
                                     test="bollen.stine", se="boot",bootstrap=1000)

summary(fit_v3b_XERw_OE_bootstrap_rev, standardized=TRUE, fit.measures=TRUE, modindices=F,rsquare=TRUE)#, modindices=T, rsquare=TRUE) #

# AIC = 4017
