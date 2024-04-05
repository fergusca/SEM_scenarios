#####################
## REVISED PATH ANALYSIS MODEL FOR PREDICTION - EPT_RICH_sc
##  XERIC wadeable
##  Model v8
##  Including habitat heterogeneity presence (richness) index
##  Phil suggested including after modifying habitat index to better capture hetergeneity - not just cover
##
##  1/31/2024
##  2/13/2024 - revising bank incision to be scaled by bankfull height
##  2/26/2024 - breaking up components of RBS - mean diameter (LSUB_DMM) and critical diameter (LDCBF_use)
##  2/20/2024 - drop critical diameter (redundant w/ specific strm power and bankfull flow) and replacing mean diameter wiht % sands/fines
##  3/6/2024  - replaced rbs with % sands fines

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
library(SEMscenarios)

###########################
## READ PROCESSED DATA
# COMPILED NRSA SURVEYS WITH SUBSET OF VARIABLES
#  INCLUDES ALL RESAMPLED SITES AND VISITS 1 & 2

# PROCESSED - all 0809 and only new sites from later surveys
#  VISITS 1 and 2 n = 4578 w/ 342 vars
dat_org <- read_csv("data_processed/Compiled/nrsa081318_nonresampled_VISIT_12.csv")

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

# CALCULATE BANK INCISION above bankfull channel
dat_proc<-dat_proc%>%
  mutate(INCISION = XINC_H-XBKF_H) # XINC_H has 71 NAs; XBKF_H has 8 NAs
# Distribution looks relatively normal so no need to transform

# MAKE % SANDS AND FINES A PROPORTION
dat_proc<-dat_proc%>%
  mutate(PCT_SAFN_sc = PCT_SAFN/100)

# TRANSFORM IMPERVIOUS SURFACE (NLCD
dat_proc<-dat_proc%>%
  mutate(asin_PCTIMP_WS = asin(sqrt(PCTIMP_WS/100)),
         asin_PCTIMP_WsRp100 = asin(sqrt(PCTIMP_WsRp100/100)),
         asin_PCTIMP_CAT = asin(sqrt(PCTIMP_CAT/100)),
         asin_PCTIMP_CATRP100 = asin(sqrt(PCTIMP_CATRP100/100)))

summary(dat_proc$asin_PCTIMP_CATRP100)


# TRANSFORM MEAN INCISION HT
dat_proc<-dat_proc%>%
  mutate(L_XINC_H = log10(XINC_H+0.01))


#########################
## SUBSET BY MANUAL ECOREGION
# XER n = 342
xer<-dat_proc%>%
  filter(AG_ECO9=="XER")%>%
  drop_na(LOE_QLow_cl)

# SCALE CUMULATIVE PRECIPITATION in XER
xer$PSUMPY_SY_WS_sc<-scale(xer$PSUMPY_SY_WS)

# WADEABLE n = 272
xer_w <- xer %>%
  filter(PROTOCOL=="WADEABLE")

#################

###########################################
## PATH MODEL v8 EPT
mymodel_v8_EPT <- '
Lpt01_XCMGW ~ asin_PCTIMP_CATRP100 + W1_HAG + W1_HNOAG + L_NABD_NrmStorWs_ratio + asin_PCTFORGRS_CATRP100 + LOE_QLow_cl+ LOE_Qbkf_cl + L_STRM_POWER + PSUMPY_SY_WS_sc
LOE_QLow_cl~ asin_PCTIMP_CATRP100 + PSUMPY_SY_WS_sc + L_NABD_NrmStorWs_ratio + asin_PCTFORGRS_CATRP100 + W1_HAG + W1_HNOAG
LOE_Qbkf_cl ~ asin_PCTIMP_CATRP100 + PSUMPY_SY_WS_sc + L_NABD_NrmStorWs_ratio + asin_PCTFORGRS_CATRP100 + W1_HAG + W1_HNOAG
evap_index_sc ~ LOE_QLow_cl + LOE_Qbkf_cl + PSUMPY_SY_WS_sc + L_NABD_NrmStorWs_ratio + asin_PCTFORGRS_CATRP100 + L_STRM_POWER
PCT_SAFN_sc ~ asin_PCTIMP_CATRP100 + L_NABD_NrmStorWs_ratio + Lpt01_XCMGW + LOE_QLow_cl + LOE_Qbkf_cl + W1_HAG + W1_HNOAG + asin_PCTFORGRS_CATRP100 + PSUMPY_SY_WS_sc + L_STRM_POWER + INCISION
L_PTL ~ asin_PCTIMP_CATRP100 + Lpt01_XCMGW + evap_index_sc + W1_HAG + W1_HNOAG + asin_PCTFORGRS_CATRP100 + LOE_QLow_cl + LOE_Qbkf_cl + PSUMPY_SY_WS_sc + drought_mean
L_NTL ~ asin_PCTIMP_CATRP100 + Lpt01_XCMGW + evap_index_sc + W1_HAG + W1_HNOAG + asin_PCTFORGRS_CATRP100 + LOE_QLow_cl + LOE_Qbkf_cl + PSUMPY_SY_WS_sc + drought_mean
L_CHLR ~ asin_PCTIMP_CATRP100 + Lpt01_XCMGW + evap_index_sc + W1_HAG + W1_HNOAG + asin_PCTFORGRS_CATRP100 + LOE_QLow_cl + LOE_Qbkf_cl + PSUMPY_SY_WS_sc + drought_mean
L_SULF ~  asin_PCTIMP_CATRP100 + Lpt01_XCMGW + evap_index_sc + W1_HAG + W1_HNOAG + asin_PCTFORGRS_CATRP100 + LOE_QLow_cl + LOE_Qbkf_cl + PSUMPY_SY_WS_sc + drought_mean
L_TURB ~ asin_PCTIMP_CATRP100 + Lpt01_XCMGW + evap_index_sc + W1_HAG + W1_HNOAG + asin_PCTFORGRS_CATRP100 + LOE_QLow_cl + LOE_Qbkf_cl + PSUMPY_SY_WS_sc + INCISION
INCISION ~ asin_PCTIMP_CATRP100 + W1_HAG + W1_HNOAG + L_NABD_NrmStorWs_ratio + asin_PCTFORGRS_CATRP100 + LOE_QLow_cl+ LOE_Qbkf_cl + L_STRM_POWER + PSUMPY_SY_WS_sc + drought_mean
SumBig_PFC ~ asin_PCTIMP_CATRP100 + L_NABD_NrmStorWs_ratio + Lpt01_XCMGW + LOE_QLow_cl + LOE_Qbkf_cl + W1_HAG + W1_HNOAG + asin_PCTFORGRS_CATRP100 + PSUMPY_SY_WS_sc + L_STRM_POWER + INCISION

EPT_RICH_sc ~ Lpt01_XCMGW + PCT_SAFN_sc + L_PTL + L_NTL + L_CHLR + L_SULF + L_TURB + LOE_QLow_cl + LOE_Qbkf_cl + asin_PCTIMP_CATRP100 + L_NABD_NrmStorWs_ratio + asin_PCTFORGRS_CATRP100 + W1_HAG + W1_HNOAG + L_STRM_POWER + SumBig_PFC + INCISION

# Covariance
LOE_QLow_cl~~LOE_Qbkf_cl
L_PTL ~~ L_NTL
L_CHLR ~~ L_SULF
L_PTL ~~ L_TURB

'

##############
## XER EPT
# ROBUST ESTIMATION MAX LIKELIHOOD METHOD
fit_v8_XERw_EPT_robust.est<- sem(mymodel_v8_EPT, data=xer_w,
                                 estimator="MLM")#estimator="MLM")

summary(fit_v8_XERw_EPT_robust.est, standardized=TRUE, fit.measures=TRUE, modindices=F,rsquare=TRUE)#, modindices=T, rsquare=TRUE) #

# request modification indices greater than 3.0 - from Grace USGS materials
mi_min <-modindices(fit_v8_XERw_EPT_robust.est)
print(mi_min[mi_min$mi >3.0,])


######################
######################
## REVISED -
mymodel_v8_EPT_rev <-'
Lpt01_XCMGW ~ 1 + i3*asin_PCTIMP_CATRP100 + ha3*W1_HAG + na4*W1_HNOAG + d2*L_NABD_NrmStorWs_ratio + sp3*L_STRM_POWER
LOE_QLow_cl~ 1 + p6*PSUMPY_SY_WS_sc + d4*L_NABD_NrmStorWs_ratio + ha4*W1_HAG + na5*W1_HNOAG
LOE_Qbkf_cl ~ 1 + p7*PSUMPY_SY_WS_sc + d5*L_NABD_NrmStorWs_ratio + ha5*W1_HAG
evap_index_sc ~ 1 + l4*LOE_QLow_cl + b2*LOE_Qbkf_cl + p5*PSUMPY_SY_WS_sc + d3*L_NABD_NrmStorWs_ratio + sp4*L_STRM_POWER
PCT_SAFN_sc ~ 1 + i2*asin_PCTIMP_CATRP100 + ha2*W1_HAG + rn2*asin_PCTFORGRS_CATRP100 + p2*PSUMPY_SY_WS_sc + sp2*L_STRM_POWER + xi2*INCISION
L_NTL ~ 1 + e3*evap_index_sc + na2*W1_HNOAG + l2*LOE_QLow_cl + p3*PSUMPY_SY_WS_sc
L_SULF ~  1 + e2*evap_index_sc + na3*W1_HNOAG + rn3*asin_PCTFORGRS_CATRP100 + l3*LOE_QLow_cl + p4*PSUMPY_SY_WS_sc
INCISION ~ 1 + i4*asin_PCTIMP_CATRP100 + b3*LOE_Qbkf_cl

EPT_RICH_sc ~ 1 + x1*Lpt01_XCMGW + sf1*PCT_SAFN_sc + n1*L_NTL + su1*L_SULF + d1*L_NABD_NrmStorWs_ratio + xi1*INCISION

# INDIRECT EFFECTS ON EPT
imp_ept:= i2*sf1 + i3*x1 + i4*xi1 + i4*xi2*sf1
dam_ept:= d2*x1 + d3*e2*su1 + d3*e3*n1 + d4*l2*n1 + d4*l3*su1 + d4*l4*e2*su1 + d4*l4*e3*n1 + d5*b2*e2*su1 + d5*b2*e3*n1 + d5*b3*xi1 + d5*b3*xi2*sf1
precip_ept:= p2*sf1 + p3*n1 + p4*su1 + p5*e2*su1 + p5*e3*n1 + p6*l2*n1 + p6*l3*su1 + p6*l4*e2*su1 + p6*l4*e3*n1 + p7*b2*e2*su1 + p7*b2*e3*n1 + p7*b3*xi1 + p7*b3*xi2*sf1
strpwr_ept:= sp2*sf1 + sp3*x1 + sp4*e2*su1 + sp4*e3*n1
hag_ept:= ha2*sf1 + ha3*x1 + ha4*l2*n1 + ha4*l3*su1 + ha4*l4*e2*su1 + ha4*l4*e3*n1 + ha5*b2*e2*su1 + ha5*b2*e3*n1 + ha5*b3*xi1 + ha5*b3*xi2*sf1
nhag_ept:= na2*n1 + na3*su1 + na4*x1 + na5*l2*n1 + na5*l3*su1 + na5*l4*e2*su1 + na5*l4*e3*n1
rnat_ept:= rn2*sf1 + rn3*su1
oelflow_ept:= l2*n1 + l3*su1 + l4*e2*su1 + l4*e3*n1
oebflow_ept:= b2*e2*su1 + b2*e3*n1 + b3*xi1 + b3*xi2*sf1
dexcess_ept:= e2*su1 + e3*n1
#xcmgw_ept:=
xinc_ept:= xi2*sf1

# TOTAL EFFECTS ON EPT
imp_tot:= imp_ept
dam_tot:= d1 + dam_ept
precip_tot:= precip_ept
strpwr_tot:= strpwr_ept
hag_tot:= hag_ept
nhag_tot:= nhag_ept
rnat_tot:= rnat_ept
oelflow_tot:= oelflow_ept
oebflow_tot:= oebflow_ept
dexcess_tot:= dexcess_ept
xcmgw_tot:= x1 #+ xcmgw_ept
xinc_tot:= xi1 + xinc_ept
sand_tot:= sf1
tn_tot:= n1
sulf_tot:= su1

# Covariance
LOE_QLow_cl~~LOE_Qbkf_cl
PCT_SAFN_sc ~~L_NTL
L_NTL ~~L_SULF

'

# ROBUST ESTIMATION MAX LIKELIHOOD METHOD
fit_v8_XERw_EPT_rev_robust.est<- sem(mymodel_v8_EPT_rev, data=xer_w,
                                     estimator="MLM")

summary(fit_v8_XERw_EPT_rev_robust.est, standardized=TRUE, fit.measures=TRUE, modindices=F,rsquare=TRUE)#, modindices=T, rsquare=TRUE) #

# request modification indices greater than 3.0 - from Grace USGS materials
mi_min <-modindices(fit_v8_XERw_EPT_rev_robust.est)
print(mi_min[mi_min$mi >3.0,])

######################
## Bollen.stine bootstrap to estimate parameters -EPT
fit_v8_XERw_EPT_bootstrap_rev  <- sem(mymodel_v8_EPT_rev, data=xer_w,
                                      #group = "ECOREG_rev",
                                      #missing="ML",
                                      test="bollen.stine", se="boot",bootstrap=1000)

summary(fit_v8_XERw_EPT_bootstrap_rev, standardized=TRUE, fit.measures=TRUE, modindices=F,rsquare=TRUE)#, modindices=T, rsquare=TRUE) #

#############
# Export R output - XERIC EPT MODEL
#https://www.r-bloggers.com/export-r-output-to-a-file/
out_fit_v8_XERw_EPT_rev<- capture.output(summary(fit_v8_XERw_EPT_bootstrap_rev, standardized=FALSE, fit.measures=TRUE,rsquare=TRUE)) #, modindices=T
write.csv(out_fit_v8_XERw_EPT_rev, "Routput/Scenario_modeling/SEM_output/XERw_v8_EPT_rev.csv" ,
          row.names=FALSE)

# Standardized estimates of bootstrap model
std_parameter_se_bootstrap_min<- standardizedSolution(fit_v8_XERw_EPT_bootstrap_rev)
write.csv(std_parameter_se_bootstrap_min, "Routput/Scenario_modeling/SEM_output/XERw_v8_EPT_CI.csv" ,
          row.names = FALSE)

# Save model fit parameters
## GRAB R2
r2_XER_EPT<-as.data.frame(lavInspect(fit_v8_XERw_EPT_bootstrap_rev, "rsquare"))

r2_XER_EPT<-cbind(variable=rownames(r2_XER_EPT),r2_XER_EPT)
rownames(r2_XER_EPT) <- 1:nrow(r2_XER_EPT) # To remove index column

colnames(r2_XER_EPT) <- c("Variable","R2")
write.csv(r2_XER_EPT, "Routput/Scenario_modeling/SEM_output/XERw_v8_EPT_R2.csv",
          row.names=FALSE)

#################
# GRAB R2 to be able to add to fit table
r2_red_EPT <- r2_XER_EPT%>%
  filter(Variable=="EPT_RICH_sc")%>%
  select(R2)%>%
  mutate(across(where(is.numeric), round,2))

# TABLE OF FIT INDICES comparing 3 RESPONSES
table_fit_EPT_xer <- matrix(NA, nrow=2, ncol=11)
colnames(table_fit_EPT_xer) = c("Model","Estimation","X2", "df","CFI","TLI", "RMSEA","SRMR","AIC","n","npar")

table_fit_EPT_xer[1,]<-c("EPT","Standard",round(fitmeasures(fit_v8_XERw_EPT_rev_robust.est,
                                                            c("chisq","df","cfi","tli",
                                                              "rmsea","srmr","aic","ntotal","npar")),2))
table_fit_EPT_xer[2,]<-c("EPT","Robust",round(fitmeasures(fit_v8_XERw_EPT_rev_robust.est,
                                                          c("chisq.scaled","df.scaled","cfi.scaled","tli.scaled",
                                                            "rmsea.robust","srmr","aic","ntotal","npar")),2))
# Add R2 for biotic response
table_fit_EPT_xer<- data.frame(table_fit_EPT_xer)%>%
  mutate(R2 = r2_red_EPT$R2)

write.csv(table_fit_EPT_xer, "Routput/Scenario_modeling/SEM_output/XERw_v8_EPT_fit.csv",
          row.names = FALSE)


################
## PATH DIAGRAM USING semPlot
#https://cran.r-project.org/web/packages/tidySEM/vignettes/Plotting_graphs.html
#chrome-extension://efaidnbmnnnibpcajpcglclefindmkaj/https://cran.r-project.org/web/packages/semPlot/semPlot.pdf
EPT_v8<-semPaths(fit_v8_XERw_EPT_rev_robust.est,what="std",intercepts=FALSE,
                 layout = "spring",exoCov = FALSE, rotation = 2,
                 residuals = FALSE) #exoCov = FALSE removes specified covariances

# Relabel nodes
#https://search.r-project.org/CRAN/refmans/semptools/html/change_node_label.html
my_label_list<-list(list(node="EPT",to="EPT"),
                    #list(node="SB_",to="Habitat"),
                    list(node="L_SU",to="SO4"),
                    list(node="L_NT",to="TN"),
                    list(node="PCT",to="Sands"),
                    list(node="e__",to="Evap"),
                    list(node="LOE_Q_",to="Bkf"),
                    list(node="LOE_QL",to="Lowf"),
                    list(node="L_ST",to="Power"),
                    list(node="L01",to="R_Cover"),
                    list(node="a_PCTI",to="Impervious"),
                    list(node="L_NA",to="DamStor"),
                    list(node="a_PCTF",to="R_Forest"),
                    list(node="PSU",to="Precip"),
                    list(node="INC",to="INC"),
                    #list(node="dr_",to="Drought"),
                    list(node="W1_HA",to="R_Ag"),
                    list(node="W1_HN",to="R_NonAg")
)

test<-change_node_label(EPT_v8,my_label_list)
plot(test)

# PRINT semPlot
tiff(filename="C:/Users/EFergus/OneDrive - Environmental Protection Agency (EPA)/a_NLA_OE_project/Project_repository/Routput/Scenario_modeling/Figures/semPlot_EPTv8.tiff",
     width=6.7, height = 9, units="in", res=300)
plot(test)
dev.off()
