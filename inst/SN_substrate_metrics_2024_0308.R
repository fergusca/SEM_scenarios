#########################
## CALCULATING SIGNAL:NOISE
##  Using script from Karen Blocksom email 2/10/23
##
## 3/8/2024
#########################

remove(list=ls())

# Getting error msg 2/9/24 "Error in initializePtr() : \n  function 'cholmod_factor_ldetA' not provided by package 'Matrix'\n"
#https://community.rstudio.com/t/error-in-initializeptr-function-cholmod-factor-ldeta-not-provided-by-package-matrix/178694
tools::package_dependencies("Matrix", which = "LinkingTo", reverse = TRUE)[[1L]]
install.packages("lme4", type = "source")

library(lme4)
library(Hmisc)
library(tidyverse)
#require(plyr)
library(dplyr)
require(gtools)
require(ggplot2)
#require(tidyr)


##########
## READ PROCESSED DATA
# COMPILED NRSA SURVEYS WITH SUBSET OF VARIABLES
#  INCLUDES ALL RESAMPLED SITES AND VISITS 1 & 2

# COMPLETE DATASET WITH ALL THREE SURVEYS COMBINED
#  VISITS 1 and 2 n = 6674
dat_org <-read.csv("data_processed/Compiled/NRSA_081318_all_O_E.csv")
#write.csv(nrsa_oe_all,"data_processed/Compiled/NRSA_081318_all_O_E.csv", row.names=FALSE)

# KEEP BOTH VISITS
dat_proc<- dat_org

###############
## PROCESS DATA DROPPING MISSING PROTOCOL
table(dat_proc$AG_ECO9)
dat_proc$PROTOCOL<-as.factor(dat_proc$PROTOCOL)
summary(dat_proc$PROTOCOL)

# n = 6644
dat_proc<- dat_proc%>%
  drop_na(PROTOCOL)%>%
  filter(PROTOCOL=="BOATABLE"|PROTOCOL=="WADEABLE")

# DROP NOPHAB class from REALM
dat_proc$PROTOCOL<-droplevels(dat_proc$PROTOCOL)
table(dat_proc$PROTOCOL)
#BOATABLE WADEABLE
#    2662     3982

#########################
## SUBSET BY AGGREGATED 9- ECOREGION
##############
## XER
xer<-dat_proc%>%
  filter(AG_ECO9=="XER")%>%
  drop_na(LOE_QLow_cl)

# SCALE CUMULATIVE PRECIPITATION in xer
xer$PSUMPY_SY_WS_sc<-scale(xer$PSUMPY_SY_WS)

# WADEABLE n = 385
xer_w <- xer %>%
  filter(PROTOCOL=="WADEABLE")

##########################
# SCRIPT FROM KAREN BLOCKSOM EMAIL 2/10/2023
# Wrote function to calculate S:N for NARS
# use UNIQUE_ID for idVars.site - if using multiple cycles of data
# use UID for idVars.samp - should be unique among all samples in the data

## Signal-to-noise Test
snTest <- function(dfIn,idVars.samp,idVars.site,year='YEAR'){
  # dfIn - This data frame is in wide format with one row per sample and assumed to contain only numeric metrics.
  #       It should also contain all visits to each site, and at least a subset of sites must have multiple visits.
  #       If there are calibration and validation subsets, this df should only contain calibration samples.
  #       Only identifying variables and metrics should be included in this data frame, and only numeric metrics
  #       can be included.
  #
  # idVars.samp - a character vector containing variables that identify individual samples.
  #
  # idVars.site - a string containing variable name that identifies sites. This cannot be the same as or a subset
  #       of variables in idVars.samp.
  #
  # year - string containing name of Year variable if sites are revisited across years (as well as within year),
  #   default is 'YEAR'. Set to NULL if no samples across years.
  #############################################################################################################################
  options(warn=2)

  # Do some error checking first
  if(idVars.site %in% idVars.samp){
    return(print("idVars.site CANNOT be a part of idVars.samp"))
  }
  # Make sure there are some repeated sites
  if(nrow(dfIn[duplicated(dfIn[,idVars.site]),])<1){
    return(print("You need multiple visits for at least 1 site"))
  }
  print(c('Number of revisits:',nrow(dfIn[duplicated(dfIn[,idVars.site]),])))
  if(!is.null(year)){
    inLong <- tidyr::pivot_longer(dfIn, cols=names(dfIn)[names(dfIn) %nin% c(idVars.samp, idVars.site, year)],
                                  names_to='variable', values_drop_na=TRUE) %>%
      dplyr::filter(!is.infinite(value)) %>%
      mutate(variable=as.character(variable),value=as.numeric(value))
    names(inLong)[names(inLong)==idVars.site] <- 'site'
    names(inLong)[names(inLong)==year] <- 'year'

  }else{
    inLong <- tidyr::pivot_longer(dfIn, cols=names(dfIn)[names(dfIn) %nin% c(idVars.samp, idVars.site)],
                                  names_to='variable', values_drop_na=TRUE) %>%
      dplyr::filter(!is.infinite(value)) %>%
      mutate(variable=as.character(variable),value=as.numeric(value))
    names(inLong)[names(inLong)==idVars.site] <- 'site'

  }

  # create vector of metric names
  parList <- unique(inLong$variable)

  ## Signal-to-noise ratio
  # Create empty data frame to accept output for each metric
  snOut <- data.frame(METRIC=character(),SIGNAL=numeric(),NOISE=numeric(),SN_RATIO=numeric(),COM=character(),stringsAsFactors=FALSE)

  # For each metric in parList, run a linear mixed-effects model with SITE_ID as a random effect
  for(i in 1:length(parList)){
    inMet <- subset(inLong,variable==parList[i])
    # print(parList[i])
    # Run model
    if(!is.null(year)){
      sn <- try(lmer(value~year + (1|year:site),inMet),
                silent=TRUE)

      if(class(sn)=='try-error'){
        sn <- try(lmer(value~year + (1|year:site),inMet, control= lmerControl(optimizer = "bobyqa",
                                                                              optCtrl = list(maxfun = 100000))),
                  silent=TRUE)
      }
    }else{
      sn <- try(lmer(value~(1|site),inMet),silent=TRUE)
    }

    # If model output is error, send to output data frame
    if(class(sn)=='try-error'){
      StoN <- data.frame(METRIC=parList[i],SIGNAL=NA,NOISE=NA,SN_RATIO=NA,COM=sn[1],stringsAsFactors=FALSE)
    }else{
      # If model output value, determine signal and noise by extracting variance components due to SITE_ID and error
      varcomp <- VarCorr(sn)
      if(!is.null(year)){
        StoN <- data.frame(METRIC=parList[i],SIGNAL=round(varcomp$'year:site'[1],2),
                           NOISE=round(attr(varcomp,"sc")^2,2)
                           ,SN_RATIO=round(varcomp$'year:site'[1]/(attr(varcomp,"sc")^2),2),COM=NA,
                           stringsAsFactors=FALSE)
      }else{
        StoN <- data.frame(METRIC=parList[i], SIGNAL= round(varcomp$site[1],2), NOISE=round(attr(varcomp,"sc")^2,2)
                           ,SN_RATIO=round(varcomp$site[1]/(attr(varcomp,"sc")^2),2),COM=NA,stringsAsFactors=FALSE)
      }
    }
    snOut <- rbind(snOut,StoN)
  }
  return(snOut)

}


########################
## SINGAL:NOISE OF BENTHIC INDICES AND METRICS
#########
############
# XERw Subset data - Substrate variables
xerw_red<-xer_w%>%
  select(UID, UNIQUE_ID, YEAR, OE_SCORE,MMI_BENT, EPT_RICH,
         PCT_SAFN, LRBS_use, LOE_RBS_use,L_NTL,L_SULF)

# CALL S:N function
xer_SN<-snTest(xerw_red,idVars.samp="UID",idVars.site="UNIQUE_ID",year='YEAR')
#[1] "Number of revisits:" "107"
xer_SN
#       METRIC SIGNAL  NOISE SN_RATIO COM
#1    OE_SCORE   0.07   0.02     2.93  NA
#2    MMI_BENT 293.08 118.81     2.47  NA
#3    EPT_RICH  23.84   5.88     4.05  NA
#4    LRBS_use   0.74   0.20     3.70  NA
#5    PCT_SAFN 810.63  60.26    13.45  NA
#6 LOE_RBS_use   0.88   0.20     4.37  NA
#7       L_NTL   0.17   0.02     7.12  NA
#8      L_SULF   0.87   0.01    88.81  NA

###################
## CALCULATE MAX r-squared using SITE_ID as random effect (Lester Yuan's approach)
##  From code Ryan shared 2/15/2023
# MaxR2 = SN/(SN+1)
##################
################
## XER wadeable - OE, MMI, EPT
# List of Unique Metric names
unqMetric <- unique(xer_SN$METRIC)

# Create empty dataframe to populate
XERmaxR2.df <- data.frame()

for(i in length(unqMetric)){
  maxR2 <- xer_SN$SN_RATIO/(xer_SN$SN_RATIO+1)
  XERmaxR2.df <- rbind(XERmaxR2.df,
                       data.frame(METRIC=unqMetric, MaxR2=maxR2))
}

XERmaxR2.df

# Create dataframe with S:N and MaxR2 values
XER_SN_maxR2<-xer_SN%>%
  select(METRIC,SN_RATIO)%>%
  left_join(XERmaxR2.df)
XER_SN_maxR2

#############
# WRITE XER S:N and R2 TABLE
write.csv(XER_SN_maxR2, "C:/Users/EFergus/OneDrive - Environmental Protection Agency (EPA)/a_NLA_OE_project/Project_repository/Routput/Scenario_modeling/XERw_SN_MaxR2.csv",
          row.names=FALSE)
