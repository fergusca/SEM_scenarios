#' Signal to Noise test SCRIPT FROM KAREN BLOCKSOM EMAIL 2/10/2023
#'
#' @param dfIn This data frame is in wide format with one row per sample and assumed to contain only numeric metrics.
#'      It should also contain all visits to each site, and at least a subset of sites must have multiple visits.
#' @param idVars.samp A character vector containing variables that identify individual samples. Use UID for idVars.samp
#'      - should be unique among all samples in the data
#' @param idVars.site A string containing variable name that identifies sites.
#'       This cannot be the same as or a subset of variables in idVars.samp. Use UNIQUE_ID for idVars.site - if using multiple cycles of data
#' @param year A string containing name of Year variable if sites are revisited across years (as well as within year),
#'      default is 'YEAR'. Set to NULL if no samples across years.
#'
#' @return dataframe of signal to noise ratios
#' @export
#'
#' @examples wmt_SN<-snTest(wmtw_red,idVars.samp="UID",idVars.site="UNIQUE_ID",year='YEAR')

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
