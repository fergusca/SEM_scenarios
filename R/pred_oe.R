# Function
#' Cross validation to evaluate predictive performance of Xeric OE candidate models
#'
#' @param index_value An index value to loop through for the k-fold cross validation
#' @param data The dataset to partition into training and test subsets and run the candidate models
#' @param log10 An argument to either log10 transform response or not
#'
#' @return A dataframe with observed O/E and predicted O/E for each of a set of candidate models using an indexed test data
#' @export
#'
#' @examples This function is to be used with purrr map to apply the function across index values
#' #index_1<-pred_oe(index_value = 1, data=xer_org,log10=FALSE)
#' #pred_oe_cv<-function(indexes,...){
#     map(indexes,m_oe,...)
#     }

pred_oe<-function(index_value,data,log10=FALSE) {
  ########
  xer_org<-data

  #################
  # Argument in the function Log transform
  if (log10) {
    xer_org$RESPONSE <- log10(xer_org$OE_SCORE + 0.01)
  } else {
    xer_org$RESPONSE <- xer_org$OE_SCORE
  }

  ###########
  # Create data subsets (test = 1; training=2 thru 5)
  test<-xer_org%>%
    filter(index==index_value)

  train<-xer_org%>%
    filter(index!=index_value)


  ###########
  # Fit candidate models using training data subset

  # REDUCE VARIABLES TO INCLUDE IN MODEL
  myvars<- c("PSUMPY_SY_WS_sc", "drought_mean",
             "asin_PCTIMP_WsRp100","L_NABD_NrmStorWs_ratio",
             "L_STRM_POWER","LOE_QLow_cl","LOE_Qbkf_cl","evap_index_sc",
             "W1_HAG","W1_HNOAG","asin_PCTFORGRS_CATRP100",
             "Lpt01_XCMGW","LRBS_use","L_NTL","L_SULF","RESPONSE")
  train_red<-train%>%
    select(myvars)

  test_red<-test%>%
    select(all_of(myvars))


  # Multiple linear regression (non-spatial)
  mlr.oe <-lm(RESPONSE~., data=train_red)# use all variables in subset except log OE
  summary(mlr.oe)


  # Predict O/E using MLR model
  pred.mlr<-predict(mlr.oe, test_red)

  #################
  # MLR w/interactions
  # Add interaction terms btw impervious surface and TN, SO4, and summer flow
  mlr.int.oe <-lm(RESPONSE~PSUMPY_SY_WS_sc+drought_mean+
                    asin_PCTIMP_WsRp100+L_NABD_NrmStorWs_ratio+
                    L_STRM_POWER+LOE_QLow_cl+LOE_Qbkf_cl+evap_index_sc+
                    W1_HAG+W1_HNOAG+asin_PCTFORGRS_CATRP100+
                    Lpt01_XCMGW+LRBS_use+L_NTL+L_SULF + asin_PCTIMP_WsRp100*LOE_QLow_cl+
                    asin_PCTIMP_WsRp100*L_NTL+asin_PCTIMP_WsRp100*L_SULF , data=train_red)# use all variables in subset except log OE
  summary(mlr.int.oe)


  # Predict O/E using MLR model
  pred.mlr.int<-predict(mlr.int.oe,test_red) #

  #################
  # SPATIAL MODEL
  # Make a simple features object
  xer_sf <- st_as_sf(x=xer_org,
                     coords=c("LON_DD83","LAT_DD83"),
                     crs=4269)%>%
    st_transform(crs=5070)# Albers equal area EPSG = 5070 - need to project to a coordinate system where x and y units are on same scale

  test_sf<-xer_sf%>%
    filter(index==index_value)

  train_sf<-xer_sf%>%
    filter(index!=index_value)

  # Modify the simple feature data so that values are greater than zero
  xer_sf_m<-xer_sf%>%
    mutate(OE_m=ifelse(OE_SCORE==0,0.00001,OE_SCORE))

  test_sf_m<-xer_sf_m%>%
    filter(index==index_value)

  train_sf_m<-xer_sf_m%>%
    filter(index!=index_value)

  ##########
  # Spatial GLM Gaussian distribution and exponential covariance structure
  spmod.oe<-splm(RESPONSE~ PSUMPY_SY_WS_sc+drought_mean+
                   asin_PCTIMP_WsRp100+L_NABD_NrmStorWs_ratio+
                   L_STRM_POWER+LOE_QLow_cl+LOE_Qbkf_cl+evap_index_sc+
                   W1_HAG+W1_HNOAG+asin_PCTFORGRS_CATRP100+
                   Lpt01_XCMGW+LRBS_use+L_NTL+L_SULF,
                 train_sf, spcov_type="exponential")# exponential is default
  ###########
  ## SP model predict
  pred.sp <- predict(spmod.oe, newdata=test_sf)


  ##########
  # Spatial GLM GAMMA Exponential Covariance structure
  spmod.gamma<-spglm(OE_m ~ PSUMPY_SY_WS_sc+drought_mean+
                       asin_PCTIMP_WsRp100+L_NABD_NrmStorWs_ratio+
                       L_STRM_POWER+LOE_QLow_cl+LOE_Qbkf_cl+evap_index_sc+
                       W1_HAG+W1_HNOAG+asin_PCTFORGRS_CATRP100+
                       Lpt01_XCMGW+LRBS_use+L_NTL+L_SULF,
                     family="Gamma",
                     train_sf_m, spcov_type="exponential")# exponential is default
  summary(spmod.gamma)


  ###########
  ## SP model predict with Gamma
  pred.sp.gamma <- predict(spmod.gamma, newdata=test_sf_m, type="response")

  ##########
  # Spatial GLM Inverse Gaussian Exponential Covariance structure
  spmod.ig<-spglm(OE_m ~ PSUMPY_SY_WS_sc+drought_mean+
                    asin_PCTIMP_WsRp100+L_NABD_NrmStorWs_ratio+
                    L_STRM_POWER+LOE_QLow_cl+LOE_Qbkf_cl+evap_index_sc+
                    W1_HAG+W1_HNOAG+asin_PCTFORGRS_CATRP100+
                    Lpt01_XCMGW+LRBS_use+L_NTL+L_SULF,
                  family="inverse.gaussian",
                  train_sf_m, spcov_type="exponential")# exponential is default
  summary(spmod.ig)


  ###########
  ## SP model predict with Gamma
  pred.sp.ig <- predict(spmod.ig, newdata=test_sf_m, type="response")


  ######################
  ## STRUCTURAL EQUATION MODEL PREDICTION
  ## USING lavPredictY
  # PATH ANLAYSIS MODEL v3
  mymodel_v3_OE <- '
Lpt01_XCMGW ~ 1 + ha3*W1_HAG + na2*W1_HNOAG + d2*L_NABD_NrmStorWs_ratio + sp3*L_STRM_POWER + p3*PSUMPY_SY_WS_sc
L_STRM_POWER ~ 1 + d3*L_NABD_NrmStorWs_ratio + ha6*W1_HAG + rn4*asin_PCTFORGRS_CATRP100
LOE_QLow_cl~ 1 + i5*asin_PCTIMP_WsRp100 + p8*PSUMPY_SY_WS_sc + d5*L_NABD_NrmStorWs_ratio + rn5*asin_PCTFORGRS_CATRP100 + ha5*W1_HAG + na4*W1_HNOAG
LOE_Qbkf_cl ~ 1 + p7*PSUMPY_SY_WS_sc + d4*L_NABD_NrmStorWs_ratio + ha4*W1_HAG + ph3*drought_mean
evap_index_sc ~ 1 + l5*LOE_QLow_cl + b3*LOE_Qbkf_cl + p6*PSUMPY_SY_WS_sc + d6*L_NABD_NrmStorWs_ratio + sp4*L_STRM_POWER
LRBS_use ~ 1 + i2*asin_PCTIMP_WsRp100 + l2*LOE_QLow_cl + b2*LOE_Qbkf_cl + ha2*W1_HAG + rn2*asin_PCTFORGRS_CATRP100 + p2*PSUMPY_SY_WS_sc + sp2*L_STRM_POWER
L_NTL ~ 1 + i3*asin_PCTIMP_WsRp100 + x2*Lpt01_XCMGW + l3*LOE_QLow_cl + p4*PSUMPY_SY_WS_sc
L_SULF ~  1 + i4*asin_PCTIMP_WsRp100 + e2*evap_index_sc + na3*W1_HNOAG + rn3*asin_PCTFORGRS_CATRP100 + l4*LOE_QLow_cl + p5*PSUMPY_SY_WS_sc + ph2*drought_mean

RESPONSE ~ 1 + r1*LRBS_use + n1*L_NTL + su1*L_SULF + boe1*LOE_Qbkf_cl + d1*L_NABD_NrmStorWs_ratio + sp1*L_STRM_POWER

# Covariance
LOE_QLow_cl~~LOE_Qbkf_cl
L_STRM_POWER~~LOE_QLow_cl
L_SULF~~L_NTL

'

  # ROBUST ESTIMATION MAX LIKELIHOOD METHOD
  fit<- sem(mymodel_v3_OE, data=train,
            estimator="MLM")

  # PREDICT Using lavPredictY - specify variables to include
  xvars<-c("PSUMPY_SY_WS_sc","drought_mean","LRBS_use","L_NTL","L_SULF","LOE_Qbkf_cl","LOE_QLow_cl","W1_HAG","W1_HNOAG","asin_PCTIMP_WsRp100","L_NABD_NrmStorWs_ratio","L_STRM_POWER","evap_index_sc","asin_PCTFORGRS_CATRP100","Lpt01_XCMGW")

  ###########
  # SEM: PREDICT Y-VALUES
  pred.sem<-lavPredictY(fit,ynames="RESPONSE",xnames=xvars,newdata = test)


  #################
  ## RANDOM FOREST
  set.seed(200)
  rf.oe <- ranger(RESPONSE~., data=train_red,
                  importance = "permutation",
                  scale.permutation.importance = TRUE)

  pred.rf <-predict(rf.oe, data = test_red)


  ############
  ## RANDOM FOREST w/SPATIAL RESIDUALS
  # Spatial model with exponential covariance structure
  rfspmod<-splmRF(RESPONSE~ PSUMPY_SY_WS_sc+drought_mean+
                    asin_PCTIMP_WsRp100+L_NABD_NrmStorWs_ratio+
                    L_STRM_POWER+LOE_QLow_cl+LOE_Qbkf_cl+evap_index_sc+
                    W1_HAG+W1_HNOAG+asin_PCTFORGRS_CATRP100+
                    Lpt01_XCMGW+LRBS_use+L_NTL+L_SULF,
                  train_sf, spcov_type="exponential")# exponential is default

  # Predict using new data
  pred.rfsp <- predict(rfspmod, newdata=test_sf)



  ##########################
  # Create dataframe of predicted and observed
  pred.df<-as.data.frame(cbind(test$RESPONSE,pred.mlr,pred.mlr.int,pred.sem,pred.sp, pred.sp.gamma, pred.sp.ig,pred.rf$predictions,pred.rfsp))%>%
    rename(obs=V1,mlr= pred.mlr,mlr_int=pred.mlr.int,sem=RESPONSE,sp=pred.sp,sp_gamma=pred.sp.gamma,sp_invgaus=pred.sp.ig,rf=V8,rfsp=pred.rfsp)

}
