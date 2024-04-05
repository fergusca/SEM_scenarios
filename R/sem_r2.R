#' Process SEM R2 output for OE,MMI, EPT relabel variable names to make into a table
#'
#' @param R2_output SEM estimated R2 for endogenous variables saved as .csv file
#'
#' @return Processed dataframe that R2 for endogenous variables in the model
#' @export
#'
#' @examples west_oe_m1_full<- r2_tab("model1") # name R2 dataframe processed from SEM output
#'
sem_r2_tab_v2 = function(R2_output){
  R2_out_proc <- R2_output %>%
    select(Variable,R2)%>%
    rename(Response = Variable, R2 = R2)%>%
    mutate(across(where(is.numeric), round,2))%>%
    mutate(Response = recode_factor(Response,
                                    L_STRM_POWER="Stream power",
                                    Lpt01_XCMGW="Riparian cover",
                                    LQLow_kmcl="Summer flow",LQbkf_kmcl="Bankfull flow",
                                    LOE_QLow_cl="OE Summer flow",LOE_Qbkf_cl="OE Bankfull flow",
                                    evap_index_sc= "Evaporation indicator",
                                    Lpt01_XFC_NAT="Instream cover",L_XINC_H="Channel incision",
                                    INCISION="Channel incision",
                                    SumBig_PFC="Habitat richness",
                                    LRBS_use="Bed stability",
                                    LSUB_DMM="Mean diameter",LDCBF_use="Critical diameter",
                                    PCT_SAFN_sc="Sands/fines",
                                    L_NTL="TN", L_SULF="Sulfate",
                                    OE_SCORE="OE", MMI_BENT_sc="MMI", EPT_RICH_sc="EPT"))
}
