#' Process SEM unstandardized parameter output for OE
#'
#' @param parameter_output SEM unstd coefficients saved as .csv file
#'
#' @return Processed dataframe with labeled model coefficients to use for scenario predictions
#' @export
#'
#' @examples
#'west_oe_m1_full<- sem_eff_tab("model1") # name dataframe loaded with raw SEM effects output
#'
unstd_coeff = function(parameter_output){
  param_proc <- parameter_output %>%
    select(lhs, op,rhs,exo,label,est,se)%>%
    filter(exo==0)%>%
    filter(!(op==":="|op=="~~"))%>% # subset data to include endogenous vars
    mutate(label=lhs)%>%
    mutate(label=recode_factor(label,
                               Lpt01_XCMGW = "b1.",
                               L_STRM_POWER = "b2.",
                               LQLow_kmcl = "b3.",
                               LQbkf_kmcl = "b4.",
                               LOE_QLow_cl = "b3.",
                               LOE_Qbkf_cl = "b4.",
                               evap_index_sc = "b5.",
                               LRBS_use = "b6.",
                               L_NTL = "b7.",
                               L_SULF = "b8.",
                               OE_SCORE = "b9.",
                               MMI_BENT_sc = "b9.",
                               EPT_RICH_sc = "b9."))%>%
    group_by(lhs)%>%
    mutate(label2=seq_along(lhs))%>%# Will sequence predictors by endogenous var
    mutate(coeff_name=str_c(label,label2))%>% # concatenate coefficient labels
    select(lhs,rhs,coeff_name,est,se)
}



# USEFUL REFERENCES
# make sequence by group:  https://stackoverflow.com/questions/11996135/create-a-sequential-number-counter-for-rows-within-each-group-of-a-dataframe
# stringr:  https://stringr.tidyverse.org/reference/str_c.html
