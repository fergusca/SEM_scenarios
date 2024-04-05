
#' Prediction function from de Rooij et al. 2023
#' #https://github.com/mjderooij/SEMpredict/blob/main/predicty.lavaan.R
#' # predict function for computing predicted values for set of response variables
#    given a set of predictor variables
#    based on the joint distribution estimated by a SEM model
#' @param object lavaan output object obtained from sem()
#' @param newdata data frame with new values for the predictors
#' @param xnames variables designated as predictors
#' @param ynames variables designated as response variables
#'
#' @return
#' @export
#'
#' @examples
predicty = function(object, newdata, xnames, ynames){
  Sxx = fitted(object)$cov[xnames , xnames]
  Sxy = fitted(object)$cov[xnames , ynames]
  mx = fitted(object)$mean[xnames]
  my = fitted(object)$mean[ynames]

  #
  Xtest = as.matrix(newdata[, xnames])
  Xtest = scale(Xtest, center = mx, scale = FALSE)
  yhat = matrix(my, nrow = nrow(Xtest), ncol = length(ynames), byrow = TRUE) + Xtest %*% solve(Sxx) %*% Sxy

  # return
  return(yhat)
}
