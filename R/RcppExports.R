# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#' @keywords internal
EstimateBetaMDESimple <- function(Y, X, D, b0, iter, critVal, type, HuberC) {
    .Call('_KoulMde_EstimateBetaMDESimple', PACKAGE = 'KoulMde', Y, X, D, b0, iter, critVal, type, HuberC)
}

#' @keywords internal
cppGet_Estimated_Img <- function(zMat, p1, p2) {
    .Call('_KoulMde_cppGet_Estimated_Img', PACKAGE = 'KoulMde', zMat, p1, p2)
}

