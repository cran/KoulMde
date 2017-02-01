// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// CLoss
double CLoss(NumericVector Y, NumericMatrix X, NumericMatrix Dstar, NumericVector beta);
RcppExport SEXP KoulMde_CLoss(SEXP YSEXP, SEXP XSEXP, SEXP DstarSEXP, SEXP betaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type Y(YSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type X(XSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type Dstar(DstarSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type beta(betaSEXP);
    rcpp_result_gen = Rcpp::wrap(CLoss(Y, X, Dstar, beta));
    return rcpp_result_gen;
END_RCPP
}
// DY
NumericMatrix DY(NumericVector Y, NumericMatrix D);
RcppExport SEXP KoulMde_DY(SEXP YSEXP, SEXP DSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type Y(YSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type D(DSEXP);
    rcpp_result_gen = Rcpp::wrap(DY(Y, D));
    return rcpp_result_gen;
END_RCPP
}
// GVec
NumericVector GVec(NumericMatrix PMat);
RcppExport SEXP KoulMde_GVec(SEXP PMatSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type PMat(PMatSEXP);
    rcpp_result_gen = Rcpp::wrap(GVec(PMat));
    return rcpp_result_gen;
END_RCPP
}
// HVec
NumericVector HVec(NumericMatrix PMat);
RcppExport SEXP KoulMde_HVec(SEXP PMatSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type PMat(PMatSEXP);
    rcpp_result_gen = Rcpp::wrap(HVec(PMat));
    return rcpp_result_gen;
END_RCPP
}
// GLI
List GLI(NumericVector E, int n, int SN, int EN);
RcppExport SEXP KoulMde_GLI(SEXP ESEXP, SEXP nSEXP, SEXP SNSEXP, SEXP ENSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type E(ESEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< int >::type SN(SNSEXP);
    Rcpp::traits::input_parameter< int >::type EN(ENSEXP);
    rcpp_result_gen = Rcpp::wrap(GLI(E, n, SN, EN));
    return rcpp_result_gen;
END_RCPP
}
// PV
NumericVector PV(NumericVector SE, int nlen);
RcppExport SEXP KoulMde_PV(SEXP SESEXP, SEXP nlenSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type SE(SESEXP);
    Rcpp::traits::input_parameter< int >::type nlen(nlenSEXP);
    rcpp_result_gen = Rcpp::wrap(PV(SE, nlen));
    return rcpp_result_gen;
END_RCPP
}
// GUI
List GUI(NumericVector E, int n, int SN, int EN);
RcppExport SEXP KoulMde_GUI(SEXP ESEXP, SEXP nSEXP, SEXP SNSEXP, SEXP ENSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type E(ESEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< int >::type SN(SNSEXP);
    Rcpp::traits::input_parameter< int >::type EN(ENSEXP);
    rcpp_result_gen = Rcpp::wrap(GUI(E, n, SN, EN));
    return rcpp_result_gen;
END_RCPP
}
// Xpm
NumericMatrix Xpm(NumericMatrix X);
RcppExport SEXP KoulMde_Xpm(SEXP XSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type X(XSEXP);
    rcpp_result_gen = Rcpp::wrap(Xpm(X));
    return rcpp_result_gen;
END_RCPP
}
// PMM
NumericMatrix PMM(int n, NumericMatrix DYM, NumericMatrix XpmM, int l, int p, NumericVector bVec);
RcppExport SEXP KoulMde_PMM(SEXP nSEXP, SEXP DYMSEXP, SEXP XpmMSEXP, SEXP lSEXP, SEXP pSEXP, SEXP bVecSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type DYM(DYMSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type XpmM(XpmMSEXP);
    Rcpp::traits::input_parameter< int >::type l(lSEXP);
    Rcpp::traits::input_parameter< int >::type p(pSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type bVec(bVecSEXP);
    rcpp_result_gen = Rcpp::wrap(PMM(n, DYM, XpmM, l, p, bVec));
    return rcpp_result_gen;
END_RCPP
}
