// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// gaston_brent
List gaston_brent(NumericVector Y, NumericMatrix X, IntegerVector p_, NumericVector Sigma, NumericMatrix U, double min_h2, double max_h2, double tol, double verbose);
RcppExport SEXP _KAML_gaston_brent(SEXP YSEXP, SEXP XSEXP, SEXP p_SEXP, SEXP SigmaSEXP, SEXP USEXP, SEXP min_h2SEXP, SEXP max_h2SEXP, SEXP tolSEXP, SEXP verboseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type Y(YSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type X(XSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type p_(p_SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type Sigma(SigmaSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type U(USEXP);
    Rcpp::traits::input_parameter< double >::type min_h2(min_h2SEXP);
    Rcpp::traits::input_parameter< double >::type max_h2(max_h2SEXP);
    Rcpp::traits::input_parameter< double >::type tol(tolSEXP);
    Rcpp::traits::input_parameter< double >::type verbose(verboseSEXP);
    rcpp_result_gen = Rcpp::wrap(gaston_brent(Y, X, p_, Sigma, U, min_h2, max_h2, tol, verbose));
    return rcpp_result_gen;
END_RCPP
}
// crossprodcpp
SEXP crossprodcpp(SEXP X);
RcppExport SEXP _KAML_crossprodcpp(SEXP XSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type X(XSEXP);
    rcpp_result_gen = Rcpp::wrap(crossprodcpp(X));
    return rcpp_result_gen;
END_RCPP
}
// geninv
SEXP geninv(SEXP GG);
RcppExport SEXP _KAML_geninv(SEXP GGSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type GG(GGSEXP);
    rcpp_result_gen = Rcpp::wrap(geninv(GG));
    return rcpp_result_gen;
END_RCPP
}
// TransData_c
void TransData_c(std::string bfile, SEXP pBigMat, long maxLine, int threads, bool verbose);
RcppExport SEXP _KAML_TransData_c(SEXP bfileSEXP, SEXP pBigMatSEXP, SEXP maxLineSEXP, SEXP threadsSEXP, SEXP verboseSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type bfile(bfileSEXP);
    Rcpp::traits::input_parameter< SEXP >::type pBigMat(pBigMatSEXP);
    Rcpp::traits::input_parameter< long >::type maxLine(maxLineSEXP);
    Rcpp::traits::input_parameter< int >::type threads(threadsSEXP);
    Rcpp::traits::input_parameter< bool >::type verbose(verboseSEXP);
    TransData_c(bfile, pBigMat, maxLine, threads, verbose);
    return R_NilValue;
END_RCPP
}
// impute_marker
void impute_marker(SEXP pBigMat, int threads, bool verbose);
RcppExport SEXP _KAML_impute_marker(SEXP pBigMatSEXP, SEXP threadsSEXP, SEXP verboseSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type pBigMat(pBigMatSEXP);
    Rcpp::traits::input_parameter< int >::type threads(threadsSEXP);
    Rcpp::traits::input_parameter< bool >::type verbose(verboseSEXP);
    impute_marker(pBigMat, threads, verbose);
    return R_NilValue;
END_RCPP
}
// hasNA
bool hasNA(SEXP pBigMat, const int threads);
RcppExport SEXP _KAML_hasNA(SEXP pBigMatSEXP, SEXP threadsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type pBigMat(pBigMatSEXP);
    Rcpp::traits::input_parameter< const int >::type threads(threadsSEXP);
    rcpp_result_gen = Rcpp::wrap(hasNA(pBigMat, threads));
    return rcpp_result_gen;
END_RCPP
}
// glm_c
SEXP glm_c(const arma::vec& y, const arma::mat& X, SEXP pBigMat, std::string barhead, const bool verbose, const int threads);
RcppExport SEXP _KAML_glm_c(SEXP ySEXP, SEXP XSEXP, SEXP pBigMatSEXP, SEXP barheadSEXP, SEXP verboseSEXP, SEXP threadsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type y(ySEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< SEXP >::type pBigMat(pBigMatSEXP);
    Rcpp::traits::input_parameter< std::string >::type barhead(barheadSEXP);
    Rcpp::traits::input_parameter< const bool >::type verbose(verboseSEXP);
    Rcpp::traits::input_parameter< const int >::type threads(threadsSEXP);
    rcpp_result_gen = Rcpp::wrap(glm_c(y, X, pBigMat, barhead, verbose, threads));
    return rcpp_result_gen;
END_RCPP
}
// mlm_c
SEXP mlm_c(const arma::vec& y, const arma::mat& X, const arma::mat& U, const double vgs, SEXP pBigMat, std::string barhead, const bool verbose, const int threads);
RcppExport SEXP _KAML_mlm_c(SEXP ySEXP, SEXP XSEXP, SEXP USEXP, SEXP vgsSEXP, SEXP pBigMatSEXP, SEXP barheadSEXP, SEXP verboseSEXP, SEXP threadsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type y(ySEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type U(USEXP);
    Rcpp::traits::input_parameter< const double >::type vgs(vgsSEXP);
    Rcpp::traits::input_parameter< SEXP >::type pBigMat(pBigMatSEXP);
    Rcpp::traits::input_parameter< std::string >::type barhead(barheadSEXP);
    Rcpp::traits::input_parameter< const bool >::type verbose(verboseSEXP);
    Rcpp::traits::input_parameter< const int >::type threads(threadsSEXP);
    rcpp_result_gen = Rcpp::wrap(mlm_c(y, X, U, vgs, pBigMat, barhead, verbose, threads));
    return rcpp_result_gen;
END_RCPP
}
// BigStat
List BigStat(SEXP pBigMat, int threads);
RcppExport SEXP _KAML_BigStat(SEXP pBigMatSEXP, SEXP threadsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type pBigMat(pBigMatSEXP);
    Rcpp::traits::input_parameter< int >::type threads(threadsSEXP);
    rcpp_result_gen = Rcpp::wrap(BigStat(pBigMat, threads));
    return rcpp_result_gen;
END_RCPP
}
// kin_cal_m
SEXP kin_cal_m(SEXP pBigMat, const Nullable<double> SUM, bool scale, const Nullable<NumericVector> wt, int threads, std::string barhead, bool verbose);
RcppExport SEXP _KAML_kin_cal_m(SEXP pBigMatSEXP, SEXP SUMSEXP, SEXP scaleSEXP, SEXP wtSEXP, SEXP threadsSEXP, SEXP barheadSEXP, SEXP verboseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type pBigMat(pBigMatSEXP);
    Rcpp::traits::input_parameter< const Nullable<double> >::type SUM(SUMSEXP);
    Rcpp::traits::input_parameter< bool >::type scale(scaleSEXP);
    Rcpp::traits::input_parameter< const Nullable<NumericVector> >::type wt(wtSEXP);
    Rcpp::traits::input_parameter< int >::type threads(threadsSEXP);
    Rcpp::traits::input_parameter< std::string >::type barhead(barheadSEXP);
    Rcpp::traits::input_parameter< bool >::type verbose(verboseSEXP);
    rcpp_result_gen = Rcpp::wrap(kin_cal_m(pBigMat, SUM, scale, wt, threads, barhead, verbose));
    return rcpp_result_gen;
END_RCPP
}
// kin_cal
SEXP kin_cal(SEXP pBigMat, const Nullable<size_t> step0, const Nullable<double> SUM, bool scale, const Nullable<NumericVector> wt, int threads, bool mkl, bool verbose);
RcppExport SEXP _KAML_kin_cal(SEXP pBigMatSEXP, SEXP step0SEXP, SEXP SUMSEXP, SEXP scaleSEXP, SEXP wtSEXP, SEXP threadsSEXP, SEXP mklSEXP, SEXP verboseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type pBigMat(pBigMatSEXP);
    Rcpp::traits::input_parameter< const Nullable<size_t> >::type step0(step0SEXP);
    Rcpp::traits::input_parameter< const Nullable<double> >::type SUM(SUMSEXP);
    Rcpp::traits::input_parameter< bool >::type scale(scaleSEXP);
    Rcpp::traits::input_parameter< const Nullable<NumericVector> >::type wt(wtSEXP);
    Rcpp::traits::input_parameter< int >::type threads(threadsSEXP);
    Rcpp::traits::input_parameter< bool >::type mkl(mklSEXP);
    Rcpp::traits::input_parameter< bool >::type verbose(verboseSEXP);
    rcpp_result_gen = Rcpp::wrap(kin_cal(pBigMat, step0, SUM, scale, wt, threads, mkl, verbose));
    return rcpp_result_gen;
END_RCPP
}
// kin_cal_s
SEXP kin_cal_s(SEXP pBigMat, const Nullable<double> SUM, bool scale, const Nullable<NumericVector> wt, int threads, bool mkl, std::string barhead, bool verbose);
RcppExport SEXP _KAML_kin_cal_s(SEXP pBigMatSEXP, SEXP SUMSEXP, SEXP scaleSEXP, SEXP wtSEXP, SEXP threadsSEXP, SEXP mklSEXP, SEXP barheadSEXP, SEXP verboseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type pBigMat(pBigMatSEXP);
    Rcpp::traits::input_parameter< const Nullable<double> >::type SUM(SUMSEXP);
    Rcpp::traits::input_parameter< bool >::type scale(scaleSEXP);
    Rcpp::traits::input_parameter< const Nullable<NumericVector> >::type wt(wtSEXP);
    Rcpp::traits::input_parameter< int >::type threads(threadsSEXP);
    Rcpp::traits::input_parameter< bool >::type mkl(mklSEXP);
    Rcpp::traits::input_parameter< std::string >::type barhead(barheadSEXP);
    Rcpp::traits::input_parameter< bool >::type verbose(verboseSEXP);
    rcpp_result_gen = Rcpp::wrap(kin_cal_s(pBigMat, SUM, scale, wt, threads, mkl, barhead, verbose));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_KAML_gaston_brent", (DL_FUNC) &_KAML_gaston_brent, 9},
    {"_KAML_crossprodcpp", (DL_FUNC) &_KAML_crossprodcpp, 1},
    {"_KAML_geninv", (DL_FUNC) &_KAML_geninv, 1},
    {"_KAML_TransData_c", (DL_FUNC) &_KAML_TransData_c, 5},
    {"_KAML_impute_marker", (DL_FUNC) &_KAML_impute_marker, 3},
    {"_KAML_hasNA", (DL_FUNC) &_KAML_hasNA, 2},
    {"_KAML_glm_c", (DL_FUNC) &_KAML_glm_c, 6},
    {"_KAML_mlm_c", (DL_FUNC) &_KAML_mlm_c, 8},
    {"_KAML_BigStat", (DL_FUNC) &_KAML_BigStat, 2},
    {"_KAML_kin_cal_m", (DL_FUNC) &_KAML_kin_cal_m, 7},
    {"_KAML_kin_cal", (DL_FUNC) &_KAML_kin_cal, 8},
    {"_KAML_kin_cal_s", (DL_FUNC) &_KAML_kin_cal_s, 8},
    {NULL, NULL, 0}
};

RcppExport void R_init_KAML(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
