// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include "../inst/include/dspCP.h"
#include <RcppArmadillo.h>
#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// draw_indicators_generic
arma::uvec draw_indicators_generic(std::vector<double> res, std::vector<int> r, int n, std::vector<double> w, std::vector<double> m, std::vector<double> s, int J);
RcppExport SEXP _dspCP_draw_indicators_generic(SEXP resSEXP, SEXP rSEXP, SEXP nSEXP, SEXP wSEXP, SEXP mSEXP, SEXP sSEXP, SEXP JSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::vector<double> >::type res(resSEXP);
    Rcpp::traits::input_parameter< std::vector<int> >::type r(rSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< std::vector<double> >::type w(wSEXP);
    Rcpp::traits::input_parameter< std::vector<double> >::type m(mSEXP);
    Rcpp::traits::input_parameter< std::vector<double> >::type s(sSEXP);
    Rcpp::traits::input_parameter< int >::type J(JSEXP);
    rcpp_result_gen = Rcpp::wrap(draw_indicators_generic(res, r, n, w, m, s, J));
    return rcpp_result_gen;
END_RCPP
}
// sample_mat
Eigen::VectorXd sample_mat(std::vector<int> row_ind, std::vector<int> col_ind, std::vector<double> val, int n, int num_entries, Eigen::VectorXd linht, Eigen::VectorXd rd, int D);
RcppExport SEXP _dspCP_sample_mat(SEXP row_indSEXP, SEXP col_indSEXP, SEXP valSEXP, SEXP nSEXP, SEXP num_entriesSEXP, SEXP linhtSEXP, SEXP rdSEXP, SEXP DSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::vector<int> >::type row_ind(row_indSEXP);
    Rcpp::traits::input_parameter< std::vector<int> >::type col_ind(col_indSEXP);
    Rcpp::traits::input_parameter< std::vector<double> >::type val(valSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< int >::type num_entries(num_entriesSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type linht(linhtSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type rd(rdSEXP);
    Rcpp::traits::input_parameter< int >::type D(DSEXP);
    rcpp_result_gen = Rcpp::wrap(sample_mat(row_ind, col_ind, val, n, num_entries, linht, rd, D));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_dspCP_draw_indicators_generic", (DL_FUNC) &_dspCP_draw_indicators_generic, 7},
    {"_dspCP_sample_mat", (DL_FUNC) &_dspCP_sample_mat, 8},
    {NULL, NULL, 0}
};

RcppExport void R_init_dspCP(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}