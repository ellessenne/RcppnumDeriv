// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// genD_cpp
Rcpp::List genD_cpp(arma::vec const& h0, Rcpp::Function const& func, arma::vec const& x, double const& d, double const& r, double const& v, double const& eps, double const& zero_tol);
RcppExport SEXP _RcppnumDeriv_genD_cpp(SEXP h0SEXP, SEXP funcSEXP, SEXP xSEXP, SEXP dSEXP, SEXP rSEXP, SEXP vSEXP, SEXP epsSEXP, SEXP zero_tolSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec const& >::type h0(h0SEXP);
    Rcpp::traits::input_parameter< Rcpp::Function const& >::type func(funcSEXP);
    Rcpp::traits::input_parameter< arma::vec const& >::type x(xSEXP);
    Rcpp::traits::input_parameter< double const& >::type d(dSEXP);
    Rcpp::traits::input_parameter< double const& >::type r(rSEXP);
    Rcpp::traits::input_parameter< double const& >::type v(vSEXP);
    Rcpp::traits::input_parameter< double const& >::type eps(epsSEXP);
    Rcpp::traits::input_parameter< double const& >::type zero_tol(zero_tolSEXP);
    rcpp_result_gen = Rcpp::wrap(genD_cpp(h0, func, x, d, r, v, eps, zero_tol));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_RcppnumDeriv_genD_cpp", (DL_FUNC) &_RcppnumDeriv_genD_cpp, 8},
    {NULL, NULL, 0}
};

RcppExport void R_init_RcppnumDeriv(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
