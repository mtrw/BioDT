// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// INTERNAL_varTableFromAlns
DataFrame INTERNAL_varTableFromAlns(DataFrame aln);
RcppExport SEXP _BioDT_INTERNAL_varTableFromAlns(SEXP alnSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< DataFrame >::type aln(alnSEXP);
    rcpp_result_gen = Rcpp::wrap(INTERNAL_varTableFromAlns(aln));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_BioDT_INTERNAL_varTableFromAlns", (DL_FUNC) &_BioDT_INTERNAL_varTableFromAlns, 1},
    {NULL, NULL, 0}
};

RcppExport void R_init_BioDT(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
