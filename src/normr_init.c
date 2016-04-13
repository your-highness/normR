#include <R.h>
#include <Rinternals.h>

#include <R_ext/Rdynload.h>

SEXP normr_logRowSum(SEXP matSEXP, SEXP nthreadsSEXP);
SEXP normr_logSumVector(SEXP vecSEXP, SEXP nthreadsSEXP);
SEXP normr_mapToOriginal(SEXP vecSEXP, SEXP m2uSEXP);
SEXP normr_mapToUniqueWithMap(SEXP vecSEXP, SEXP m2uSEXP);
SEXP normr_computeEnrichment(SEXP rSEXP, SEXP sSEXP, SEXP lnPostSEXP, SEXP thetaSEXP, SEXP FSEXP, SEXP BSEXP, SEXP diffCallSEXP, SEXP nthreadsSEXP);
SEXP normr_normr_core(SEXP rSEXP, SEXP sSEXP, SEXP modelsSEXP, SEXP epsSEXP, SEXP iterationsSEXP, SEXP bgIdxSEXP, SEXP diffCallSEXP, SEXP verboseSEXP, SEXP nthreadsSEXP);

R_CallMethodDef callMethods[]  = {
  {"normr_logRowSum", (DL_FUNC) &normr_logRowSum, 2},
  {"normr_logSumVector", (DL_FUNC) &normr_logSumVector, 2},
  {"normr_mapToOriginal", (DL_FUNC) &normr_mapToOriginal, 2},
  {"normr_mapToUniqueWithMap", (DL_FUNC) &normr_mapToUniqueWithMap, 2},
  {"normr_computeEnrichment", (DL_FUNC) &normr_computeEnrichment, 8},
  {"normr_normr_core", (DL_FUNC) &normr_normr_core, 9},
  {NULL, NULL, 0}
};

void R_init_normr(DllInfo *info) {
  R_registerRoutines(info, NULL, callMethods, NULL, NULL);
  R_useDynamicSymbols(info, FALSE);
}

