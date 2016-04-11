#include <R.h>
#include <Rinternals.h>

#include <R_ext/Rdynload.h>

SEXP normr_logRowSum(SEXP matSEXP, SEXP nthreadsSEXP);
SEXP normr_logSumVector(SEXP vecSEXP, SEXP nthreadsSEXP);
SEXP normr_mapToUniquePairs(SEXP rSEXP, SEXP sSEXP);
SEXP normr_mapToOriginal(SEXP vecSEXP, SEXP m2uSEXP);
SEXP normr_mapToUniqueWithMap(SEXP vecSEXP, SEXP m2uSEXP);
SEXP normr_computeEnrichmentWithMap(SEXP lnPostSEXP, SEXP m2uSEXP, SEXP thetaSEXP, SEXP FSEXP, SEXP BSEXP, SEXP diffCallSEXP, SEXP nthreadsSEXP);
SEXP normr_computeEnrichment(SEXP rSEXP, SEXP sSEXP, SEXP lnPostSEXP, SEXP thetaSEXP, SEXP FSEXP, SEXP BSEXP, SEXP diffCallSEXP, SEXP nthreadsSEXP);
SEXP normr_getLnP(SEXP sSEXP, SEXP rSEXP, SEXP pSEXP, SEXP twoTailedSEXP, SEXP epsSEXP);
SEXP normr_tthreshold(SEXP pSEXP, SEXP epsSEXP, SEXP diffCallSEXP);
SEXP normr_normr_core(SEXP rSEXP, SEXP sSEXP, SEXP modelsSEXP, SEXP epsSEXP, SEXP iterationsSEXP, SEXP bgIdxSEXP, SEXP diffCallSEXP, SEXP verboseSEXP, SEXP nthreadsSEXP);

R_CallMethodDef callMethods[]  = {
  {"normr_logRowSum", (DL_FUNC) &normr_logRowSum, 2},
  {"normr_logSumVector", (DL_FUNC) &normr_logSumVector, 2},
  {"normr_mapToUniquePairs", (DL_FUNC) &normr_mapToUniquePairs, 2},
  {"normr_mapToOriginal", (DL_FUNC) &normr_mapToOriginal, 2},
  {"normr_mapToUniqueWithMap", (DL_FUNC) &normr_mapToUniqueWithMap, 2},
  {"normr_computeEnrichmentWithMap", (DL_FUNC) &normr_computeEnrichmentWithMap, 7},
  {"normr_computeEnrichment", (DL_FUNC) &normr_computeEnrichment, 8},
  {"normr_getLnP", (DL_FUNC) &normr_getLnP, 5},
  {"normr_tthreshold", (DL_FUNC) &normr_tthreshold, 3},
  {"normr_normr_core", (DL_FUNC) &normr_normr_core, 9},
  {NULL, NULL, 0}
};

void R_init_normr(DllInfo *info) {
  R_registerRoutines(info, NULL, callMethods, NULL, NULL);
  R_useDynamicSymbols(info, FALSE);
}

