#include <R.h>
#include <Rinternals.h>

#include <R_ext/Rdynload.h>

SEXP normr_logRowSum(SEXP matSEXP, SEXP nthreadsSEXP);
SEXP normr_logSumVector(SEXP vecSEXP, SEXP nthreadsSEXP);
SEXP normr_mapToOriginal(SEXP vecSEXP, SEXP m2uSEXP);
SEXP normr_mapToUniqueWithMap(SEXP vecSEXP, SEXP m2uSEXP);
SEXP normr_em(SEXP m2u_subSEXP, SEXP modelsSEXP, SEXP epsSEXP, SEXP verboseSEXP, SEXP nthreadsSEXP);
SEXP normr_computeEnrichmentWithMap(SEXP lnPostSEXP, SEXP m2uSEXP, SEXP thetaSEXP, SEXP fgSEXP, SEXP bgSEXP, SEXP diffCallSEXP, SEXP standardizedSEXP, SEXP nthreadsSEXP);
SEXP normr_normr_core(SEXP rSEXP, SEXP sSEXP, SEXP modelsSEXP, SEXP epsSEXP, SEXP iterationsSEXP, SEXP minPThreshSEXP, SEXP bgIdxSEXP, SEXP diffCallSEXP, SEXP verboseSEXP, SEXP nthreadsSEXP, SEXP binFilterSEXP);

R_CallMethodDef callMethods[]  = {
  {"normr_logRowSum", (DL_FUNC) &normr_logRowSum, 2},
  {"normr_logSumVector", (DL_FUNC) &normr_logSumVector, 2},
  {"normr_mapToOriginal", (DL_FUNC) &normr_mapToOriginal, 2},
  {"normr_mapToUniqueWithMap", (DL_FUNC) &normr_mapToUniqueWithMap, 2},
  {"normr_em", (DL_FUNC) &normr_em, 5},
  {"normr_computeEnrichmentWithMap", (DL_FUNC) &normr_computeEnrichmentWithMap, 8},
  {"normr_normr_core", (DL_FUNC) &normr_normr_core, 11},
  {NULL, NULL, 0}
};

void R_init_normr(DllInfo *info) {
  R_registerRoutines(info, NULL, callMethods, NULL, NULL);
  R_useDynamicSymbols(info, FALSE);
}

