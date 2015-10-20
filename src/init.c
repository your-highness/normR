#include <R.h>
#include <Rinternals.h>

#include <R_ext/Rdynload.h>

SEXP diffr_logical2Int(SEXP idxSEXP);
SEXP diffr_logical2Count(SEXP vecSEXP, SEXP nthreadsSEXP);
SEXP diffr_logRowSum(SEXP matSEXP, SEXP nthreadsSEXP);
SEXP diffr_logSumVector(SEXP vecSEXP, SEXP nthreadsSEXP);
SEXP diffr_mapToUniquePairs(SEXP rSEXP, SEXP sSEXP);
SEXP diffr_mapToOriginal(SEXP vecSEXP, SEXP mapSEXP);
SEXP diffr_em(SEXP rSEXP, SEXP sSEXP, SEXP modelsSEXP, SEXP epsSEXP, SEXP verboseSEXP, SEXP nthreadsSEXP);

R_CallMethodDef callMethods[]  = {
  {"Cdiffr_logical2Int", (DL_FUNC) &diffr_logical2Int, 1},
  {"Cdiffr_logical2Count", (DL_FUNC) &diffr_logical2Count, 2},
  {"Cdiffr_logRowSum", (DL_FUNC) &diffr_logRowSum, 2},
  {"Cdiffr_logSumVector", (DL_FUNC) &diffr_logSumVector, 2},
  {"Cdiffr_mapToUniquePairs", (DL_FUNC) &diffr_mapToUniquePairs, 2},
  {"Cdiffr_mapToOriginal", (DL_FUNC) &diffr_mapToOriginal, 2},
  {"Cdiffr_em", (DL_FUNC) &diffr_em, 6},
  {NULL, NULL, 0}
};

void R_init_(DllInfo *info) {
  R_registerRoutines(info, NULL, callMethods, NULL, NULL);
  R_useDynamicSymbols(info, FALSE);
}

