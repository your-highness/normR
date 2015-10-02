#include <R.h>
#include <Rinternals.h>

#include <R_ext/Rdynload.h>

SEXP diffr_logSum(SEXP matSEXP);
SEXP diffr_logSumVector(SEXP vecSEXP);
SEXP diffr_em(SEXP sSEXP, SEXP rSEXP, SEXP modelsSEXP, SEXP epsSEXP, SEXP verboseSEXP);

R_CallMethodDef callMethods[]  = {
  {"Cdiffr_logSum", (DL_FUNC) &diffr_logSum, 1},
  {"Cdiffr_logSumVector", (DL_FUNC) &diffr_logSumVector, 1},
  {"Cdiffr_em", (DL_FUNC) &diffr_em, 5},
  {NULL, NULL, 0}
};

void R_init_(DllInfo *info) {
  R_registerRoutines(info, NULL, callMethods, NULL, NULL);
  R_useDynamicSymbols(info, FALSE);
}

