/* Copyright (C) 2017 Johannes Helmuth & Ho-Ryun Chung
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 */
/**
 * Expectation Maximization Implementation for normr R-package
 * em.cpp
 * Functions defined herein fit a binomial mixture model to count data
 *
 * @author Johannes Helmuth
 * @version 1.3.1 2016-06-28
 */
#include <algorithm>
#include <Rcpp.h>

#ifndef MY_OMP_H_
#define MY_OMP_H_

#ifdef _OPENMP
#include <omp.h>
#if _OPENMP >= 201307
#define OMP_VER_4
#endif
#endif

#define STRINGIFY(a) #a

// Insert SIMD pragma if supported
#ifdef OMP_VER_4
#define SAFE_SIMD(s)  _Pragma(STRINGIFY(omp simd s))
#define SAFE_FOR_SIMD(s)  _Pragma(STRINGIFY(omp for simd s))
#define SAFE_PARALLEL_SIMD(s)  _Pragma(STRINGIFY(omp parallel for simd s))
#else
#define SAFE_SIMD(s)  _Pragma(STRINGIFY(omp s))
#define SAFE_FOR_SIMD(s)  _Pragma(STRINGIFY(omp for s))
#define SAFE_PARALLEL_SIMD(s)  _Pragma(STRINGIFY(omp parallel for s))
#endif
#endif


using namespace Rcpp;


/*
 * USE C++ TYPE FOR ACCESSING R DATA
 */
///convert a Rcpp::IntegerVector to standard vector<int>
std::vector<int> asVector(const IntegerVector& x) {
    return as<std::vector<int> >(x);
}
///convert a Rcpp::NumericVector to standard vector<double>
std::vector<double> asVector(const NumericVector& x) {
    return as<std::vector<double> >(x);
}

/*
 * Propagate logging to R
 */
void message(std::string s) {
  //Rcpp::Function msg("message");
  //msg(s);
}

/*
 * Rcpp helper functions
 */
///convert a Rcpp::LogicalVector to a Rcpp::IntegerVector containing indizes
IntegerVector logical2Int(const LogicalVector& idx) {
    std::vector<int> index(idx.length());//allocate enough space
    int counter = 0;
    for (int i = 0; i < idx.length(); ++i) {
        if (idx[i]) {
            index[counter] = idx[i];
            counter++;
        }
    }
    index.resize(counter+1);
    index.shrink_to_fit();//shrink vector to content (C++11)
    return IntegerVector(index.begin(), index.end());
}
///count number of TRUE elements in a Rcpp::LogicalVector
int logical2Count(const LogicalVector& vec, int nthreads=1) {
    int count = 0;
    SAFE_PARALLEL_SIMD(reduction(+:count) num_threads(nthreads))
    for (int i = 0; i < vec.size(); ++i) {
        if (vec[i] == TRUE) {
            count++;
        }
    }
    return count;
}


/*
 * MATRIX OPERATIONS
 */
///returns a Rcpp::NumericMatrix with <col> columns removed
NumericMatrix subMatrix(const NumericMatrix& mat, int col) {
    NumericMatrix matslice (mat.nrow(), mat.ncol()-1);
    for (int i=0, j=0; i < mat.ncol()-1; i++) {
        if ( i != col ) {
            matslice(_,j) = mat(_, i);
            ++j;
        }
    }
    return(matslice);
}
///computes logarhythmic row sums for a Rcpp::NumericMatrix
// [[Rcpp::export]]
NumericVector logRowSum(const NumericMatrix& mat, int nthreads=1){
    const int rows = mat.nrow();
    const int cols = mat.ncol();
    NumericVector rowSum(rows);
    // iterate over each row: 1. find max value 2. calculate ln sum
    SAFE_PARALLEL_SIMD(schedule(static) num_threads(nthreads))
    for (int i = 0; i < rows; ++i) {
        double max = -DBL_MAX;
        double tmp = 0;
        for (int j = 0; j < cols; ++j) {
            if (mat(i,j) > max) {
                max = mat(i, j);
            }
        }
        for (int j = 0; j < cols; ++j) {
            tmp += exp(mat(i,j) - max);
        }
        rowSum(i) = max + log(tmp);
    }
    return rowSum;
}
///computes maximum value of a Rcpp::NumericVector/Rcpp::NumericMatrix
double max_parallel(const NumericVector& vec, int nthreads=1) {
    const int len = vec.size();
    double max = -DBL_MAX;
    #pragma omp parallel num_threads(nthreads)
    {
        int priv_max = 0;
        SAFE_FOR_SIMD(schedule(static) nowait)
        for (int i = 0; i < len; ++i) {
            if (vec(i) > priv_max) {
                max = vec(i);
            }
        }
        #pragma omp flush(max)
        if (priv_max > max) { //check before entering critical (save time)
        #pragma omp critical
            {
                if (priv_max > max) max = priv_max;
            }
        }
    }
    return max;
}
///computes logarithmic sum for Rcpp::NumericVector
// [[Rcpp::export]]
double logSumVector(const NumericVector& vec, int nthreads=1) {
    const int len = vec.size();
    //better to use unparallized here
    double max = -DBL_MAX;
    for (int i = 0; i < len; ++i) {
        if (vec(i) > max) {
            max = vec(i);
        }
    }

    double sum = 0.0, comp = 0.0;
    SAFE_PARALLEL_SIMD(reduction(+:sum,comp) num_threads(nthreads))
    for (int i = 0; i < len; ++i) {
        double y = exp(vec(i) - max) - comp;
        double t = sum + y;
        comp = (t-sum) - y;
        sum = t;
    }
    return max + log(sum-comp);
}
///computes sum of a Rcpp::NumericVector
double sumVector(const NumericVector& vec, int nthreads=1) {
    const int len = vec.size();

    double sum = 0.0, comp = 0.0;
    SAFE_PARALLEL_SIMD(reduction(+:sum,comp) num_threads(nthreads))
    for (int i = 0; i < len; ++i) {
        double y = vec(i) - comp;
        double t = sum + y;
        comp = (t-sum) - y;
        sum = t;
    }
    return sum;
}

/*
 * GIVE ORDERING OF AN ARRAY
 */
///A naive wrapper for storing values of an vector and doing an index sort
struct CompareIndex {
    const std::vector<double> values;
    CompareIndex(const std::vector<double> vec) : values(vec) {};
    bool operator() (unsigned int l, unsigned int r)  {
        return values[l] < values[r];
    }
};
///retrieve indices of sorted data vector
std::vector<unsigned int> indexSort(std::vector<double> data) {
    std::vector<unsigned int> index(data.size());
    for (unsigned int i = 0; i < index.size(); ++i) index[i] = i;
    //struct takes values and sorts a given vector based on internal values
    CompareIndex compix(data);
    std::sort(index.begin(), index.end(), compix);
    return index;
}

/*
 *
 * MAP2UNIQUEPAIRS ROUTINES
 *
 *
 */
///A naive wrapper for tupel
struct Pair {
    int r;
    int s;
    int pos;
    Pair(int _r, int _s, int _pos): r(_r), s(_s), pos(_pos){}
};
///Comparison of Pair based on underlying tupel
static inline bool pairCompare(const Pair& a, const Pair& b){
    if (a.r < b.r) {
        return true;
    } else if (a.r > b.r) {
        return false;
    } else {
        if (a.s < b.s) {
            return true;
        } else {
            return false;
        }
    }
}
//Extracts unique tupels from (r,s) and keep a map to retrieve original vectors
static inline void map2uniquePairs_core(std::vector<int> r, std::vector<int> s,
  std::vector<int>& map, std::vector<int>& idx, std::vector<int>& amount,
  std::vector<int>& ur, std::vector<int>& us)  {
    if ((r.size() <= 0) | (s.size() <= 0)) {
        return;
    }
    if (r.size() != s.size()) {
        return; //Dunno what's best practise in this case
    }

    //Construct Pair vector
    Pair empty(0,0,0);
    std::vector<Pair> pairs(r.size(), empty);
    for (unsigned int i = 0;  i < r.size(); ++i) {
        pairs[i].r = r[i];
        pairs[i].s = s[i];
        pairs[i].pos = i;
    }
    std::sort(pairs.begin(), pairs.end(), pairCompare);

    //fill in unique values and map
    int lastR = pairs[0].r;
    ur.push_back(lastR);
    int lastS = pairs[0].s;
    us.push_back(lastS);
    amount.push_back(0);
    //map[pairs[0].pos] = 0;
    for (unsigned int i = 0; i < pairs.size(); ++i) {
        Pair& p = pairs[i];
        if (p.r != lastR || p.s != lastS) {
            lastR = p.r;
            lastS = p.s;
            ur.push_back(lastR);
            us.push_back(lastS);
            amount.push_back(1);
        } else {
            amount.back()++;
        }
        map[p.pos] = ur.size() - 1;
    }
}

// Group unique tupels from two integer vectors r and s
//
// @param r An \code{integer()}-vector If elements are not integers, they
//  will be casted to integers.
// @param s An \code{integer()}-vector. If elements are not integers, they
//  will be casted to integers.
// @return a list with the following items:
//        \item{values}{unique and sorted values of \code{r} and \code{s}}
//        \item{map}{a vector such that
//                   \code{cbind(r,s)[i,] = values[,map[i]]} for every i}
List mapToUniquePairs(const IntegerVector& r, const IntegerVector& s){
  if (r.length() != s.length()) {
    stop("Lengths differ.");
  }
  std::vector<int> map(r.length());
  std::vector<int> idx;
  std::vector<int> amount;
  std::vector<int> ur;
  std::vector<int> us;

  map2uniquePairs_core(asVector(r), asVector(s), map, idx, amount, ur, us);

  NumericMatrix values(2,ur.size());
  for (int i = 0; i < values.ncol(); ++i) {
    values(0,i) = ur[i];
    values(1,i) = us[i];
  }
  return List::create(Named("values")=values,
      Named("amount")=NumericVector(amount.begin(), amount.end()),
      Named("map")=IntegerVector(map.begin(), map.end())+1
      );
}

// Retrieve original vector of values with a supplied map
//
// When computing only on unique values, it is often desired to retrieve the
// original ordering for output. This function does this by creating a new
// numeric().
//
// @param vec The \code{numeric()}-vector to be mapped back.
// @param map A map computed by \code{\link{mapToUniquePairs}}.
// @return a \code{numeric()}-vector of mapped back values.
//
// @seealso \code{\link{mapToUniquePairs}} for generation of map
// [[Rcpp::export]]
NumericVector mapToOriginal(const NumericVector& vec, const List& m2u) {
  IntegerVector map = as<IntegerVector>(m2u["map"]);
  NumericVector out(map.size());
  for (int i = 0; i < out.size(); ++i) {
    out[i] = vec[map[i]-1];
  }
  return out;
}

// [[Rcpp::export]]
NumericVector mapToUniqueWithMap(const NumericVector& vec, const List& m2u) {
  int l = as<NumericMatrix>(m2u["values"]).ncol();
  NumericVector out(l);
  std::fill(out.begin(), out.end(), NumericVector::get_na());

  IntegerVector map = as<IntegerVector>(m2u["map"]);
  for (int i = 0; i < map.size(); i++) {
    if (NumericVector::is_na(out[map[i]-1])) out[map[i]-1] = vec[i];
  }

  return out;
}


/*
 * normr_core() HELPERS
 */
///calculates the posteriors based on a normr model fit
static inline void calculatePost(NumericMatrix& lnPost, NumericVector& lnZ,
    const NumericVector& r, const NumericVector& s, const NumericVector&
    lnprior, const NumericVector& lntheta, const NumericVector& lnftheta,
  const int nthreads=1) {
    int models = lnprior.length();
    for (int k = 0; k < models; ++k) {
      lnPost(_,k) = r * lnftheta[k] + s * lntheta[k] + lnprior[k];
    }
    lnZ = logRowSum(lnPost, nthreads);
    for (int k = 0; k < models; ++k) {
      lnPost(_,k) = lnPost(_, k) - lnZ;
    }
}

///EXPECTATION MAXIMIZATION
//[[Rcpp::export]]
List em(const List& m2u_sub, const int models=2, const double eps=1e-5,
    const bool verbose=false, const int nthreads=1) {
  //Get values from mapToUniquePairs structure as NumericVector
  NumericVector ur_sub = as<NumericMatrix>(m2u_sub["values"]).row(0);
  NumericVector us_sub = as<NumericMatrix>(m2u_sub["values"]).row(1);
  NumericVector uamount_sub = as<NumericVector>(m2u_sub["amount"]);
  IntegerVector map_sub = as<IntegerVector>(m2u_sub["map"]);

  /*
   * Initialization of Mixture Parameters
   */
  NumericVector lnprior(models), lntheta(models), lnftheta(models);
  lnprior = runif(models, 0, 1);
  lnprior = log(lnprior / sum(lnprior));
  double qstar = sum(us_sub * uamount_sub) /
    sum(ur_sub * uamount_sub + us_sub * uamount_sub);
  lntheta = log(rep(qstar, models) - runif(models, 0, (qstar - eps)));
  lnftheta = log(1 - exp(lntheta));
  if (verbose) {
    message("  ...initiatilizing prior and theta (theta*=" +
        std::to_string(qstar) + ")");
  }

  //working objects
  int lsub = ur_sub.size();
  NumericVector lnL(1, -DBL_MAX), lnZ(lsub);
  NumericMatrix lnPost(lsub, models);
  NumericVector lnPost_us_sub_log(lsub), lnPost_un_sub_log(lsub);
  //Compute log counts, s.t. we don't have to do it again all the time
  NumericVector ur_sub_log = log(ur_sub);
  NumericVector us_sub_log = log(us_sub);
  NumericVector un_sub_log = log(ur_sub + us_sub);
  NumericVector uamount_sub_log = log(uamount_sub);

  //loop
  int runs = 0;
  bool notConverged = true;
  double lnPostSum = 0;
  while ((runs < 25) | notConverged) { // run at least 25 times for burn in
    /*
     * Expectation. Works only on unique values and uses indexes of logRowSum
     * result to restore original matrix
     */
    calculatePost(lnPost, lnZ, ur_sub, us_sub, lnprior, lntheta, lnftheta,
        nthreads);

    /*
     * Maximization (done on map reduced data)
     */
    for (int k = 0; k < models; ++k) {
      //add quantities of each i from mapping (to get overall model posterior)
      lnPost(_,k) = lnPost(_, k) + uamount_sub_log;
    }
    //Calculation of likelihood and model parameters needs amount of unqiues
    lnPostSum = logSumVector(lnPost, nthreads);
    for (int k = 0; k < models; ++k) {
      lnPost_us_sub_log = lnPost(_, k) + us_sub_log;
      lnPost_un_sub_log = lnPost(_, k) + un_sub_log;
      lntheta[k] = logSumVector(lnPost_us_sub_log, nthreads)
                   - logSumVector(lnPost_un_sub_log, nthreads);
      lnftheta[k] = log(1 - exp(lntheta[k]));
      lnprior[k] = logSumVector(lnPost(_, k), nthreads) - lnPostSum;
    }

    /*
     * Convergence
     */
    lnZ = lnZ * uamount_sub;
    double lnLNew = sum(lnZ); //multiply by amount of unique
    if (runs > 25 && lnLNew - lnL[lnL.size() - 1] < eps) {
      notConverged = false;
    }
    lnL.push_back(lnLNew);

    /*
     * Logging
     */
    runs++;
    if (verbose && (!notConverged || (runs % 10) == 0)) {
      std::string priorstr = ""; std::string thetastr = "";
      for (int k = 0; k < models; k++) { //loop over model components
        priorstr += " "+std::to_string(exp(lnprior[k]));
        thetastr += " "+std::to_string(exp(lntheta[k]));
      }
      message("  Run " + std::to_string(runs) + ": lnL=" +
          std::to_string(lnLNew) + ", step=" +
          std::to_string(lnLNew - lnL[lnL.size() - 2]) + ", prior=(" +
          priorstr + ") theta=(" + thetastr + ")");
    }
  }

  if (verbose) {
    message("  ...converged.");
  }
  return List::create( Named("qstar")=qstar, Named("lnprior")=lnprior,
      Named("lntheta")=lntheta, Named("lnL")=lnL);
}

///compute posteriors with a map and a log posterior matrix on the unique
//values
// [[Rcpp::export]]
NumericVector computeEnrichmentWithMap(const NumericMatrix& lnPost,
     const List& m2u, const NumericVector& theta, const int fg=1,
     const int bg=0, const bool diffCall=false, const bool standardized=true,
     const int nthreads=1) {
  if (bg < 0 || bg >= lnPost.ncol() || bg >= theta.size()) {
    stop("invalid bg argument");
  }
  if (lnPost.ncol() != theta.size()) stop("lnPost and theta not matching");

  NumericVector ur_log = log(as<NumericMatrix>(m2u["values"]).row(0));
  NumericVector us_log = log(as<NumericMatrix>(m2u["values"]).row(1));
  NumericVector uamount_log = log(as<NumericVector>(m2u["amount"]));

  if (lnPost.nrow() != ur_log.size()) {
    stop("lnPost and m2u do not match in size");
  }

  //calculate pseudo-counts from fit bg
  NumericVector lnP(lnPost.nrow());
  lnP = lnPost(_,bg) + uamount_log;
  double lnP_sum = logSumVector(lnP, nthreads);
  lnP = lnPost(_,bg) + uamount_log + ur_log;
  double lnPseu_r = logSumVector(lnP, nthreads) - lnP_sum;
  lnP = lnPost(_,bg) + uamount_log + us_log;
  double lnPseu_s = logSumVector(lnP, nthreads) - lnP_sum;
  double rglrz = lnPseu_r - lnPseu_s;

  //compute regularized, standardized fold change
  NumericVector out(ur_log.size());
  SAFE_PARALLEL_SIMD(schedule(static) num_threads(nthreads))
  for (int i = 0; i < out.size(); ++i) {
    ur_log[i] = exp(ur_log[i]) + exp(lnPseu_r);
    us_log[i] = exp(us_log[i]) + exp(lnPseu_s);
    out[i] = log(us_log[i]/ur_log[i]) + rglrz;
  }
  if (standardized) {
    if (diffCall) {//standardization dependent on algebraic sign of fc
      double stdrzC = -log(theta[0]/(1-theta[0])*(1-theta[bg])/theta[bg]);
      double stdrzT = log(theta[2]/(1-theta[2])*(1-theta[bg])/theta[bg]);
      SAFE_PARALLEL_SIMD(schedule(static) num_threads(nthreads))
      for (int i = 0; i < out.size(); ++i) {
        if (out[i] < 0) {
          out[i] /= stdrzC;
        } else {
          out[i] /= stdrzT;
         }
      }
    } else {
      double stdrz = log(theta[fg]/(1-theta[fg])*(1-theta[bg])/theta[bg]);
      SAFE_PARALLEL_SIMD(schedule(static) num_threads(nthreads))
      for (int i = 0; i < out.size(); ++i) {
        out[i] /= stdrz;
      }
    }
  }

  return out;
}

double getLnP(const int s, const int r, const double p,
    const bool twoTailed=false, const double eps = .0000001) {
  int n = r+s;
  if (twoTailed) {
    double m = n * p;
    double off = m-(double)s;
    if (std::abs(off) < eps) return 0;
    double d = R::dbinom(s, n, p, 1);
    int y = 0;
    double err = 1 + eps;
    if (off > eps) {
      for (int i=ceil(m); i <= n; ++i) {
        if (R::dbinom(i, n, p, 1) <= d) y+=err;
      }
      double up = R::pbinom(n-y, n, p, 0, 1);
      double dw = R::pbinom(s, n, p, 1, 1);
      if (up < -DBL_MAX) return dw;
      else return up + log(1 + exp(dw - up));
    } else {
      for (int i=0; i <= floor(m); ++i) {
        if (R::dbinom(i, n, p, 1) <= d) y+=err;
      }
      double up = R::pbinom(s-1, n, p, 0, 1);
      return up + log(1 + exp(R::pbinom(y-1, n, p, 1, 1) - up));
    }
  } else {
    if (s == 0) return 0;
    return R::pbinom(s-1, n, p, 0, 1);
  }
}

NumericVector getPWithMap(const List& m2u, const double theta,
    const bool diffCall=false, const int nthreads=1) {
  NumericVector ur = as<NumericMatrix>(m2u["values"]).row(0);
  NumericVector us = as<NumericMatrix>(m2u["values"]).row(1);

  NumericVector out(ur.size());
  SAFE_PARALLEL_SIMD(schedule(static) num_threads(nthreads))
  for (int i = 0; i < out.size(); ++i) {
    out[i] = getLnP(us[i], ur[i], theta, diffCall);
  }
  return out;
}

//a T filter implementation: What is the margin were significance can be
//achieved?
int tthreshold(const double p, const double minPThresh=0.05,
    const bool diffCall=false) {
  if (p < 0 || p > 1) stop("invalid p");

  int marg = 0;
  double thresh = log(minPThresh);
  bool run = true;
  while (run) {
    if (getLnP(0, marg, p, diffCall) <= thresh) break;
    if ((marg-1) > 0) {
      for (int i=(marg-1); i >= 1; --i) {
        if (getLnP((marg-i), i, p, diffCall) <= thresh) {
          run = false;
          break;
        }
        if ((marg-i) != i) {
          if(getLnP(i, (marg-i), p, diffCall) <= thresh) {
            run = false;
            break;
          }
        }
      }
    }
    marg++;
  }

  return marg;
}

//gives indices in R notation, i.e. >=1
IntegerVector filterIdx(const List& m2u, const double theta,
    int& threshold, const double minPThresh=0.05, const bool diffCall=false) {
  if (theta < 0 || theta > 1) stop("invalid theta");
  if (minPThresh < 0 || minPThresh > 1) stop("invalid minPThresh");

  threshold = tthreshold(theta, minPThresh, diffCall);

  NumericVector n = as<NumericMatrix>(m2u["values"]).row(0) +
    as<NumericMatrix>(m2u["values"]).row(1);
  std::vector<int> idx;
  for (int i=0; i < n.size(); ++i) {
    if (n[i] >= threshold) {
      idx.push_back(i+1);
    }
  }
  return as<IntegerVector>(wrap(idx));
}

// Deconvolute bivariate count data in multiple enrichment regimes. Bivariate
// data is modeled as a mixture of binomial distributions. Fitting is done
// with Expectation Maximization (EM) on data points were \code{r > 0 & s > 0}.
//
// In a first step, a map of unique non-zero (r,s) values is generated. This
// allows for faster runtime. Second, EM is run with the given number of
// components on these reduced data representation. If \code{iterations > 1},
// the fit with the highest likelihood is selected. In a third step, a map for
// all (r,s) values is generated. Based on this map the complete posterior
// matrix, P-values and enrichment is calculated. The returned values are
// holding results exclusively for unique (r,s) values and go directly into a
// \code{\link{NormRFit-class}} object.
//
// @param r \code{integer}-vector of counts (e.g. control counts in enrichment
// calls). If elements are not integers, they will be casted to integers.
// @param s \code{integer}-vector of counts (e.g. treatment counts in
// enrichment calls). If elements are not integers, they will be casted to
// integers.
// @param models \code{integer} specifying number of mixture components which
// should be >= 2 (default=2).
// @param eps \code{double} specifying termination criterion for EM fit
// (default=0.001).
// @param iterations \code{integer} specifying the number of individual EM runs
// with differential initial parameters to be done. Adjust this argument to
// find global maxima (default=5).
// @param minPThresh \code{double} specifying the minimum achieveable P-value
// for the T method (default=0.05).
// @param bgIdx \code{integer} giving the index of the background component. In
// enrichment and regime calls this should be 0. In difference calls, this
// value can be > 0 (default=0).
// @param diffCall \code{logical} specifying if difference calling is done such
// that a two-sided significance test will be conducted (default=FALSE).
// @param verbose \code{logical} specifying if logging should be performed
// (default=FALSE).
// @param nthreads \code{integer} specifying number of cores to use
// (default=1).
// @return a list with the following items:
//  \item{qstar}{naive enrichment ratio \code{s/(r + s)}. Basis for EM fit.}
//  \item{map}{a map of unique (r,s) values. See
//   \code{map2uniquePairs()}}
//  \item{lntheta}{ln parametrization of mixture binomials}
//  \item{lnprior}{ln mixture proportions of mixture binomials}
//  \item{lnL}{log likelihood trace of EM}
//  \item{lnposterior}{ln posteriors for unique (r,s) according to map}
//  \item{lnenrichment}{ln enrichment for unique (r,s) according to map}
//  \item{lnpvals}{ln P-values for mixture component \code{bgIdx} for each
//   unique (r,s)}
//  \item{filtered}{unique (r,s) tupels passing the T filter with \code{eps}}
// [[Rcpp::export]]
List normr_core(const IntegerVector& r, const IntegerVector& s,
    const int models=2, const double eps=1e-5, const int iterations=5,
    const double minPThresh = 0.05, const int bgIdx=0,
    const bool diffCall=false, const bool verbose=false,
    const int nthreads=1, const String binFilter="zeroSum") {
  if (models < 2) {
      stop("Error: need at least 2 models (i.e. background and foreground)");
  }
  if (eps < 0 || eps > 1) stop("invalid eps specified");

  if (verbose) {
      message("[Started] normR mixture modeling");
  }

  /*
   * Expectation Maximization on (r,s) with non-zero entries
   *
   * Subset data for excluding zero regions. Uses Rcpp sugar. Therefore ursub
   * and ssub have to be NumericVectors. Recreating a new map for subset is
   * much faster than reorganizing the old one.
   */
  LogicalVector idx(r.size(), true);
  if (binFilter == "zero")
    idx = ((r > 0) & (s > 0));
  else if (binFilter == "zeroSum") {
    idx = ((r + s) > 0);
  }
  if (verbose) {
    message("\t... removing (r == 0) or (s == 0) regions [" +
      std::to_string(logical2Count(idx)) + " of " +
      std::to_string(r.length()) +" regions kept].");
  }
  List m2u_sub = mapToUniquePairs(r[idx], s[idx]);
  List fit;
  for (int i = 0; i < iterations; ++i) {
    if (verbose) message("\n*** Iteration " + std::to_string(i+1) + ":");
    List fit_new = em(m2u_sub, models, eps, verbose, nthreads);
    if (fit.size() == 0) {
      fit = fit_new;
    } else {
      NumericVector lnLOld = as<NumericVector>(fit["lnL"]);
      NumericVector lnLNew = as<NumericVector>(fit_new["lnL"]);
      if (lnLNew[lnLNew.size()-1] > lnLOld[lnLOld.size()-1]) {
        fit = fit_new;
        if (verbose) message("+++ Iteration " + std::to_string(i+1) +
            " best so far. Fit updated.");
      }
    }
  }

  /*
   * Map2UniquePairs on complete data for model calculations.
   */
  List m2u = mapToUniquePairs(r, s);
  NumericVector ur = as<NumericMatrix>(m2u["values"]).row(0);
  NumericVector us = as<NumericMatrix>(m2u["values"]).row(1);
  IntegerVector umap = as<IntegerVector>(m2u["map"]);

  //sort theta and calculate posterior for all data
  if (verbose) message("...computing posterior for all data.");
  std::vector<unsigned int> o =
    indexSort(as<std::vector<double> >(as<NumericVector>(fit["lntheta"])));
  NumericVector lnprior(models), lntheta(models);
  for (int k = 0; k < models; k++) { //loop over model components
      lnprior[k] = as<NumericVector>(fit["lnprior"])(o[k]);
      lntheta[k] = as<NumericVector>(fit["lntheta"])(o[k]);
  }
  NumericVector lnftheta = log(1 - exp(lntheta));
  NumericMatrix lnPost(ur.length(), models);
  NumericVector lnZ(ur.length());
  calculatePost(lnPost, lnZ, ur, us, lnprior, lntheta, lnftheta, nthreads);

  //calculate enrichment on map
  if (verbose) message("...computing enrichment.");
  NumericVector enr = computeEnrichmentWithMap(lnPost, m2u, exp(lntheta),
      (models-1), bgIdx, diffCall, true, nthreads);

  //calculate p-values on map
  if (verbose) message("...computing P-values.");
  NumericVector pvals = getPWithMap(m2u, exp(lntheta[bgIdx]), diffCall,
      nthreads);

  //indizes of p-values passing T filter
  if (verbose) {
    message("...applying T filter with threshold eps=" + std::to_string(eps) +
        ".");
  }
  int threshold = 0;
  IntegerVector filteredT = filterIdx(m2u, exp(lntheta[bgIdx]), threshold,
    minPThresh, diffCall);

  if (verbose) message("[Finished] normR mixture modeling");
  return List::create(Named("qstar")=fit["qstar"], Named("map")=m2u,
      Named("lntheta")=lntheta, Named("lnprior")=lnprior,
      Named("lnposterior")=lnPost, Named("lnL")=as<NumericVector>(fit["lnL"]),
      Named("lnenrichment")=enr, Named("lnpvals")=pvals,
      Named("filteredT")=filteredT, Named("Tthreshold")=threshold);
}
