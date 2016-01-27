/* Copyright (C) 2016 Johannes Helmuth
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
 * @version 0.1.0 1/26/2016
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

// Insert SIMD pragma if supported
#ifdef OMP_VER_4
#define SAFE_SIMD _Pragma("omp simd")
#define SAFE_FOR_SIMD _Pragma("omp for simd")
#else
#define SAFE_SIMD
#define SAFE_FOR_SIMD
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
    #pragma omp parallel for reduction(+:count) num_threads(nthreads)
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
NumericVector logRowSum(const NumericMatrix& mat, int nthreads=1){
    const int rows = mat.nrow();
    const int cols = mat.ncol();
    NumericVector rowSum(rows);
    // iterate over each row: 1. find max value 2. calculate ln sum
    #pragma omp parallel for schedule(static) num_threads(nthreads)
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
        int priv_max;
        #pragma omp for schedule(static) nowait
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
    #pragma omp parallel for reduction(+:sum,comp) num_threads(nthreads)
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
    #pragma omp parallel for reduction(+:sum,comp) num_threads(nthreads)
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
  std::vector<int>& map, std::vector<int>& amount, std::vector<int>& ur, 
  std::vector<int>& us)  {
    if (r.size() <= 0 | s.size() <= 0) {
        return;
    }
    if (r.size() != s.size()) {
        return; //Dunno what's best practise in this case
    }

    //Construct Pair vector
    Pair empty(0,0,0);
    std::vector<Pair> pairs(r.size(), empty);
    for (int i = 0;  i < r.size(); ++i) {
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
    for (int i = 0; i < pairs.size(); ++i) {
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

//' Group unique tupels from two integer vectors r and s
//'
//' @param r An \code{integer()}-vector If elements are not integers, they
//  will be casted to integers.
//' @param s An \code{integer()}-vector. If elements are not integers, they 
//  will be casted to integers.
//' @return a list with the following items:
//'        \item{values}{unique and sorted values of \code{r} and \code{s}}
//'        \item{map}{a vector such that 
//'                   \code{cbind(r,s)[i,] = values[,map[i]]} for every i}
//' @export
//[[Rcpp::export]]
List mapToUniquePairs(const IntegerVector& r, const IntegerVector& s){
     if (r.length() != s.length()) {
        stop("Lengths differ.");
    }
    std::vector<int> map(r.length());
    std::vector<int> amount;
    std::vector<int> ur;
    std::vector<int> us;

    map2uniquePairs_core(asVector(r), asVector(s), map, amount, ur, us);

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

//' Retrieve original vector of values with a supplied map
//'
//' When computing only on unique values, it is often desired to retrieve the
//' original ordering for output. This function does this by creating a new
//' numeric().
//'
//' @param vec The \code{numeric()}-vector to be mapped back
//' @param map The \code{integer()}-vector containg the maps
//' @return a \code{numeric()}-vector of mapped back values
//'
//' @seealso \code{\link{mapToUniquePairs}} for generation of map
//
//' @export
//[[Rcpp::export]]
NumericVector mapToOriginal(const NumericVector& vec, const IntegerVector& map)
{
    NumericVector out(map.length());
    for (int i = 0; i < map.length(); ++i) {
        out[i] = vec[map[i]-1];
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

//' Get normalized enrichment from a diffR fit
//'
//' @param r vector of counts in control/condition1
//' @param s vector of counts in treatment/condition2
//' @param posteriors posterior matrix as computed by diffR
//' @param B column index of background component in posteriors (DEFAULT=1)
//' @return a numeric with enrichment values in log space
//'
//' @export
//[[Rcpp::export]]
NumericVector getEnrichment(const IntegerVector& r, const IntegerVector& s, 
    const NumericMatrix& posteriors, const int B=1, const int nthreads=1) { 
  NumericVector p = posteriors(_,B);
  double p_sum = sum(p);
  double pseu_r = sum(p * r) / p_sum;
  double pseu_s = sum(p * s) / p_sum;
  return log((s + pseu_s)/(r + pseu_r) ) / log(pseu_r / pseu_s)
}

///compute posteriors with a map
Rcpp::NumericVector getEnrichmentWithMap(const Rcpp::NumericMatrix& lnPost,
    const Rcpp::List& m2u, const int B=1) {
  NumericVector lnP = lnPost(_,B) + m2u["amount"];
  double lnP_sum = logSumVector(p);
  double lnPseu_r = logSumVector(p * r) / p_sum;
  double lnPseu_s = logSumVector(p * s) / p_sum;

  return log((s + pseu_s)/(r + pseu_r) ) / log(pseu_r / pseu_s)
}


///EXPECTATION MAXIMIZATION
List em(const List& m2u_sub, const int models=2, const double eps=0.0001, 
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
  double jitter = qstar-.01;
  lntheta = log(rep(qstar, models) + runif(models, -jitter, +jitter));
  lnftheta = log(1 - exp(lntheta));
  if (verbose) {
    message("  ...initiatilizing prior and theta (theta*=" + qstar + ")");
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
        priorstr += exp(lnprior[k]);
        thetastr += exp(lntheta[k]);
      }
      message("  Run " + runs + ": lnL=" + lnLNew + ", step=" + 
          (lnLNew - lnL[lnL.size() - 2]) + ", prior=(" + priorstr + 
          ") theta=(" + thetastr + ")");
    }
  }

  if (verbose) {
    message("  ...converged.");
  }
  return List::create( Named("qstar")=qstar, Named("lnprior")=lnprior,
      Named("lntheta")=lntheta, Named("lnL")=lnL);
}

//TODO write doc
//' Deconvolute bivariate count data in multiple enrichment regimes. Bivariate
//'  data is modeled as a mixture of binomial distributions. Fitting is done
//'  with Expectation Maximization (EM).
//'
//' @param s \code{integer}-vector of counts (e.g. treatment counts in enrichment calls).
//'     If elements are not integers, they will be casted to integers.
//' @param r \code{integer}-vector of counts (e.g. control counts in enrichment calls).
//'     If elements are not integers, they will be casted to integers.
//' @param models \code{integer} specifying number of mixture components which should be
//'     >= 2 (default=2).
//' @param eps \code{double} specifying termination criterion for EM fit (default=0.001).
//' @param verbose \code{logical} specifying if logging should be performed (default=FALSE).
//' @param nthreads \code{integer} specifying number of cores to use (default=1).
//' @return a list with the following items:
//'        \item{qstar}{naive enrichment ratio \code{s/(r + s)}. Basis for EM fit.}
//'        \item{theta}{Parametrization for \code{models} binomial distributiions}
//'        \item{prior}{Mixture proportions for \code{models} binomial distributiions}
//'        \item{posterior}{Posteriormatrix over all bins for \code{models} binomial distributiions}
//'        \item{lnL}{log likelihood trace}
//[[Rcpp::export]]
List normr_core(const IntegerVector& r, const IntegerVector& s, 
    const int models=2, const double eps=.0001, const int iterations=5, 
    const int bgIdx=1, const bool verbose=false, const int nthreads=1) {
    if (models < 2) {
        stop("Error: need at least 2 models (i.e. background and foreground)");
    }
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
    LogicalVector idx = (r > 0 & s > 0);
    if (verbose) {
       message("\t... removing (r == 0) or (s == 0) regions [" +
           logical2Count(idx) + " of " + r.length() +" regions kept].");
    }
    List m2u_sub = mapToUniquePairs(r[idx], s[idx]);
    List fit();
    for (int i = 0; i < iterations; ++i) {
      message("\n*** Iteration " + i + ":");
      List fit_new = em(m2u_sub, models, eps, verbose, nthreads);
      if (fit.size() == 0) {
        fit = fit_new;
      } else {
        if (as<NumericVector>(fit_new["lnL"]) > as<NumericVector>(fit["lnL"])) {
          fit = fit_new;
          message("+++ Iteration " + i + " best so far. Fit updated.")
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

    //posterior matrix
    if (verbose) {
      message("...computing posterior for all data.");
    }
    NumericMatrix lnPost(ur.length(), models);
    NumericVector lnZ(ur.length());
    NumericVector lnprior = as<NumericVector>(fit["lnprior"]);
    NumericVector lntheta = as<NumericVector>(fit["lntheta"]);
    NumericVector lnftheta = log(1 - exp(lntheta));
    calculatePost(lnPost, lnZ, ur, us, lnprior, lntheta, lnftheta, nthreads);

    //calculate enrichment
    if ( verbose ) {
        message("...computing enrichment.");
    }
    std::vector<unsigned int> o = indexSort(as<std::vector<double> >(lntheta));
    NumericMatrix fc( r.length(), models-1 );

    for (int k = 1; k < models; ++k) {
        //foldchange for k against all other components
        if (models == 2) {
            fc(_,0) = lnPost2(_, o[1] ) - lnPost2(_, o[0]);
        } else {
            fc(_,k-1) = lnPost2(_, o[k]) - logSum(subMatrix(lnPost2, k));
        }
    }

    // P value calculation
    //TODO

    //sort data by theta and compute exponential
    std::vector<unsigned int> o = indexSort(as<std::vector<double> >(lntheta));
    NumericMatrix post(r.length(), models);
    NumericVector prior(models), theta(models);
    for (int k = 0; k < models; k++) { //loop over model components
        NumericVector postk = exp(lnPost.column(o[k]));
        post(_,k) = mapToOriginal(postk, umap);
        prior[k] = exp(lnprior(o[k]));
        theta[k] = exp(lntheta(o[k]));
    }

    //output of re-mapped results
    if (verbose) {
        message("[Finished] normR mixture modeling");
    }
    return List::create( Named("qstar")=as<NumericVector>(fit["qstar"]),
        Named("theta")=theta, Named("prior")=prior, Named("posterior")=post,
        Named("lnL")=lnL);
}


/*
 * FOLLOWING ROUTINES HAVE TO BE REFACTORED
 *
//computation of the binomial coefficient
double binomialCoeff(double n, double k) {
    double res = 1;
        // since C(n, k) = C(n, n-k)
        if ( k > n - k ) {
            k = n - k;
        }
    // calculate value of [n * (n-1) *---* (n-k+1)] / [k * (k-1) *----* 1]
    for (int i = 0; i < k; ++i) {
        res *= (n - i);
        res /= (i + 1);
    }
    return res;
}
*/

