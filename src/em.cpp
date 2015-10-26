/* Copyright (C) 2015 Johannes Helmuth
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

#include <algorithm>
#include <Rcpp.h>

/*
 * USE C++ TYPE FOR ACCESSING R DATA
 */
std::vector<int> asVector(const Rcpp::IntegerVector& x) {
	return Rcpp::as<std::vector<int> >(x);
}
std::vector<double> asVector(const Rcpp::NumericVector& x) {
	return Rcpp::as<std::vector<double> >(x);
}


/*
 * Rcpp helper functions
 */
//[[Rcpp::export]]
Rcpp::IntegerVector logical2Int(const Rcpp::LogicalVector& idx) {
	std::vector<int> index(idx.length());//allocate enough space
	for (int i = 0; i < idx.length(); ++i) {
		if (idx[i])
			index.push_back(idx[i]);
	}
	index.shrink_to_fit();//shrink vector to content (C++11)
	return Rcpp::IntegerVector(index.begin(), index.end());
}
//[[Rcpp::export]]
int logical2Count(const Rcpp::LogicalVector& vec, int nthreads=1) {
	int count = 0;
	#pragma omp parallel for reduction(+:count) num_threads(nthreads)
	for (int i = 0; i < vec.size(); ++i) {
		if (vec[i] == TRUE) {
			count++;
		}
	}
	return count;
}


//CORE candidates
/*
 * MATRIX OPERATIONS
 */
//computes logSum for a Rcpp::NumericMatrix
//[[Rcpp::export]]
Rcpp::NumericVector logRowSum(const Rcpp::NumericMatrix& mat, int nthreads=1){
	const int rows = mat.nrow();
	const int cols = mat.ncol();
	Rcpp::NumericVector rowSum(rows);
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
//computes logSum for Rcpp::NumericVector
//implements Kahan summation algorithm
//[[Rcpp::export]]
double logSumVector(const Rcpp::NumericVector& vec, int nthreads=1) {
	const int len = vec.size();
	double max = -DBL_MAX;
	for (int i = 0; i < len; ++i) {
		if (vec(i) > max) {
			max = vec(i);
		}
	}

	//bug: per added core results differs by -4.050094e-13
	/*
	 * Implementation of Kahan summation algorithm
	 *
	 * The precison is not full but close. Otherwise use http://code.activestate.com/recipes/393090/
	 */
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
//returns a Rcpp::NumericMatrix with col columns removed
Rcpp::NumericMatrix subMatrix(const Rcpp::NumericMatrix& mat, int col) {
	Rcpp::NumericMatrix matslice (mat.nrow(), mat.ncol()-1);
	for (int i=0, j=0; i < mat.ncol()-1; i++) {
		if ( i != col ) {
			matslice(Rcpp::_,j) = mat(Rcpp::_, i);
			++j;
		}
	}
	return(matslice);
}


/*
 * GIVE ORDERING OF AN ARRAY
 */
struct CompareIndex {
	const std::vector<double> values;
	CompareIndex(const std::vector<double> vec) : values(vec) {};
	bool operator() (unsigned int l, unsigned int r)  {
		return values[l] < values[r];
	} 
};
std::vector<unsigned int> indexSort(std::vector<double> data) {
	std::vector<unsigned int> index(data.size());
	for (unsigned int i = 0; i < index.size(); ++i) index[i] = i;
	// struct that takes values and sorts a given vector based on internal values
	CompareIndex compix(data);
	std::sort(index.begin(), index.end(), compix);
	return index;
}


/*
 * MAP2UNIQUEPAIRS ROUTINES
 */
//TODO use pointers to original R workspace data here
struct Pair {
	int r;
	int s;
	int pos;
	Pair(int _r, int _s, int _pos): r(_r), s(_s), pos(_pos){}
};
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
static inline void map2uniquePairs_core(std::vector<int> r, std::vector<int> s, std::vector<int>& map, std::vector<int>& amount, std::vector<int>& ur, std::vector<int>& us)  {
	if (r.size() <= 0 | s.size() <= 0) {
		return;
	}

	if (r.size() != s.size()) {
		return; //Dunno what's best practise in this case
	}

	Pair empty(0,0,0);
	std::vector<Pair> pairs(r.size(), empty);
	for (int i = 0;  i < r.size(); ++i) {
		pairs[i].r = r[i];
		pairs[i].s = s[i];
		pairs[i].pos = i;
	}

	//potentially parallizeable (TODO)
	std::sort(pairs.begin(), pairs.end(), pairCompare);

	//fill in uniquevalues and map
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
//' Group unique pairs from two vectors
//'
//' @param r first vector of integers. If elements are not integers, they will be
//'     casted to integers.
//' @param s second vector of integers. If elements are not integers, they will be
//'     casted to integers.
//' @return a list with the following items:
//'        \item{values}{unique and sorted values of \code{r} and \code{s}}
//'        \item{map}{a vector such that \code{cbind(r,s)[i,] = values[,map[i]]} for every i}
//' @export
//[[Rcpp::export]]
Rcpp::List mapToUniquePairs(const Rcpp::IntegerVector& r, const Rcpp::IntegerVector& s){
	 if (r.length() != s.length()) {
		Rcpp::stop("Lengths differ.");
	}
	std::vector<int> map(r.length());
	std::vector<int> amount;
	std::vector<int> ur;
	std::vector<int> us;

	map2uniquePairs_core(asVector(r), asVector(s), map, amount, ur, us);

	Rcpp::IntegerMatrix values(2,ur.size());
	for (int i = 0; i < values.ncol(); ++i) {
		values(0,i) = ur[i];
		values(1,i) = us[i];
	}
	return Rcpp::List::create(
			Rcpp::Named("values")=values,
			Rcpp::Named("amount")=Rcpp::IntegerVector(amount.begin(), amount.end()),
			Rcpp::Named("map")=Rcpp::IntegerVector(map.begin(), map.end())+1
 			);
}
//[[Rcpp::export]]
Rcpp::NumericVector mapToOriginal(const Rcpp::NumericVector& vec, const Rcpp::IntegerVector& map) {
	Rcpp::NumericVector out(map.length());
	for (int i = 0; i < map.length(); ++i) {
		out[i] = vec[map[i]-1];
	}
	return out;
}

/*
 * diffr_core() AND em_core() HELPERS
 */
static inline void calculatePost(Rcpp::NumericMatrix& lnPost, Rcpp::NumericVector& lnZ, const Rcpp::NumericVector& r, const Rcpp::NumericVector& s, const Rcpp::NumericVector& lnprior, const Rcpp::NumericVector& lntheta, const Rcpp::NumericVector& lnftheta, const int nthreads=1) {
		int models = lnprior.length();
		for (int k = 0; k < models; ++k) { 
			lnPost(Rcpp::_,k) = r * lnftheta[k] + s * lntheta[k] + lnprior[k];
		}
		lnZ = logRowSum(lnPost, nthreads);
		for (int k = 0; k < models; ++k) {
			lnPost(Rcpp::_,k) = lnPost(Rcpp::_, k) - lnZ;
		}
}

/*
 * diffr_core() EXPECTATION MAXIMIZATION
 */
Rcpp::List em_core(const Rcpp::List m2u_sub, const int models=2, const double eps=.0001, const bool verbose=false, const int nthreads=1) {
	/* 
	 * Get values from mapToUniquePairs structure
	 */
	Rcpp::IntegerVector ur_sub = Rcpp::as<Rcpp::IntegerMatrix>(m2u_sub["values"]).row(0);
	Rcpp::IntegerVector us_sub = Rcpp::as<Rcpp::IntegerMatrix>(m2u_sub["values"]).row(1);
	Rcpp::IntegerVector uamount_sub = Rcpp::as<Rcpp::IntegerVector>(m2u_sub["amount"]); 
	Rcpp::IntegerVector map_sub = Rcpp::as<Rcpp::IntegerVector>(m2u_sub["map"]); 

	/*
	 * Initialization of Mixture Parameters
	 */
	Rcpp::NumericVector lnprior(models), lntheta(models), lnftheta(models);
	if (verbose) {
		Rcpp::Rcout << "\t... initiatilizing prior and theta" << std::endl;
	}
	lnprior = Rcpp::runif(models, 0, 1);
	lnprior = log(lnprior / sum(lnprior));
	//Casting necessary otherwise integer division is performed
	double qstar = ((double) sum(us_sub * uamount_sub)) / ((double) sum(ur_sub * uamount_sub + us_sub * uamount_sub));
	lntheta = log(Rcpp::rep(qstar, models) + Rcpp::runif(models, 0, .01));
	lnftheta = log(1 - exp(lntheta));
	if (verbose) {
		Rcpp::Rcout << "\t... q* = " << qstar << std::endl;
	}

	if (verbose) {
		Rcpp::Rcout << "\t... starting Expectation Maximization" << std::endl;
	}
	int runs = 0, lsub = ur_sub.size();

	Rcpp::NumericVector lnL(1, -DBL_MAX), lnZ(lsub);

	/*
	 * Variables only locally allocated
	 */
	bool notConverged = true;
	double lnPostSum = 0;
	Rcpp::NumericMatrix lnPost(lsub, models);
	Rcpp::NumericVector lnPost_us_sub_log(lsub), lnPost_un_sub_log(lsub);
	//for correct computation I need to convert IntegerVector to NumericVector
	Rcpp::NumericVector ur_sub_numeric(ur_sub); 
	Rcpp::NumericVector us_sub_numeric(us_sub); 
	Rcpp::NumericVector uamount_sub_numeric(uamount_sub); 
	//Compute log counts, s.t. we don't have to do it again all the time
	Rcpp::NumericVector ur_sub_log = log(ur_sub_numeric);
	Rcpp::NumericVector us_sub_log = log(us_sub_numeric);
	Rcpp::NumericVector un_sub_log = log(ur_sub_numeric + us_sub_numeric);
	Rcpp::NumericVector uamount_sub_log = log(uamount_sub_numeric);

	while ((runs < 25) | notConverged) { // run at least 50 times for burn in
		/*
		 * Expectation. Works only on unique values and uses indexes of logRowSum result to restore original matrix
		 */
		calculatePost(lnPost, lnZ, ur_sub_numeric, us_sub_numeric, lnprior, lntheta, lnftheta, nthreads);

		/*
		 * Maximization (done on map reduced data)
		 */
		for (int k = 0; k < models; ++k) {
			//add quantities of each i from mapping (to get overall model posterior)
			lnPost(Rcpp::_,k) = lnPost(Rcpp::_, k) + uamount_sub_log;
		}
		//Calculation of likelihood and model parameters needs amount information from mapping
		lnPostSum = logSumVector(lnPost, nthreads);
		for (int k = 0; k < models; ++k) {
			lnPost_us_sub_log = lnPost(Rcpp::_, k) + us_sub_log;
			lnPost_un_sub_log = lnPost(Rcpp::_, k) + un_sub_log;
			lntheta[k] = logSumVector(lnPost_us_sub_log, nthreads) - logSumVector(lnPost_un_sub_log, nthreads);
			lnftheta[k] = log(1 - exp(lntheta[k]));
			lnprior[k] = logSumVector(lnPost(Rcpp::_, k), nthreads) - lnPostSum;
		}

		/*
		 * Convergence
		 */
		lnZ = lnZ * uamount_sub_numeric;
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
			Rcpp::Rcout << "Run " << runs << 
				": lnL=" << lnLNew << 
				", step=" << (lnLNew - lnL[ lnL.size() - 2]) << 
				", prior=( ";
			for (int k = 0; k < models; k++) { //loop over model components
				Rcpp::Rcout << exp(lnprior[k]) << " ";
			}
			Rcpp::Rcout << ") theta=( ";
			for (int k = 0; k < models; k++) { //loop over model components
				Rcpp::Rcout << exp(lntheta[k]) << " ";
			}
			Rcpp::Rcout << ")" << std::endl;
		}
	}

	if (verbose) {
		Rcpp::Rcout << "\t... finished EM." << std::endl;
	}

	return Rcpp::List::create(
			Rcpp::Named("qstar")=qstar,
			Rcpp::Named("lnprior")=lnprior,
			Rcpp::Named("lntheta")=lntheta,
			Rcpp::Named("lnL")=lnL
 			);
}

/*
 * DIFFR DECONVOLUTION FUNCTION
 */
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
//' @export
//[[Rcpp::export]]
Rcpp::List diffr_core(const Rcpp::IntegerVector& r, const Rcpp::IntegerVector& s, const int models=2, const double eps=.0001, const bool verbose=false, const int nthreads=1) {
	if (verbose) {
		Rcpp::Rcout << "Starting diffR mixture deconvolution" << std::endl;
	}

	if (models < 2) {
		Rcpp::stop("Error: need at least 2 models (i.e. background and foreground)");
	}

	/*
	 * Expectation Maximization on subsetdata
	 *
	 * Subset data for excluding zero regions. 
	 * Uses Rcpp sugar. Therefore ursub and ssub have to be NumericVectors. 
	 * Recreating a new map for subset is much faster than reorganizing the old one.
	 */
	Rcpp::LogicalVector idx = (r > 0 & s > 0);
	if (verbose) {
		Rcpp::Rcout << "\t... removing (r == 0) or (s == 0) regions [" << logical2Count(idx) << " of " << r.length() << " regions kept]." << std::endl;
	}
	Rcpp::List m2u_sub = mapToUniquePairs(r[idx], s[idx]);
	Rcpp::List fit = em_core(m2u_sub, models, eps, verbose, nthreads);
	Rcpp::NumericVector qstar = Rcpp::as<Rcpp::NumericVector>(fit["qstar"]);
	Rcpp::NumericVector lnprior = Rcpp::as<Rcpp::NumericVector>(fit["lnprior"]);
	Rcpp::NumericVector lntheta = Rcpp::as<Rcpp::NumericVector>(fit["lntheta"]);
	Rcpp::NumericVector lnL = Rcpp::as<Rcpp::NumericVector>(fit["lnL"]);

	/*
	 * Map2UniquePairs. This is the data to do computation on and reference computed values by m2u$map.
	 */
	Rcpp::List m2u = mapToUniquePairs(r, s);
	Rcpp::NumericVector ur = Rcpp::as<Rcpp::NumericMatrix>(m2u["values"]).row(0);
	Rcpp::NumericVector us = Rcpp::as<Rcpp::NumericMatrix>(m2u["values"]).row(1);
	Rcpp::IntegerVector umap = Rcpp::as<Rcpp::IntegerVector>(m2u["map"]);

	/*
	 * Compute whole posterior matrix
	 */ 
	if (verbose) {
		Rcpp::Rcout << "\t... computing Posterior for all bins." << std::endl;
	}
	Rcpp::NumericMatrix lnPost(ur.length(), models);
	Rcpp::NumericVector lnZ(ur.length());
	Rcpp::NumericVector lnftheta = log(1 - exp(lntheta));
	calculatePost(lnPost, lnZ, ur, us, lnprior, lntheta, lnftheta, nthreads);

	/*
	 * Sort theta by value and also corresponding prior and posteriors.
	 */
	std::vector<unsigned int> o = indexSort(Rcpp::as<std::vector<double> >(lntheta));
	Rcpp::NumericMatrix post(r.length(), models);
	Rcpp::NumericVector prior(models), theta(models);
	for (int k = 0; k < models; k++) { //loop over model components
		Rcpp::NumericVector postk = exp(lnPost.column(o[k]));
		post(Rcpp::_,k) = mapToOriginal(postk, umap); 
		prior[k] = exp(lnprior(o[k]));
		theta[k] = exp(lntheta(o[k]));
	}

	/*
	 * Return result
	 */
	if (verbose) {
		Rcpp::Rcout << "Finished Expectation maximization em()." << std::endl;
	}
	return Rcpp::List::create(
				Rcpp::Named("qstar")=qstar,
				Rcpp::Named("theta")=theta,
				Rcpp::Named("prior")=prior,
				Rcpp::Named("posterior")=post,
				Rcpp::Named("lnL")=lnL
	);
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
