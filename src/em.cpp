/**
 * Binomial mixture deconvolution with Expectation Minimization
 *
 * uses Rcpp sugar with represents vectorized functions in Rcpp http://adv-r.had.co.nz/Rcpp.html#rcpp-sugar
 * 
 * To avoid numerical under- and overflows, the EM is done in logspace.
 *
 * 2014-11-18
 *
 * helmuth <helmuth@molgen.mpg.de> 
 *
 * TODO refactor code for packaging
 *
 */
#include <Rcpp.h>

using namespace Rcpp;

// [[Rcpp::export]]
NumericVector logSum( const NumericMatrix& mat ){
 	NumericVector rowSum(mat.nrow());
	NumericVector::iterator rowSum_iter = rowSum.begin();
	double tmp, max;

	// iterate over each row: 1. find max value 2. calculate ln sum
	for (int i = 0; i < mat.nrow(); ++i, ++rowSum_iter) {
		max = -DBL_MAX;
		for (int j = 0; j < mat.ncol(); ++j) {
			if ( mat(i,j) > max ) {
				max = mat(i, j);
			}
		}
		for (int j = 0; j < mat.ncol(); ++j) {
			tmp += exp( mat(i,j) - max );  
		}
		*rowSum_iter = max + log( tmp );
		tmp = 0;
	}
				
	return rowSum;
}

// [[Rcpp::export]]
double logSumVector( const NumericVector& vec ) {
	return logSum( NumericMatrix( 1, vec.length(), vec.begin() ) )[0];
}

struct CompareIndex {
	const std::vector< double > values;

	CompareIndex(const std::vector< double > vec) : values( vec ) {};

	bool operator() ( unsigned int l, unsigned int r)  {
		return values[l] < values[r];
	} 
};
std::vector< unsigned int > indexSort( std::vector< double > data ) {
	std::vector< unsigned int > index ( data.size() );
	for ( unsigned int i = 0; i < index.size(); ++i ) index[i] = i;

	// struct that takes values and sorts a given vector based on internal values
	CompareIndex compix( data );
	std::sort(index.begin(), index.end(), compix);

	return index;
}

double binomialCoeff( double n, double k) {
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

// returns a NumericMatrix mat[ , -col]
NumericMatrix subMatrix( const NumericMatrix& mat, int col ) {
	NumericMatrix matslice (mat.nrow(), mat.ncol()-1);

	for (int i=0, j=0; i < mat.ncol()-1; i++) {
		if ( i != col ) {
			matslice(_,j) = mat(_, i);
			++j;
		}
	}

	return(matslice);
}

// [[Rcpp::export]]
List em( const NumericVector& r, const NumericVector& s, int models = 2, double eps = .001, bool verbose = false ) {
	if ( verbose ) {
		Rcout << "Starting Expectation maximization em()" << std::endl;
	}

	if (models < 2) {
		Rcout << "Error: need at least 2 models (i.e. background and foreground)" << std::endl;
		return List::create(); 
	}

	//subset data (uses Rcpp sugar)
	LogicalVector idx = (r + s) > 0;
	NumericVector rsub = r[ idx ];
	NumericVector ssub = s[ idx ];
	double qstar = sum( ssub ) / sum( rsub + ssub );
	if ( verbose ) {
		Rcout << "\t... removing (r + s) == 0 [ " << ( r.length() - rsub.length() ) << " bins ]" << std::endl;
	}

	//initialize data structures
	if ( verbose ) {
		Rcout << "\t... initiatilizing prior and theta" << std::endl;
	}
	NumericVector lnprior(models), lntheta(models), lnftheta(models);
	lnprior = runif(models, 0, 1);
	lnprior = log( lnprior / sum(lnprior) );
	lntheta = rep( qstar, models ) + runif(models, 0, .01);
	lnftheta = log( 1 - lntheta );
	lntheta = log( lntheta );

	//EM algorithm in log space to avoid over- and underflows
	if ( verbose ) {
		Rcout << "\t... starting EM" << std::endl;
	}
	int runs = 0, lsub = rsub.size();
	bool notConverged = true;
	double lnPostSum = 0;
	NumericVector lnL(1, -DBL_MAX), lnZ(lsub);
	NumericMatrix lnPost(lsub, models);
	while ( (runs < 10) | notConverged ) {

		//Expectation
		for (int k = 0; k < models; ++k) { //loop over model components and Rcpp sugar
			lnPost(_,k) = rsub * lnftheta( k ) + ssub * lntheta( k ) + lnprior( k );
		}
		lnZ = logSum( lnPost );

		//Maximization
		for (int k = 0; k < models; ++k) { //loop over model components and Rcpp sugar
			lnPost(_,k) = lnPost(_,k) - lnZ;
		}
		lnPostSum = logSumVector( lnPost );
		for (int k = 0; k < models; ++k) { //loop over model components and Rcpp sugar
			lntheta( k ) = log( sum( exp( lnPost(_,k) ) * ssub ) / sum( exp( lnPost(_,k) ) * (rsub + ssub) ) ); //TODO Can we do better here?
			lnftheta( k ) = log( 1 - exp( lntheta( k ) ) );
			lnprior( k ) = logSumVector( lnPost(_,k) ) - lnPostSum;
		}

		//Convergence
		double lnLNew = sum( lnZ );
		if ( runs > 10 && lnLNew - lnL[ lnL.size() - 1] < eps ) {
			notConverged = false;
		}
		lnL.push_back( lnLNew );

		//Logging
		runs++;
		if ( verbose && ( !notConverged || (runs % 10) == 0 ) ) {
			Rcout << "Run " << runs << ": lnL=" << lnL[ lnL.size()-1 ] << ", step=" << (lnL[lnL.size() - 1] - lnL[ lnL.size() - 2]) << ", prior=";
			for (int k = 0; k < models; k++) { //loop over model components
				Rcout << exp( lnprior[k] ) << " ";
			}
			Rcout << ", theta=";
			for (int k = 0; k < models; k++) { //loop over model components
				Rcout << exp( lntheta[k] ) << " ";
			}
			Rcout << std::endl;
		}
	}
	if ( verbose ) {
		Rcout << "\t... finished EM." << std::endl;
	}

	//binomial computation of posteriors to get posterios that add up to one
	if ( verbose ) {
		Rcout << "\t... computing overall Posterior" << std::endl;
	}
	NumericMatrix lnPost2( r.length(), models);
	NumericVector lnBinomCoeff( r.length() ); 
	NumericVector::iterator r_it = r.begin(), s_it = s.begin(), lnBinomCoeff_it = lnBinomCoeff.begin();
	for (; r_it != r.end(); ++r_it, ++s_it, ++lnBinomCoeff_it) {
		*lnBinomCoeff_it = binomialCoeff( (*r_it + *s_it), *s_it );
	}
	lnBinomCoeff = log( lnBinomCoeff );
	for (int k = 0; k < models; ++k) { //loop over model components
		lnPost2(_,k) = lnBinomCoeff + r * lnftheta( k ) + s * lntheta( k ) + lnprior( k );
	}
	lnZ = logSum( lnPost2 );
	for (int k = 0; k < models; k++) { //loop over model components
		lnPost2(_,k) = lnPost2(_,k) - lnZ; //Rcpp sugar 
	}

	//folchange calculation
	if ( verbose ) {
		Rcout << "\t... computing Foldchange" << std::endl;
	}
	std::vector< unsigned int > o = indexSort( as<std::vector<double> >( lntheta ) );
	NumericMatrix fc( r.length(), models-1 );
	for (int k = 1; k < models; ++k) {

		//foldchange for k against all other components
		if (models == 2) {
			fc(_,0) = lnPost2(_, o[1] ) - lnPost2(_, o[0]);
		} else {
			fc(_,k-1) = lnPost2(_, o[ k ] ) - logSum( subMatrix( lnPost2, k ) );
		}
	}

	// Posteriors and sorted
	NumericMatrix post( r.length(), models);
	NumericVector prior(models), theta(models);
	for (int k = 0; k < models; k++) { //loop over model components
		post(_,k) = exp( lnPost2(_,o[k]) ); //Rcpp sugar 
		prior(k) = exp( lnprior( o[k] ) );
		theta(k) = exp( lntheta( o[k] ) );
	}

	if ( verbose ) {
		Rcout << "Finished Expectation maximization em()." << std::endl;
	}
	return List::create(_("posterior")=post, _("foldchange")=fc,_("fit")=List::create(_("qstar")=qstar,_("theta")=theta, _("prior")=prior, _("lnL")=lnL));
}
