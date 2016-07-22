
/*
 * Code originally adapted from:
 *   http://systematicinvestor.github.io/Run-Correlation-Rcpp 
 *
 * version: 0.9
 */

#include <Rcpp.h>
using namespace Rcpp;
using namespace std;

inline NumericMatrix c_cov_helper(const NumericMatrix& mat, const int rstart, const int rend
	, NumericVector means, NumericVector vars, NumericVector sds
	, const int cstart, const int cend) {
	int nc = mat.ncol();
	int nperiod = rend - rstart;
	int ksize = cend - cstart + 1;
	NumericMatrix routl(ksize, 3);

	// pre-compute outliers
	for (int row = 0; row < nc; row++) {
		for (int col = 0; col < ksize; col++) {
			double sXY = 0;
			double aux;

			for (int r = rstart; r < rend; r++)
	    			sXY += mat(r, row) * mat(r, col + cstart);
		   	aux = sXY / (nperiod - 1);
			routl(col, 0) += (aux / sds[row] - routl(col, 0)) / (row+1); //aggregated mean
			routl(col, 1) += (aux * (means[row] / vars[row]) - routl(col, 1)) / (row+1); //aggregated mean
			//routl(col, 2) += (aux / vars[row] - routl(col, 0)) / (row+1);
			routl(col, 2) += aux / vars[row];
		}
	}
	// perform mean over routl
	for (int col = 0; col < ksize; col++) {
		//routl(col, 0) /= nc;
		//routl(col, 1) /= nc;
		routl(col, 2) /= nc;
	}
	return routl;
}

// [[Rcpp::export]]
NumericMatrix c_cov(NumericMatrix mat
	, NumericVector means, NumericVector vars, NumericVector sds
	, int c_ini, int c_end) {
	return c_cov_helper(mat, 0, mat.nrow(), means, vars, sds, c_ini-1, c_end-1);
}

