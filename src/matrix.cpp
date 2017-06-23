#include "matrix.h"

using Rcpp::NumericVector;
using Rcpp::IntegerVector;

CoinPackedMatrix* from_dgTMatrix(Rcpp::S4 mat)
{
    IntegerVector dim = mat.slot("Dim");
    IntegerVector mi  = mat.slot("i");
    IntegerVector mj  = mat.slot("j");
    NumericVector mx  = mat.slot("x");

    CoinPackedMatrix* cmat = new CoinPackedMatrix(false, mi.begin(), mj.begin(), mx.begin(), mx.length());
    cmat->setDimensions(dim[0], dim[1]);

    return cmat;
}

CoinPackedMatrix* from_dgCMatrix(Rcpp::S4 mat)
{
    IntegerVector dim  = mat.slot("Dim");
    IntegerVector mi   = mat.slot("i");
    IntegerVector mp   = mat.slot("p");
    IntegerVector mlen = Rcpp::diff(mp);
    NumericVector mx   = mat.slot("x");

    CoinPackedMatrix* cmat = new CoinPackedMatrix(true, dim[0], dim[1], mx.length(),
                                     mx.begin(), mi.begin(), mp.begin(), mlen.begin());
    cmat->setDimensions(dim[0], dim[1]);

    return cmat;
}

CoinPackedMatrix* from_dgRMatrix(Rcpp::S4 mat)
{
    IntegerVector dim  = mat.slot("Dim");
    IntegerVector mj   = mat.slot("j");
    IntegerVector mp   = mat.slot("p");
    IntegerVector mlen = Rcpp::diff(mp);
    NumericVector mx   = mat.slot("x");

    CoinPackedMatrix* cmat = new CoinPackedMatrix(false, dim[1], dim[0], mx.length(),
                                     mx.begin(), mj.begin(), mp.begin(), mlen.begin());
    cmat->setDimensions(dim[0], dim[1]);

    return cmat;
}
