#include "matrix.h"

using Rcpp::NumericVector;
using Rcpp::IntegerVector;
using Rcpp::List;

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

Rcpp::S4 from_CoinPackedMatrix(const CoinPackedMatrix* mat)
{
    const bool col_major = mat->isColOrdered();
    std::string klass = col_major ? "dgCMatrix" : "dgRMatrix";
    const int nr           = mat->getNumRows();
    const int nc           = mat->getNumCols();
    const int nnz          = mat->getNumElements();
    const double* mat_data = mat->getElements();
    const int* i_or_j      = mat->getIndices();
    const int* p           = mat->getVectorStarts();

    Rcpp::S4 res(klass);
    res.slot(col_major ? "i" : "j") = IntegerVector(i_or_j, i_or_j + nnz);
    res.slot("p")        = IntegerVector(p, p + mat->getSizeVectorStarts());
    res.slot("Dim")      = IntegerVector::create(nr, nc);
    res.slot("Dimnames") = List::create(R_NilValue, R_NilValue);
    res.slot("x")        = NumericVector(mat_data, mat_data + nnz);
    res.slot("factors")  = List();

    return res;
}
