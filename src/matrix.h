#ifndef CLPLITE_MATRIX_H
#define CLPLITE_MATRIX_H

#include <Rcpp.h>
#include <CoinPackedMatrix.hpp>

CoinPackedMatrix* from_dgTMatrix(Rcpp::S4 mat);
CoinPackedMatrix* from_dgCMatrix(Rcpp::S4 mat);
CoinPackedMatrix* from_dgRMatrix(Rcpp::S4 mat);

Rcpp::S4 from_CoinPackedMatrix(const CoinPackedMatrix* mat);


#endif // CLPLITE_MATRIX_H
