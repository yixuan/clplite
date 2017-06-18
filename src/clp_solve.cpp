#include <Rcpp.h>
#include <ClpSimplex.hpp>
#include <CoinPackedMatrix.hpp>
#include "matrix.h"

using Rcpp::NumericVector;
using Rcpp::IntegerVector;
using Rcpp::S4;
using Rcpp::List;

// [[Rcpp::export]]
NumericVector clp_solve_(
    NumericVector obj, S4 mat,
    NumericVector constr_lb, NumericVector constr_ub,
    NumericVector var_lb, NumericVector var_ub,
    bool obj_max,
    List control
)
{
    // Create model
    ClpSimplex model;

    // -1 for maximization, 1 for minimization
    model.setOptimizationDirection(obj_max ? -1 : 1);

    // Logging level
    model.setLogLevel(Rcpp::as<int>(control["verbose"]));

    // Constraint matrix
    std::string type = mat.attr("class");
    CoinPackedMatrix* constr;
    if(type == "dgTMatrix")
    {
        constr = from_dgTMatrix(mat);
    } else if(type == "dgCMatrix") {
        constr = from_dgCMatrix(mat);
    } else if(type == "dgRMatrix") {
        constr = from_dgRMatrix(mat);
    } else {
        Rcpp::stop("unsupported matrix type");
    }

    // Add constraints
    model.loadProblem(*constr, var_lb.begin(), var_ub.begin(), obj.begin(), constr_lb.begin(), constr_ub.begin());
    delete constr;

    // Solve
    model.primal();

    // Get pointer to solution
    const int nc = model.numberColumns();
    const double* solution = model.primalColumnSolution();

    // Wrap result
    NumericVector sol(solution, solution + nc);

    return sol;
}
