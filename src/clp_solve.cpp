#include <Rcpp.h>
#include <ClpSimplex.hpp>
#include <CoinPackedMatrix.hpp>

using Rcpp::NumericVector;
using Rcpp::IntegerVector;

// [[Rcpp::export]]
NumericVector clp_solve_(
    NumericVector obj,
    int nrow, int ncol, IntegerVector mi, IntegerVector mj, NumericVector mx,
    NumericVector constr_lb, NumericVector constr_ub,
    NumericVector var_lb, NumericVector var_ub,
    bool obj_max
)
{
    // Create model
    ClpSimplex model;

    // -1 for maximization, 1 for minimization
    model.setOptimizationDirection(obj_max ? -1 : 1);

    // Constraint matrix
    CoinPackedMatrix constr(false, mi.begin(), mj.begin(), mx.begin(), mx.length());
    constr.setDimensions(nrow, ncol);

    // Add constraints
    model.loadProblem(constr, var_lb.begin(), var_ub.begin(), obj.begin(), constr_lb.begin(), constr_ub.begin());

    // Solve
    model.primal();

    // Get pointer to solution
    const int nc = model.numberColumns();
    const double* solution = model.primalColumnSolution();

    // Wrap result
    NumericVector sol(solution, solution + nc);

    return sol;
}
