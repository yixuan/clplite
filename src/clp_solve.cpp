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

/*
 * minimize   -1 x0 - 1 x1
 * s.t         1 x0 + 2 x1 <= 3
 *             2 x0 + 1 x1 <= 3
 *             x0 >= 0
 *             x1 >= 0
 */

// [[Rcpp::export]]
void test()
{
    ClpSimplex model;

    IntegerVector mi = IntegerVector::create(0, 1, 0, 1);
    IntegerVector mj = IntegerVector::create(0, 0, 1, 1);
    NumericVector mx = NumericVector::create(1, 2, 2, 1);
    CoinPackedMatrix constr(false, mi.begin(), mj.begin(), mx.begin(), mx.length());
    constr.setDimensions(2, 2);

    NumericVector lhs = NumericVector::create(-std::numeric_limits<double>::max(), -std::numeric_limits<double>::max());
    NumericVector rhs = NumericVector::create(3.0, 3.0);
    NumericVector obj = NumericVector::create(-1.0, -1.0);
    NumericVector lb  = NumericVector::create(0.0, 0.0);
    NumericVector ub  = NumericVector::create(std::numeric_limits<double>::max(), std::numeric_limits<double>::max());

    model.loadProblem(constr, lb.begin(), ub.begin(),
                      obj.begin(), lhs.begin(), rhs.begin());
    model.primal();

    int numberColumns = model.numberColumns();
    double * columnPrimal = model.primalColumnSolution();
    double * columnDual = model.dualColumnSolution();

    int iColumn;

    for (iColumn=0;iColumn<numberColumns;iColumn++)
            Rprintf("Column %d, primal %g, dual %g\n",iColumn,
               columnPrimal[iColumn],columnDual[iColumn]);
}
