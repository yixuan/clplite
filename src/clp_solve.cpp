#include <Rcpp.h>
#include <ClpSimplex.hpp>
#include <CoinPackedMatrix.hpp>
#include "matrix.h"

using Rcpp::NumericMatrix;
using Rcpp::NumericVector;
using Rcpp::IntegerVector;
using Rcpp::S4;
using Rcpp::List;
using Rcpp::Named;

// [[Rcpp::export]]
List clp_solve_(
    NumericMatrix obj, S4 mat,
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
    const int nproblem = obj.ncol();
    NumericMatrix sol(nc, nproblem);
    NumericVector objval(nproblem);
    IntegerVector status(nproblem);

    // Obtain the solution for the first problem
    std::copy(solution, solution + nc, sol.begin());
    objval[0] = model.getObjValue();
    status[0] = model.status();

    // If we have more than one problems
    for(int i = 1; i < nproblem; i++)
    {
        model.chgObjCoefficients(&obj(0, i));
        model.primal();
        const int nc = model.numberColumns();
        const double* solution = model.primalColumnSolution();

        std::copy(solution, solution + nc, &sol(0, i));
        objval[i] = model.getObjValue();
        status[i] = model.status();
    }

    return List::create(
        Named("solution") = sol,
        Named("objval")   = objval,
        Named("status")   = status
    );
}


// [[Rcpp::export]]
List clp_solve_parallel_(
    NumericMatrix obj, S4 mat,
    NumericVector constr_lb, NumericVector constr_ub,
    NumericVector var_lb, NumericVector var_ub,
    bool obj_max,
    List control
)
{
    // Number of independent tasks
    const int nproblem = obj.ncol();

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

    // Results
    const int nvar = obj.nrow();
    NumericMatrix sol(nvar, nproblem);
    NumericVector objval(nproblem);
    IntegerVector status(nproblem);

    // Pointers to R objects
    double* sol_ptr = sol.begin();
    double* objval_ptr = objval.begin();
    int* status_ptr = status.begin();

    // Number of threads for OpenMP
    const int nthread = Rcpp::as<int>(control["verbose"]);

    // Solving
#ifdef _OPENMP
#pragma omp parallel for num_threads(nthread) schedule(static)
#endif
    for(int i = 0; i < nproblem; i++)
    {
        // Create model
        ClpSimplex model;

        // -1 for maximization, 1 for minimization
        model.setOptimizationDirection(obj_max ? -1 : 1);

        // Logging level -- we need to turn off logging in parallel mode
        model.setLogLevel(0);

        // Add constraints
        model.loadProblem(*constr, var_lb.begin(), var_ub.begin(), obj.begin() + nvar * i, constr_lb.begin(), constr_ub.begin());

        // Solve
        model.primal();

        // Get pointer to solution
        const double* solution = model.primalColumnSolution();

        // Copy solution to R objects
        std::copy(solution, solution + nvar, sol_ptr + nvar * i);
        objval_ptr[i] = model.getObjValue();
        status_ptr[i] = model.status();
    }

    // Release resource
    delete constr;

    return List::create(
        Named("solution") = sol,
        Named("objval")   = objval,
        Named("status")   = status
    );
}
