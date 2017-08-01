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
List read_mps_(std::string mps_file)
{
    // Create model
    ClpSimplex model;

    // Read model from MPS file
    int status = model.readMps(mps_file.c_str());
    if(status)
        Rcpp::stop("reading MPS file '" + mps_file + "' failed");

    // Objective coefficients
    const int nr = model.getNumRows();
    const int nc = model.getNumCols();
    const double* obj_data = model.getObjCoefficients();
    NumericVector obj(obj_data, obj_data + nc);

    // Constraint matrix
    CoinPackedMatrix* constr_data = model.matrix();
    Rcpp::S4 constr = from_CoinPackedMatrix(constr_data);

    // Bounds
    const double* constr_lb_data = model.getRowLower();
    const double* constr_ub_data = model.getRowUpper();
    const double* var_lb_data = model.getColLower();
    const double* var_ub_data = model.getColUpper();
    NumericVector constr_lb(constr_lb_data, constr_lb_data + nr);
    NumericVector constr_ub(constr_ub_data, constr_ub_data + nr);
    NumericVector var_lb(var_lb_data, var_lb_data + nc);
    NumericVector var_ub(var_ub_data, var_ub_data + nc);

    return List::create(
        Named("obj")       = obj,
        Named("constr")    = constr,
        Named("constr_lb") = constr_lb,
        Named("constr_ub") = constr_ub,
        Named("var_lb")    = var_lb,
        Named("var_ub")    = var_ub
    );
}

// [[Rcpp::export]]
List clp_mps_solve_(
    std::string mps_file,
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

    // Read model from MPS file
    int status = model.readMps(mps_file.c_str());
    if(status)
        Rcpp::stop("reading MPS file '" + mps_file + "' failed");

    // Solve
    model.primal();

    // Get pointer to solution
    const int nc = model.getNumCols();
    const double* solution = model.getColSolution();

    // Wrap result -- MPS file only has one objective
    NumericVector sol(solution, solution + nc);

    const double objval = model.getObjValue();
    status = model.status();

    return List::create(
        Named("solution") = sol,
        Named("objval")   = objval,
        Named("status")   = status
    );
}
