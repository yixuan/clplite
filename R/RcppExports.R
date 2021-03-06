# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

clp_solve_ <- function(obj, mat, constr_lb, constr_ub, var_lb, var_ub, obj_max, control) {
    .Call('_clplite_clp_solve_', PACKAGE = 'clplite', obj, mat, constr_lb, constr_ub, var_lb, var_ub, obj_max, control)
}

clp_solve_parallel_ <- function(obj, mat, constr_lb, constr_ub, var_lb, var_ub, obj_max, control) {
    .Call('_clplite_clp_solve_parallel_', PACKAGE = 'clplite', obj, mat, constr_lb, constr_ub, var_lb, var_ub, obj_max, control)
}

read_mps_ <- function(mps_file) {
    .Call('_clplite_read_mps_', PACKAGE = 'clplite', mps_file)
}

clp_mps_solve_ <- function(mps_file, obj_max, control) {
    .Call('_clplite_clp_mps_solve_', PACKAGE = 'clplite', mps_file, obj_max, control)
}

