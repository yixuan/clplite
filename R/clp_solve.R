#' Linear Programming Solver Using Clp
#'
#' High level R interface to the COIN-OR linear programming (Clp) library for
#' solving large-scale linear programming problems. Sparse constraints are fully
#' supported.
#'
#' This function solves linear programming problems of the following form:
#' \preformatted{
#' min     c' * x
#' s.t.    constr_lb <= A * x <= constr_ub
#'            var_lb <= x     <= var_ub
#' }
#'
#' @param obj A vector representing the objective coefficients \eqn{c}.
#' @param A   A matrix representing the constraint matrix \eqn{A}. Supported
#'            matrix types include \code{matrix}, \code{dgTMatrix}, \code{dgCMatrix},
#'            \code{dgRMatrix}, and other ones that can be coerced to these
#'            types.
#' @param constr_lb Lower bounds of the constraints. Infinite values are allowed.
#' @param constr_ub Upper bounds of the constraints. Infinite values are allowed.
#' @param var_lb Lower bounds of the objective variables. Infinite values are allowed.
#' @param var_ub Upper bounds of the objective variables. Infinite values are allowed.
#' @param max Logical value. \code{TRUE} for maximization and \code{FALSE} for
#'            minimization.
clp_solve = function(obj, A,
                     constr_lb = rep(-Inf, nrow(A)), constr_ub = rep(Inf, nrow(A)),
                     var_lb = rep(0, ncol(A)), var_ub = rep(Inf, ncol(A)),
                     max = FALSE, control = list())
{
    supported_classes = c("dgTMatrix", "dgCMatrix", "dgRMatrix")
    if(!(class(A) %in% supported_classes))
    {
        ## dgRMatrix has the highest priority since Clp use row-oriented sparse matrix
        if(canCoerce(A, "dgRMatrix"))
        {
            A = as(A, "dgRMatrix")
        } else if(canCoerce(A, "dgTMatrix")) {
            A = as(A, "dgTMatrix")
        } else if(canCoerce(A, "dgCMatrix")) {
            A = as(A, "dgCMatrix")
        } else if(canCoerce(A, "matrix")) {
            A = as(as(A, "matrix"), "dgRMatrix")
        } else {
            stop("unsupported matrix type")
        }
    }

    ## Control parameters
    control_param = list(
        verbose = 0
    )
    pars = intersect(names(control_param), names(control))
    control_param[pars] = control[pars]

    clp_solve_(obj, A, constr_lb, constr_ub, var_lb, var_ub, max, control_param)
}
