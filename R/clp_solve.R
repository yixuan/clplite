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
#' @param obj Typically a vector representing the objective coefficients \eqn{c}.
#'            If \code{obj} is a matrix of \eqn{N} columns, then \eqn{N} linear
#'            programming problems that have the same constraints but different
#'            objective coefficients will be solved.
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
#'
#' @return A list with the following components:
#'
#' \describe{
#'   \item{\code{solution}}{If the input \code{obj} is a vector, then \code{solution}
#'                          contains the vector of optimal coefficients. If \code{obj}
#'                          is a matrix, then the solutions are combined horizontally.}
#'   \item{\code{objval}}{The value of the objective function at the solution.}
#'   \item{\code{status}}{Status code for each problem solved.}
#' }
#'
#' @author Yixuan Qiu \email{yixuanq@@gmail.com}
#'
#' @export
#' @keywords optimize
#' @examples
#' ## Example from the Rglpk pakcage
#' ##
#' ## maximize:   2 x1 + 4 x2 + 3 x3
#' ## subject to: 3 x1 + 4 x2 + 2 x3 <= 60
#' ##             2 x1 +   x2 + 2 x3 <= 40
#' ##               x1 + 3 x2 + 2 x3 <= 80
#' ## x1, x2, x3 are non-negative real numbers
#'
#' obj = c(2, 4, 3)
#' mat = matrix(c(3, 2, 1, 4, 1, 3, 2, 2, 2), nrow = 3)
#' constr_ub = c(60, 40, 80)
#' var_lb = c(0, 0, 0)
#'
#' clp_solve(obj, mat, constr_ub = constr_ub, var_lb = var_lb, max = TRUE)
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
