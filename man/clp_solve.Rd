% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/clp_solve.R
\name{clp_solve}
\alias{clp_solve}
\title{Linear Programming Solver Using Clp}
\usage{
clp_solve(obj, A, constr_lb = rep(-Inf, nrow(A)), constr_ub = rep(Inf,
  nrow(A)), var_lb = rep(0, ncol(A)), var_ub = rep(Inf, ncol(A)),
  max = FALSE, control = list())
}
\arguments{
\item{obj}{Typically a vector representing the objective coefficients \eqn{c}.
If \code{obj} is a matrix of \eqn{N} columns, then \eqn{N} linear
programming problems that have the same constraints but different
objective coefficients will be solved.}

\item{A}{A matrix representing the constraint matrix \eqn{A}. Supported
matrix types include \code{matrix}, \code{dgTMatrix}, \code{dgCMatrix},
\code{dgRMatrix}, and other ones that can be coerced to these
types.}

\item{constr_lb}{Lower bounds of the constraints. Infinite values are allowed.}

\item{constr_ub}{Upper bounds of the constraints. Infinite values are allowed.}

\item{var_lb}{Lower bounds of the objective variables. Infinite values are allowed.}

\item{var_ub}{Upper bounds of the objective variables. Infinite values are allowed.}

\item{max}{Logical value. \code{TRUE} for maximization and \code{FALSE} for
minimization.}
}
\value{
A list with the following components:

\describe{
  \item{\code{solution}}{If the input \code{obj} is a vector, then \code{solution}
                         contains the vector of optimal coefficients. If \code{obj}
                         is a matrix, then the solutions are combined horizontally.}
  \item{\code{objval}}{The value of the objective function at the solution.}
  \item{\code{status}}{Status code for each problem solved.}
}
}
\description{
High level R interface to the COIN-OR linear programming (Clp) library for
solving large-scale linear programming problems. Sparse constraints are fully
supported.
}
\details{
This function solves linear programming problems of the following form:
\preformatted{
min     c' * x
s.t.    constr_lb <= A * x <= constr_ub
           var_lb <= x     <= var_ub
}
}
\examples{
## Example from the Rglpk pakcage
##
## maximize:   2 x1 + 4 x2 + 3 x3
## subject to: 3 x1 + 4 x2 + 2 x3 <= 60
##             2 x1 +   x2 + 2 x3 <= 40
##               x1 + 3 x2 + 2 x3 <= 80
## x1, x2, x3 are non-negative real numbers

obj = c(2, 4, 3)
mat = matrix(c(3, 2, 1, 4, 1, 3, 2, 2, 2), nrow = 3)
constr_ub = c(60, 40, 80)
var_lb = c(0, 0, 0)

clp_solve(obj, mat, constr_ub = constr_ub, var_lb = var_lb, max = TRUE)

## Handling sparse matrices
if(require("Matrix")) {
    spmat = as(mat, "sparseMatrix")
    clp_solve(obj, spmat, constr_ub = constr_ub, var_lb = var_lb, max = TRUE)
}
}
\author{
Yixuan Qiu <\email{yixuanq@gmail.com}>
}
\keyword{optimize}
