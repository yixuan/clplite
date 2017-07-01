## High Performance Linear Programming Solver Based on Clp

This package provides an R interface to the high performance linear
programming library [Clp](https://projects.coin-or.org/Clp). Relevant
source code of **Clp** is included in the R package for easy installation.
A C++11 compiler is required if the package is built from source.

Below is an example problem taken from
[IBM Optimization Library Guide and Reference](https://www.cenapad.unicamp.br/parque/manuais/OSL/oslweb/features/feat24DT.htm).

```
  Minimize or maximize Z = x1 + 2x5 - x8

  Subject to:

  2.5 <=   3x1 +  x2         -  2x4 - x5               -    x8
                 2x2 + 1.1x3                                   <=  2.1
                          x3              +  x6                ==  4.0
  1.8 <=                      2.8x4             -1.2x7         <=  5.0
  3.0 <= 5.6x1                      + x5               + 1.9x8 <= 15.0

  where:

  2.5 <= x1
    0 <= x2 <= 4.1
    0 <= x3
    0 <= x4
  0.5 <= x5 <= 4.0
    0 <= x6
    0 <= x7
    0 <= x8 <= 4.3
```

Using **clplite**, we can solve this problem by the following R code:

```r
library(clplite)
obj = c(1, 0, 0, 0, 2, 0, 0, -1)
A = matrix(c(
      3,   1,   0,  -2,  -1,   0,    0,   -1,
      0,   2, 1.1,   0,   0,   0,    0,    0,
      0,   0,   1,   0,   0,   1,    0,    0,
      0,   0,   0, 2.8,   0,   0, -1.2,    0,
    5.6,   0,   0,   0,   1,   0,    0,  1.9
), nrow = 5, byrow = TRUE)
constr_lb = c(2.5, -Inf,  4, 1.8,  3)
constr_ub = c(Inf,  2.1,  4,   5, 15)
var_lb    = c(2.5,   0,   0,   0, 0.5,   0,   0,   0)
var_ub    = c(Inf, 4.1, Inf, Inf,   4, Inf, Inf, 4.3)

clp_solve(obj, A, constr_lb, constr_ub, var_lb, var_ub, max = TRUE)
```

```
$solution
[1] 2.5000000 0.0000000 0.0000000 0.6428571 1.0000000 4.0000000 0.0000000 0.0000000

$objval
[1] 4.5

$status
[1] 0
```

```r
clp_solve(obj, A, constr_lb, constr_ub, var_lb, var_ub, max = FALSE)
```

```
$solution
[1] 2.5000000 0.0000000 0.0000000 0.6428571 0.5000000 4.0000000 0.0000000 0.2631579

$objval
[1] 3.236842

$status
[1] 0
```

A significant feature of **Clp** and **clplite** is the ability to handle
sparse constraint matrices, which is critical in handling large-scale
optimization problems. Currently **clplite** supports matrices of types
`matrix`, `dgTMatrix`, `dgCMatrix`, `dgRMatrix`, and types that can be
coerced to them via the `as()` function.

```r
library(Matrix)
As = as(A, "sparseMatrix")
clp_solve(obj, As, constr_lb, constr_ub, var_lb, var_ub, max = TRUE)
clp_solve(obj, As, constr_lb, constr_ub, var_lb, var_ub, max = FALSE)
```
