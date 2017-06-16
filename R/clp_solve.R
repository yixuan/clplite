clp_solve = function(obj, A,
                     constr_lb = rep(-Inf, nrow(A)), constr_ub = rep(Inf, nrow(A)),
                     var_lb = rep(0, ncol(A)), var_ub = rep(Inf, ncol(A)),
                     max = FALSE)
{
    A = as(A, "dgTMatrix")

    clp_solve_(obj, nrow(A), ncol(A), A@i, A@j, A@x,
               constr_lb, constr_ub,
               var_lb, var_ub, max)
}
