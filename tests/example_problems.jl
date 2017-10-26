 workspace()

# # include pure julia osqp solver
include("../osqp_julia.jl")
using OSQPSolver

# # include python interface for native osqp solver
# using PyCall
# @pyimport osqp as osqp
# @pyimport scipy.sparse as spar
# @pyimport numpy as np

# Simple problem
# m = 50;
# n = 100;
# A  = sparse(randn(m,n));
# l = -rand(m,1) * 2;
# u = +rand(m,1) * 2;
# P = Symmetric(sprand(n,n,0.1));
# P = P + n*speye(n)
# P = sparse(P)
# q = randn(n,1);
# settings = qpSettings(rho=1.0,verbose=true)
# res = solveOSQP(P,q,A,l,u,settings);


# Primal infeasible problem
n = 50;
m = 500;
Pt = sprandn(n, n, 0.6);
P = Pt' * Pt;
q = randn(n, 1);
A = sprandn(m, n, 0.8);
u = 3 + randn(m, 1);
l = -3 + randn(m, 1);

# # Make random problem primal infeasible
nhalf = Int64(floor(n/2));
A[nhalf, :] = A[nhalf + 1, :];
l[nhalf] = u[nhalf + 1] + 10 * rand();
u[nhalf] = l[nhalf] + 0.5;

settings = qpSettings(rho=1.0,verbose=true)
res = solveOSQP(P,q,A,l,u,settings);

# # Dual infeasible problem
# P = sparse(diagm([4; 0]));
# q = [0.0; 2];
# A = sparse([1.0 1; -1 1]);
# l = [-Inf; -Inf];
# u = [2.0; 3];
# settings = qpSettings(rho=1.0,verbose=true)
# res = solveOSQP(P,q,A,l,u,settings);

# # Setup and solve the problem
# # Setup settings
# settings.alpha = 1.6;
# settings.rho = 0.1;
# settings.sigma = 0.1;
# settings.eps_prim_inf = 1e-5;
# settings.eps_dual_inf = 1e-5;
# settings.eps_rel = 1e-5;
# settings.eps_abs = 1e-5;
# settings.max_iter = 2500;
# settings.verbose = 1;
# settings.scaling = 0; % Disable scaling. Pure Matlab implementation does not support it yet



# #  Solve with osqp
# solver = osqp;
# solver.setup(problem.P, problem.q, problem.A, problem.l, problem.u, settings);
# resOSQP = solver.solve();


# # Solve with osqpmatlab
# [xmatlab, ymatlab, costmatlab, statusmatlab, itermatlab] = osqpmatlab(problem,[], settings);


# # Check solution if the problem is solved
# if resOSQP.info.status_val == 1
#     assert(norm(xmatlab - resOSQP.x) < settings.eps_abs*10, 'Error primal solution')
#     assert(norm(ymatlab - resOSQP.y) < settings.eps_abs*10, 'Error dual solution')
#     assert(norm(costmatlab - resOSQP.info.obj_val) < settings.eps_abs*10, 'Error cost function')
# end