from scipy import sparse
import numpy as np
from multigrid.solver import Solver
import multigrid
from numerical_solver import NumericalSolver
import matplotlib.pyplot as plt

# Load the linear system that we want to solve
A = sparse.load_npz('data/axb_10k/A.npz')
b = np.load('data/axb_10k/b.npy')


multigrid.init()
# Initialize multigrid settings:
# Defailt values are:
#       multigrid.pre_smoother = 'jacobi'
#       multigrid.pre_smooth_iters = 2
#       multigrid.post_smoother = 'jacobi'
#       multigrid.post_smooth_iters = 2

# Note: If you have prolongators available, load them into
#       multigrid.prolongators as scipy.sparse matrix.
#       E.g: If the prolongator for the first level of v-cycle is p1, then:
#            multigrid.prolongators = {1: p1}


# Phase I: Just a couple of iterations, so we could see multigrid process
#          more clearly.
nm_solver = NumericalSolver(A, b, n=50)
phase1_x, phase1_norms = nm_solver.jacobi()

# Phase II: Multigrid cycle
solver = Solver(A, b, initial_guess=phase1_x)
x = solver.solve(number_of_levels=5)

# Phase III: More jacobi iterations, so we could see multigrid process
#            more clearly.
nm_solver = NumericalSolver(A, b, n=100, x=x)
x, norms = nm_solver.jacobi()

# Plotting the result of multigrid solver
norms = [*phase1_norms, *norms]
plt.plot(range(len(norms)), norms, label='Multigrid on iteration 50')


# Phase IV: A simple jacobi solver with the same number of total iterations
#           to compare with multigrid solver
nm_solver = NumericalSolver(A, b, n=100, x=phase1_x)
x, norms = nm_solver.jacobi()
norms = [*phase1_norms, *norms]

# Plotting jacobi solver results
plt.plot(range(len(norms)), norms, label='No Multigrid')

# Setting of matplotlob
plt.legend()
plt.xlabel('Iteration')
plt.ylabel('L2 norm of residual')
plt.title('Multigrid V-Cycle vs No Multigrid')

# Saving the graph
plt.savefig('result.png')
