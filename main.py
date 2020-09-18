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
# Default values are:
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
# One V-Cycle
solver = Solver(A, b, initial_guess=phase1_x)
x_single_cycle = solver.solve(number_of_levels=10)
# Multiple V-Cycles
x_multiple_cycles = phase1_x
for _ in range(10):
    solver = Solver(A, b, initial_guess=x_multiple_cycles)
    x_multiple_cycles = solver.solve(number_of_levels=10)
    x_multiple_cycles, _ = NumericalSolver(
        A, b, n=5, x=x_multiple_cycles).jacobi()

# Phase III: More jacobi iterations, so we could see multigrid process
#            more clearly.
nm_solver = NumericalSolver(A, b, n=100, x=x_single_cycle)
x, norms_single_v_cycle = nm_solver.jacobi()

nm_solver = NumericalSolver(A, b, n=100, x=x_multiple_cycles)
x, norms_multiple_v_cycles = nm_solver.jacobi()

# Plotting the result of multigrid solver
norms_single_cycle = [*phase1_norms, *norms_single_v_cycle]
plt.plot(range(len(norms_single_cycle)), norms_single_cycle,
         label='Multigrid on iteration 50 - One V-Cycle')

norms_single_cycle = [*phase1_norms, *norms_multiple_v_cycles]
plt.plot(range(len(norms_single_cycle)), norms_single_cycle,
         label='Multigrid on iteration 50 - Ten V-Cycles')


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
