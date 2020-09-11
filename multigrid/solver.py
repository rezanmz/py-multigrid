import numpy as np


class Solver:
    def __init__(self, A, b, initial_guess=None):
        self.A = A
        self.b = b
        if initial_guess is not None:
            self.initial_guess = initial_guess
        else:
            self.initial_guess = np.zeros(A.shape[0])

    def solve(self,
              coarsening_method='beck',
              number_of_levels=1):

        from numerical_solver import NumericalSolver
        from multigrid.multigrid import Multigrid
        import multigrid as ml

        # Preparing x for multigrid cycle
        print('Pre-Smoothing...')
        smoother = NumericalSolver(
            A=self.A,
            b=self.b,
            n=ml.pre_smooth_iters,
            x=self.initial_guess,
            verbosity=True
        )
        pre_smooth_x, _ = smoother.jacobi()

        # Calculating residual
        residual = self.b - self.A.dot(pre_smooth_x)

        # Multigrid cycle
        multigrid = Multigrid(residual=residual,
                              adjacency_matrix=self.A)
        correction, prolongator = multigrid.v_cycle(
            coarsening_method=coarsening_method,
            number_of_levels=number_of_levels,
            p=ml.prolongators[multigrid.id] if multigrid.id in ml.prolongators else None
        )

        # Saving prolongator obtained in multigrid cycle for future use
        ml.prolongators[multigrid.id] = prolongator

        # Post smoothing the correction obtained from multigrid cycle
        print('Post-Smoothing...')
        smoother = NumericalSolver(
            A=self.A,
            b=self.b,
            n=ml.post_smooth_iters,
            x=correction,
            verbosity=True
        )
        post_smooth_x, _ = smoother.jacobi()

        # Adding the smoothed correction to the approximate answers in
        # fine grid and returning the result.
        return self.initial_guess + post_smooth_x
