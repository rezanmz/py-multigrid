class Multigrid:
    _total_levels = 0

    # Constructor
    def __init__(self, residual, adjacency_matrix):
        # With each creation of a Multigrid instance, _total_levels is
        # increased in the constructor of the class.
        Multigrid._total_levels += 1
        # We keep the current level number in 'id' attribute of the object
        # As we go deeper in the cycle, the we have a higher 'id'
        self.id = Multigrid._total_levels
        # In a basic multigrid solver, the residual and
        # the adjacency matrix will be passed to multigrid
        # cycle from a finer level
        self.residual = residual
        self.adjacency_matrix = adjacency_matrix
        self.number_of_equations = len(residual)

    # Destructor
    def __del__(self):
        """
        In the destructor, as we are done with the current corse level,
        we will go to a finer level and decrease the total level by 1
        """
        Multigrid._total_levels -= 1

    def v_cycle(self, coarsening_method='beck', number_of_levels=1, p=None):
        from multigrid.fine2coarse import Fine2Coarse
        from multigrid.prolongator import Prolongator
        from scipy.sparse import linalg as scipy_linalg
        from multigrid.solver import Solver

        print('Current Level:', self.id)
        print('Number of unknowns:', self.number_of_equations)
        # Check if prolongator is already provided in input
        if p is None:
            if coarsening_method == 'beck':
                fine2coarse = Fine2Coarse(
                    adjacency_matrix=self.adjacency_matrix)
                fine, coarse = fine2coarse.beck()
                prolongator = Prolongator(
                    adjacency_matrix=self.adjacency_matrix,
                    fine=fine,
                    coarse=coarse
                )
                p = prolongator.beck()
            else:
                raise SyntaxError('Coarsening method not found.')

        # Restriction operator is the transpose of prolongation operator
        r = p.T
        restricted_adjacency = r.dot(self.adjacency_matrix).dot(p)
        restricted_residual = r.dot(self.residual)

        # If the current level is the same as maximum multigrid levels
        # specified by user, solve the coarse correction directly.
        if self.id == number_of_levels:
            coarse_correction = scipy_linalg.spsolve(
                restricted_adjacency, restricted_residual)

        # If we haven't reached the desired depth(level) in V-cycle,
        # perform another multigrid V-cycle to find the coarse correction.
        else:
            solver = Solver(
                A=restricted_adjacency,
                b=restricted_residual,
            )
            coarse_correction = solver.solve(
                coarsening_method=coarsening_method,
                number_of_levels=number_of_levels
            )

        # Finally, interpolate coarse correction to the fine grid using
        # the prolongator operator ---> P.CoarseCorrection = FineCorrection
        fine_correction = p.dot(coarse_correction)
        return fine_correction, p
