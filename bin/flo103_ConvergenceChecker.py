

class flo103_ConvergenceChecker:

    def __init__(self, input):
        self.num_cycles = 0
        self.n_runs = 5

    # monitor convergence
    def is_converged(self, residuals):
        self.num_cycles += 1
        print(self.num_cycles)
        return self.num_cycles >= self.n_runs