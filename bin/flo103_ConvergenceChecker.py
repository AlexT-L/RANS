"""This module runs a set number of cycles. 
    Libraries/Modules:    
        Would use: numpy
        """

class flo103_ConvergenceChecker:
    """ Current implementation is to run a pre-determined 
        number of cycles, and decide that convergence is when that number
        of cycles has been reached. 
        A future improvement would be to actually implement a way to check
        the convergence of the solution, based on residuals or errors, and
        continue runnning until that convergence creiteria is below a 
        certain threshold. 

    Attributes:
        self.num_cycles: Cycle number
        self.n_runs: Number of runs
    Notes: 
        Later implementation would have more sophisticated convergence."""
    def __init__(self, input):
        """
        Sets the number of cycles and number of runs. 
        """
        self.num_cycles = 0
        self.n_runs = 5

    # monitor convergence
    def is_converged(self, residuals):
        """
        Cycle number is increased by 1 on each call. 
        Compared the cycle number to number of runs. 
        Returns true if all cycles have completed. 
        Ideally, would have a way to actually monitor convergence. 
        Then would return true when a certain criterion is met. 
        """
        self.num_cycles += 1
        print(self.num_cycles)
        return self.num_cycles >= self.n_runs