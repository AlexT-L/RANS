import Input
import SqrtGrid
import Multicycle

if __name__ == '__main__':
    # Comment later
    filename = 'data.dat'
    input = Input(filename)
    grid = SqrtGrid(input)
    solution, convergence_info = Multicycle(grid, input)
    
    # Take solution and plot and save info