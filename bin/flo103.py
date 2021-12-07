import Input, SqrtGrid, NavierStokes, Workspace, MultiGrid

if __name__ == '__main__':
    # Comment later
    filename = 'data.dat'
    TOLERANCE = 0.001
    
    input = Input(filename) # Will actually take all command line inputs
    grid = SqrtGrid(input)
    model = NavierStokes(input.flo_params)
    workspace = Workspace(grid, model, input)
    mg = MultiGrid(model, workspace, input)
    
    while mg.res < TOLERANCE:
        mg.step()
    
    sol = mg.workspace
    
    # Take solution and plot and save info