import Input, AirfoilMap, NavierStokes, NS_AirfoilBC, MultiGrid, CellCenterWS
import PostProcessor, ConvergenceChecker
from bin.ImplicitEuler import ImplicitEuler

if __name__ == '__main__':
    # Comment later
    filename = 'data.dat'
    TOLERANCE = 0.001
    physicsUpdateFrequency = 1
    
    # Command line inputs: Cycle type, Integrator type
    input = Input(filename) # Will actually take all command line inputs
    grid = AirfoilMap(input.grid)
    workspace = CellCenterWS(grid)
    bcmodel = NS_AirfoilBC
    model = NavierStokes(input.model)
    state = model.init_state(workspace)
    integrator = ImplicitEuler(model, input.integrator)
    mg = MultiGrid(workspace, model, integrator, state, input.multigrid)
    watcher = ConvergenceChecker(input)
    post = PostProcessor(input)

    # Update model viscosity with current state
    
    # initialize trackers
    CONVERGED = False
    num_iterations = 0

    while not CONVERGED:
        # update rev and rlv at specified interval
        updatePhysics = True
        if num_iterations != 0:
            if (num_iterations % physicsUpdateFrequency) != 0:
                updatePhysics = False
        if updatePhysics: 
            model.update_physics(workspace, state)
        
        # perform an interation of the multigrid cycle
        mg.performCycle()

        # get the residuals
        resid = mg.residuals()

        # output the convergence
        post.print_convergence(resid)

        # update convergence checker
        CONVERGED = watcher.is_converged(resid)
    

    # retrieve final solution and print results
    state = mg.solution()
    post.print_solution(state)
    
    # Take solution and plot and save info