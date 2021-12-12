from numpy.core.numeric import Infinity
import Input, AirfoilMap, NavierStokes, MultiGrid, CellCenterWS
from bin.Field import Field
import flo103_PostProcessor, flo103_ConvergenceChecker
from bin.ImplicitEuler import ImplicitEuler
from bin.NS_AirfoilBC import NS_AirfoilBC

if __name__ == '__main__':
    # Comment later
    filename = 'data.dat'
    physicsUpdateFrequency = 1
    
    # Command line inputs: Cycle type, Integrator type

    # read in input
    input = Input(filename) # Will actually take all command line inputs

    # create geometry objects
    grid = AirfoilMap(input.grid)
    workspace = CellCenterWS(grid)

    # create physics objects
    modelInput = input.add_dicts(input.flo_param, input.solv_param)
    bcmodel = NS_AirfoilBC(input)
    model = NavierStokes(bcmodel, modelInput)
    integrator = ImplicitEuler(model, input.solv_param)

    # create multigrid cycle objects
    mg = MultiGrid(workspace, model, integrator, input.solv_param)
    watcher = flo103_ConvergenceChecker(input)
    post = flo103_PostProcessor(input)

    # initialize trackers
    CONVERGED = False
    num_iterations = 0
    
    # create fields for tracking state and residuals
    field_size = workspace.field_size()
    stateDim = model.dim()
    state = Field(field_size, stateDim)
    resid = Field(field_size, stateDim)

    # get initial state
    mg.get_solution(state)

    # enforce cfl < 10 on first few cycles
    integrator.update_cfl_limit(10)

    while not CONVERGED:
        # update rev and rlv at specified interval
        updatePhysics = True
        if num_iterations != 0:
            if (num_iterations % physicsUpdateFrequency) != 0:
                updatePhysics = False
        if updatePhysics: 
            model.update_physics(workspace, state)

        # after the first few cycles, relax restriction on cfl
        integrator.update_cfl_limit(Infinity)
        
        # perform an interation of the multigrid cycle
        mg.performCycle()

        # get the new state and residuals
        mg.get_solution(state)
        mg.get_residuals(resid)

        # output the convergence
#        post.print_convergence(resid)

        # update convergence checker
        CONVERGED = watcher.is_converged(resid)
        
    

    # print results
#    post.print_solution(state)
    
    # Take solution and plot and save info