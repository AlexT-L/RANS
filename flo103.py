"""Flo103

    Solves the Euler equations for an airfoil using a multigrid cycle.
    Method of lines integration is used to solve the Partial Differential
    Euler equations. To speed up convergence, the solution is calculated on
    the desired mesh size, and then a new solution is found on successively 
    smaller meshes using the solution at the previous mesh refinement as a 
    guess at the state. The solutions found on coarser meshes are then used
    as a correction to the state on the coarser mesh, and a new solution is found
    on the fine mesh after applying the corrections.

    Libraries/Modules:
        Input\n
        Field\n
        AirfoilMap\n
        CellCenterWS\n
        NavierStokes\n
        ImplicitEuler\n
        MultiGrid\n

    Notes:
        Currently in development

    Authors:
        Satya Butler, Nick Conlin, Vedin Dewan, Andy Rothstein, Alex Taylor-Lash, and Brian Wynne. \n
        """
from bin.Field import Field, max, mean
from bin.Input import Input
from bin.flo103_PostProcessor import flo103_PostProcessor
from bin.flo103_ConvergenceChecker import flo103_ConvergenceChecker
from bin.ImplicitEuler import ImplicitEuler
from bin.NS_Airfoil import NS_Airfoil
from bin.AirfoilMap import AirfoilMap
from bin.CellCenterWS import CellCenterWS
from bin.NavierStokes import NavierStokes
from bin.MultiGrid import MultiGrid
from time import time

if __name__ == '__main__':

    # Comment later
    filename = 'rae9-s1.data'
    filename = 'rae9e-s3.data'
    physicsUpdateFrequency = 1
    
    # Command line inputs: Cycle type, Integrator type

    # read in input
    input = Input(filename) # Will actually take all command line inputs

    # format input
    input.geo_param["inflation_layer"] = (input.flo_param["kvis"] != 0)
    gridInput = input.add_dicts(input.geo_param, input.in_var)
    grid_dim = [input.dims['nx'], input.dims['ny']]
    modelInput = input.add_dicts(input.flo_param, input.solv_param)

    # create geometry objects
    grid = AirfoilMap.from_file(grid_dim, gridInput)
    workspace = CellCenterWS(grid)

    # create physics objects
    bcmodel = NS_Airfoil(modelInput)
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
    [nx, ny] = workspace.field_size()
    stateDim = model.dim()
    shape = (nx, ny, stateDim)
    state = Field(shape)
    resid = Field(shape)

    # get initial state
    mg.solution(state)

    # enforce cfl < 10 on first few cycles
    model.update_cfl_limit(10.0)

    start = time()
    while not CONVERGED:
        # update ev and lv at specified interval
        updatePhysics = True
        if num_iterations != 0:
            if (num_iterations % physicsUpdateFrequency) != 0:
                updatePhysics = False
        if updatePhysics: 
            model.update_physics(workspace, state)

        # after the first few cycles, relax restriction on cfl
        model.update_cfl_limit()
        
        # perform an interation of the multigrid cycle
        mg.performCycle()

        # get the new state and residuals
        mg.solution(state)
        mg.residuals(resid)

        # output the convergence
        #post.print_convergence(resid)

        # update convergence checker
        CONVERGED = watcher.is_converged(resid)
    stop = time()
    
    # print results
    #post.print_solution(state)

    rho = state[:,:,0]
    print(max(rho))
    print(mean(rho))
    print(resid)
    
    print("Total time: ", stop-start, " s")

    # Take solution and plot and save info