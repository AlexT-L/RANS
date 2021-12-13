import numpy as np
from numpy.core.numeric import Infinity
from Field import Field
from Input import Input
from flo103_PostProcessor import flo103_PostProcessor
from flo103_ConvergenceChecker import flo103_ConvergenceChecker
from ImplicitEuler import ImplicitEuler
from NS_AirfoilBC import NS_AirfoilBC
from AirfoilMap import AirfoilMap
from CellCenterWS import CellCenterWS
from NavierStokes import NavierStokes
from MultiGrid import MultiGrid
import Contractinator as con

if __name__ == '__main__':

    # # Testing fields
    # field = Field([12,8],4)
    # field2 = Field([6,4],4)
    # weights = Field([12,8],4)
    # print(field.vals)

    # print(field.shape())

    # print("sum")
    # print(sum(sum(field[0:2,0:2,0])))

    # for z in range(4):
    #     field[0:12,0:8,z] = 3
    #     weights[0:12,0:8,z] = 2
    
    # print(field.vals)

    # print(field.shape())

    # con.conservative4way(field, field2, weights)

    # print(field2.vals)

    # # print(field.vals)

    # exit()

    




    # Comment later
    filename = 'rae9-s1.data'
    physicsUpdateFrequency = 1
    
    # Command line inputs: Cycle type, Integrator type

    # read in input
    input = Input(filename) # Will actually take all command line inputs

    # create geometry objects
    grid = AirfoilMap(input)
    workspace = CellCenterWS(grid)

    # create physics objects
    modelInput = input.add_dicts(input.flo_param, input.solv_param)
    bcmodel = NS_AirfoilBC(modelInput)
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
    mg.solution(state)

    # enforce cfl < 10 on first few cycles
    model.update_cfl_limit(10.0)

    while not CONVERGED:
        # update rev and rlv at specified interval
        updatePhysics = True
        if num_iterations != 0:
            if (num_iterations % physicsUpdateFrequency) != 0:
                updatePhysics = False
        if updatePhysics: 
            model.update_physics(workspace, state)

        # after the first few cycles, relax restriction on cfl
        model.update_cfl_limit(Infinity)
        
        # perform an interation of the multigrid cycle
        mg.performCycle()

        # get the new state and residuals
        mg.solution(state)
        mg.residuals(resid)

        # output the convergence
#        post.print_convergence(resid)

        # update convergence checker
        CONVERGED = watcher.is_converged(resid)
        
    
    # print results
#    post.print_solution(state)
    
    # Take solution and plot and save info