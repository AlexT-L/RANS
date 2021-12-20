"""
Tests the NavierStokes object to make sure it behaves as expected

Libraries/Modules:
    pytest 
    Input
    NavierStokes
    AirfoilMap
    CellCenterWS
    NS_Airfoil

Notes:
    Runs the following tests:
    1. Checks that the constructor works


"""
from bin.Input import Input
from bin.NavierStokes import NavierStokes
from bin.AirfoilMap import AirfoilMap
from bin.CellCenterWS import CellCenterWS
from bin.NS_Airfoil import NS_Airfoil

# read in input
filename = 'rae9-s1.data'
input = Input(filename)

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


def test_constructor():
    """
    Asserts that we can create a NavierStokes object
    
    """
    # Assert correct type
    model = NavierStokes(bcmodel, modelInput)
    assert type(model) is NavierStokes


