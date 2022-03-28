"""This module tests the python version of eflux against the fortran version

    Libraries/Modules:
        pytest\n
        Field\n
        NavierStokes\n
    
        """

import sys
sys.path.append("../../RANS/bin")

import pytest
from bin.Field import Field, max, min, mean, abs, save, load
from bin.NavierStokes import UPDATE_FORTRAN_DATA
from tests.validation.validation import get_environment


def test_eflux_validation():    
    # set up environment
    [model, ws, state, [nx, ny]] = get_environment()
    
    # create faux fields
    dw = Field.create((nx, ny, model.dim()))
    
    # perform test of python script
    model.test(ws,state,dw,'eflux')

    # compare with fortran
    if UPDATE_FORTRAN_DATA:
        dw_fortran = Field.create(dw.shape)
        model.test(ws,state,dw_fortran,'eflux','fortran')
        save('tests/validation/eflux', dw_fortran)
    TOL = 1e-5
    dw_fortan = load('tests/validation/eflux')

    print ("max(dw_fortran) = "+str(max(dw_fortan)))
    print ("max(dw) = "+str(max(dw)))
    print ("mean(dw_fortran-dw) = "+str(mean(dw_fortan-dw)))
    print ("min(dw_fortran-dw) = "+str(min(abs(dw_fortan-dw))))
    assert max(abs(dw_fortan - dw)) < TOL
    