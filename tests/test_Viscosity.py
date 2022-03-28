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


def test_ev_validation():    
    # set up environment
    [model, ws, state, [nx, ny]] = get_environment()
    
    # create faux fields
    ev = Field.create((nx,ny))
    
    # perform test of python script
    model.test(ws,state,ev,'ev')

    # compare with fortran
    if UPDATE_FORTRAN_DATA:
        ev_fortran = Field.create(ev.shape)
        model.test(ws,state,ev_fortran,'ev','fortran')
        save('tests/validation/ev', ev_fortran)
    TOL = 1e-8
    ev_fortan = load('tests/validation/ev')

    print ("max(ev_fortran) = "+str(max(abs(ev_fortan))))
    print ("max(ev) = "+str(max(abs(ev))))
    print ("mean(ev_fortran-ev) = "+str(mean(ev_fortan-ev)))
    print ("min(ev_fortran-ev) = "+str(min(abs(ev_fortan-ev))))
    assert max(abs(ev_fortan - ev)) < TOL
    