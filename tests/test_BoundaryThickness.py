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


def test_ynot_validation():    
    # set up environment
    [model, ws, state, [nx, ny]] = get_environment()
    
    # create faux fields
    ynot = Field.create(nx)
    
    # perform test of python script
    model.test(ws,state,ynot,'ynot')

    # compare with fortran
    if UPDATE_FORTRAN_DATA:
        ynot_fortran = Field.create(ynot.shape)
        model.test(ws,state,ynot_fortran,'ynot','fortran')
        save('tests/validation/ynot', ynot_fortran)
    TOL = 1e-5
    ynot_fortan = load('tests/validation/ynot')

    print ("max(ynot_fortran) = "+str(max(ynot_fortan)))
    print ("max(ynot) = "+str(max(ynot)))
    print ("mean(ynot_fortran-ynot) = "+str(mean(ynot_fortan-ynot)))
    print ("min(ynot_fortran-ynot) = "+str(min(abs(ynot_fortan-ynot))))
    assert max(abs(ynot_fortan - ynot)) < TOL


def test_dsti_validation():
    # set up environment
    [model, ws, state, [nx, ny]] = get_environment()
    
    # create faux fields
    dsti = Field.create(nx)
    
    # perform test of python script
    model.test(ws,state,dsti,'dsti')

    # compare with fortran
    if UPDATE_FORTRAN_DATA:
        dsti_fortran = Field.create(dsti.shape)
        model.test(ws,state,dsti_fortran,'dsti','fortran')
        save('tests/validation/dsti', dsti_fortran)
    TOL = 1e-8
    dsti_fortan = load('tests/validation/dsti')

    print ("max(dsti_fortran) = "+str(max(abs(dsti_fortan))))
    print ("max(dsti) = "+str(max(abs(dsti))))
    print ("mean(dsti_fortran-dsti) = "+str(mean(dsti_fortan-dsti)))
    print ("min(dsti_fortran-dsti) = "+str(min(abs(dsti_fortan-dsti))))
    assert max(abs(dsti_fortan - dsti)) < TOL
    