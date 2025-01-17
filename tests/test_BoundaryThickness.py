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
from tests.validation.validation import run_test


def test_ynot_validation():    
    run_test("ynot")


def test_dsti_validation():
    run_test("dsti")