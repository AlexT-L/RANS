"""This module tests the python version of eflux against the fortran version

    Libraries/Modules:
        pytest\n
        Field\n
        Validation\n
    
        """

import sys
sys.path.append("../../RANS/bin")

import pytest
from tests.validation.validation import run_test


def test_ev_validation():    
    run_test("ev")

def test_lv_validation():    
    run_test("lv")