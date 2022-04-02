"""This module tests the python version of viscf against the fortran version

    Libraries/Modules:
        pytest\n
        validation\n
    
        """

import sys
sys.path.append("../../RANS/bin")

import pytest
from tests.validation.validation import run_test


def test_radi_validation():    
    run_test("radi")

'''

def test_radj_validation():    
    run_test("radj")
    
def test_rfli_validation():    
    run_test("rfli")

def test_rflj_validation():    
    run_test("rflj")

def test_dtl_validation():    
    run_test("dtl")

def test_dtl_validation():    
    run_test("dtlc")

'''