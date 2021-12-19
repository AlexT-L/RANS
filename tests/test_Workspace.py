"""
Description
-----------
Tests the Workspace object to see if it works as expected

Libraries/Modules
-----------------
-pytest \n
-Workspace

Notes
-----
Runs the following tests:\n
Verify x, xc, vol can be retrieved
Verify x, xc, vol are non-zero
Verify that init_vars works as expected
then check that init_vars? exists return True
check that has_dict? exists and has_dict return false when they should
check is_finest()
"""
from numpy.core.numeric import NaN
import pytest
from bin.Workspace import Workspace
import numpy as np


def test_x():
    """
    Asserts that we can retrieve x and that it is non zero
    """
    x = Workspace.get_field('x')
    # Asserts non zero, as np.any is true when any are non-zero
    assert np.any(x)

def test_xc():
    """
    Asserts that we can retrieve xc and that it is non zero
    """
    xc = Workspace.get_field('xc')
    # Asserts non zero, as np.any is true when any are non-zero
    assert np.any(xc)

def test_vol():
    """
    Asserts that we can retrieve vol and that it is non zero
    """
    vol = Workspace.get_field('vol')
    # Asserts non zero, as np.any is true when any are non-zero
    assert np.any(vol)

def test_init_vars():
    '''
    Verify that init_vars works as expected
    then check that init_vars? exists return True
    '''
    # not sure about this yet
    a=True
    assert a

def test_finest():
    finest_check= Workspace.isFinest
    '''
    Asserts that isFinest 
    '''
    # but not always? cross check with mode?
    assert finest_check