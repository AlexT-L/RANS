import pytest
from bin.model_funcs.dummy import dummy_function

def test_dumb():
    assert dummy_function() == 'We are all huge dummies'