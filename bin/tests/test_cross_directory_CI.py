import pytest
from .. import cross_directory_CI

def test_my_func():
    the_string = cross_directory_CI.my_func()
    assert the_string == "It works!"
