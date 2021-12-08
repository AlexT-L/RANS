import pytest
from .. import cross_directory_testing

def test_my_func():
    the_string = cross_directory_testing.my_func()
    assert the_string == "It works!"
