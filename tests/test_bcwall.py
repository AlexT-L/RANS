"""This module tests the python version of bcwall against the fortran version

    Libraries/Modules:
        pytest\n
        validation\n
    
        """

import sys
sys.path.append("../../RANS/bin")
import pytest
from tests.validation.validation import run_test


def test_bcwall_validation():    
    run_test('bcwall')
    