from model_funcs.eflux_wrap import eflux
import numpy as np
from Field import Field

def test_flux():
    assert eflux() == 0

