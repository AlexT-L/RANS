"""
Example test program.

Description
-----------
Example pytest program with proper comment styles.

Libraries/Modules
-----------------
-pytest \n
-fruit (local)

Notes
-----
Here's a note.

Author(s)
---------
Satya Butler, Nick Conlin, Vedin Dewan, Andy Rothstein, Alex Taylor-Lash, and Brian Wynne. \n

"""
import pytest
from bin.fruit import bananas

def test_bananas():
    """
    Description of 'test_bananas' method goes here.

    Parameters
    ----------
    
    Returns
    -------
    :
        Nothing, but it asserts if 'this_is_bananas' is in fact the string 'bananas'.
    
    """
    this_is_bananas = bananas("cantaloupe")
    assert this_is_bananas == "bananas"
