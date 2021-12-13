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
Satya, Joey Salads, V, Alex, My sexuality is rich people/dank meme Andy, and Brian. \n
Created on 12/11/2021. \n
Last modified on 12/12/2021.

"""
import pytest
from .. import fruit

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
    this_is_bananas = fruit.bananas("cantaloupe")
    assert this_is_bananas == "bananas"
