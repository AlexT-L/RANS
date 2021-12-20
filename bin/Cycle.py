"""
Cycle.py 

Description
-----------
Class describing multigrid cycle. \n
Contains information on cycle shape and depth.

Libraries/Modules
-----------------
None.

Notes
-----
Could be expanded to include standard options like "V" and "W" cycle

Author(s)
---------
Satya Butler, Nick Conlin, Vedin Dewan, Andy Rothstein, Alex Taylor-Lash, and Brian Wynne. \n

"""

from abc import ABC
import numpy as np

class Cycle(ABC):
    
    # Constructor
    def __init__(self, input):
        self.grids = []
        self.pattern = [-1,-1,1,-1,1,1]
        self.levels = 3

    # array of directions to follow
    def path(self):
        return np.copy(self.pattern)

    # depth of cycle
    def depth(self):
        return int(self.levels)