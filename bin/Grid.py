
""" This module contains an abstract base class Grid"""
from abc import ABC, abstractmethod
import numpy as np
from bin.Field import Field

class Grid(ABC):
    """ Abstract base class, never directly instantiated

        AirfoilMap is a child class of this ABC
    """
    # get number of vertices
    @abstractmethod
    def get_size(self):
        pass

    @abstractmethod
    def get_geometry(self):
        pass